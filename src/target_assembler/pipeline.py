from __future__ import annotations

import copy
import csv
import errno
import json
import logging
import os
import re
import shutil
import tempfile
from collections import defaultdict
from pathlib import Path

import yaml

from .config import choose_validation_profile, infer_layout, resolve_validation
from .errors import ExternalToolError, TargetAssemblerError
from .external import Runner, check_stream
from .io import (
    base_read_id,
    count_units,
    exact_dereplicate,
    extract_library_by_ids,
    fastq_to_fasta,
    iter_fasta,
    iter_fastq,
    normalize_library,
    parse_nhmmer_tbl,
    reverse_complement,
    sample_library,
    write_fasta,
)
from .models import Candidate, FastaRecord, HMMHit, ValidationMetrics
from .proportion import hit_nucleotides, summarise_replicates
from .validation import alignment_mode_flag, evaluate_metrics, measure_candidates, validate_candidates


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def internal_clip(pos1: int, cigar: str, ref_length: int) -> tuple[int, int]:
    """Soft/hard-clipped bases that had REFERENCE to align to, left and right.

    A clip is "external" where the read runs off an end of the reference and
    "internal" where the reference continues. Only the internal part is returned:
    min(left clip, bases of reference before the alignment) and the same on the
    right. An unknown reference length (0) makes every right clip internal.
    """
    ops = [(int(n), op) for n, op in _CIGAR.findall(cigar)]
    if not ops:
        return 0, 0
    left = 0
    for n, op in ops:
        if op in "SH":
            left += n
        else:
            break
    right = 0
    for n, op in reversed(ops):
        if op in "SH":
            right += n
        else:
            break
    ref_span = sum(n for n, op in ops if op in "MDN=X")
    before = max(0, pos1 - 1)
    after = max(0, ref_length - (pos1 - 1 + ref_span)) if ref_length else right
    return min(left, before), min(right, after)


class Pipeline:
    def __init__(self, config: dict, dry_run: bool = False):
        self.config = config
        self.layout = infer_layout(config)
        self.paired = self.layout == "paired"
        self.outdir = Path(config["project"]["outdir"])
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.logger = _make_logger(self.outdir / "pipeline.log")
        self.runner = Runner(self.logger, dry_run=dry_run)
        self.dry_run = dry_run
        temp_parent = config["project"].get("temp_dir")
        if temp_parent:
            Path(temp_parent).mkdir(parents=True, exist_ok=True)
        self.workdir = Path(tempfile.mkdtemp(prefix="target_assembler_", dir=temp_parent))
        self.mira_dir = self.workdir / "mira_work"
        self.mira_dir.mkdir(parents=True, exist_ok=True)
        self.hmm = Path(config["input"]["hmm"]) if config["input"].get("hmm") else None
        self.summary: dict = {
            "status": "RUNNING",
            "layout": self.layout,
            "config": config["_config_path"],
            "events": [],
        }
        self.clean_reads: tuple[str, str | None] | None = None
        self.sample_reads: tuple[str, str | None] | None = None
        self.target_reads: tuple[str, str | None] | None = None
        self.expanded_reads: tuple[str, str | None] | None = None
        self.target_units = 0
        self.round_used = 0

    def event(self, name: str, value=None) -> None:
        self.summary["events"].append({"name": name, "value": value})
        self.logger.info("%s: %s", name, value)

    def dependencies(self) -> list[str]:
        mode = self.config["project"]["mode"]
        programs: set[str] = set()
        inp = self.config["input"]
        if not inp.get("clean_reads1") and not inp.get("validation_reads1") and self.config["preprocessing"]["enabled"]:
            programs.add("fastp")
        if mode in {"all", "variants"}:
            programs.update({"nhmmer", "hmmstat", "makeblastdb", "blastn", "mira", "bowtie2", "bowtie2-build"})
        if mode == "validate":
            programs.update({"bowtie2", "bowtie2-build"})
        if mode in {"all", "proportion"} and self.config["proportion"]["enabled"]:
            programs.update({"nhmmer", "mafft", "hmmbuild"})
        return sorted(programs)

    def preflight(self, require_programs: bool = True) -> None:
        paths = [self.hmm] if self.hmm else []
        inp = self.config["input"]
        for key in (
            "reads1", "reads2", "clean_reads1", "clean_reads2", "sample_reads1", "sample_reads2",
            "variants_fasta", "validation_reads1", "validation_reads2",
        ):
            if inp.get(key):
                paths.append(Path(inp[key]))
        missing = [str(path) for path in paths if not path.exists()]
        if missing:
            raise TargetAssemblerError("Input files not found: " + ", ".join(missing))
        if require_programs:
            self.runner.require(self.dependencies())

    def run(self) -> dict:
        try:
            self.preflight(require_programs=not self.dry_run)
            self._write_effective_config()
            if self.dry_run:
                self.summary["status"] = "DRY_RUN_OK"
                self.summary["required_programs"] = self.dependencies()
                return self.summary

            mode = self.config["project"]["mode"]
            if mode == "validate" and self.config["input"].get("validation_reads1"):
                self.clean_reads = (
                    self.config["input"]["validation_reads1"],
                    self.config["input"].get("validation_reads2"),
                )
                self.sample_reads = None
                self.event("validation_library", "provided")
            else:
                self.clean_reads, self.sample_reads = self.prepare_libraries(need_sample=mode in {"all", "variants"})
            accepted_fasta: Path | None = None
            if mode in {"all", "variants"}:
                accepted_fasta = self.run_variants()
            elif mode == "validate":
                self.run_validation_comparison()
            if mode in {"all", "proportion"} and self.config["proportion"]["enabled"]:
                if accepted_fasta is None:
                    accepted_fasta = self.find_existing_verified()
                if accepted_fasta:
                    self.calculate_proportion(accepted_fasta)
                elif mode == "all":
                    self.event("proportion_skipped", "no candidates passed validation")
                else:
                    raise TargetAssemblerError("No verified FASTA found for proportion-only mode")
            self.summary["status"] = "SUCCESS"
            return self.summary
        except OSError as exc:
            if exc.errno == errno.ENOSPC:
                free_gb = shutil.disk_usage(self.workdir).free / (1024 ** 3)
                message = (
                    f"Temporary filesystem ran out of space in {self.workdir.parent} "
                    f"({free_gb:.2f} GiB free). Set project.temp_dir to a directory on a larger "
                    "filesystem and/or cap sampling.maximum_units."
                )
                self.summary["status"] = "FAILED"
                self.summary["error"] = message
                raise TargetAssemblerError(message) from exc
            self.summary["status"] = "FAILED"
            self.summary["error"] = str(exc)
            raise
        except Exception as exc:
            self.summary["status"] = "FAILED"
            self.summary["error"] = str(exc)
            raise
        finally:
            if self.config["runtime"]["keep_work"]:
                self.summary["workdir"] = str(self.workdir)
            self._write_summary()
            if not self.config["runtime"]["keep_work"]:
                shutil.rmtree(self.workdir, ignore_errors=True)

    def _write_effective_config(self) -> None:
        public = {key: value for key, value in self.config.items() if not key.startswith("_")}
        (self.outdir / "effective_config.yaml").write_text(
            yaml.safe_dump(public, sort_keys=False, allow_unicode=True), encoding="utf-8"
        )

    def _write_summary(self) -> None:
        try:
            (self.outdir / "pipeline_summary.json").write_text(
                json.dumps(self.summary, indent=2, ensure_ascii=False), encoding="utf-8"
            )
        except OSError:
            pass

    def prepare_libraries(self, need_sample: bool) -> tuple[tuple[str, str | None], tuple[str, str | None] | None]:
        inp = self.config["input"]
        if inp.get("clean_reads1"):
            clean = (inp["clean_reads1"], inp.get("clean_reads2"))
            self.event("clean_library", "reused")
        elif self.config["preprocessing"]["enabled"]:
            clean = self.run_fastp(inp["reads1"], inp.get("reads2"))
        else:
            clean = (inp["reads1"], inp.get("reads2"))
            self.event("clean_library", "raw input used without preprocessing")

        if self.config["outputs"]["save_intermediate_libraries"]:
            self._save_library(clean, "clean")
        if not need_sample:
            return clean, None
        if inp.get("sample_reads1"):
            sample = self.normalized_copy(
                inp["sample_reads1"], inp.get("sample_reads2"), "provided_sample"
            )
            self.event("assembly_sample", "reused")
        else:
            sample = self.make_assembly_sample(clean)
        if self.config["outputs"]["save_intermediate_libraries"]:
            self._save_library(sample, "sample")
        return clean, sample

    def run_fastp(self, reads1: str, reads2: str | None) -> tuple[str, str | None]:
        pre = self.config["preprocessing"]
        out1 = self.workdir / "clean_R1.fastq.gz"
        out2 = self.workdir / "clean_R2.fastq.gz" if reads2 else None
        command = [
            "fastp", "-i", reads1, "-o", str(out1), "--length_required", str(pre["min_length"]),
            "--thread", str(self.config["runtime"]["threads"]), "--json", str(self.outdir / "fastp.json"),
            "--html", str(self.outdir / "fastp.html"),
        ]
        if reads2 and out2:
            command.extend(["-I", reads2, "-O", str(out2)])
            if pre["detect_adapters"]:
                command.append("--detect_adapter_for_pe")
        if pre["cut_tail"]:
            command.append("--cut_tail")
        command.extend(str(value) for value in pre["extra_args"])
        self.runner.run(command, stderr_path=self.outdir / "fastp.stderr.log")
        self.event("clean_library", "fastp")
        return str(out1), str(out2) if out2 else None

    def normalized_copy(self, reads1: str, reads2: str | None, prefix: str) -> tuple[str, str | None]:
        out1 = self.workdir / f"{prefix}_R1.fastq.gz"
        out2 = self.workdir / f"{prefix}_R2.fastq.gz" if reads2 else None
        normalize_library(reads1, reads2, out1, out2)
        return str(out1), str(out2) if out2 else None

    def _save_library(self, reads: tuple[str, str | None], prefix: str) -> None:
        destination = self.outdir / "intermediate_libraries"
        destination.mkdir(parents=True, exist_ok=True)
        saved: list[str] = []
        reused = 0
        for index, source in enumerate((reads[0], reads[1]), start=1):
            if not source:
                continue
            suffix = ".fastq.gz" if str(source).endswith(".gz") else ".fastq"
            target = destination / f"{prefix}_R{index}{suffix}"
            saved.append(str(target))
            # A library supplied through input.clean_reads1/2 or
            # input.sample_reads1/2 may BE the file this would copy it to: that is
            # what happens when a previous run's intermediate_libraries/ is fed
            # back in to skip preprocessing. shutil.copyfile then raises
            #   SameFileError: '<path>' and '<path>' are the same file
            # and the run dies after "clean_library: reused" with nothing to fix,
            # because the file is already exactly where it is wanted.
            try:
                if target.exists() and os.path.samefile(source, target):
                    reused += 1
                    continue
            except OSError:
                pass
            shutil.copyfile(source, target)
        if reused:
            self.event(f"saved_{prefix}_library_already_in_place", reused)
        self.event(f"saved_{prefix}_library", saved)

    def make_assembly_sample(self, clean: tuple[str, str | None]) -> tuple[str, str | None]:
        sampling = self.config["sampling"]
        total = count_units(*clean)
        if total == 0:
            raise TargetAssemblerError("The cleaned library is empty")
        self.event("clean_library_units", total)
        mode = sampling["mode"]
        space_probe: tuple[int, int, int] | None = None
        if mode == "all":
            requested = total
        elif mode == "fixed":
            requested = int(sampling["fixed_units"])
            if requested <= 0:
                raise TargetAssemblerError("sampling.fixed_units must be positive in fixed mode")
        else:
            micro_n = min(total, int(sampling["micro_units"]))
            micro = self.sample_to_work(clean, micro_n, total, "micro", self.config["runtime"]["seed"])
            micro_fasta = self.workdir / "micro.fasta"
            fastq_to_fasta(*micro, micro_fasta)
            table = self.workdir / "micro.tbl"
            self.run_nhmmer(micro_fasta, table)
            found = len({base_read_id(hit.target) for hit in self.filter_hmm_hits(parse_nhmmer_tbl(table))})
            self.event("micro_hmm_hit_units", found)
            if found:
                requested = round(micro_n * int(sampling["target_hmm_hits"]) / found)
            else:
                requested = int(sampling["zero_hit_fallback_units"])
            requested = max(int(sampling["minimum_units"]), requested)
            micro_fastq_bytes = sum(Path(path).stat().st_size for path in micro if path)
            space_probe = (micro_n, micro_fastq_bytes, micro_fasta.stat().st_size)
        if sampling.get("maximum_units"):
            requested = min(requested, int(sampling["maximum_units"]))
        requested = min(total, requested)
        if space_probe:
            self.check_projected_temp_space(requested, *space_probe)
            if not self.config["runtime"]["keep_work"]:
                self.remove_paths(
                    *(Path(path) for path in micro if path),
                    micro_fasta,
                    table,
                    Path(str(table) + ".stdout"),
                    Path(str(table) + ".stderr"),
                )
        sample = self.sample_to_work(clean, requested, total, "sample", self.config["runtime"]["seed"])
        self.event("assembly_sample_units", requested)
        return sample

    def sample_to_work(
        self,
        reads: tuple[str, str | None],
        number: int,
        total: int,
        prefix: str,
        seed: int,
    ) -> tuple[str, str | None]:
        out1 = self.workdir / f"{prefix}_R1.fastq.gz"
        out2 = self.workdir / f"{prefix}_R2.fastq.gz" if reads[1] else None
        sample_library(reads[0], reads[1], out1, out2, number, total, seed, normalize_names=True)
        return str(out1), str(out2) if out2 else None

    def check_projected_temp_space(
        self,
        requested_units: int,
        probe_units: int,
        probe_fastq_bytes: int,
        probe_fasta_bytes: int,
    ) -> None:
        if not self.config["runtime"].get("temp_space_check", True) or probe_units <= 0:
            return
        scale = requested_units / probe_units
        # Peak additional data: sampled compressed FASTQ, search FASTA, and BLAST database.
        base_estimate = scale * (probe_fastq_bytes + (1.5 * probe_fasta_bytes))
        safety = float(self.config["runtime"].get("temp_space_safety_factor", 1.20))
        projected = int(base_estimate * safety)
        margin = int(float(self.config["runtime"].get("temp_space_margin_gb", 1.0)) * (1024 ** 3))
        free = shutil.disk_usage(self.workdir).free
        self.event(
            "temp_space_estimate",
            {
                "filesystem": str(self.workdir.parent),
                "projected_additional_gib": round(projected / (1024 ** 3), 2),
                "free_gib": round(free / (1024 ** 3), 2),
                "required_margin_gib": round(margin / (1024 ** 3), 2),
            },
        )
        if free < projected + margin:
            raise TargetAssemblerError(
                f"Insufficient temporary space in {self.workdir.parent}: approximately "
                f"{projected / (1024 ** 3):.2f} GiB additional workspace plus "
                f"{margin / (1024 ** 3):.2f} GiB safety margin is required, but only "
                f"{free / (1024 ** 3):.2f} GiB is free. Set project.temp_dir to a larger "
                "filesystem or reduce sampling.maximum_units."
            )

    @staticmethod
    def remove_paths(*paths: Path) -> None:
        for path in paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    def hmm_thresholds(self) -> tuple[float | None, float | None]:
        """selection.hmm_evalue (-E) and selection.hmm_score_min (-T), or (None, None)."""
        selection = self.config.get("selection", {})
        evalue = selection.get("hmm_evalue")
        score = selection.get("hmm_score_min")
        return (None if evalue is None else float(evalue),
                None if score is None else float(score))

    def filter_hmm_hits(self, hits: list[HMMHit]) -> list[HMMHit]:
        """Re-apply the configured thresholds after parsing.

        nhmmer already applied them, so this normally changes nothing. It matters
        when a table on disk was written by an earlier run with different
        settings, and it makes the enforcement point explicit rather than
        implicit in a command line built elsewhere.
        """
        evalue, score = self.hmm_thresholds()
        if evalue is None and score is None:
            return hits
        kept = [
            hit for hit in hits
            if (evalue is None or hit.evalue <= evalue)
            and (score is None or hit.score >= score)
        ]
        if len(kept) != len(hits):
            self.event("hmm_hits_below_threshold", len(hits) - len(kept))
        return kept

    def run_nhmmer(self, fasta: str | Path, table: str | Path, hmm: str | Path | None = None) -> None:
        chosen_hmm = hmm or self.hmm
        if chosen_hmm is None:
            raise TargetAssemblerError("No HMM was supplied for nhmmer")
        command = [
            "nhmmer", "--tblout", str(table), "--cpu", str(self.config["runtime"]["threads"]),
        ]
        # -T and -E are mutually exclusive in HMMER: a score threshold replaces
        # the E-value model entirely, so -T wins when both are set and nothing
        # silently applies two incompatible cutoffs.
        evalue, score = self.hmm_thresholds()
        if score is not None:
            command += ["-T", repr(score), "--incT", repr(score)]
        elif evalue is not None:
            command += ["-E", repr(evalue), "--incE", repr(evalue)]
        command += [str(chosen_hmm), str(fasta)]
        self.runner.run(
            command,
            stdout_path=str(table) + ".stdout",
            stderr_path=str(table) + ".stderr",
        )

    def run_variants(self) -> Path | None:
        assert self.clean_reads and self.sample_reads
        self.target_reads, self.target_units = self.select_target_reads(self.sample_reads)
        if self.config["outputs"]["save_intermediate_libraries"]:
            self._save_library(self.target_reads, "target")
        if not self.config["runtime"]["keep_work"]:
            self.remove_paths(*(Path(path) for path in self.sample_reads if path))
        if self.target_units == 0:
            raise TargetAssemblerError("No target-associated reads were selected")
        round1 = self.run_mira(self.target_reads, 1)
        contigs1 = self.filter_contigs(round1, 1)
        final_contigs = contigs1
        self.round_used = 1

        if self.config["assembly"]["rounds"] == 2:
            self.expanded_reads = self.expand_reads(contigs1, self.clean_reads)
            if self.expanded_reads and self.config["outputs"]["save_intermediate_libraries"]:
                self._save_library(self.expanded_reads, "expanded")
            round2_input = self.expanded_reads or self.target_reads
            try:
                round2 = self.run_mira(round2_input, 2)
                final_contigs = self.filter_contigs(round2, 2)
                self.round_used = 2
            except (ExternalToolError, TargetAssemblerError) as exc:
                if not self.config["assembly"]["fallback_to_round1"]:
                    raise
                self.event("round2_fallback", str(exc))

        candidates, candidates_fasta = self.extract_candidates(final_contigs, self.round_used)
        if not candidates:
            raise TargetAssemblerError("No complete repeat candidates were extracted")
        return self.validate(candidates, candidates_fasta)

    def select_target_reads(self, sample: tuple[str, str | None]) -> tuple[tuple[str, str | None], int]:
        # Internally sampled/provided libraries are already normalized. Avoid making a
        # second full uncompressed copy before the FASTA/BLAST stage.
        normalized = sample
        all_fasta = self.workdir / "sample_all.fasta"
        fastq_to_fasta(*normalized, all_fasta)
        table = self.workdir / "sample_hits.tbl"
        self.run_nhmmer(all_fasta, table)
        hits = self.filter_hmm_hits(parse_nhmmer_tbl(table))
        hit_names = {hit.target for hit in hits}
        selected_bases = {base_read_id(name) for name in hit_names}
        self.event("nhmmer_hit_reads", len(hit_names))

        hit_sequences: dict[str, str] = {}
        for path in normalized:
            if path:
                for record in iter_fastq(path):
                    token = record.name.split()[0]
                    if token in hit_names:
                        hit_sequences[token] = record.sequence
        selection = self.config["selection"]
        flank_records: list[FastaRecord] = []
        flank_size = int(selection["flank_bp"])
        minimum_flank = int(selection["minimum_flank_length"])
        for index, hit in enumerate(hits, start=1):
            sequence = hit_sequences.get(hit.target)
            if not sequence:
                continue
            left = sequence[max(0, hit.start0 - flank_size) : hit.start0]
            right = sequence[hit.end0 : min(len(sequence), hit.end0 + flank_size)]
            if len(left) >= minimum_flank:
                flank_records.append(FastaRecord(f"flank_{index}_L", left))
            if len(right) >= minimum_flank:
                flank_records.append(FastaRecord(f"flank_{index}_R", right))
        flank_records = [item[0] for item in exact_dereplicate(flank_records)]
        self.event("unique_flank_queries", len(flank_records))
        if flank_records:
            query = self.workdir / "flank_queries.fasta"
            write_fasta(flank_records, query)
            database = self.workdir / "sample_db"
            self.runner.run(
                ["makeblastdb", "-in", str(all_fasta), "-dbtype", "nucl", "-out", str(database)],
                stdout_path=self.workdir / "makeblastdb.out",
                stderr_path=self.workdir / "makeblastdb.err",
            )
            blast_out = self.workdir / "flank_blast.tsv"
            self.runner.run(
                [
                    "blastn", "-query", str(query), "-db", str(database),
                    "-outfmt", "6 qseqid sseqid pident length qlen evalue bitscore",
                    "-max_target_seqs", str(selection["blast_max_targets"]),
                    "-evalue", str(selection["blast_evalue"]),
                    "-perc_identity", str(selection["blast_identity"]),
                    "-num_threads", str(self.config["runtime"]["threads"]),
                    "-out", str(blast_out),
                ],
                stderr_path=self.workdir / "blastn.err",
            )
            blast_names: set[str] = set()
            with open(blast_out, encoding="utf-8") as handle:
                for line in handle:
                    columns = line.rstrip("\n").split("\t")
                    if len(columns) < 7:
                        continue
                    query_coverage = 100.0 * int(columns[3]) / max(1, int(columns[4]))
                    if query_coverage >= float(selection["blast_min_query_coverage"]):
                        blast_names.add(columns[1])
            selected_bases.update(base_read_id(name) for name in blast_names)
            self.event("blast_associated_reads", len(blast_names))

        out1 = self.workdir / "target_R1.fastq.gz"
        out2 = self.workdir / "target_R2.fastq.gz" if normalized[1] else None
        units = extract_library_by_ids(
            normalized[0], normalized[1], selected_bases, out1, out2, normalize_names=True
        )
        self.event("target_units", units)
        if not self.config["runtime"]["keep_work"]:
            blast_files = list(self.workdir.glob("sample_db*"))
            self.remove_paths(
                all_fasta,
                self.workdir / "sample_hits.tbl",
                self.workdir / "sample_hits.tbl.stdout",
                self.workdir / "sample_hits.tbl.stderr",
                *blast_files,
            )
        return (str(out1), str(out2) if out2 else None), units

    def mira_threads(self) -> int:
        """MIRA's thread count: assembly.mira_threads if set, else runtime.threads.
        Appended LAST, so a -GE:not inside assembly.parameters never takes effect."""
        value = self.config["assembly"].get("mira_threads")
        return int(value) if value is not None else int(self.config["runtime"]["threads"])

    def run_mira(self, reads: tuple[str, str | None], round_number: int) -> Path:
        total = count_units(*reads)
        limit = int(self.config["assembly"]["mira_max_units"] or 0)
        if limit and total > limit:
            reads = self.sample_for_mira(reads, limit, total, round_number)
            self.event(f"mira_round{round_number}_downsample", {"from": total, "to": limit})
        project = _safe_project_name(self.config["project"]["name"]) + f"_r{round_number}"
        r1 = self.mira_dir / f"{project}_R1.fastq"
        r2 = self.mira_dir / f"{project}_R2.fastq" if reads[1] else None
        normalize_library(reads[0], reads[1], r1, r2)
        manifest = self.mira_dir / f"{project}.manifest.conf"
        assembly = self.config["assembly"]
        parameters = [str(item) for item in assembly["parameters"]]
        parameters.append(f"-GE:not={self.mira_threads()}")
        lines = [
            f"project = {project}",
            f"job = {assembly['job']}",
            f"parameters = {' '.join(parameters)}",
            f"readgroup = reads_r{round_number}",
            f"data = {r1.name}" + (f" {r2.name}" if r2 else ""),
        ]
        if r2:
            lines.extend(
                [
                    f"template_size = {assembly['insert_size_min']} {assembly['insert_size_max']}",
                    f"segment_placement = {assembly['segment_placement']}",
                ]
            )
        lines.append(f"technology = {assembly['technology']}")
        manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")

        assembly_dir = self.mira_dir / f"{project}_assembly"
        if assembly_dir.exists():
            shutil.rmtree(assembly_dir)
        result = self.runner.run(
            ["mira", manifest.name], cwd=self.mira_dir,
            stdout_path=self.mira_dir / f"{project}.stdout.log",
            stderr_path=self.mira_dir / f"{project}.stderr.log",
            check=False,
        )
        fasta = assembly_dir / f"{project}_d_results" / f"{project}_out.unpadded.fasta"
        if not fasta.exists():
            raise ExternalToolError(
                f"MIRA round {round_number} produced no contig FASTA (exit {result.returncode}); "
                f"see {project}.stdout.log and {project}.stderr.log"
            )
        self.event(f"mira_round{round_number}", "success")
        return fasta

    def sample_for_mira(
        self, reads: tuple[str, str | None], number: int, total: int, round_number: int
    ) -> tuple[str, str | None]:
        out1 = self.mira_dir / f"mira_r{round_number}_sample_R1.fastq"
        out2 = self.mira_dir / f"mira_r{round_number}_sample_R2.fastq" if reads[1] else None
        sample_library(
            reads[0], reads[1], out1, out2, number, total,
            self.config["runtime"]["seed"] + round_number, normalize_names=True,
        )
        return str(out1), str(out2) if out2 else None

    def filter_contigs(self, contigs: Path, round_number: int) -> Path:
        table = self.workdir / f"contigs_r{round_number}.tbl"
        self.run_nhmmer(contigs, table)
        wanted = {hit.target for hit in self.filter_hmm_hits(parse_nhmmer_tbl(table))}
        minimum = int(self.config["assembly"]["contig_min_length"])
        records = [record for record in iter_fasta(contigs) if record.name in wanted and len(record.sequence) >= minimum]
        if not records:
            raise TargetAssemblerError(f"No round-{round_number} contigs passed HMM/length filtering")
        output = self.workdir / f"contigs_r{round_number}_filtered.fasta"
        write_fasta(records, output)
        self.event(f"contigs_round{round_number}_filtered", len(records))
        return output

    def expand_reads(
        self, reference: Path, clean: tuple[str, str | None]
    ) -> tuple[str, str | None] | None:
        if not self.config["expansion"]["enabled"]:
            return None
        expansion = self.config["expansion"]
        prefix = str(self.workdir / "expansion_index")
        self.runner.run(
            ["bowtie2-build", str(reference), prefix],
            stdout_path=self.workdir / "expansion_build.out",
            stderr_path=self.workdir / "expansion_build.err",
        )
        command = [
            "bowtie2", "-x", prefix,
            alignment_mode_flag(expansion["sensitivity"], expansion["alignment"]),
            "-k", str(expansion["max_alignments"]), "-p", str(self.config["runtime"]["threads"]),
            "--no-unal",
        ]
        if clean[1]:
            command.extend(["-1", clean[0], "-2", clean[1]])
        else:
            command.extend(["-U", clean[0]])
        process, stderr_handle = self.runner.popen_stdout(
            command, stderr_path=self.workdir / "expansion_bowtie2.log"
        )
        wanted: set[str] = set()
        # expansion.max_internal_clip: a local alignment may soft-clip a read only
        # where the read runs off the END of the reference contig. Clipping where
        # the contig still continues means the read shares a stretch with the
        # contig and then diverges - foreign rDNA aligned to the conserved 5.8S
        # block, for instance. null keeps the upstream behaviour (no test).
        clip_limit = expansion.get("max_internal_clip")
        ref_lengths: dict[str, int] = {}
        if clip_limit is not None:
            ref_lengths = {record.name: len(record.sequence) for record in iter_fasta(reference)}
        clip_rejected = 0
        assert process.stdout is not None
        try:
            for line in process.stdout:
                if line.startswith("@"):
                    continue
                columns = line.rstrip("\n").split("\t")
                if len(columns) >= 5 and not (int(columns[1]) & 0x4):
                    if int(columns[4]) >= int(expansion["minimum_mapq"]):
                        if clip_limit is not None and len(columns) >= 6:
                            left, right = internal_clip(
                                int(columns[3]), columns[5], ref_lengths.get(columns[2], 0)
                            )
                            if max(left, right) > int(clip_limit):
                                clip_rejected += 1
                                continue
                        wanted.add(base_read_id(columns[0]))
        finally:
            process.stdout.close()
        check_stream(process, stderr_handle, "Bowtie2 library expansion")
        if clip_limit is not None:
            self.event("expansion_internal_clip_rejected", {"limit": int(clip_limit), "alignments": clip_rejected})
        if not wanted:
            self.event("expanded_units", 0)
            return None
        out1 = self.workdir / "expanded_R1.fastq.gz"
        out2 = self.workdir / "expanded_R2.fastq.gz" if clean[1] else None
        count = extract_library_by_ids(clean[0], clean[1], wanted, out1, out2, normalize_names=True)
        self.event("expanded_units", count)
        return (str(out1), str(out2) if out2 else None) if count else None

    def hmm_model_length(self) -> int:
        if self.hmm is None:
            raise TargetAssemblerError("input.hmm is required for repeat extraction")
        result_path = self.workdir / "hmmstat.txt"
        self.runner.run(["hmmstat", str(self.hmm)], stdout_path=result_path, stderr_path=self.workdir / "hmmstat.err")
        for line in result_path.read_text(encoding="utf-8").splitlines():
            if line.strip() and not line.startswith("#"):
                columns = line.split()
                if len(columns) >= 6:
                    return int(columns[5])
        raise TargetAssemblerError("Could not parse HMM model length from hmmstat")

    @staticmethod
    def tandem_downstream_hit(first: HMMHit, second: HMMHit) -> HMMHit:
        """Return the downstream hit in repeat orientation for a genomic-order pair."""
        return second if first.strand == "+" else first

    def extract_candidates(self, contigs: Path, assembly_round: int) -> tuple[list[Candidate], Path]:
        renamed: list[FastaRecord] = []
        source_names: dict[str, str] = {}
        for index, record in enumerate(iter_fasta(contigs), start=1):
            name = f"contig_r{assembly_round}_{index:05d}"
            renamed.append(FastaRecord(name, record.sequence))
            source_names[name] = record.name
        renamed_fasta = self.workdir / "contigs_renamed.fasta"
        write_fasta(renamed, renamed_fasta)
        first_table = self.workdir / "orientation_hits.tbl"
        self.run_nhmmer(renamed_fasta, first_table)
        first_hits = self.filter_hmm_hits(parse_nhmmer_tbl(first_table))
        best_hit: dict[str, HMMHit] = {}
        for hit in first_hits:
            if hit.target not in best_hit or hit.score > best_hit[hit.target].score:
                best_hit[hit.target] = hit
        oriented: list[FastaRecord] = []
        for record in renamed:
            hit = best_hit.get(record.name)
            if not hit:
                continue
            sequence = reverse_complement(record.sequence) if hit.strand == "-" else record.sequence
            oriented.append(FastaRecord(record.name, sequence))
        oriented_fasta = self.workdir / "contigs_oriented.fasta"
        write_fasta(oriented, oriented_fasta)
        table = self.workdir / "oriented_hits.tbl"
        self.run_nhmmer(oriented_fasta, table)
        hits_by_contig: dict[str, list[HMMHit]] = defaultdict(list)
        for hit in self.filter_hmm_hits(parse_nhmmer_tbl(table)):
            hits_by_contig[hit.target].append(hit)
        for hits in hits_by_contig.values():
            hits.sort(key=lambda item: item.start0)
        sequences = {record.name: record.sequence for record in oriented}
        model_length = self.hmm_model_length()
        extraction = self.config["extraction"]
        minimum_fraction = float(extraction["minimum_model_fraction"])
        raw_records: list[FastaRecord] = []
        raw_meta: dict[str, dict] = {}
        tandem_pairs = 0
        skipped_incomplete = 0
        skipped_mixed_strand = 0
        for contig, hits in hits_by_contig.items():
            sequence = sequences[contig]
            if extraction["mode"] == "tandem":
                iterator = ((index, hits[index], hits[index + 1]) for index in range(len(hits) - 1))
                for index, first, second in iterator:
                    tandem_pairs += 1
                    if first.strand != second.strand:
                        skipped_mixed_strand += 1
                        continue
                    downstream = self.tandem_downstream_hit(first, second)
                    # A tandem unit is cut from the end of one HMM hit to the end
                    # of the next.  The first (boundary) hit may therefore be
                    # partial; only the downstream hit must reach the required
                    # fraction of the model's 3-prime end.  This preserves the
                    # extraction rule used by the original pipeline.
                    if downstream.model_end_fraction(model_length) < minimum_fraction:
                        skipped_incomplete += 1
                        continue
                    if first.strand == "+":
                        start = max(0, first.end0 - int(extraction["flank5"]))
                        end = min(len(sequence), second.end0 + int(extraction["flank3"]))
                        piece = sequence[start:end]
                    else:
                        start = max(0, first.start0 - int(extraction["flank3"]))
                        end = min(len(sequence), second.start0 + int(extraction["flank5"]))
                        piece = reverse_complement(sequence[start:end])
                    name = f"{contig}_tandem_{index + 1}"
                    if end > start:
                        raw_records.append(FastaRecord(name, piece))
                        raw_meta[name] = {"contig": contig, "kind": "tandem", "strand": first.strand, "start": start, "end": end}
            else:
                for index, hit in enumerate(hits, start=1):
                    if hit.model_fraction(model_length) < minimum_fraction:
                        continue
                    start = max(0, hit.start0 - int(extraction["flank5"]))
                    end = min(len(sequence), hit.end0 + int(extraction["flank3"]))
                    piece = sequence[start:end]
                    if hit.strand == "-":
                        piece = reverse_complement(piece)
                    name = f"{contig}_site_{index}"
                    if end > start:
                        raw_records.append(FastaRecord(name, piece))
                        raw_meta[name] = {"contig": contig, "kind": "site", "strand": hit.strand, "start": start, "end": end}
        if extraction["mode"] == "tandem":
            self.event("tandem_hit_pairs", tandem_pairs)
            self.event("skipped_incomplete_coding", skipped_incomplete)
            self.event("skipped_mixed_strand_pairs", skipped_mixed_strand)
        self.event("extracted_sequences_count", len(raw_records))
        write_fasta(raw_records, self.outdir / "Repeats_extracted_raw.fasta")
        dereplicated = exact_dereplicate(raw_records)
        candidates: list[Candidate] = []
        for index, (representative, members) in enumerate(dereplicated, start=1):
            meta = raw_meta[representative.name]
            variant_id = f"variant_{index:05d}"
            candidates.append(
                Candidate(
                    variant_id=variant_id,
                    sequence=representative.sequence,
                    source_contig=source_names.get(meta["contig"], meta["contig"]),
                    extraction=meta["kind"],
                    assembly_round=assembly_round,
                    strand=meta["strand"],
                    start0=meta["start"],
                    end0=meta["end"],
                    derep_size=len(members),
                    member_ids=members,
                )
            )
        candidates_fasta = self.outdir / "Repeats_candidates.fasta"
        write_fasta(
            [
                FastaRecord(
                    item.variant_id, item.sequence,
                    f"source={item.source_contig} round={item.assembly_round} derep_size={item.derep_size}",
                )
                for item in candidates
            ],
            candidates_fasta,
        )
        self.event("repeat_candidates", len(candidates))
        return candidates, candidates_fasta

    def validate(self, candidates: list[Candidate], candidates_fasta: Path) -> Path | None:
        profile = choose_validation_profile(self.config, self.layout, self.target_units)
        validation = resolve_validation(self.config, profile)
        validation.setdefault("insert_size_min", self.config["assembly"]["insert_size_min"])
        validation.setdefault("insert_size_max", self.config["assembly"]["insert_size_max"])
        self.summary["validation_profile"] = profile
        self.summary["validation_parameters"] = validation
        source_name = validation["source"]
        sources = {
            "clean": self.clean_reads,
            "expanded": self.expanded_reads,
            "assembly": self.target_reads,
        }
        source = sources.get(source_name)
        if source is None:
            self.event("validation_source_fallback", {"requested": source_name, "used": "assembly"})
            source = self.target_reads
            source_name = "assembly"
        assert source is not None
        total = count_units(*source)
        limit = validation.get("max_units")
        if limit and total > int(limit):
            source = self.sample_to_work(
                source, int(limit), total, "validation",
                self.config["runtime"]["seed"] + 500,
            )
            used = int(limit)
        else:
            source = self.normalized_copy(*source, "validation")
            used = total
        self.event("validation_library", {"source": source_name, "units": used, "available": total})
        metrics = validate_candidates(
            self.runner, str(candidates_fasta), candidates, source[0], source[1], validation,
            int(self.config["runtime"]["threads"]), self.workdir,
        )
        self.write_validation_outputs(candidates, metrics)
        accepted = {item.variant_id for item in metrics if item.decision == "PASS"}
        verified = self.outdir / "Repeats_verified.fasta"
        write_fasta(
            [
                FastaRecord(candidate.variant_id, candidate.sequence, f"profile={profile} source={candidate.source_contig}")
                for candidate in candidates if candidate.variant_id in accepted
            ],
            verified,
        )
        rejected = self.outdir / "Repeats_rejected.fasta"
        write_fasta(
            [FastaRecord(candidate.variant_id, candidate.sequence) for candidate in candidates if candidate.variant_id not in accepted],
            rejected,
        )
        self.event("verified_candidates", len(accepted))
        return verified if accepted else None

    def run_validation_comparison(self) -> None:
        """Measure one candidate/library combination once and apply several profiles."""
        assert self.clean_reads
        configured = self.config["input"].get("variants_fasta")
        candidates_fasta = Path(configured) if configured else self.outdir / "Repeats_candidates.fasta"
        if not candidates_fasta.exists():
            raise TargetAssemblerError(
                "No candidate FASTA found for validate mode; set input.variants_fasta or "
                "place Repeats_candidates.fasta in project.outdir"
            )
        candidates = self.load_external_candidates(candidates_fasta)
        if not candidates:
            raise TargetAssemblerError(f"Candidate FASTA is empty: {candidates_fasta}")

        total = count_units(*self.clean_reads)
        limit = self.config["validation"].get("comparison_max_units")
        if limit and total > int(limit):
            source = self.sample_to_work(
                self.clean_reads,
                int(limit),
                total,
                "validation_comparison",
                self.config["runtime"]["seed"] + 750,
            )
            used = int(limit)
        else:
            source = self.clean_reads
            used = total
        self.event("validation_comparison_library", {"units": used, "available": total})

        measurement = resolve_validation(self.config, "balanced")
        measurement["insert_size_min"] = self.config["assembly"]["insert_size_min"]
        measurement["insert_size_max"] = self.config["assembly"]["insert_size_max"]
        measured = measure_candidates(
            self.runner,
            str(candidates_fasta),
            candidates,
            source[0],
            source[1],
            measurement,
            int(self.config["runtime"]["threads"]),
            self.workdir,
        )

        profiles = list(dict.fromkeys(self.config["validation"]["compare_profiles"]))
        evaluated: dict[str, list[ValidationMetrics]] = {}
        for profile in profiles:
            thresholds = resolve_validation(self.config, profile)
            rows = [
                evaluate_metrics(copy.deepcopy(item), thresholds, paired=source[1] is not None)
                for item in measured
            ]
            evaluated[profile] = rows
            self.write_validation_outputs(candidates, rows, f"variants_{profile}.tsv")
            accepted = {item.variant_id for item in rows if item.decision == "PASS"}
            write_fasta(
                [
                    FastaRecord(candidate.variant_id, candidate.sequence, f"profile={profile}")
                    for candidate in candidates
                    if candidate.variant_id in accepted
                ],
                self.outdir / f"Repeats_verified_{profile}.fasta",
            )
            write_fasta(
                [
                    FastaRecord(candidate.variant_id, candidate.sequence, f"profile={profile}")
                    for candidate in candidates
                    if candidate.variant_id not in accepted
                ],
                self.outdir / f"Repeats_rejected_{profile}.fasta",
            )
        self.write_validation_comparison(candidates, measured, profiles, evaluated)
        self.summary["validation_comparison_profiles"] = profiles
        self.summary["validation_comparison_passed"] = {
            profile: sum(item.decision == "PASS" for item in rows)
            for profile, rows in evaluated.items()
        }

    def load_external_candidates(self, path: Path) -> list[Candidate]:
        candidates: list[Candidate] = []
        for record in iter_fasta(path):
            tags = dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", record.description))
            try:
                assembly_round = int(tags.get("round", "0"))
            except ValueError:
                assembly_round = 0
            try:
                derep_size = int(tags.get("derep_size", "1"))
            except ValueError:
                derep_size = 1
            candidates.append(
                Candidate(
                    variant_id=record.name,
                    sequence=record.sequence,
                    source_contig=tags.get("source", "external"),
                    extraction=self.config["extraction"]["mode"],
                    assembly_round=assembly_round,
                    strand="+",
                    start0=0,
                    end0=len(record.sequence),
                    derep_size=derep_size,
                    member_ids=[record.name],
                )
            )
        return candidates

    def write_validation_comparison(
        self,
        candidates: list[Candidate],
        measured: list[ValidationMetrics],
        profiles: list[str],
        evaluated: dict[str, list[ValidationMetrics]],
    ) -> None:
        measured_by_id = {item.variant_id: item for item in measured}
        evaluated_by_profile = {
            profile: {item.variant_id: item for item in rows}
            for profile, rows in evaluated.items()
        }
        fields = [
            "variant_id", "length", "mapped_fragments", "primary_mapped_fragments",
            "single_variant_mapped_fragments", "cross_variant_multimapped_fragments",
            "proper_pair_fragments", "average_coverage", "covered_percent",
            "longest_zero_run", "max_adjacent_window_fold", "softclip_percent",
            "unique_kmers_total", "unique_kmers_supported", "unique_kmer_fraction",
        ]
        for profile in profiles:
            fields.extend([f"{profile}_decision", f"{profile}_failure_reasons"])
        with open(self.outdir / "validation_comparison.tsv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for candidate in candidates:
                item = measured_by_id[candidate.variant_id]
                row = {
                    "variant_id": candidate.variant_id,
                    "length": len(candidate.sequence),
                    "mapped_fragments": item.mapped_fragments,
                    "primary_mapped_fragments": item.primary_mapped_fragments,
                    "single_variant_mapped_fragments": item.single_variant_mapped_fragments,
                    "cross_variant_multimapped_fragments": item.cross_variant_multimapped_fragments,
                    "proper_pair_fragments": item.proper_pair_fragments,
                    "average_coverage": f"{item.average_coverage:.4f}",
                    "covered_percent": f"{item.covered_percent:.4f}",
                    "longest_zero_run": item.longest_zero_run,
                    "max_adjacent_window_fold": f"{item.max_adjacent_window_fold:.4f}",
                    "softclip_percent": f"{item.softclip_percent:.4f}",
                    "unique_kmers_total": item.unique_kmers_total,
                    "unique_kmers_supported": item.unique_kmers_supported,
                    "unique_kmer_fraction": "NA" if item.unique_kmer_fraction is None else f"{item.unique_kmer_fraction:.4f}",
                }
                for profile in profiles:
                    decision = evaluated_by_profile[profile][candidate.variant_id]
                    row[f"{profile}_decision"] = decision.decision
                    row[f"{profile}_failure_reasons"] = ",".join(decision.reasons)
                writer.writerow(row)

    def write_validation_outputs(
        self,
        candidates: list[Candidate],
        metrics: list[ValidationMetrics],
        filename: str = "variants.tsv",
    ) -> None:
        metric_by_id = {item.variant_id: item for item in metrics}
        path = self.outdir / filename
        fields = [
            "variant_id", "source_contig", "assembly_round", "extraction", "strand", "start0", "end0",
            "length", "derep_size", "members", "mapped_fragments", "primary_mapped_fragments",
            "single_variant_mapped_fragments", "cross_variant_multimapped_fragments",
            "proper_pair_fragments",
            "average_coverage", "covered_percent", "longest_zero_run", "max_adjacent_window_fold",
            "softclip_percent", "unique_kmers_total", "unique_kmers_supported", "unique_kmer_fraction",
            "decision", "failure_reasons",
        ]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for candidate in candidates:
                item = metric_by_id[candidate.variant_id]
                writer.writerow(
                    {
                        "variant_id": candidate.variant_id,
                        "source_contig": candidate.source_contig,
                        "assembly_round": candidate.assembly_round,
                        "extraction": candidate.extraction,
                        "strand": candidate.strand,
                        "start0": candidate.start0,
                        "end0": candidate.end0,
                        "length": len(candidate.sequence),
                        "derep_size": candidate.derep_size,
                        "members": ",".join(candidate.member_ids),
                        "mapped_fragments": item.mapped_fragments,
                        "primary_mapped_fragments": item.primary_mapped_fragments,
                        "single_variant_mapped_fragments": item.single_variant_mapped_fragments,
                        "cross_variant_multimapped_fragments": item.cross_variant_multimapped_fragments,
                        "proper_pair_fragments": item.proper_pair_fragments,
                        "average_coverage": f"{item.average_coverage:.4f}",
                        "covered_percent": f"{item.covered_percent:.4f}",
                        "longest_zero_run": item.longest_zero_run,
                        "max_adjacent_window_fold": f"{item.max_adjacent_window_fold:.4f}",
                        "softclip_percent": f"{item.softclip_percent:.4f}",
                        "unique_kmers_total": item.unique_kmers_total,
                        "unique_kmers_supported": item.unique_kmers_supported,
                        "unique_kmer_fraction": "NA" if item.unique_kmer_fraction is None else f"{item.unique_kmer_fraction:.4f}",
                        "decision": item.decision,
                        "failure_reasons": ",".join(item.reasons),
                    }
                )

    def find_existing_verified(self) -> Path | None:
        configured = self.config["input"].get("variants_fasta")
        if configured:
            path = Path(configured)
            return path if path.exists() and any(True for _ in iter_fasta(path)) else None
        path = self.outdir / "Repeats_verified.fasta"
        return path if path.exists() and any(True for _ in iter_fasta(path)) else None

    def calculate_proportion(self, verified_fasta: Path) -> None:
        assert self.clean_reads
        config = self.config["proportion"]
        alignment = self.workdir / "verified.aligned.fasta"
        profile_hmm = self.workdir / "verified_variants.hmm"
        self.runner.run(
            ["mafft", "--auto", str(verified_fasta)], stdout_path=alignment,
            stderr_path=self.workdir / "mafft.err",
        )
        self.runner.run(
            ["hmmbuild", str(profile_hmm), str(alignment)],
            stdout_path=self.workdir / "hmmbuild.out", stderr_path=self.workdir / "hmmbuild.err",
        )
        total_available = count_units(*self.clean_reads)
        if total_available == 0:
            raise TargetAssemblerError("The cleaned library is empty; genome proportion cannot be estimated")
        sample_n = min(total_available, int(config["sample_units"]))
        rows: list[dict] = []
        for replicate in range(1, int(config["replicates"]) + 1):
            sample = self.sample_to_work(
                self.clean_reads, sample_n, total_available, f"proportion_r{replicate}",
                self.config["runtime"]["seed"] + 1000 + replicate,
            )
            fasta = self.workdir / f"proportion_r{replicate}.fasta"
            _, total_nt = fastq_to_fasta(*sample, fasta)
            table = self.workdir / f"proportion_r{replicate}.tbl"
            self.run_nhmmer(fasta, table, profile_hmm)
            hits = self.filter_hmm_hits(parse_nhmmer_tbl(table))
            aligned_nt = hit_nucleotides(hits, "aligned_nt")
            envelope_nt = hit_nucleotides(hits, "envelope_nt")
            rows.append(
                {
                    "replicate": replicate,
                    "sample_units": sample_n,
                    "sample_total_nt": total_nt,
                    "aligned_repeat_nt": aligned_nt,
                    "aligned_nt_percent": 100.0 * aligned_nt / max(1, total_nt),
                    "envelope_repeat_nt": envelope_nt,
                    "envelope_nt_percent": 100.0 * envelope_nt / max(1, total_nt),
                }
            )
        output = self.outdir / "genome_proportion.tsv"
        keep_envelope = bool(config["include_envelope_metric"]) or config["report_metric"] == "envelope_nt"
        fields = list(rows[0])
        if not keep_envelope:
            fields = [key for key in fields if not key.startswith("envelope_")]
        with open(output, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {key: (f"{row[key]:.8f}" if isinstance(row[key], float) else row[key]) for key in fields}
                )
        metric = config["report_metric"] + "_percent"
        summary = summarise_replicates([row[metric] for row in rows])
        summary["reported_metric"] = config["report_metric"]
        summary["unit"] = "percent"
        summary["layout"] = self.layout
        (self.outdir / "genome_proportion_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        self.summary["genome_proportion"] = summary

    def cleanup(self) -> None:
        if self.workdir.exists() and not self.config["runtime"]["keep_work"]:
            shutil.rmtree(self.workdir, ignore_errors=True)


def _safe_project_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_")
    return cleaned or "target_assembly"


def _make_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger(f"target_assembler.{path.resolve()}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    file_handler = logging.FileHandler(path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger
