from __future__ import annotations

import io
import logging
import tempfile
import unittest
from pathlib import Path

import yaml

from target_assembler.config import choose_validation_profile, infer_layout, load_config, resolve_validation
from target_assembler.io import (
    count_units,
    exact_dereplicate,
    iter_fastq,
    normalize_library,
    parse_nhmmer_tbl,
    sample_library,
)
from target_assembler.models import Candidate, FastaRecord, ValidationMetrics
from target_assembler.models import HMMHit
from target_assembler.pipeline import Pipeline
from target_assembler.proportion import hit_nucleotides, merge_intervals, summarise_replicates
from target_assembler.validation import (
    build_validation_reference,
    candidate_unique_kmers,
    collect_alignment_metrics,
    evaluate_metrics,
    longest_zero_run,
)


def write_fastq(path: Path, names: list[str], mate: int | None = None) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for name in names:
            suffix = f"/{mate}" if mate else ""
            sequence = "ACGT" * 5
            handle.write(f"@{name}{suffix} description\n{sequence}\n+\n{'I' * len(sequence)}\n")


def minimal_config(directory: Path, paired: bool = True) -> Path:
    hmm = directory / "target.hmm"
    hmm.write_text("HMMER3/f\n", encoding="utf-8")
    r1 = directory / "R1.fq"
    write_fastq(r1, ["r1", "r2"], 1 if paired else None)
    data = {
        "project": {"mode": "variants", "outdir": str(directory / "out")},
        "input": {"reads1": str(r1), "layout": "paired" if paired else "single", "hmm": str(hmm)},
    }
    if paired:
        r2 = directory / "R2.fq"
        write_fastq(r2, ["r1", "r2"], 2)
        data["input"]["reads2"] = str(r2)
    path = directory / "config.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


class CoreTests(unittest.TestCase):
    def test_config_infers_and_resolves_profiles(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = load_config(minimal_config(Path(temporary)))
            self.assertEqual(infer_layout(config), "paired")
            self.assertEqual(choose_validation_profile(config, "paired", 99), "lenient")
            self.assertEqual(choose_validation_profile(config, "paired", 100), "balanced")
            self.assertEqual(choose_validation_profile(config, "paired", 1000), "strict")
            resolved = resolve_validation(config, "strict")
            self.assertEqual(resolved["source"], "clean")
            self.assertEqual(resolved["minimum_mapq"], 0)
            self.assertEqual(resolve_validation(config, "lenient")["source"], "expanded")
            self.assertEqual(resolved["reference_mode"], "tandem")

    def test_single_layout_selects_lenient(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = load_config(minimal_config(Path(temporary), paired=False))
            self.assertEqual(infer_layout(config), "single")
            self.assertEqual(choose_validation_profile(config, "single", 100000), "lenient")

    def test_proportion_only_does_not_require_input_hmm(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            r1, r2 = directory / "clean_R1.fq", directory / "clean_R2.fq"
            write_fastq(r1, ["r1", "r2"], 1)
            write_fastq(r2, ["r1", "r2"], 2)
            variants = directory / "accepted.fasta"
            variants.write_text(">v1\nACGTACGT\n", encoding="utf-8")
            config_path = directory / "proportion.yaml"
            config_path.write_text(
                yaml.safe_dump(
                    {
                        "project": {"mode": "proportion", "outdir": str(directory / "out")},
                        "input": {
                            "clean_reads1": str(r1),
                            "clean_reads2": str(r2),
                            "layout": "paired",
                            "variants_fasta": str(variants),
                        },
                    }
                ),
                encoding="utf-8",
            )
            config = load_config(config_path)
            self.assertIsNone(config["input"]["hmm"])
            self.assertEqual(config["input"]["variants_fasta"], str(variants))

    def test_paired_normalization_and_exact_sampling(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            r1, r2 = directory / "a.fq", directory / "b.fq"
            write_fastq(r1, [f"read{i}" for i in range(10)], 1)
            write_fastq(r2, [f"read{i}" for i in range(10)], 2)
            n1, n2 = directory / "n1.fq", directory / "n2.fq"
            self.assertEqual(normalize_library(r1, r2, n1, n2), 10)
            self.assertEqual(next(iter_fastq(n1)).name, "read0/1")
            s1, s2 = directory / "s1.fq", directory / "s2.fq"
            self.assertEqual(sample_library(n1, n2, s1, s2, 4, 10, 42, True), 4)
            self.assertEqual(count_units(s1, s2), 4)

    def test_hmmer_parser_and_interval_metrics(self):
        with tempfile.TemporaryDirectory() as temporary:
            table = Path(temporary) / "hits.tbl"
            table.write_text(
                "# header\nread/1 - profile - 1 100 5 104 4 106 150 + 1e-20 50.0 0.0 desc\n"
                "read/1 - profile - 1 100 100 140 98 142 150 + 1e-10 40.0 0.0 desc\n",
                encoding="utf-8",
            )
            hits = parse_nhmmer_tbl(table)
            self.assertEqual(hits[0].model_fraction(100), 1.0)
            self.assertEqual(hit_nucleotides(hits, "aligned_nt"), 136)
            self.assertEqual(hit_nucleotides(hits, "envelope_nt"), 139)
            self.assertEqual(merge_intervals([(1, 5), (5, 7), (10, 10)]), 8)

    def test_tandem_extraction_allows_partial_boundary_hit(self):
        first = HMMHit("ctg", 61, 120, 1, 60, 1, 60, 500, "+", 1e-20, 50.0)
        second = HMMHit("ctg", 12, 120, 250, 358, 250, 358, 500, "+", 1e-20, 50.0)
        model_length = 120

        self.assertLess(first.model_fraction(model_length), 0.90)
        self.assertGreaterEqual(second.model_end_fraction(model_length), 0.90)
        self.assertIs(Pipeline.tandem_downstream_hit(first, second), second)

        reverse_first = HMMHit("ctg", 1, 120, 100, 219, 100, 219, 500, "-", 1e-20, 50.0)
        reverse_second = HMMHit("ctg", 31, 120, 300, 389, 300, 389, 500, "-", 1e-20, 50.0)
        self.assertIs(Pipeline.tandem_downstream_hit(reverse_first, reverse_second), reverse_first)

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            contigs = directory / "contigs.fasta"
            contigs.write_text(">mira_contig\n" + "ACGT" * 100 + "\n", encoding="utf-8")

            pipeline = object.__new__(Pipeline)
            pipeline.workdir = directory / "work"
            pipeline.workdir.mkdir()
            pipeline.outdir = directory / "out"
            pipeline.outdir.mkdir()
            pipeline.config = {
                "extraction": {
                    "mode": "tandem",
                    "flank5": 0,
                    "flank3": 0,
                    "minimum_model_fraction": 0.90,
                }
            }
            pipeline.summary = {"events": []}
            pipeline.logger = logging.getLogger("target_assembler_test_extraction")
            pipeline.hmm_model_length = lambda: model_length

            def fake_nhmmer(_fasta, table, _hmm=None):
                table = Path(table)
                if table.name == "orientation_hits.tbl":
                    rows = "contig_r2_00001 - profile - 61 120 1 60 1 60 400 + 1e-20 50.0 0.0 first\n"
                else:
                    rows = (
                        "contig_r2_00001 - profile - 61 120 1 60 1 60 400 + 1e-20 50.0 0.0 boundary\n"
                        "contig_r2_00001 - profile - 12 120 200 308 200 308 400 + 1e-20 50.0 0.0 downstream\n"
                    )
                table.write_text(rows, encoding="utf-8")

            pipeline.run_nhmmer = fake_nhmmer
            candidates, _ = pipeline.extract_candidates(contigs, 2)

            self.assertEqual(len(candidates), 1)
            events = {item["name"]: item["value"] for item in pipeline.summary["events"]}
            self.assertEqual(events["tandem_hit_pairs"], 1)
            self.assertEqual(events["skipped_incomplete_coding"], 0)
            self.assertEqual(events["extracted_sequences_count"], 1)

    def test_derep_and_unique_kmers(self):
        records = [FastaRecord("a", "AAAACCCC"), FastaRecord("b", "AAAACCCC"), FastaRecord("c", "AAAAGGGG")]
        derep = exact_dereplicate(records)
        self.assertEqual([len(members) for _, members in derep], [2, 1])
        candidates = [
            Candidate("v1", "AAAACCCC", "c1", "site", 1, "+", 0, 8),
            Candidate("v2", "AAAAGGGG", "c2", "site", 1, "+", 0, 8),
        ]
        unique = candidate_unique_kmers(candidates, 4)
        self.assertIn("AACC", unique["v1"])
        self.assertNotIn("AACC", unique["v2"])

    def test_tandem_validation_reference_is_longer_than_insert_context(self):
        with tempfile.TemporaryDirectory() as temporary:
            candidate = Candidate("v1", "ACGT" * 50, "c1", "tandem", 1, "+", 0, 200)
            path = Path(temporary) / "validation.fasta"
            tandem_ids = build_validation_reference(
                [candidate], path, {"reference_mode": "tandem", "insert_size_max": 450, "tandem_min_copies": 3}
            )
            sequence = "".join(line.strip() for line in path.read_text().splitlines() if not line.startswith(">"))
            self.assertIn("v1", tandem_ids)
            self.assertGreaterEqual(len(sequence), 1100)
            self.assertEqual(len(sequence) % 200, 0)

    def test_tandem_coverage_wraps_and_multimapping_is_reported_separately(self):
        class FakeProcess:
            def __init__(self, sam: str):
                self.stdout = io.StringIO(sam)

            def wait(self):
                return 0

        class FakeRunner:
            def __init__(self, sam: str):
                self.sam = sam

            def popen_stdout(self, command, stderr_path):
                return FakeProcess(self.sam), open(stderr_path, "w", encoding="utf-8")

        sam = (
            "read1/1\t99\tv1\t191\t0\t20M\t=\t221\t50\t" + "A" * 20 + "\t" + "I" * 20 + "\n"
            "read1/2\t147\tv1\t221\t0\t20M\t=\t191\t-50\t" + "A" * 20 + "\t" + "I" * 20 + "\n"
            "read2/1\t0\tv1\t1\t0\t20M\t*\t0\t0\t" + "A" * 20 + "\t" + "I" * 20 + "\n"
            "read2/1\t256\tv2\t1\t0\t20M\t*\t0\t0\t" + "A" * 20 + "\t" + "I" * 20 + "\n"
        )
        candidates = [
            Candidate("v1", "A" * 200, "c1", "tandem", 1, "+", 0, 200),
            Candidate("v2", "C" * 200, "c2", "tandem", 1, "+", 0, 200),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            coverage, mapped, proper, _, primary, single, cross = collect_alignment_metrics(
                FakeRunner(sam),
                "unused",
                "R1.fq",
                "R2.fq",
                candidates,
                {
                    "sensitivity": "very-sensitive",
                    "alignment": "local",
                    "max_alignments": 10,
                    "insert_size_min": 1,
                    "insert_size_max": 500,
                    "minimum_mapq": 0,
                    "primary_only": True,
                    "softclip_minimum_bases": 12,
                },
                1,
                Path(temporary),
                {"v1", "v2"},
            )
        self.assertGreater(coverage["v1"][0], 0)
        self.assertGreater(coverage["v1"][199], 0)
        self.assertIn("read1", proper["v1"])
        self.assertIn("read1", single["v1"])
        self.assertIn("read2", cross["v1"])
        self.assertIn("read2", cross["v2"])
        self.assertNotIn("read2", primary["v2"])
        self.assertEqual(mapped["v1"], {"read1", "read2"})

    def test_validate_mode_accepts_dedicated_library_without_hmm(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            r1, r2 = directory / "validation_R1.fq", directory / "validation_R2.fq"
            write_fastq(r1, ["r1", "r2"], 1)
            write_fastq(r2, ["r1", "r2"], 2)
            candidates = directory / "candidates.fasta"
            candidates.write_text(">v1\nACGTACGT\n", encoding="utf-8")
            config_path = directory / "validate.yaml"
            config_path.write_text(
                yaml.safe_dump(
                    {
                        "project": {"mode": "validate", "outdir": str(directory / "out")},
                        "input": {
                            "layout": "paired",
                            "validation_reads1": str(r1),
                            "validation_reads2": str(r2),
                            "variants_fasta": str(candidates),
                        },
                    }
                ),
                encoding="utf-8",
            )
            config = load_config(config_path)
            self.assertIsNone(config["input"]["hmm"])
            self.assertEqual(config["project"]["mode"], "validate")

    def test_se_validation_ignores_pair_thresholds(self):
        config = {
            "minimum_average_coverage": 2,
            "minimum_covered_percent": 50,
            "minimum_mapped_fragments": 3,
            "maximum_zero_run": 10,
            "maximum_adjacent_window_fold": 10,
            "maximum_softclip_percent": 50,
            "pair_rule": "all",
            "minimum_proper_pairs": 100,
            "minimum_proper_pair_percent": 100,
            "minimum_unique_kmer_fraction": 0,
            "require_unique_kmers": False,
        }
        metrics = ValidationMetrics("v1", 100, 5, 0, 3, 80, 2, 2, 0, 0, 0, None)
        self.assertEqual(evaluate_metrics(metrics, config, paired=False).decision, "PASS")
        self.assertEqual(longest_zero_run([0, 0, 1, 0]), 2)

    def test_replicate_summary(self):
        summary = summarise_replicates([1.0, 2.0, 3.0])
        self.assertEqual(summary["mean"], 2.0)
        self.assertLess(summary["ci95_low"], 2.0)
        self.assertGreater(summary["ci95_high"], 2.0)


if __name__ == "__main__":
    unittest.main()


# --- selection.hmm_evalue / selection.hmm_score_min --------------------------
# Added with the threshold keys: nhmmer's own default is -E 10, which for a
# conserved rRNA profile recruits across kingdoms. These tests pin both the
# command line and the post-parse filter.

def test_hmm_threshold_defaults_are_none():
    from target_assembler.config import DEFAULTS
    assert DEFAULTS["selection"]["hmm_evalue"] is None
    assert DEFAULTS["selection"]["hmm_score_min"] is None


def test_hmm_threshold_validation_rejects_nonsense():
    import pytest
    from target_assembler.config import DEFAULTS, validate_config
    from target_assembler.errors import ConfigError
    import copy

    def cfg(**sel):
        c = copy.deepcopy(DEFAULTS)
        c["project"]["mode"] = "proportion"
        c["input"]["reads1"] = "r.fq"
        c["input"]["layout"] = "single"
        c["selection"].update(sel)
        return c

    validate_config(cfg())                          # both null: fine
    validate_config(cfg(hmm_score_min=98.4))        # our calibrated value
    validate_config(cfg(hmm_evalue=1e-10))
    with pytest.raises(ConfigError):
        validate_config(cfg(hmm_evalue=0))
    with pytest.raises(ConfigError):
        validate_config(cfg(hmm_score_min=-1))


def test_filter_hmm_hits_applies_score_and_evalue():
    from target_assembler.models import HMMHit

    class Fake:
        """Only what filter_hmm_hits touches."""
        def __init__(self, **sel):
            self.config = {"selection": dict(sel)}
            self.summary = {"events": []}
            import logging
            self.logger = logging.getLogger("t")
        filter_hmm_hits = None

    from target_assembler.pipeline import Pipeline
    Fake.filter_hmm_hits = Pipeline.filter_hmm_hits
    Fake.hmm_thresholds = Pipeline.hmm_thresholds
    Fake.event = Pipeline.event

    def hit(score, evalue):
        return HMMHit(target="r", hmm_from=1, hmm_to=10, ali_from=1, ali_to=10,
                      env_from=1, env_to=10, target_length=150, strand="+",
                      evalue=evalue, score=score)

    hits = [hit(169.8, 1e-40), hit(98.4, 1e-20), hit(40.0, 0.5)]

    # no thresholds: unchanged, and the object identity is preserved
    assert Fake().filter_hmm_hits(hits) is hits

    # score threshold at the calibrated separation point
    kept = Fake(hmm_score_min=98.4).filter_hmm_hits(hits)
    assert [h.score for h in kept] == [169.8, 98.4]

    # E-value threshold
    kept = Fake(hmm_evalue=1e-30).filter_hmm_hits(hits)
    assert [h.score for h in kept] == [169.8]

    # both: each is applied
    kept = Fake(hmm_score_min=50.0, hmm_evalue=1e-30).filter_hmm_hits(hits)
    assert [h.score for h in kept] == [169.8]

    # the drop is recorded as an event, not silent
    f = Fake(hmm_score_min=1000.0)
    assert f.filter_hmm_hits(hits) == []
    assert f.summary["events"][-1]["name"] == "hmm_hits_below_threshold"
    assert f.summary["events"][-1]["value"] == 3


def test_run_nhmmer_command_carries_the_threshold(tmp_path):
    from target_assembler.pipeline import Pipeline

    recorded = {}

    class FakeRunner:
        def run(self, command, stdout_path=None, stderr_path=None):
            recorded["command"] = list(command)

    class Fake:
        hmm = "t.hmm"
        def __init__(self, **sel):
            self.config = {"selection": dict(sel), "runtime": {"threads": 4}}
            self.runner = FakeRunner()
    Fake.hmm_thresholds = Pipeline.hmm_thresholds
    Fake.run_nhmmer = Pipeline.run_nhmmer

    Fake().run_nhmmer("in.fa", tmp_path / "t.tbl")
    assert "-T" not in recorded["command"] and "-E" not in recorded["command"]

    Fake(hmm_score_min=98.4).run_nhmmer("in.fa", tmp_path / "t.tbl")
    assert "-T" in recorded["command"] and "--incT" in recorded["command"]
    assert "-E" not in recorded["command"]

    Fake(hmm_evalue=1e-10).run_nhmmer("in.fa", tmp_path / "t.tbl")
    assert "-E" in recorded["command"] and "--incE" in recorded["command"]
    assert "-T" not in recorded["command"]

    # -T and -E are mutually exclusive in HMMER; -T must win, alone
    Fake(hmm_score_min=98.4, hmm_evalue=1e-10).run_nhmmer("in.fa", tmp_path / "t.tbl")
    assert "-T" in recorded["command"] and "-E" not in recorded["command"]
