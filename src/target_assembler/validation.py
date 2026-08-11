from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from .external import Runner, check_stream
from .io import base_read_id, iter_fastq, reverse_complement, write_fasta
from .models import Candidate, ValidationMetrics
from .models import FastaRecord


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def canonical_kmer(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return sequence if sequence <= reverse else reverse


def candidate_unique_kmers(candidates: Iterable[Candidate], k: int) -> dict[str, set[str]]:
    owners: dict[str, set[str]] = defaultdict(set)
    for candidate in candidates:
        sequence = candidate.sequence.upper()
        for start in range(0, len(sequence) - k + 1):
            kmer = sequence[start : start + k]
            if set(kmer) <= {"A", "C", "G", "T"}:
                owners[canonical_kmer(kmer)].add(candidate.variant_id)
    result: dict[str, set[str]] = defaultdict(set)
    for kmer, variant_ids in owners.items():
        if len(variant_ids) == 1:
            result[next(iter(variant_ids))].add(kmer)
    return result


def observed_unique_kmers(
    reads1: str,
    reads2: str | None,
    unique_by_variant: dict[str, set[str]],
    k: int,
) -> dict[str, set[str]]:
    lookup: dict[str, str] = {}
    for variant_id, kmers in unique_by_variant.items():
        for kmer in kmers:
            lookup[kmer] = variant_id
    observed: dict[str, set[str]] = defaultdict(set)
    if not lookup:
        return observed
    for path in (reads1, reads2):
        if not path:
            continue
        for record in iter_fastq(path):
            sequence = record.sequence.upper()
            for start in range(0, len(sequence) - k + 1):
                kmer = sequence[start : start + k]
                if set(kmer) <= {"A", "C", "G", "T"}:
                    canonical = canonical_kmer(kmer)
                    owner = lookup.get(canonical)
                    if owner:
                        observed[owner].add(canonical)
    return observed


def alignment_mode_flag(sensitivity: str, alignment: str) -> str:
    sensitivity = sensitivity.removeprefix("--")
    if alignment == "local" and not sensitivity.endswith("-local"):
        sensitivity += "-local"
    return "--" + sensitivity


def build_validation_reference(
    candidates: list[Candidate], destination: Path, config: dict
) -> set[str]:
    """Write linear or tandemized references and return IDs requiring modulo projection."""
    mode = config.get("reference_mode", "tandem")
    tandem_ids: set[str] = set()
    records: list[FastaRecord] = []
    insert_max = int(config.get("insert_size_max", 0) or 0)
    minimum_copies = int(config.get("tandem_min_copies", 3))
    for candidate in candidates:
        sequence = candidate.sequence.upper()
        if mode == "tandem":
            # The synthetic molecule must contain internal placements long enough for a
            # complete PE fragment, with repeat copies on both sides of the scored unit.
            minimum_length = (2 * insert_max) + len(sequence)
            copies = max(minimum_copies, math.ceil(minimum_length / max(1, len(sequence))))
            sequence = sequence * copies
            tandem_ids.add(candidate.variant_id)
        records.append(FastaRecord(candidate.variant_id, sequence))
    write_fasta(records, destination)
    return tandem_ids


def collect_alignment_metrics(
    runner: Runner,
    index_prefix: str,
    reads1: str,
    reads2: str | None,
    candidates: list[Candidate],
    config: dict,
    threads: int,
    workdir: Path,
    tandem_ids: set[str],
) -> tuple[
    dict[str, list[int]],
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, list[int]],
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, set[str]],
]:
    command = [
        "bowtie2",
        "-x",
        index_prefix,
        alignment_mode_flag(config["sensitivity"], config["alignment"]),
        "-k",
        str(config["max_alignments"]),
        "-p",
        str(threads),
        "--no-unal",
    ]
    if reads2:
        command.extend(
            [
                "-1", reads1, "-2", reads2,
                "-I", str(config["insert_size_min"]),
                "-X", str(config["insert_size_max"]),
            ]
        )
    else:
        command.extend(["-U", reads1])

    coverage = {candidate.variant_id: [0] * len(candidate.sequence) for candidate in candidates}
    mapped: dict[str, set[str]] = defaultdict(set)
    proper: dict[str, set[str]] = defaultdict(set)
    clipping: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    primary_mapped: dict[str, set[str]] = defaultdict(set)
    fragment_variants: dict[str, set[str]] = defaultdict(set)
    process, stderr_handle = runner.popen_stdout(
        command, stderr_path=workdir / "validation_bowtie2.log"
    )
    assert process.stdout is not None
    try:
        for line in process.stdout:
            if line.startswith("@"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 11:
                continue
            flag = int(columns[1])
            if flag & 0x4:
                continue
            mapq = int(columns[4])
            if mapq < int(config["minimum_mapq"]):
                continue
            variant_id = columns[2]
            if variant_id not in coverage:
                continue
            fragment = base_read_id(columns[0])
            mapped[variant_id].add(fragment)
            fragment_variants[fragment].add(variant_id)
            if flag & 0x2:
                proper[variant_id].add(fragment)

            is_secondary = bool(flag & (0x100 | 0x800))
            if config["primary_only"] and is_secondary:
                continue
            primary_mapped[variant_id].add(fragment)

            reference_position = int(columns[3]) - 1
            clipped = False
            for length_text, operation in CIGAR_PATTERN.findall(columns[5]):
                length = int(length_text)
                if operation in {"M", "=", "X"}:
                    for position in range(reference_position, reference_position + length):
                        if variant_id in tandem_ids:
                            coverage[variant_id][position % len(coverage[variant_id])] += 1
                        elif 0 <= position < len(coverage[variant_id]):
                            coverage[variant_id][position] += 1
                    reference_position += length
                elif operation in {"D", "N"}:
                    reference_position += length
                elif operation == "S" and length >= int(config["softclip_minimum_bases"]):
                    clipped = True
            clipping[variant_id][1] += 1
            if clipped:
                clipping[variant_id][0] += 1
    finally:
        process.stdout.close()
    check_stream(process, stderr_handle, "Bowtie2 validation mapping")
    single_variant: dict[str, set[str]] = defaultdict(set)
    cross_variant: dict[str, set[str]] = defaultdict(set)
    for fragment, variant_ids in fragment_variants.items():
        destination = single_variant if len(variant_ids) == 1 else cross_variant
        for variant_id in variant_ids:
            destination[variant_id].add(fragment)
    return coverage, mapped, proper, clipping, primary_mapped, single_variant, cross_variant


def longest_zero_run(values: list[int]) -> int:
    longest = current = 0
    for value in values:
        if value == 0:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def adjacent_window_fold(values: list[int], window_size: int) -> float:
    if not values:
        return math.inf
    means = [
        sum(values[start : start + window_size]) / len(values[start : start + window_size])
        for start in range(0, len(values), window_size)
    ]
    if len(means) < 2:
        return 1.0
    return max((max(a, b) + 1.0) / (min(a, b) + 1.0) for a, b in zip(means, means[1:]))


def evaluate_metrics(metrics: ValidationMetrics, config: dict, paired: bool) -> ValidationMetrics:
    failures: list[str] = []
    checks = (
        (metrics.average_coverage >= float(config["minimum_average_coverage"]), "low_average_coverage"),
        (metrics.covered_percent >= float(config["minimum_covered_percent"]), "insufficient_breadth"),
        (metrics.mapped_fragments >= int(config["minimum_mapped_fragments"]), "too_few_mapped_fragments"),
        (metrics.longest_zero_run <= int(config["maximum_zero_run"]), "long_uncovered_run"),
        (
            metrics.max_adjacent_window_fold <= float(config["maximum_adjacent_window_fold"]),
            "abrupt_coverage_change",
        ),
        (metrics.softclip_percent <= float(config["maximum_softclip_percent"]), "excess_soft_clipping"),
    )
    failures.extend(reason for passed, reason in checks if not passed)

    if paired and config["pair_rule"] != "disabled":
        pair_percent = 100.0 * metrics.proper_pair_fragments / max(1, metrics.mapped_fragments)
        absolute_ok = metrics.proper_pair_fragments >= int(config["minimum_proper_pairs"])
        percent_ok = pair_percent >= float(config["minimum_proper_pair_percent"])
        pair_ok = (absolute_ok or percent_ok) if config["pair_rule"] == "any" else (absolute_ok and percent_ok)
        if not pair_ok:
            failures.append("insufficient_proper_pair_support")

    if metrics.unique_kmers_total:
        if (metrics.unique_kmer_fraction or 0.0) < float(config["minimum_unique_kmer_fraction"]):
            failures.append("insufficient_variant_specific_kmer_support")
    elif config.get("require_unique_kmers", False):
        failures.append("no_informative_unique_kmers")

    metrics.reasons = failures
    metrics.decision = "PASS" if not failures else "FAIL"
    return metrics


def measure_candidates(
    runner: Runner,
    candidates_fasta: str,
    candidates: list[Candidate],
    reads1: str,
    reads2: str | None,
    config: dict,
    threads: int,
    workdir: Path,
) -> list[ValidationMetrics]:
    reference = workdir / "validation_reference.fasta"
    tandem_ids = build_validation_reference(candidates, reference, config)
    index_prefix = str(workdir / "candidate_index")
    runner.run(
        ["bowtie2-build", str(reference), index_prefix],
        stdout_path=workdir / "bowtie2_build.out",
        stderr_path=workdir / "bowtie2_build.err",
    )
    coverage, mapped, proper, clipping, primary_mapped, single_variant, cross_variant = collect_alignment_metrics(
        runner, index_prefix, reads1, reads2, candidates, config, threads, workdir, tandem_ids
    )
    k = int(config["kmer_size"])
    unique = candidate_unique_kmers(candidates, k)
    observed = observed_unique_kmers(reads1, reads2, unique, k)
    metrics: list[ValidationMetrics] = []
    for candidate in candidates:
        values = coverage[candidate.variant_id]
        clipped, alignments = clipping[candidate.variant_id]
        unique_total = len(unique.get(candidate.variant_id, set()))
        unique_supported = len(observed.get(candidate.variant_id, set()))
        item = ValidationMetrics(
            variant_id=candidate.variant_id,
            length=len(candidate.sequence),
            mapped_fragments=len(mapped.get(candidate.variant_id, set())),
            proper_pair_fragments=len(proper.get(candidate.variant_id, set())),
            average_coverage=sum(values) / max(1, len(values)),
            covered_percent=100.0 * sum(value > 0 for value in values) / max(1, len(values)),
            longest_zero_run=longest_zero_run(values),
            max_adjacent_window_fold=adjacent_window_fold(values, int(config["window_size"])),
            softclip_percent=100.0 * clipped / max(1, alignments),
            unique_kmers_total=unique_total,
            unique_kmers_supported=unique_supported,
            unique_kmer_fraction=(unique_supported / unique_total) if unique_total else None,
            primary_mapped_fragments=len(primary_mapped.get(candidate.variant_id, set())),
            single_variant_mapped_fragments=len(single_variant.get(candidate.variant_id, set())),
            cross_variant_multimapped_fragments=len(cross_variant.get(candidate.variant_id, set())),
        )
        metrics.append(item)
    return metrics


def validate_candidates(
    runner: Runner,
    candidates_fasta: str,
    candidates: list[Candidate],
    reads1: str,
    reads2: str | None,
    config: dict,
    threads: int,
    workdir: Path,
) -> list[ValidationMetrics]:
    measured = measure_candidates(
        runner, candidates_fasta, candidates, reads1, reads2, config, threads, workdir
    )
    return [evaluate_metrics(item, config, paired=reads2 is not None) for item in measured]
