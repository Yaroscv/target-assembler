from __future__ import annotations

import gzip
import random
import re
from contextlib import contextmanager
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Iterator, TextIO

from .errors import TargetAssemblerError
from .models import FastaRecord, FastqRecord, HMMHit


def open_text(path: str | Path, mode: str = "rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def base_read_id(name: str) -> str:
    token = name.lstrip("@>").split(maxsplit=1)[0]
    return re.sub(r"/[12]$", "", token)


def iter_fastq(path: str | Path) -> Iterator[FastqRecord]:
    with open_text(path, "rt") as handle:
        line_no = 0
        while True:
            header = handle.readline()
            if not header:
                return
            line_no += 1
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            line_no += 3
            if not sequence or not plus or not quality:
                raise TargetAssemblerError(f"Truncated FASTQ record near line {line_no} in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise TargetAssemblerError(f"Malformed FASTQ record near line {line_no - 3} in {path}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if len(sequence) != len(quality):
                raise TargetAssemblerError(
                    f"Sequence/quality length mismatch for {header.strip()} in {path}"
                )
            yield FastqRecord(header[1:].strip(), sequence, quality)


def iter_library(reads1: str | Path, reads2: str | Path | None = None):
    if reads2 is None:
        for record in iter_fastq(reads1):
            yield (record,)
        return
    for index, pair in enumerate(zip_longest(iter_fastq(reads1), iter_fastq(reads2)), start=1):
        r1, r2 = pair
        if r1 is None or r2 is None:
            raise TargetAssemblerError("Paired FASTQ files contain different numbers of records")
        if base_read_id(r1.name) != base_read_id(r2.name):
            raise TargetAssemblerError(
                f"Unsynchronised pair {index}: {r1.name!r} versus {r2.name!r}"
            )
        yield r1, r2


def write_fastq_record(handle: TextIO, record: FastqRecord, name: str | None = None) -> None:
    handle.write(f"@{name or record.name}\n{record.sequence}\n+\n{record.quality}\n")


def count_units(reads1: str | Path, reads2: str | Path | None = None) -> int:
    return sum(1 for _ in iter_library(reads1, reads2))


def sample_library(
    reads1: str | Path,
    reads2: str | Path | None,
    out1: str | Path,
    out2: str | Path | None,
    number: int,
    total: int,
    seed: int,
    normalize_names: bool = False,
) -> int:
    """Select exactly ``number`` units without holding the library in memory."""
    number = min(max(0, int(number)), total)
    rng = random.Random(seed)
    selected = 0
    remaining_total = total
    remaining_needed = number
    with open_text(out1, "wt") as h1:
        h2_cm = open_text(out2, "wt") if out2 else _null_context()
        with h2_cm as h2:
            for records in iter_library(reads1, reads2):
                take = remaining_needed > 0 and rng.random() < (remaining_needed / remaining_total)
                remaining_total -= 1
                if not take:
                    continue
                base = base_read_id(records[0].name)
                write_fastq_record(h1, records[0], f"{base}/1" if normalize_names and reads2 else (base if normalize_names else None))
                if reads2 and h2 is not None:
                    write_fastq_record(h2, records[1], f"{base}/2" if normalize_names else None)
                selected += 1
                remaining_needed -= 1
    if selected != number:
        raise TargetAssemblerError(f"Sampling selected {selected} units; expected {number}")
    return selected


@contextmanager
def _null_context():
    yield None


def normalize_library(
    reads1: str | Path,
    reads2: str | Path | None,
    out1: str | Path,
    out2: str | Path | None,
) -> int:
    count = 0
    with open_text(out1, "wt") as h1:
        h2_cm = open_text(out2, "wt") if out2 else _null_context()
        with h2_cm as h2:
            for records in iter_library(reads1, reads2):
                base = base_read_id(records[0].name)
                write_fastq_record(h1, records[0], f"{base}/1" if reads2 else base)
                if reads2 and h2 is not None:
                    write_fastq_record(h2, records[1], f"{base}/2")
                count += 1
    return count


def extract_library_by_ids(
    reads1: str | Path,
    reads2: str | Path | None,
    wanted: set[str],
    out1: str | Path,
    out2: str | Path | None,
    normalize_names: bool = True,
) -> int:
    count = 0
    with open_text(out1, "wt") as h1:
        h2_cm = open_text(out2, "wt") if out2 else _null_context()
        with h2_cm as h2:
            for records in iter_library(reads1, reads2):
                base = base_read_id(records[0].name)
                if base not in wanted:
                    continue
                write_fastq_record(h1, records[0], f"{base}/1" if reads2 and normalize_names else base)
                if reads2 and h2 is not None:
                    write_fastq_record(h2, records[1], f"{base}/2" if normalize_names else records[1].name)
                count += 1
    return count


def fastq_to_fasta(
    reads1: str | Path,
    reads2: str | Path | None,
    destination: str | Path,
) -> tuple[int, int]:
    units = 0
    total_nt = 0
    with open(destination, "w", encoding="utf-8") as handle:
        for records in iter_library(reads1, reads2):
            for record in records:
                handle.write(f">{record.name.split()[0]}\n{record.sequence}\n")
                total_nt += len(record.sequence)
            units += 1
    return units, total_nt


def total_fastq_nt(reads1: str | Path, reads2: str | Path | None) -> int:
    return sum(len(record.sequence) for records in iter_library(reads1, reads2) for record in records)


def iter_fasta(path: str | Path) -> Iterator[FastaRecord]:
    name: str | None = None
    description = ""
    chunks: list[str] = []
    with open_text(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield FastaRecord(name, "".join(chunks), description)
                full = line[1:]
                parts = full.split(maxsplit=1)
                name = parts[0]
                description = parts[1] if len(parts) == 2 else ""
                chunks = []
            else:
                if name is None:
                    raise TargetAssemblerError(f"FASTA sequence before header in {path}")
                chunks.append(line)
    if name is not None:
        yield FastaRecord(name, "".join(chunks), description)


def write_fasta(records: Iterable[FastaRecord], path: str | Path, width: int = 80) -> int:
    count = 0
    with open(path, "w", encoding="utf-8") as handle:
        for record in records:
            description = f" {record.description}" if record.description else ""
            handle.write(f">{record.name}{description}\n")
            sequence = str(record.sequence)
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")
            count += 1
    return count


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")
    return sequence.translate(table)[::-1]


def parse_nhmmer_tbl(path: str | Path) -> list[HMMHit]:
    hits: list[HMMHit] = []
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) < 15:
                raise TargetAssemblerError(f"Unexpected nhmmer tblout at {path}:{line_no}")
            try:
                hits.append(
                    HMMHit(
                        target=columns[0],
                        hmm_from=int(columns[4]),
                        hmm_to=int(columns[5]),
                        ali_from=int(columns[6]),
                        ali_to=int(columns[7]),
                        env_from=int(columns[8]),
                        env_to=int(columns[9]),
                        target_length=int(columns[10]),
                        strand=columns[11],
                        evalue=float(columns[12]),
                        score=float(columns[13]),
                    )
                )
            except ValueError as exc:
                raise TargetAssemblerError(f"Invalid nhmmer numeric field at {path}:{line_no}") from exc
    return hits


def exact_dereplicate(records: Iterable[FastaRecord]) -> list[tuple[FastaRecord, list[str]]]:
    by_sequence: dict[str, tuple[FastaRecord, list[str]]] = {}
    order: list[str] = []
    for record in records:
        key = record.sequence.upper()
        if key not in by_sequence:
            by_sequence[key] = (record, [record.name])
            order.append(key)
        else:
            by_sequence[key][1].append(record.name)
    return [by_sequence[key] for key in order]

