from __future__ import annotations

from collections import defaultdict
from statistics import mean, stdev
from typing import Iterable

from .models import HMMHit


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    if not ordered:
        return 0
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end
    return total + end - start + 1


def hit_nucleotides(hits: Iterable[HMMHit], metric: str) -> int:
    by_read: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for hit in hits:
        if metric == "aligned_nt":
            by_read[hit.target].append((hit.ali_from, hit.ali_to))
        elif metric == "envelope_nt":
            by_read[hit.target].append((hit.env_from, hit.env_to))
        else:
            raise ValueError(f"Unknown interval metric: {metric}")
    return sum(merge_intervals(intervals) for intervals in by_read.values())


def summarise_replicates(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"replicates": 0, "mean": None, "sd": None, "ci95_low": None, "ci95_high": None}
    average = mean(values)
    if len(values) == 1:
        return {"replicates": 1, "mean": average, "sd": None, "ci95_low": None, "ci95_high": None}
    deviation = stdev(values)
    margin = 1.96 * deviation / (len(values) ** 0.5)
    return {
        "replicates": len(values),
        "mean": average,
        "sd": deviation,
        "ci95_low": average - margin,
        "ci95_high": average + margin,
    }

