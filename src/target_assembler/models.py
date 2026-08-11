from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class FastqRecord:
    name: str
    sequence: str
    quality: str


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str
    description: str = ""


@dataclass(frozen=True)
class HMMHit:
    target: str
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    target_length: int
    strand: str
    evalue: float
    score: float

    @property
    def start0(self) -> int:
        return min(self.ali_from, self.ali_to) - 1

    @property
    def end0(self) -> int:
        return max(self.ali_from, self.ali_to)

    def model_fraction(self, model_length: int) -> float:
        return (abs(self.hmm_to - self.hmm_from) + 1) / model_length

    def model_end_fraction(self, model_length: int) -> float:
        """Fraction of the HMM model reached by the alignment's 3-prime end."""
        return max(self.hmm_from, self.hmm_to) / model_length


@dataclass
class Candidate:
    variant_id: str
    sequence: str
    source_contig: str
    extraction: str
    assembly_round: int
    strand: str
    start0: int
    end0: int
    derep_size: int = 1
    member_ids: list[str] = field(default_factory=list)


@dataclass
class ValidationMetrics:
    variant_id: str
    length: int
    mapped_fragments: int
    proper_pair_fragments: int
    average_coverage: float
    covered_percent: float
    longest_zero_run: int
    max_adjacent_window_fold: float
    softclip_percent: float
    unique_kmers_total: int
    unique_kmers_supported: int
    unique_kmer_fraction: float | None
    primary_mapped_fragments: int = 0
    single_variant_mapped_fragments: int = 0
    cross_variant_multimapped_fragments: int = 0
    decision: str = "PENDING"
    reasons: list[str] = field(default_factory=list)
