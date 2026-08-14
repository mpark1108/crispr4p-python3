"""CRISPR4P result types."""

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping

from .coordinates import cut_from_pam
from .guides import SUFFIX_LENGTHS


@dataclass(frozen=True)
class GuideCandidate:
    """One guide from a design result."""

    chromosome: str
    seed: str
    forward_cloning_oligo: str
    reverse_cloning_oligo: str
    pam_coordinates: tuple[int, int]
    cut_coordinates: tuple[int, int]
    strand: int
    pam: str
    match_counts: Mapping[int, int]
    _legacy_row: Any = field(repr=False, compare=False)

    @classmethod
    def from_legacy_row(cls, legacy_row, chromosome):
        """Read a nine-item guide row."""
        if not isinstance(legacy_row, (list, tuple)):
            raise TypeError("legacy guide row must be a list or tuple")
        if len(legacy_row) != 9:
            raise ValueError("legacy guide row must contain nine items")

        primer = legacy_row[1]
        if not isinstance(primer, (list, tuple)):
            raise TypeError("legacy guide primer must be a list or tuple")
        if len(primer) != 6:
            raise ValueError("legacy guide primer must contain six items")

        pam_coordinates = tuple(primer[3])
        strand = primer[4]
        return cls(
            chromosome=chromosome,
            seed=legacy_row[0],
            forward_cloning_oligo=primer[1],
            reverse_cloning_oligo=primer[2],
            pam_coordinates=pam_coordinates,
            cut_coordinates=cut_from_pam(
                pam_coordinates[0],
                pam_coordinates[1],
                strand,
            ),
            strand=strand,
            pam=primer[5],
            match_counts={
                length: count
                for length, count in zip(SUFFIX_LENGTHS, legacy_row[2:])
            },
            _legacy_row=legacy_row,
        )

    def __post_init__(self):
        object.__setattr__(self, "pam_coordinates", tuple(self.pam_coordinates))
        object.__setattr__(self, "cut_coordinates", tuple(self.cut_coordinates))
        object.__setattr__(
            self,
            "match_counts",
            MappingProxyType(dict(self.match_counts)),
        )

    def to_legacy(self):
        """Return the original guide row."""
        return self._legacy_row


@dataclass(frozen=True)
class DesignResult:
    """Result of a gene or region design."""

    guide_table: Any
    hr_dna: Any
    checking_primers: Any
    name: Any
    chromosome: Any
    start: Any
    end: Any
    _legacy_result: tuple = field(repr=False, compare=False)
    _guide_candidates: tuple[GuideCandidate, ...] | None = field(
        default=None,
        repr=False,
        compare=False,
    )

    @classmethod
    def from_legacy(cls, legacy_result):
        """Read a seven-item ``runWeb`` result."""
        if not isinstance(legacy_result, tuple):
            raise TypeError("legacy design result must be a tuple")
        if len(legacy_result) != 7:
            raise ValueError("legacy design result must contain seven items")

        try:
            guide_candidates = tuple(
                GuideCandidate.from_legacy_row(row, legacy_result[4])
                for row in legacy_result[0]
            )
        except (IndexError, TypeError, ValueError):
            # Defer validation of old nested rows until guides are requested.
            guide_candidates = None

        return cls(
            guide_table=legacy_result[0],
            hr_dna=legacy_result[1],
            checking_primers=legacy_result[2],
            name=legacy_result[3],
            chromosome=legacy_result[4],
            start=legacy_result[5],
            end=legacy_result[6],
            _legacy_result=legacy_result,
            _guide_candidates=guide_candidates,
        )

    @property
    def guides(self):
        """Return named guide candidates."""
        if self._guide_candidates is not None:
            return self._guide_candidates
        return tuple(
            GuideCandidate.from_legacy_row(row, self.chromosome)
            for row in self.guide_table
        )

    def to_legacy(self):
        """Return the original seven-item tuple."""
        return self._legacy_result


@dataclass(frozen=True)
class OligoMatch:
    """One full-length genomic match for an oligo query."""

    chromosome: str
    pam_coordinates: tuple[int, int]
    cut_coordinates: tuple[int, int]
    strand: int
    seed: str
    pam: str

    def __post_init__(self):
        object.__setattr__(self, "pam_coordinates", tuple(self.pam_coordinates))
        object.__setattr__(self, "cut_coordinates", tuple(self.cut_coordinates))


@dataclass(frozen=True)
class OligoAnalysisResult:
    """Genome-wide oligo analysis."""

    oligo_sequence: str
    seed: str
    pam: str
    n_mismatch: int
    spedit_forward: str
    spedit_reverse: str
    has_internal_bsai: bool
    match_counts: Mapping[int, int]
    full_matches: tuple[OligoMatch, ...]

    def __post_init__(self):
        object.__setattr__(
            self,
            "match_counts",
            MappingProxyType(dict(self.match_counts)),
        )
        object.__setattr__(self, "full_matches", tuple(self.full_matches))
