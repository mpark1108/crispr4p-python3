"""Structured application results with explicit legacy compatibility."""

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class DesignResult:
    """Named top-level fields from a CRISPR4P region-design operation.

    Nested values deliberately remain in their legacy forms during this
    checkpoint. Keeping the original tuple lets existing consumers recover
    the exact positional layout without copying or reordering its contents.
    """

    guide_table: Any
    hr_dna: Any
    checking_primers: Any
    name: Any
    chromosome: Any
    start: Any
    end: Any
    _legacy_result: tuple = field(repr=False, compare=False)

    @classmethod
    def from_legacy(cls, legacy_result):
        """Build a named result from the seven-item ``runWeb`` tuple."""
        if not isinstance(legacy_result, tuple):
            raise TypeError("legacy design result must be a tuple")
        if len(legacy_result) != 7:
            raise ValueError("legacy design result must contain seven items")

        return cls(
            guide_table=legacy_result[0],
            hr_dna=legacy_result[1],
            checking_primers=legacy_result[2],
            name=legacy_result[3],
            chromosome=legacy_result[4],
            start=legacy_result[5],
            end=legacy_result[6],
            _legacy_result=legacy_result,
        )

    def to_legacy(self):
        """Return the original seven-item tuple used by legacy consumers."""
        return self._legacy_result
