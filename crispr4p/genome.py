"""Genome-wide PAM indexing independent of web and CLI concerns."""

import re
from collections.abc import Mapping
from types import MappingProxyType


class GenomePamIndex(Mapping):
    """Read-only suffix lookup for PAM-adjacent genomic hits.

    Bucket order and hit order intentionally match the legacy scanner. Hit
    payloads are supplied by the caller so CRISPR4P can retain its historical
    ``NGG`` class and cached-pickle compatibility during this extraction.
    """

    PAM_SUFFIXES = ("GG", "AG")
    COMPLEMENTS = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }

    def __init__(self, buckets):
        frozen_buckets = {
            suffix: tuple(hits)
            for suffix, hits in buckets.items()
        }
        self._buckets = MappingProxyType(frozen_buckets)
        self.hit_count = sum(len(hits) for hits in frozen_buckets.values())

    @classmethod
    def build(cls, chromosome_sequences, hit_factory, seed_length=20):
        """Build an index using the legacy NGG/NAG scanning algorithm."""
        buckets = {}
        for chromosome, forward_sequence in chromosome_sequences.items():
            strands = (
                (1, forward_sequence),
                (-1, cls.reverse_complement(forward_sequence)),
            )
            for strand, sequence in strands:
                for pam_suffix in cls.PAM_SUFFIXES:
                    # Lookahead preserves overlapping PAMs such as both GG
                    # starts in GGG, matching the corrected legacy scanner.
                    for match in re.finditer(f"(?={pam_suffix})", sequence):
                        position = match.start()
                        pam = sequence[position - 1:position + 2]
                        seed = sequence[
                            position - seed_length - 1:position - 1
                        ]
                        suffix = seed[-8:]
                        buckets.setdefault(suffix, []).append(
                            hit_factory(
                                chromosome,
                                position,
                                strand,
                                seed,
                                pam,
                            )
                        )

        return cls(buckets)

    @classmethod
    def reverse_complement(cls, sequence):
        return "".join(cls.COMPLEMENTS[base] for base in sequence)[::-1]

    @property
    def by_suffix(self):
        """Expose the read-only mapping for legacy ``PrimerDesign.NGGs``."""
        return self._buckets

    def __getitem__(self, suffix):
        return self._buckets[suffix]

    def __iter__(self):
        return iter(self._buckets)

    def __len__(self):
        return len(self._buckets)
