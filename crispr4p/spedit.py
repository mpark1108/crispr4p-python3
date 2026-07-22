SPEDIT_FORWARD_PREFIX = "CTAGAGGTCTCGGACT"
SPEDIT_FORWARD_SUFFIX = "GTTTCGAGACCCTTCC"

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_VALID_BASES = frozenset("ACGT")


def _normalize_guide(guide: str) -> str:
    """Normalize and validate a 20-nt DNA guide sequence."""
    if not isinstance(guide, str):
        raise TypeError("Guide sequence must be a string.")

    guide = guide.strip().upper()

    if len(guide) != 20:
        raise ValueError("SpEDIT requires a 20-nt guide sequence.")

    invalid = set(guide) - _VALID_BASES
    if invalid:
        raise ValueError(
            f"Guide contains invalid nucleotide characters: {sorted(invalid)}"
        )

    return guide


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an uppercase DNA sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def make_spedit_oligos(guide: str) -> tuple[str, str]:
    """
    Generate the 52-nt forward and reverse SpEDIT/pLSB
    BsaI Golden Gate oligonucleotides.
    """
    guide = _normalize_guide(guide)

    forward = SPEDIT_FORWARD_PREFIX + guide + SPEDIT_FORWARD_SUFFIX
    reverse = reverse_complement(forward)

    if len(forward) != 52 or len(reverse) != 52:
        raise RuntimeError("Internal error: SpEDIT oligos must be 52 nt.")

    return forward, reverse


def has_internal_bsai_site(guide: str) -> bool:
    """Return True if the 20-nt guide contains a BsaI recognition site."""
    guide = _normalize_guide(guide)
    return "GGTCTC" in guide or "GAGACC" in guide
