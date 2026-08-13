"""Guide discovery and matching independent of query/presentation state."""

import re


SUFFIX_LENGTHS = (8, 10, 12, 14, 16, 18, 20)


def discover_region_guides(
    chromosome_sequence,
    chromosome_name,
    start,
    end,
    hit_factory,
    reverse_complement,
    seed_length=20,
):
    """Find unique NGG-adjacent guide seeds in one genomic interval.

    Ordering, strand-relative positions, and duplicate-seed filtering match
    the historical ``PrimerDesign._getUserNGGs`` implementation.
    """
    interval_sequence = chromosome_sequence[start:end + 1]
    hits_by_seed = {}
    strands = (
        (1, interval_sequence),
        (-1, reverse_complement(interval_sequence)),
    )

    for strand, sequence in strands:
        # Lookahead retains both PAMs when the interval contains GGG.
        for match in re.finditer(r"(?=GG)", sequence):
            position = match.start()
            pam = sequence[position - 1:position + 2]
            seed = sequence[
                position - seed_length - 1:position - 1
            ]
            if seed:
                hits_by_seed.setdefault(seed, []).append(
                    hit_factory(
                        chromosome_name,
                        position,
                        strand,
                        seed,
                        pam,
                    )
                )

    if not hits_by_seed:
        raise AssertionError("No nGG found in your input")

    return [hits[0] for hits in hits_by_seed.values() if len(hits) == 1]


def build_guide_primer_tuple(
    chromosome_sequence,
    start,
    end,
    guide,
    reverse_complement,
):
    """Build the legacy guide oligos and normalized PAM tuple."""
    if guide.strand == 1:
        start_index = start + 1 + guide.pos
        guide_sequence = chromosome_sequence[
            start_index - 22:start_index - 2
        ]
        pam = chromosome_sequence[start_index - 2:start_index + 1]
        pam_start = start_index - 2
    else:
        start_index = end - 1 - guide.pos
        pam = reverse_complement(
            chromosome_sequence[start_index:start_index + 3]
        )
        guide_sequence = reverse_complement(
            chromosome_sequence[start_index + 3:start_index + 23]
        )
        pam_start = start_index

    forward_oligo = (
        guide_sequence[-10:] + "gttttagagctagaaatagcaagttaaaataa"
    )
    reverse_oligo = (
        reverse_complement(guide_sequence[:10])
        + "ttcttcggtacaggttatgttttttggcaaca"
    )

    # Input slices use zero-based indexes; user-visible PAM coordinates are
    # one-based and inclusive.
    pam_coordinates = (pam_start + 1, pam_start + 3)
    return (
        guide_sequence,
        forward_oligo,
        reverse_oligo,
        pam_coordinates,
        guide.strand,
        pam,
    )


def sequences_match(first_sequence, second_sequence, allowed_mismatches):
    """Return the legacy positional mismatch comparison result."""
    if allowed_mismatches == 0:
        return first_sequence == second_sequence

    mismatch_count = len(
        [
            position
            for position in range(len(first_sequence))
            if first_sequence[position] != second_sequence[position]
        ]
    )
    return allowed_mismatches >= mismatch_count


def build_suffix_match_table(guide, genome_hits_by_suffix, n_mismatch):
    """Build CRISPR4P's cumulative 8-to-20-nt suffix match table.

    The 8-nt bucket is an independent mutable list for compatibility. Later
    buckets preserve the legacy filtering and ``set``-based deduplication.
    """
    suffix = guide.seed[-8:]
    genome_hits = list(genome_hits_by_suffix.get(suffix, ()))
    match_table = {8: genome_hits}

    for suffix_length in SUFFIX_LENGTHS[1:]:
        remaining_hits = []
        for genome_hit in genome_hits:
            if sequences_match(
                guide.seed[-suffix_length:],
                genome_hit.seed[-suffix_length:],
                n_mismatch,
            ):
                remaining_hits.append(genome_hit)
        genome_hits = list(set(remaining_hits))
        match_table[suffix_length] = genome_hits

    return guide, match_table
