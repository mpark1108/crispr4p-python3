"""Guide matching computations independent of query and presentation state."""


SUFFIX_LENGTHS = (8, 10, 12, 14, 16, 18, 20)


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
