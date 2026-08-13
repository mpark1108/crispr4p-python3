"""Coordinate conversions independent of web, CLI, and legacy engine state."""


def cas9_cut_boundary_from_pam(pam_start, pam_end, strand):
    """Return the 1-based reference bases immediately around the Cas9 cut.

    ``pam_start`` and ``pam_end`` are a 1-based inclusive PAM interval. The
    returned ``(cut_left, cut_right)`` identifies the two adjacent reference
    bases separated by the cut.
    """
    if strand == 1:
        cut_left = pam_start - 4
    elif strand == -1:
        cut_left = pam_end + 3
    else:
        raise ValueError("Oligo hit strand must be 1 or -1")
    return cut_left, cut_left + 1


def indexed_hit_coordinates(
    hit_position,
    hit_strand,
    hit_pam,
    chromosome_sequence,
    reverse_complement,
):
    """Normalize one legacy genome-index hit to reference PAM/cut coordinates.

    ``hit_position`` is the zero-based first G/A position in whichever strand
    sequence was indexed. The resulting PAM is checked against the supplied
    forward-reference chromosome sequence before coordinates are returned.
    """
    if hit_strand == 1:
        pam_start = hit_position
        pam_end = hit_position + 2
    elif hit_strand == -1:
        chromosome_length = len(chromosome_sequence)
        pam_start = chromosome_length - hit_position - 1
        pam_end = chromosome_length - hit_position + 1
    else:
        raise ValueError("Oligo hit strand must be 1 or -1")

    reference_pam = chromosome_sequence[pam_start - 1:pam_end]
    if hit_strand == -1:
        reference_pam = reverse_complement(reference_pam)
    if reference_pam != hit_pam:
        raise ValueError("Normalized PAM coordinates do not match the FASTA")

    pam_coordinates = (pam_start, pam_end)
    cut_coordinates = cas9_cut_boundary_from_pam(
        pam_start,
        pam_end,
        hit_strand,
    )
    return pam_coordinates, cut_coordinates
