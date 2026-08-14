"""Homology sequences and primer3 design."""

from primer3 import bindings as primer3


def design_primers(sequence_args, global_args):
    """Call primer3-py, including its older API spelling."""
    design = getattr(primer3, "design_primers", None)
    if design is None:
        design = primer3.designPrimers
    return design(sequence_args, global_args)


def build_hr_dna(
    chromosome_sequence,
    start,
    end,
    sequence_complement,
):
    """Build the original HR oligos and joined flanks."""
    previous_250 = chromosome_sequence[start - 250:start]
    next_250 = chromosome_sequence[end:end + 250]
    forward = previous_250[-80:] + next_250[:20]
    reverse = (
        "".join(reversed(next_250[:80]))
        + "".join(reversed(previous_250[-20:]))
    )
    reverse = sequence_complement(reverse)
    return forward, reverse, previous_250 + next_250


def checking_primers(
    chromosome_sequence,
    start,
    end,
    width,
    number_of_alternatives,
    primer_designer=design_primers,
):
    """Design checking primers around a region."""
    previous_sequence = chromosome_sequence[start - width:start]
    next_sequence = chromosome_sequence[end:end + width]
    sequence_args = {
        "SEQUENCE_ID": "MH1000",
        "SEQUENCE_TEMPLATE": previous_sequence + next_sequence,
        "SEQUENCE_INCLUDED_REGION": [0, 2 * width],
        "SEQUENCE_EXCLUDED_REGION": [[width - 80, 160]],
    }
    global_args = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_PICK_INTERNAL_OLIGO": 1,
        "PRIMER_INTERNAL_MAX_SELF_END": 8,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_MAX_POLY_X": 100,
        "PRIMER_INTERNAL_MAX_POLY_X": 100,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 12,
        "PRIMER_MAX_SELF_END": 8,
        "PRIMER_PAIR_MAX_COMPL_ANY": 12,
        "PRIMER_PAIR_MAX_COMPL_END": 8,
        "PRIMER_PRODUCT_SIZE_RANGE": [[width - 75, 2 * width]],
    }

    primer3_result = primer_designer(sequence_args, global_args)
    result_keys = (
        "PRIMER_LEFT_%s_SEQUENCE",
        "PRIMER_LEFT_%s_SEQUENCE",
        "PRIMER_RIGHT_%s_SEQUENCE",
        "PRIMER_LEFT_%s_TM",
        "PRIMER_RIGHT_%s_TM",
        "PRIMER_LEFT_%s_GC_PERCENT",
        "PRIMER_RIGHT_%s_GC_PERCENT",
        "PRIMER_PAIR_%s_PRODUCT_SIZE",
        "PRIMER_LEFT_%s_TM",
        "PRIMER_RIGHT_%s_TM",
    )

    alternatives = []
    for index in range(number_of_alternatives):
        alternative = {}
        for key_template in result_keys:
            key = key_template % index
            alternative[key] = primer3_result[key]
        alternative["negative_result"] = (
            primer3_result["PRIMER_PAIR_%s_PRODUCT_SIZE" % index]
            + (end - start)
        )
        alternatives.append(alternative)
    return alternatives
