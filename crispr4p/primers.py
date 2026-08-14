"""Homology sequences and primer3 design."""

from dataclasses import dataclass

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


def _primer_settings(product_range):
    return {
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
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_range)],
    }


class PrimerNotFoundError(RuntimeError):
    """Raised when Primer3 cannot design a checking pair."""


@dataclass(frozen=True, slots=True)
class InsertionPrimerPair:
    """One PCR pair spanning a cassette insertion."""

    forward: str
    reverse: str
    forward_tm: float
    reverse_tm: float
    wt_product_size: int
    insert_length: int

    @property
    def disrupted_product_size(self):
        return self.wt_product_size + self.insert_length


def insertion_primers(
    reference,
    cut,
    arm_length=80,
    window=300,
    insert_length=23,
    primer_designer=design_primers,
):
    """Design one checking pair outside the insertion donor arms."""
    try:
        cut_left, cut_right = cut
    except (TypeError, ValueError):
        raise ValueError("cut coordinates must contain two bases") from None

    integers = (cut_left, cut_right, arm_length, window, insert_length)
    if any(not isinstance(value, int) or isinstance(value, bool)
           for value in integers):
        raise ValueError("primer coordinates and lengths must be integers")
    if cut_right != cut_left + 1:
        raise ValueError("cut coordinates must describe adjacent bases")
    if arm_length < 1 or window <= arm_length or insert_length < 1:
        raise ValueError("primer lengths and window are invalid")

    reference = reference.upper()
    if not 0 < cut_left < len(reference):
        raise ValueError("cut is outside the reference sequence")

    left = reference[max(0, cut_left - window):cut_left]
    right = reference[cut_left:min(len(reference), cut_left + window)]
    if len(left) < arm_length or len(right) < arm_length:
        raise ValueError("reference does not contain complete homology arms")

    template = left + right
    junction = len(left)
    right_start = junction + arm_length
    left_region = junction - arm_length
    right_region = len(template) - right_start
    sequence_args = {
        "SEQUENCE_ID": "insertion_check",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_INCLUDED_REGION": [0, len(template)],
        "SEQUENCE_EXCLUDED_REGION": [[left_region, 2 * arm_length]],
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
            [0, left_region, right_start, right_region]
        ],
    }
    global_args = _primer_settings((200, 300))
    global_args["PRIMER_PICK_INTERNAL_OLIGO"] = 0
    global_args["PRIMER_NUM_RETURN"] = 1

    answer = primer_designer(sequence_args, global_args)
    if answer.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        raise PrimerNotFoundError("Primer3 could not find a checking pair")

    return InsertionPrimerPair(
        forward=answer["PRIMER_LEFT_0_SEQUENCE"],
        reverse=answer["PRIMER_RIGHT_0_SEQUENCE"],
        forward_tm=answer["PRIMER_LEFT_0_TM"],
        reverse_tm=answer["PRIMER_RIGHT_0_TM"],
        wt_product_size=answer["PRIMER_PAIR_0_PRODUCT_SIZE"],
        insert_length=insert_length,
    )


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
    global_args = _primer_settings((width - 75, 2 * width))

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
