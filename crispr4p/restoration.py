"""Wild-type restoration donors."""

from dataclasses import dataclass

from .spedit import reverse_complement


@dataclass(frozen=True, slots=True)
class RestorationDonor:
    arm_length: int
    left_arm: str
    right_arm: str

    @property
    def sequence(self):
        return self.left_arm + self.right_arm

    @property
    def reverse(self):
        return reverse_complement(self.sequence)

    @property
    def total_length(self):
        return len(self.sequence)


def build_donor(reference, cut, arm_length):
    """Join the wild-type sequence on each side of a cut."""
    cut_left, cut_right = cut
    if cut_right != cut_left + 1:
        raise ValueError("cut coordinates must describe adjacent bases")
    if (
        not isinstance(arm_length, int)
        or isinstance(arm_length, bool)
        or arm_length < 1
    ):
        raise ValueError("arm length must be a positive integer")

    reference = reference.upper()
    if not 0 < cut_left < len(reference):
        raise ValueError("cut is outside the reference sequence")
    if cut_left < arm_length or len(reference) - cut_left < arm_length:
        raise ValueError("reference does not contain complete homology arms")

    return RestorationDonor(
        arm_length=arm_length,
        left_arm=reference[cut_left - arm_length:cut_left],
        right_arm=reference[cut_left:cut_left + arm_length],
    )
