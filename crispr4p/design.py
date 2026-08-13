"""Single-query CRISPR design orchestration and legacy table construction."""

if __package__:
    from .guides import SUFFIX_LENGTHS
else:  # Direct ``python crispr4p/crispr4p.py`` compatibility.
    from guides import SUFFIX_LENGTHS


class TableSorting:
    """Legacy stable bubble sort over selected row positions."""

    def __init__(self, posList, reversed_sort):
        self.reversed = reversed_sort
        self.posList = posList

    def bubbleSort(self, alist):
        for passnum in range(len(alist) - 1, 0, -1):
            for index in range(passnum):
                if self._biggerThanTuple(alist[index], alist[index + 1]):
                    temporary = alist[index]
                    alist[index] = alist[index + 1]
                    alist[index + 1] = temporary

    def _biggerThanTuple(self, first, second):
        positions = list(range(self.posList[0], self.posList[1] + 1))
        if self.reversed:
            positions = list(reversed(positions))

        for position in positions:
            if first[position] > second[position]:
                return True
            if first[position] < second[position]:
                return False

    def sortByPosCriteria(self, table):
        self.bubbleSort(table)
        return table


def build_guide_result_table(guide_matches):
    """Build and sort the historical browser/CLI guide rows."""
    table = []
    for guide, matches_by_length in guide_matches.items():
        row = [guide.seed, guide.primer]
        row.extend(
            len(matches_by_length[length])
            for length in SUFFIX_LENGTHS
        )
        table.append(row)

    # Deliberately retains the legacy IndexError when no guide rows exist.
    return TableSorting(
        (2, len(table[0]) - 1),
        reversed_sort=True,
    ).sortByPosCriteria(table)


def finalize_design_query(
    chromosome,
    start,
    end,
    guide_matches,
    build_guide_primer,
    build_hr_dna,
    build_checking_primers,
):
    """Attach guide primers and construct the historical four-item result."""
    for guide in guide_matches:
        guide.primer = build_guide_primer(chromosome, start, end, guide)

    table = build_guide_result_table(guide_matches)
    hr_dna = build_hr_dna(chromosome, start, end)
    checking_primers = build_checking_primers(chromosome, start, end)
    return table, hr_dna, checking_primers, guide_matches


def run_design_query(
    chromosome_name,
    start,
    end,
    n_mismatch,
    regression,
    chromosomes,
    guide_matches,
    ensure_genome_index,
    discover_guides,
    match_guides,
    build_guide_primer,
    build_hr_dna,
    build_checking_primers,
):
    """Execute one legacy design query through injected core operations."""
    if not regression:
        ensure_genome_index()
    chromosome = chromosomes.get(chromosome_name, None)
    discover_guides(chromosome, start, end)
    match_guides(n_mismatch)
    return finalize_design_query(
        chromosome,
        start,
        end,
        guide_matches,
        build_guide_primer,
        build_hr_dna,
        build_checking_primers,
    )
