"""Cut-site annotation CLI."""

import argparse
import re
from pathlib import Path

from .annotations import GenomeAnnotations, VIABILITY_LABELS
from .coordinates import cut_from_pam as pam_cut


PROJECT_DATA_DIRECTORY = Path(__file__).resolve().parent.parent / "data"
DEFAULT_GFF = (
    PROJECT_DATA_DIRECTORY / "Schizosaccharomyces_pombe_all_chromosomes.gff3"
)
DEFAULT_VIABILITY = PROJECT_DATA_DIRECTORY / "gene_viability.tsv"

REGION_LABELS = {
    "CDS": "exon (CDS)",
    "five_prime_UTR": "5' UTR",
    "three_prime_UTR": "3' UTR",
    "intron": "intron",
    "exon": "non-coding exon",
    "gene": "gene (no child feature at this coordinate)",
}

def cut_from_pam(text, strand):
    """Parse a copied PAM and return its cut boundary."""
    numbers = [int(number) for number in re.findall(r"\d+", text)]
    if len(numbers) != 2:
        raise ValueError("PAM input must contain exactly two coordinates")
    pam_start, pam_end = numbers
    if pam_end != pam_start + 2:
        raise ValueError("a 1-based inclusive PAM must span exactly three bases")

    strand = str(strand).strip()
    if strand in {"+", "+1", "1"}:
        normalized_strand = "+"
        numeric_strand = 1
    elif strand in {"-", "-1"}:
        normalized_strand = "-"
        numeric_strand = -1
    else:
        raise ValueError("strand must be 1/+1/+ or -1/-")

    cut_after, _ = pam_cut(
        pam_start,
        pam_end,
        numeric_strand,
    )
    return pam_start, pam_end, cut_after, normalized_strand


def format_block(block):
    """Format a region block."""
    return (
        f"{REGION_LABELS[block.feature_type]} "
        f"{block.start}-{block.end} ({block.length} bp)"
    )


def format_neighbor(label, direction, neighbor):
    """Format a neighboring region."""
    prefix = f"next {label} region ({direction} genomic coordinates): "
    if neighbor is None:
        return prefix + "none within transcript"
    return (
        prefix
        + f"{format_block(neighbor.block)}; {neighbor.distance} bp from cut"
    )


def format_gene(label, nearby_gene):
    """Format a gene beside an intergenic cut."""
    prefix = f"nearest {label}-coordinate gene: "
    if nearby_gene is None:
        return prefix + "none"
    gene = nearby_gene.gene
    name = f" ({gene.name})" if gene.name else ""
    return (
        f"{prefix}{gene.gene_id}{name} {gene.start}-{gene.end}; "
        f"{nearby_gene.distance} bp from cut; gene viability (PomBase): "
        f"{VIABILITY_LABELS[gene.viability]}"
    )


def format_cut(annotation):
    """Format a cut-site report."""
    chromosome = annotation.chromosome
    cut_left, cut_right = annotation.cut_coordinates
    lines = [
        f"Cas9 cut boundary: {chromosome}:{cut_left} | "
        f"{chromosome}:{cut_right}"
    ]

    if annotation.is_intergenic:
        lines.extend(
            [
                "region: intergenic on both sides of cut",
                format_gene("lower", annotation.lower_gene),
                format_gene("higher", annotation.higher_gene),
            ]
        )
        return "\n".join(lines)

    gene_ids = sorted(gene.gene_id for gene in annotation.genes)
    if len(gene_ids) > 1:
        lines.append(
            f"overlapping genes at cut boundary: {len(gene_ids)} "
            f"({', '.join(gene_ids)})"
        )

    for index, context in enumerate(annotation.transcripts):
        if index or len(gene_ids) > 1:
            lines.append("")
        gene = context.gene
        name = f" ({gene.name})" if gene.name else ""
        lines.extend(
            [
                f"gene: {gene.gene_id}{name}",
                "gene viability (PomBase): "
                f"{VIABILITY_LABELS[gene.viability]}",
                f"transcript: {context.transcript_id}  "
                f"strand: {context.strand}",
            ]
        )

        if context.relation == "within":
            lines.extend(
                [
                    f"region block: {format_block(context.block)}",
                    "cut position within block: "
                    f"{context.lower_bases} bp toward lower coordinates; "
                    f"{context.higher_bases} bp toward higher coordinates",
                ]
            )
            if context.cds_position is not None:
                cds = context.cds_position
                lines.append(
                    f"CDS: base {cds.base}/{cds.total} "
                    f"({cds.percent:.1f}%)"
                )
        else:
            left_label = (
                format_block(context.left.block)
                if context.left is not None
                else "intergenic"
            )
            right_label = (
                format_block(context.right.block)
                if context.right is not None
                else "intergenic"
            )
            lines.append(
                f"cut crosses region boundary: {left_label} -> {right_label}"
            )

        lines.extend(
            [
                format_neighbor(
                    "upstream/5'",
                    context.upstream_direction,
                    context.upstream,
                ),
                format_neighbor(
                    "downstream/3'",
                    context.downstream_direction,
                    context.downstream,
                ),
            ]
        )

    lines.extend(
        [
            "distance 0 bp means that the named region touches the cut boundary",
        ]
    )
    return "\n".join(lines)


def build_parser():
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description="Annotate a 1-based genomic Cas9 cut using PomBase GFF3."
    )
    parser.add_argument("chromosome", nargs="?", help="for example: I, II, or III")
    parser.add_argument(
        "coordinate",
        nargs="?",
        type=int,
        help="1-based reference base immediately before the cut",
    )
    parser.add_argument(
        "--pam",
        "--crispr4p-pam",
        dest="pam",
        metavar='"START - END"',
        help=(
            "corrected CRISPR4P PAM copied as one value, "
            'e.g. "1316795 - 1316797"'
        ),
    )
    parser.add_argument(
        "--strand",
        help="guide strand shown by CRISPR4P: 1 or -1",
    )
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument(
        "--viability",
        type=Path,
        default=DEFAULT_VIABILITY,
        help="PomBase two-column gene viability TSV",
    )
    return parser


def main(argv=None):
    """Run the cut-site report."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.chromosome is None:
        args.chromosome = input("chromosome (example: III): ").strip()

    if args.coordinate is None and args.pam is None:
        raw = input(
            "1-based base immediately before cut, or PAM range (START - END): "
        ).strip()
        numbers = re.findall(r"\d+", raw)
        if len(numbers) == 1:
            args.coordinate = int(numbers[0])
        elif len(numbers) == 2:
            args.pam = raw
        else:
            parser.error("enter one coordinate or a three-base PAM range")

    if args.pam:
        if args.coordinate is not None:
            parser.error("use either a coordinate or --pam, not both")
        if args.strand is None:
            args.strand = input(
                "guide strand shown by CRISPR4P (1 or -1): "
            ).strip()
        try:
            pam_start, pam_end, cut_after, strand = cut_from_pam(
                args.pam,
                args.strand,
            )
        except ValueError as error:
            parser.error(str(error))
        print(
            f"true PAM: {args.chromosome}:{pam_start}-{pam_end}  "
            f"guide strand: {strand}"
        )
        args.coordinate = cut_after
    elif args.strand is not None:
        parser.error("--strand is only used with --pam")

    if args.coordinate < 1:
        parser.error("coordinate must be at least 1")

    annotations = GenomeAnnotations.from_files(args.gff, args.viability)
    result = annotations.annotate_cut(
        args.chromosome,
        (args.coordinate, args.coordinate + 1),
    )
    print(format_cut(result))
if __name__ == "__main__":
    main()
