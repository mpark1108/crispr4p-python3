# PomBase annotation snapshots

CRISPR4P ships fixed copies of two PomBase datasets so cut-site annotation is
reproducible and works from a fresh clone without network access.

## Genome feature annotations

- File: `Schizosaccharomyces_pombe_all_chromosomes.gff3`
- Purpose: gene, transcript, CDS, intron, UTR, and non-coding-exon coordinates
- PomBase source directory:
  <https://www.pombase.org/data/genome_sequence_and_features/gff3/>
- Local source-file timestamp: 2026-07-30
- SHA-256:
  `88e4a26c16762c6d97e7f6a1600cd7dac193bdaeff4753cc1dd2123a79f6a025`

## Gene viability summary

- File: `gene_viability.tsv`
- Purpose: PomBase null/deletion viability summary by systematic gene ID
- PomBase documentation:
  <https://www.pombase.org/downloads/phenotype-annotations>
- Local source-file timestamp: 2026-08-07
- SHA-256:
  `e9399024327be0a2a6618c8fda1dfaef6ef72c493eb1050c7a6d86d41a3d3d09`

The viability file contains the raw states `viable`, `inviable`,
`condition-dependent`, and `unknown`.

The upstream files do not embed an exact PomBase release identifier. The
timestamps above are those of the reviewed local copies, and the hashes are
the authoritative identifiers for the snapshots used by this project.
