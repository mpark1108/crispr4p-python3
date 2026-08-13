"""Public application boundary for CRISPR4P design operations."""

import os
from pathlib import Path

from .crispr4p import NGG, PrimerDesign
from .models import DesignResult, OligoAnalysisResult, OligoMatch
from .spedit import has_internal_bsai_site, make_spedit_oligos


PROJECT_DATA_DIRECTORY = Path(__file__).resolve().parent.parent / "data"


class OligoLengthError(ValueError):
    """Raised when an oligo is neither a 20-nt seed nor seed plus PAM."""

    def __init__(self, sequence_length):
        self.sequence_length = sequence_length
        super().__init__(
            "Oligo sequence must be 20 bp (seed only) or 23 bp "
            f"(seed + PAM). Received length: {sequence_length}"
        )


def normalize_oligo_query(oligo_sequence):
    """Normalize a 20-nt seed or a 23-nt seed-plus-PAM query."""
    normalized_sequence = oligo_sequence.upper().strip()
    if len(normalized_sequence) == 20:
        return normalized_sequence, normalized_sequence, "NGG"
    if len(normalized_sequence) == 23:
        return normalized_sequence, normalized_sequence[:20], normalized_sequence[20:]
    raise OligoLengthError(len(normalized_sequence))


class Crispr4pService:
    """Delegate design requests to a fresh legacy engine instance.

    The service exposes named structured results while preserving the exact
    ``runWeb`` seven-item tuple through ``DesignResult.to_legacy()``.
    """

    def __init__(
        self,
        sequence_file,
        coordinates_file,
        synonyms_file,
        precomputed_folder="precomputed",
        designer_factory=PrimerDesign,
        designer_verbose=False,
        genome_index=None,
        reference_resources=None,
    ):
        self.sequence_file = os.fspath(sequence_file)
        self.coordinates_file = os.fspath(coordinates_file)
        self.synonyms_file = os.fspath(synonyms_file)
        self.precomputed_folder = os.fspath(precomputed_folder)
        self._designer_factory = designer_factory
        self.designer_verbose = designer_verbose
        self._genome_index = genome_index
        self._reference_resources = reference_resources

    @classmethod
    def from_project_data(
        cls,
        precomputed_folder="precomputed",
        designer_factory=PrimerDesign,
        designer_verbose=False,
        genome_index=None,
        reference_resources=None,
    ):
        """Create a service using the reference files shipped with CRISPR4P."""
        return cls(
            PROJECT_DATA_DIRECTORY
            / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            PROJECT_DATA_DIRECTORY / "COORDINATES.txt",
            PROJECT_DATA_DIRECTORY / "SYNONIMS.txt",
            precomputed_folder=precomputed_folder,
            designer_factory=designer_factory,
            designer_verbose=designer_verbose,
            genome_index=genome_index,
            reference_resources=reference_resources,
        )

    def design_gene(self, name, n_mismatch=0):
        """Design guides and primers for a gene name or systematic ID."""
        return self._run_design(name=name, n_mismatch=n_mismatch)

    def design_region(
        self,
        chromosome,
        start,
        end,
        strand=None,
        n_mismatch=0,
    ):
        """Design guides and primers for an explicit chromosome interval."""
        return self._run_design(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=strand,
            n_mismatch=n_mismatch,
        )

    def analyze_oligo(self, oligo_sequence, n_mismatch=0):
        """Analyze a 20-nt seed or 23-nt seed-plus-PAM genome-wide."""
        normalized_sequence, seed, pam = normalize_oligo_query(oligo_sequence)

        spedit_forward, spedit_reverse = make_spedit_oligos(seed)
        designer = self._new_designer()
        designer.getNGGsFromGenome()
        self._remember_genome_index(designer)

        query = NGG(
            chro="query",
            pos=0,
            strand=1,
            seed=seed,
            pam=pam,
        )
        _, legacy_matches = designer._single_table_worker(
            query,
            n_mismatch,
        )

        match_counts = {
            length: len(legacy_matches.get(length, []))
            for length in (8, 10, 12, 14, 16, 18, 20)
        }
        full_matches = []
        for match in legacy_matches.get(20, []):
            pam_coordinates, cut_coordinates = (
                designer.getOligoHitCoordinates(match)
            )
            full_matches.append(
                OligoMatch(
                    chromosome=match.chromosome,
                    pam_coordinates=pam_coordinates,
                    cut_coordinates=cut_coordinates,
                    strand=match.strand,
                    seed=match.seed,
                    pam=match.pam,
                )
            )

        return OligoAnalysisResult(
            oligo_sequence=normalized_sequence,
            seed=seed,
            pam=pam,
            n_mismatch=n_mismatch,
            spedit_forward=spedit_forward,
            spedit_reverse=spedit_reverse,
            has_internal_bsai=has_internal_bsai_site(seed),
            match_counts=match_counts,
            full_matches=full_matches,
        )

    def _run_design(
        self,
        name=None,
        chromosome=None,
        start=None,
        end=None,
        strand=None,
        n_mismatch=0,
    ):
        designer = self._new_designer()
        result = DesignResult.from_legacy(
            designer.runWeb(
                name=name,
                cr=chromosome,
                start=start,
                end=end,
                strand=strand,
                nMismatch=n_mismatch,
            )
        )
        self._remember_genome_index(designer)
        return result

    def _new_designer(self):
        # PrimerDesign stores query-specific mutable state, so do not reuse an
        # instance across service calls or concurrent HTTP requests.
        options = {
            "precomputed_folder": self.precomputed_folder,
            "verbose": self.designer_verbose,
        }
        if self._genome_index is not None:
            options["genome_index"] = self._genome_index
        if self._reference_resources is not None:
            options["reference_resources"] = self._reference_resources

        designer = self._designer_factory(
            self.sequence_file,
            self.coordinates_file,
            self.synonyms_file,
            **options,
        )
        self._remember_reference_resources(designer)
        return designer

    @property
    def genome_index(self):
        return self._genome_index

    @property
    def reference_resources(self):
        return self._reference_resources

    def _remember_genome_index(self, designer):
        genome_index = getattr(designer, "genome_index", None)
        if genome_index is not None:
            self._genome_index = genome_index

    def _remember_reference_resources(self, designer):
        resources = getattr(designer, "reference_resources", None)
        if resources is not None:
            self._reference_resources = resources
