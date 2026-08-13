"""Public application boundary for CRISPR4P design operations."""

import os
from pathlib import Path

from .crispr4p import PrimerDesign
from .models import DesignResult


PROJECT_DATA_DIRECTORY = Path(__file__).resolve().parent.parent / "data"


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
    ):
        self.sequence_file = os.fspath(sequence_file)
        self.coordinates_file = os.fspath(coordinates_file)
        self.synonyms_file = os.fspath(synonyms_file)
        self.precomputed_folder = os.fspath(precomputed_folder)
        self._designer_factory = designer_factory

    @classmethod
    def from_project_data(
        cls,
        precomputed_folder="precomputed",
        designer_factory=PrimerDesign,
    ):
        """Create a service using the reference files shipped with CRISPR4P."""
        return cls(
            PROJECT_DATA_DIRECTORY
            / "Schizosaccharomyces_pombe.ASM294v2.26.dna.toplevel.fa",
            PROJECT_DATA_DIRECTORY / "COORDINATES.txt",
            PROJECT_DATA_DIRECTORY / "SYNONIMS.txt",
            precomputed_folder=precomputed_folder,
            designer_factory=designer_factory,
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

    def _run_design(
        self,
        name=None,
        chromosome=None,
        start=None,
        end=None,
        strand=None,
        n_mismatch=0,
    ):
        # PrimerDesign stores query-specific mutable state, so do not reuse an
        # instance across service calls or concurrent HTTP requests.
        designer = self._designer_factory(
            self.sequence_file,
            self.coordinates_file,
            self.synonyms_file,
            precomputed_folder=self.precomputed_folder,
        )
        return DesignResult.from_legacy(
            designer.runWeb(
                name=name,
                cr=chromosome,
                start=start,
                end=end,
                strand=strand,
                nMismatch=n_mismatch,
            )
        )
