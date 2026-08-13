"""Read-only reference data loading independent of design queries."""

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType


@dataclass(frozen=True, slots=True, init=False, eq=False)
class chromosomeFasta:
    """One parsed legacy FASTA record.

    The historical lowercase class name and raw-record constructor are kept
    so imports from ``crispr4p.crispr4p`` remain compatible.
    """

    header: str
    sequence: str
    name: str

    def __init__(self, data):
        lines = data.split("\n")
        header = lines[0]
        object.__setattr__(self, "header", header)
        object.__setattr__(self, "sequence", "".join(lines[1:]))
        object.__setattr__(self, "name", header[:header.index(" ")])

    def __str__(self):
        return " ".join(
            [
                "chromosome:",
                self.name,
                "Length:",
                str(len(self.sequence)),
                "Header:",
                self.header,
            ]
        )


def read_fasta(sequence_file):
    """Parse a FASTA file with the exact record ordering used historically."""
    with open(sequence_file, "r", encoding="utf-8") as fasta_file:
        data = fasta_file.read()

    chromosomes = {}
    for raw_record in data.split(">"):
        if raw_record:
            chromosome = chromosomeFasta(raw_record)
            chromosomes[chromosome.name] = chromosome
    return MappingProxyType(chromosomes)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class AnnotationParser:
    """Read-only coordinate and synonym records with the legacy lookup API."""

    coordinates_: tuple
    synonims_: tuple

    def __init__(self, coordinates_txt, synonims_txt):
        coordinates = self.readCoordinates_(coordinates_txt)
        synonyms = self.readSynonims_(synonims_txt)
        object.__setattr__(
            self,
            "coordinates_",
            tuple(tuple(row) for row in coordinates),
        )
        object.__setattr__(
            self,
            "synonims_",
            tuple(tuple(row) for row in synonyms),
        )

    def readCoordinates_(self, coordinates_txt):
        with open(coordinates_txt, "r", encoding="utf-8") as coordinates_file:
            data = [line.rstrip() for line in coordinates_file.readlines()][1:]
        return [line.split("\t") for line in data]

    def readSynonims_(self, synonims_txt):
        with open(synonims_txt, "r", encoding="utf-8") as synonyms_file:
            data = [line.rstrip() for line in synonyms_file.readlines()][2:]
        rows = [line.split("\t") for line in data]
        return [[value for value in row if value] for row in rows]

    def normalize_name(self, name):
        name = name.upper().strip()
        if not any(row[0] == name for row in self.coordinates_):
            name = name.lower()
        return name

    def getCoordsFromName(self, name):
        input_name = name
        name = self.normalize_name(name)

        try:
            found = next(row for row in self.synonims_ if name in row)[0]
        except StopIteration:
            raise Exception(
                '"%s" name not found, check the name is correct' % input_name
            )

        coordinates = next(
            row for row in self.coordinates_ if row[0] == found
        )[1:]
        # The legacy API returned a new list for every successful lookup.
        return list(coordinates)


@dataclass(frozen=True, slots=True, eq=False)
class ReferenceResources:
    """Query-independent FASTA and name/coordinate annotation data."""

    chromosomes: Mapping
    annotations: AnnotationParser

    def __post_init__(self):
        # Snapshot even an incoming proxy so no external backing dictionary
        # can mutate a shared resource after construction.
        object.__setattr__(
            self,
            "chromosomes",
            MappingProxyType(dict(self.chromosomes)),
        )

    @classmethod
    def from_files(cls, sequence_file, coordinates_file, synonyms_file):
        return cls(
            chromosomes=read_fasta(sequence_file),
            annotations=AnnotationParser(coordinates_file, synonyms_file),
        )
