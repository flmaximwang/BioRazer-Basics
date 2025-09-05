from dataclasses import dataclass
from abc import abstractmethod
from pathlib import Path
from collections import OrderedDict
from .basic import Entry, Library


@dataclass
class SequenceEntry(Entry):
    """
    A sequence entry represents an entry based on a biological sequence.
    This type of entry is useful for designing proteins or peptides from scratch

    When constructing a SequenceEntry, from samples, the samples DataFrame should contain a column "sequence",
    which can either contain a single string (for single-sequence design) or a dictionary of strings (for multi-sequence design).
    """

    sequence: dict[str] = None

    def get_basic_info(self):
        sequence = self.get_samples().iloc[0]["sequence"]
        if isinstance(sequence, str):
            self.sequence = {"default": sequence}
        else:
            self.sequence = self.get_samples().iloc[0]["sequence"]

    def get_sequence(self) -> dict[str]:
        self.sequence: dict[str]
        return self.sequence


class SequenceLibrary(Library):
    """
    A sequence library represents a library of sequence-based entries.
    This type of library is useful for designing multiple proteins or peptides from scratch
    """

    entry_type = SequenceEntry

    def get_sequences(self):
        """
        Write the sequences of all entries in the library to a FASTA file.

        Parameters
        ----------
        fasta: str or Path
            The path to the output FASTA file.
        """
        self.entries: list[SequenceEntry]
        res = []
        for entry in self.entries:
            seq_dict = entry.get_sequence()
            for name, seq in seq_dict.items():
                res.append((f"{entry.marker}_{name}", seq))
        res = OrderedDict(res)
        return res
