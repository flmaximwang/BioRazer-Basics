from biotite.sequence.io import fasta
from biorazer.io import Converter


class SEQ2FASTA(Converter):
    def write(self, tmp: dict):
        """
        Write dictionary of sequences to a FASTA file.

        Parameters
        ----------
        tmp: dict[str, str]
            A dictionary with sequence IDs as keys and sequences as values.
        """

        fasta_file = fasta.FastaFile()
        for key, value in tmp.items():
            fasta_file[key] = value
        fasta_file.write(self.output_file)
