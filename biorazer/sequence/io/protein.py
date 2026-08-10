from biorazer.io import Converter
from biorazer.sequence.io import StrDict_Fasta
from biotite.sequence.io import fasta
from biotite.sequence import ProteinSequence

class Fasta_ProteinSequence(Converter):
    """
    Support a single-sequence fasta
    """

    def read(self):
        sequence_dict: dict = fasta.get_sequences(fasta_file=fasta.FastaFile.read(self.input_io), seq_type=ProteinSequence)
        if len(sequence_dict) > 1:
            raise ValueError("Use Fasta_ProteinSequenceDict for multi-sequence fasta")
        return list(sequence_dict.values())[0]

class Fasta_ProteinSequenceDict(Converter):
    """
    Support a multi-sequence fasta
    """

    def read(self):
        sequence_dict = fasta.get_sequences(fasta_file=fasta.FastaFile.read(self.input_io), seq_type=ProteinSequence)
        return sequence_dict

class ProteinSequenceDict_Fasta(Converter):

    def write(self, tmp: dict[ProteinSequence]):
        """
        Write dictionary of sequences to a FASTA file.

        Parameters
        ----------
        tmp: dict[str, str]
            A dictionary with sequence IDs as keys and sequences as values.
        """

        fasta_file = fasta.FastaFile()
        for key, value in tmp.items():
            fasta_file[key] = str(value)
        fasta_file.write(self.output_io)