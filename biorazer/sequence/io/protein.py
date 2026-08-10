from biorazer.io import Converter
from biorazer.sequence.io import StrDict_Fasta
from biotite.sequence.io import fasta
from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment

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


class A3m_Alignment(Converter):
    """Read an A3M multiple sequence alignment into a biotite Alignment."""

    def read(self):
        sequences = []
        with self._text_io(self.input_io) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                cleaned_seq = "".join(
                    [c for c in line.strip() if not c.islower()]
                )  # lowercase letters indicate insertions, which will cause inconsistencies in alignment length
                sequences.append(cleaned_seq)
        traces = Alignment.trace_from_strings(sequences)
        protein_seqs = [
            ProteinSequence(seq.replace("-", "")) for seq in sequences
        ]  # Clear all "-" gaps for biotite
        alignment = Alignment(protein_seqs, traces)
        return alignment

    def read_labels(self) -> list[str]:
        labels = []
        with self._text_io(self.input_io) as f:
            for line in f:
                if line.startswith(">"):
                    labels.append(line[1:].strip())
        return labels


class Alignment_A3m(Converter):
    """Write a biotite Alignment to A3M format."""

    def write(self, alignment: Alignment, labels: list[str] = []):
        sequences = alignment.get_gapped_sequences()
        with self._text_io(self.output_io, "w") as f:
            if len(labels) == 0:
                for i, seq in enumerate(sequences):
                    f.write(f">{i}\n")
                    f.write(str(seq) + "\n")
            else:
                for label, seq in zip(labels, sequences):
                    f.write(f">{label}\n")
                    f.write(str(seq) + "\n")