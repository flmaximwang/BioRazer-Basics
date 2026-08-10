from biotite.sequence import ProteinSequence
from biorazer.io import Converter
from .util import Alignment


class A3M2ALIGN(Converter):

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


class ALIGN2A3M(Converter):

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
