import numpy as np
from biotite.sequence import SequenceProfile
from biotite.sequence.align import Alignment

class AlignmentHelper:

    @staticmethod
    def get(alignment: Alignment, key):
        return Alignment(alignment.sequences[key], alignment.trace[:, key])

    @staticmethod
    def concat_alignments(alignments: list[Alignment]):
        """
        Merge multiple alignments into a single alignment.

        The alignments are merged by concatenating the sequences
        in the order they are provided.

        Parameters
        ----------
        alignments : list of Alignment
            A list of `Alignment` objects to be merged.

        Returns
        -------
        Alignment
            A new `Alignment` object containing the merged sequences.
        """
        if not alignments:
            raise ValueError("The list of alignments is empty.")

        # Ensure all alignments have traces of the same length
        trace_length = alignments[0].trace.shape[0]
        for aln in alignments[1:]:
            if aln.trace.shape[0] != trace_length:
                raise ValueError("All alignments must have the same trace length.")

        trace_concat = np.concatenate([aln.trace for aln in alignments], axis=1)
        sequences_concat = []
        for aln in alignments:
            sequences_concat.extend(aln.sequences)
        alignment_concat = Alignment(sequences_concat, trace_concat)

        return alignment_concat

def parse_matrix_like(matrix_like: np.ndarray | SequenceProfile | Alignment):
    """
    Calculate position-wise entropy for a 2D-matrix like object
    """

    if isinstance(matrix_like, np.ndarray):
        matrix = matrix_like
        if matrix_like.shape.__len__() != 2:
            raise ValueError("Support only 2D matrix")
    elif isinstance(matrix_like, SequenceProfile):
        matrix = matrix_like.symbols
    elif isinstance(matrix_like, Alignment):
        profile = SequenceProfile.from_alignment(matrix_like)
        matrix = profile.symbols
    else:
        raise TypeError(f"Unsopported matrix type {type(alphabet_like)}")

    return matrix

def parse_matrix_like_with_alphabet_like(
    matrix_like: np.ndarray | SequenceProfile | Alignment,
    alphabet_like: str | list[str] | None = None,
    entropy_cutoff=1.0,
):
    """
    Calculate position-wise entropy for a 2D-matrix like object
    """

    if isinstance(matrix_like, np.ndarray):
        matrix = matrix_like
        if matrix_like.shape.__len__() != 2:
            raise ValueError("Support only 2D matrix")
        sequence_len, alphabet_len = matrix.shape
        if alphabet_like is None:
            raise ValueError("Alphabet must be provided for an array input")
        if alphabet_like.__len__() != alphabet_len:
            raise ValueError(
                f"Alphabet length {alphabet_len.__len__()} doesn't match the input matrix column number {alphabet_len}"
            )
        if isinstance(alphabet_like, str):
            alphabet = list(alphabet_like)
        elif isinstance(alphabet_like, list):
            for i in alphabet_like:
                if not isinstance(i, str):
                    raise ValueError("Alphabet list should only contain strings")
                if i.__len__() > 1:
                    raise ValueError("Alphabet should only contain 1-char symbols")
            alphabet = alphabet_like
        else:
            raise TypeError(f"Unsupported alphabet type {type(alphabet_like)}")
    elif isinstance(matrix_like, SequenceProfile):
        matrix = matrix_like.symbols
        alphabet = list(matrix_like.alphabet)
    elif isinstance(matrix_like, Alignment):
        profile = SequenceProfile.from_alignment(matrix_like)
        matrix = matrix_like.symbols
        alphabet = list(matrix_like.alphabet)
    else:
        raise TypeError(f"Unsopported matrix type {type(alphabet_like)}")

    return matrix, alphabet