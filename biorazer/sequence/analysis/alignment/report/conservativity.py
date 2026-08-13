import numpy as np
from biotite.sequence import SequenceProfile
from biotite.sequence.align import Alignment
from biorazer.sequence.analysis.alignment.calculation import (
    calculate_entropy_position_wise,
)
from biorazer.sequence.analysis.alignment.util import (
    parse_matrix_like,
    parse_matrix_like_with_alphabet_like,
)


def report_conserved_positions(
    matrix_like: np.ndarray | SequenceProfile | Alignment,
    entropy_cutoff=1.0,
):

    matrix = parse_matrix_like(matrix_like=matrix_like)
    entropy_position_wise = calculate_entropy_position_wise(matrix_like)
    sequence_len, alphabet_len = matrix.shape

    result = []
    for i in range(sequence_len):
        if entropy_position_wise[i] <= entropy_cutoff:
            result.append(i)

    return result


def report_conserved_sequence(
    matrix_like: np.ndarray | SequenceProfile | Alignment,
    alphabet_like: str | list[str] | None = None,
    entropy_cutoff=1.0,
):
    """
    Calculate position-wise entropy for a 2D-matrix like object
    """

    matrix, alphabet = parse_matrix_like_with_alphabet_like(
        matrix_like=matrix_like, alphabet_like=alphabet_like
    )
    sequence_len, alphabet_len = matrix.shape

    entropy_position_wise = calculate_entropy_position_wise(matrix)
    result_string = ""
    for i in range(sequence_len):
        if entropy_position_wise[i] <= entropy_cutoff:
            max_idx = np.argmax(matrix[i])
            conserved_res = alphabet[max_idx]
        else:
            conserved_res = "-"
        result_string += conserved_res

    return result_string


def report_different_positions(alignment: Alignment):
    """
    Report alignment columns where any sequence differs from the first sequence.

    The first sequence (``alignment.sequences[0]``) serves as the reference.
    For every alignment column, the residue of each remaining sequence is
    compared with the reference residue in the same column; a column is
    reported as soon as one sequence differs from the reference.

    .. note::
        The reported values are 1-based alignment column positions, not
        residue IDs and not sequence indices. The comparison is based on
        residue types (amino acids), not on any residue numbering.

        The implementation indexes the sequences directly with the column
        index (``alignment.sequences[j][i]``), so it only works correctly
        for gapless alignments, where the number of columns equals the
        sequence length.

    Parameters
    ----------
    alignment : Alignment
        A Biotite Alignment object containing the aligned sequences.

    Returns
    -------
    list of int
        1-based indices of the alignment columns where at least one sequence
        differs from the first sequence.

    Examples
    --------
    >>> import numpy as np
    >>> from biotite.sequence import ProteinSequence
    >>> from biotite.sequence.align import Alignment
    >>> from biorazer.sequence.analysis.alignment.report import report_different_res_ids
    >>> seq1 = ProteinSequence("ACDEF")
    >>> seq2 = ProteinSequence("ACDEF")
    >>> seq3 = ProteinSequence("ACDDF")
    >>> trace = np.arange(5).reshape(-1, 1).repeat(3, axis=1)
    >>> alignment = Alignment([seq1, seq2, seq3], trace)
    >>> # Only column 4 (1-based) holds a residue different from seq1 ('E' vs 'D')
    >>> report_different_res_ids(alignment)
    [4]
    """
    different_res_id_indices = []
    # trace[:, 0]: positions of the reference sequence (sequence 0) at each column
    reference_trace = alignment.trace[:, 0]
    # Iterate over the alignment columns
    for i in range(alignment.trace.shape[0]):
        # Residue of the reference sequence at this column
        seq_0_aa = alignment.sequences[0][i]
        # Compare it with the residue of every other sequence in the same column
        for seq_id in range(1, alignment.trace.shape[1]):
            seq_i_aa = alignment.sequences[seq_id][i]
            if seq_0_aa != seq_i_aa:
                # Record the 1-based column index and move on to the next column
                different_res_id_indices.append(i)
                break
    return different_res_id_indices
