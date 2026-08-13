import numpy as np
from numbers import Number
from biotite.sequence import SequenceProfile
from biotite.sequence.align import Alignment
from .scaler import calculate_entropy
from biorazer.sequence.analysis.alignment.util import parse_matrix_like

def calculate_entropy_position_wise(matrix_like: np.ndarray | SequenceProfile | Alignment, base=2):
    """
    Calculate position-wise entropy for a 2D-matrix like object
    """

    matrix = parse_matrix_like(matrix_like=matrix_like)
    sequence_len, alphabet_len = matrix.shape

    result = np.zeros(sequence_len)
    for i in range(sequence_len):
        result[i] = calculate_entropy(matrix[i], base=base)

    return result