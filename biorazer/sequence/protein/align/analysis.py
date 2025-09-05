import numpy as np
from biotite.sequence.align import Alignment
from .alignment import AlignmentHelper


def report_different_res_ids(alignment: Alignment):
    """
    Report the indices of sequences in the alignment that have different residue IDs.

    Parameters
    ----------
    alignment : Alignment
        A Biotite Alignment object containing the aligned sequences.

    Returns
    -------
    list of int
        A list of indices of sequences that have different residue IDs.
    """
    different_res_id_indices = []
    reference_trace = alignment.trace[:, 0]
    for i in range(alignment.trace.shape[0]):
        seq_0_aa = alignment.sequences[0][i]
        for seq_id in range(1, alignment.trace.shape[1]):
            seq_i_aa = alignment.sequences[seq_id][i]
            if seq_0_aa != seq_i_aa:
                different_res_id_indices.append(i + 1)
                break
    return different_res_id_indices
