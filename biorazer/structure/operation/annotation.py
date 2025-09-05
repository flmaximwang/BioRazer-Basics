import numpy as np
import biotite.structure as bio_struct


def get_renumbered_res_ids(
    atom_array: bio_struct.AtomArray, res_id_start=1
) -> np.ndarray:
    """
    Get renumbered residue IDs from an AtomArray.

    This function returns a list of residue IDs that are renumbered
    based on the unique residues present in the AtomArray.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The AtomArray from which to extract residue IDs.

    Returns
    -------
    np.array
        An array of renumbered residue IDs, starting from `res_id_start`.
        The IDs are unique and ordered based on the appearance of residues
        in the AtomArray.
    """
    old_res_ids = atom_array.get_annotation("res_id")
    res_id_shift = old_res_ids[0] - res_id_start
    new_res_ids = old_res_ids - res_id_shift
    return new_res_ids
