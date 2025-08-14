import numpy as np
import biotite.structure as bio_struct


def remove_side_chains(atom_array: bio_struct.AtomArray):
    """
    Remove side chains from a biotite structure array.
    This function modifies the input array in place.
    """
    backbone_masks = np.isin(
        atom_array.get_annotation("atom_name"), ["N", "CA", "C", "O"]
    )
    backbone_struct = atom_array[backbone_masks]
    backbone_len = len(backbone_struct)
    if backbone_len == 0:
        raise ValueError("No backbone atoms found in the provided structure.")
    backbone_struct.set_annotation("res_name", ["GLY"] * backbone_len)
    return backbone_struct
