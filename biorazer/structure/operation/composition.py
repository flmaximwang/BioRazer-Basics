import numpy as np
import biotite.structure as bio_struct
import hydride


def add_hydrogens(atom_array: bio_struct.AtomArray):
    if not hasattr(atom_array, "bonds") or not atom_array.bonds:
        bond_list = bio_struct.connect_via_residue_names(atom_array)
        atom_array.bonds = bond_list
    if not hasattr(atom_array, "charge") or not atom_array.charge:
        atom_array.set_annotation("charge", np.zeros(len(atom_array)))
    atom_array, _ = hydride.add_hydrogen(atom_array)
    atom_array.coord = hydride.relax_hydrogen(atom_array)
    return atom_array


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
