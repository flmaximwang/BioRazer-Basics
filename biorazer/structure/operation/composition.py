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


def remove_side_chains(
    atom_array: bio_struct.AtomArray,
    mask: np.ndarray | None = None,
):
    """
    Remove side chains from residues selected by ``mask``.

    Parameters
    ----------
    atom_array
        Input atom array.
    mask
        Optional per-atom boolean mask that selects residues to be converted
        to backbone-only GLY. If ``None``, all residues are selected
        (backward-compatible behavior).

    Returns
    -------
    biotite.structure.AtomArray
        A new atom array where selected residues keep only backbone atoms
        (N/CA/C/O) and are renamed to GLY. Unselected residues are unchanged.
    """
    if mask is None:
        target_mask = np.ones(len(atom_array), dtype=bool)
    else:
        target_mask = np.asarray(mask, dtype=bool)
        if target_mask.shape != (len(atom_array),):
            raise ValueError(
                "mask must be a 1D boolean array with the same length as atom_array"
            )

    backbone_mask = np.isin(
        atom_array.get_annotation("atom_name"), ["N", "CA", "C", "O"]
    )
    keep_mask = (~target_mask) | (target_mask & backbone_mask)
    out = atom_array[keep_mask]

    if len(out) == 0:
        raise ValueError("No atoms left after side-chain removal.")

    if np.any(target_mask):
        target_chain_ids = atom_array.chain_id[target_mask]
        target_res_ids = atom_array.res_id[target_mask]
        target_residues = set(zip(target_chain_ids.tolist(), target_res_ids.tolist()))

        mutate_mask = np.array(
            [(c, r) in target_residues for c, r in zip(out.chain_id, out.res_id)],
            dtype=bool,
        )
        out.res_name[mutate_mask] = "GLY"

    return out
