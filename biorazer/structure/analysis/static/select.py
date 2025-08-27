import warnings
import numpy as np
from scipy.spatial import KDTree
import biotite.structure as bio_struct
import hydride
from . import check, report


def revert_mask(atom_array: bio_struct.AtomArray, mask1: np.ndarray, mask2: np.ndarray):
    """
    Return
    -------
    mask: np.ndarray
        The returned mask makes atom_array[mask1][mask2] = atom_array[mask]
    """

    # mask1 选出一部分原子，mask2 在这部分再选
    # 需要返回原始 atom_array 的 mask，使得 atom_array[mask] == atom_array[mask1][mask2]
    idx1 = np.where(mask1)[0]
    idx2 = idx1[mask2]
    mask = np.zeros(atom_array.shape, dtype=bool)
    mask[idx2] = True
    return mask


def by_res_id(atom_array: bio_struct.AtomArray, mask: np.ndarray):
    """
    This function extends the original mask to include all atoms in the same residue
    """

    res_ids_in_mask = bio_struct.get_residues(atom_array[mask])[0]
    return np.isin(atom_array.res_id, res_ids_in_mask)


def mask_atoms_within_distance(
    atom_array_center: bio_struct.AtomArray,
    atom_array_query: bio_struct.AtomArray,
    distance: float,
):
    """
    Find atoms in atom_array_query within a certain distance from any atom in atom_array_center.
    """

    kdtree_center = KDTree(atom_array_center.coord)
    kdtree_around = KDTree(atom_array_query.coord)
    neighbors = kdtree_center.query_ball_tree(kdtree_around, distance)
    indices = set()
    for neighbor_indices in neighbors:
        for index in neighbor_indices:
            indices.add(index)
    indices = list(indices)
    indices.sort()
    mask = np.zeros(atom_array_query.shape, dtype=bool)
    mask[indices] = True
    return mask


def mask_interface_atoms(
    atom_array: bio_struct.AtomArray,
    selection1: np.ndarray,
    selection2: np.ndarray,
    distance_cutoff: float = 3.5,
):
    """
    This function selects interface atoms from a biotite structure array.

    Parameters
    ----------
    bio_structure: bio_struct.AtomArray
        The biotite structure array from which to select interface atoms.
    selection_1: np.ndarray
        The atom mask for the 1st selection.
    selection_2: np.ndarray
        The atom mask for the 2nd selection.

    Returns
    -------
    interface_mask_1: bio_struct.AtomArray
        The interface atom mask from the 1st selection.
    interface_mask_2: bio_struct.AtomArray
        The interface atom mask from the 2nd selection.
    atom_array: bio_struct.AtomArray
        The biotite structure array.
    """

    heavy_atom_mask = atom_array.element != "H"
    heavy_atom_indices_1 = np.where(heavy_atom_mask & selection1)[0]
    heavy_atom_indices_2 = np.where(heavy_atom_mask & selection2)[0]
    heavy_atom_coords_1 = atom_array.coord[heavy_atom_mask & selection1]
    heavy_atom_coords_2 = atom_array.coord[heavy_atom_mask & selection2]
    tree_1 = KDTree(heavy_atom_coords_1)
    tree_2 = KDTree(heavy_atom_coords_2)
    neighbors = tree_1.query_ball_tree(tree_2, distance_cutoff)
    interface_1_indices = []
    interface_2_indices = []
    for tree_1_index, tree_2_indices in enumerate(neighbors):
        if len(tree_2_indices) == 0:
            continue
        interface_1_indices.append(heavy_atom_indices_1[tree_1_index])
        for tree_2_index in tree_2_indices:
            interface_2_indices.append(heavy_atom_indices_2[tree_2_index])
    interface_mask_1 = np.zeros(atom_array.shape, dtype=bool)
    interface_mask_2 = np.zeros(atom_array.shape, dtype=bool)
    interface_mask_1[interface_1_indices] = True
    interface_mask_2[interface_2_indices] = True

    return interface_mask_1, interface_mask_2, atom_array


def mask_hbond_atoms(
    atom_array: bio_struct.AtomArray,
    selection_1: np.ndarray,
    selection_2: np.ndarray,
    **kwargs
):
    """
    This function selects hydrogen bond atoms from a biotite structure array.

    Parameters
    ----------
    atom_array: bio_struct.AtomArray
        The biotite structure array from which to select hydrogen bond atoms.
    selection_1: np.ndarray
        The atom mask for the 1st selection.
    selection_2: np.ndarray
        The atom mask for the 2nd selection.

    Returns
    -------
    hbond_mask_1: bio_struct.AtomArray
        The hydrogen bond atom mask from the 1st selection.
    hbond_mask_2: bio_struct.AtomArray
        The hydrogen bond atom mask from the 2nd selection.
    atom_array: bio_struct.AtomArray
        The biotite structure array with hydrogen atoms added if not present.
    """

    if not check.is_hydrided(atom_array):
        atom_array, _ = hydride.add_hydrogen(atom_array)
        atom_array.coord = hydride.relax_hydrogen(atom_array)

    hbonds = bio_struct.hbond(
        atom_array,
        selection1=selection_1,
        selection2=selection_2,
    )

    hbond_indices = np.unique(np.concatenate(hbonds))
    hbond_indices.sort()

    hbond_mask = np.zeros(atom_array.shape, dtype=bool)
    hbond_mask[hbond_indices] = True
    hbond_mask_1 = hbond_mask & selection_1
    hbond_mask_2 = hbond_mask & selection_2

    return hbond_mask_1, hbond_mask_2, atom_array


def mask_buried_unsat_hbond_atoms(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    sasa_kwargs=dict(sasa_cutoff=0.5, probe_radius=1.4),
    interface_kwargs=dict(distance_cutoff=3.5),
    hbond_kwargs=dict(),
):
    """
    Mask unsatisfied hydrogen bonds in the interface between two selections of atoms.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The biotite structure array containing the atoms.
    selection1 : np.ndarray
        The atom mask for the 1st selection.
    selection2 : np.ndarray
        The atom mask for the 2nd selection.
    format : str, optional
        The output format, by default "pymol".

    Warnings
    --------
    - Unsaturated heavy atoms exposed to inner cavities will not be found by this method
    """

    if not check.is_hydrided(atom_array):
        atom_array, _ = hydride.add_hydrogen(atom_array)
        atom_array.coord = hydride.relax_hydrogen(atom_array)

    interface_mask_1, interface_mask_2, _ = mask_interface_atoms(
        atom_array, selection1, selection2, **interface_kwargs
    )
    interface_mask = interface_mask_1 | interface_mask_2
    interface_res_mask = by_res_id(atom_array, interface_mask)

    hbond_mask, _, _ = mask_hbond_atoms(
        atom_array, interface_res_mask, np.ones(atom_array.shape, bool), **hbond_kwargs
    )

    NOS_mask = np.isin(atom_array.element, ["N", "O", "S"])
    heavy_mask = atom_array.element != "H"

    # 接下来, 使用无 H 的 SASA 计算 buried atoms
    sasa_kwargs_copy = sasa_kwargs.copy()
    sasa_cutoff = sasa_kwargs_copy.pop("sasa_cutoff")
    atom_array_without_h = atom_array[heavy_mask]
    atom_sasa = bio_struct.sasa(
        atom_array_without_h, vdw_radii="Single", **sasa_kwargs_copy
    )
    sasa_mask = atom_sasa < sasa_cutoff
    unsat_hbonds_mask_1 = revert_mask(
        atom_array,
        heavy_mask,
        (interface_mask_1 & NOS_mask & ~hbond_mask)[heavy_mask] & sasa_mask,
    )
    unsat_hbonds_mask_2 = revert_mask(
        atom_array,
        heavy_mask,
        (interface_mask_2 & NOS_mask & ~hbond_mask)[heavy_mask] & sasa_mask,
    )

    return unsat_hbonds_mask_1, unsat_hbonds_mask_2, atom_array
