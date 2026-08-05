import numpy as np
import biotite.structure as bio_struct
import hydride
from ...selector.mask import by_res_id as _by_res_id, revert_mask as _revert_mask
from ...selector.spatial import interface
from . import check


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

    interface_mask_1, interface_mask_2 = interface(
        atom_array, selection1, selection2, **interface_kwargs
    )
    interface_mask = interface_mask_1 | interface_mask_2
    interface_res_mask = _by_res_id(atom_array, interface_mask)

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
    unsat_hbonds_mask_1 = _revert_mask(
        atom_array,
        heavy_mask,
        (interface_mask_1 & NOS_mask & ~hbond_mask)[heavy_mask] & sasa_mask,
    )
    unsat_hbonds_mask_2 = _revert_mask(
        atom_array,
        heavy_mask,
        (interface_mask_2 & NOS_mask & ~hbond_mask)[heavy_mask] & sasa_mask,
    )

    return unsat_hbonds_mask_1, unsat_hbonds_mask_2, atom_array
