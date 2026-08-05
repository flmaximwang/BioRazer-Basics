"""
Spatial (geometry-based) mask functions for biotite AtomArray structures.

The functions in this module derive per-atom boolean masks from atomic
coordinates and distance cutoffs. They are aligned with the input atom
array, like the mask transforms in :mod:`biorazer.structure.selector.mask`.
"""
import numpy as np
from scipy.spatial import KDTree


def within_distance(atom_array_query, atom_array_center, distance):
    """
    Mask the atoms of ``atom_array_query`` that lie within ``distance``
    of any atom of ``atom_array_center``.

    Parameters
    ----------
    atom_array_query : biotite.structure.AtomArray
        The atom array the returned mask is aligned with.
    atom_array_center : biotite.structure.AtomArray
        The atom array whose atoms act as the distance centers.
    distance : float
        Distance cutoff in Å.

    Returns
    -------
    numpy.ndarray
        Boolean mask on ``atom_array_query``, ``True`` for atoms within
        ``distance`` of at least one atom of ``atom_array_center``.
    """
    kdtree_center = KDTree(atom_array_center.coord)
    kdtree_query = KDTree(atom_array_query.coord)
    neighbors = kdtree_center.query_ball_tree(kdtree_query, distance)
    indices = sorted({i for neighbor_indices in neighbors for i in neighbor_indices})
    mask = np.zeros(atom_array_query.shape, dtype=bool)
    mask[indices] = True
    return mask


def interface(atom_array, selection1, selection2, distance_cutoff=3.5):
    """
    Mask the atoms of two selections that form a contact interface.

    An atom belongs to the interface if it lies within
    ``distance_cutoff`` of at least one heavy atom of the other
    selection. Hydrogen atoms are ignored for interface detection.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array both selections are defined on.
    selection1, selection2 : numpy.ndarray
        Boolean masks selecting the two sides of the interface.
    distance_cutoff : float, optional
        Distance cutoff in Å.

    Returns
    -------
    (numpy.ndarray, numpy.ndarray)
        The interface masks of ``selection1`` and ``selection2``
        respectively, both aligned with ``atom_array``.
    """
    heavy_atom_mask = atom_array.element != "H"
    heavy_atom_indices_1 = np.where(heavy_atom_mask & selection1)[0]
    heavy_atom_indices_2 = np.where(heavy_atom_mask & selection2)[0]
    heavy_atom_coords_1 = atom_array.coord[heavy_atom_indices_1]
    heavy_atom_coords_2 = atom_array.coord[heavy_atom_indices_2]
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

    return interface_mask_1, interface_mask_2
