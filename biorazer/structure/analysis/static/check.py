from scipy.spatial import KDTree
import numpy as np


def _to_coordinate_array(atom_coords) -> np.ndarray:
    """Normalize AtomArray-like input or raw coordinates into shape (n, 3)."""
    coords = atom_coords.coord if hasattr(atom_coords, "coord") else atom_coords
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("atom_coords must have shape (n, 3).")
    return coords


def is_intra_steric_clashed(atom_coords, cutoff: float = 1.2):
    """
    Check if any atoms within one structure are sterically clashed.

    Parameters
    ----------
    - atom_coords: shape (n, 3) or AtomArray-like
        Coordinates of atoms to check for clashes, or any object exposing
        a .coord attribute (e.g. biotite AtomArray).
    - cutoff: float, optional
        Distance below which atoms are considered to be in a steric clash.

    Returns
    ----------
    True if there is a steric clash, False otherwise.
    """

    atom_coords = _to_coordinate_array(atom_coords)

    if atom_coords.shape[0] < 2:
        raise ValueError("At least two atoms are required to check for steric clashes.")
    tree = KDTree(atom_coords)
    distances, _ = tree.query(atom_coords, k=2)

    nearest_distances = distances[:, 1]  # Exclude self-distance

    for distance in nearest_distances:
        if distance < cutoff:
            return True
    return False


def is_inter_steric_clashed(atom_coords_a, atom_coords_b, cutoff: float = 1.2):
    """
    Check if two structures have steric clashes between each other.

    Parameters
    ----------
    - atom_coords_a: shape (n, 3) or AtomArray-like
        Coordinates of the first structure, or any object exposing
        a .coord attribute (e.g. biotite AtomArray).
    - atom_coords_b: shape (m, 3) or AtomArray-like
        Coordinates of the second structure, or any object exposing
        a .coord attribute (e.g. biotite AtomArray).
    - cutoff: float, optional
        Distance below which two atoms are considered to be in a steric clash.

    Returns
    ----------
    True if there is an inter-structure steric clash, False otherwise.
    """

    atom_coords_a = _to_coordinate_array(atom_coords_a)
    atom_coords_b = _to_coordinate_array(atom_coords_b)

    if atom_coords_a.shape[0] == 0 or atom_coords_b.shape[0] == 0:
        return False

    tree_b = KDTree(atom_coords_b)
    nearest_distances, _ = tree_b.query(atom_coords_a, k=1)
    return bool(np.any(nearest_distances < cutoff))


def is_hydrided(atom_array):
    """
    Check if the atom array contains hydrogen atoms.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The biotite structure array containing the atoms.

    Returns
    -------
    bool
        True if the atom array contains hydrogen atoms, False otherwise.
    """
    return np.any(atom_array.element == "H")
