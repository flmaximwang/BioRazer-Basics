from scipy.spatial import KDTree
import numpy as np


def is_steric_clashed(atom_coords: np.ndarray, cutoff: float = 1.2):
    """
    Check if any atoms in the provided coordinates are sterically clashed.

    Parameters
    ----------
    - atom_coords: shape (n, 3)
        Coordinates of atoms to check for clashes.
    - cutoff: float, optional
        Distance below which atoms are considered to be in a steric clash.

    Returns
    ----------
    True if there is a steric clash, False otherwise.
    """

    if len(atom_coords) < 2:
        raise ValueError("At least two atoms are required to check for steric clashes.")
    tree = KDTree(atom_coords)
    distances, _ = tree.query(atom_coords, k=2)

    nearest_distances = distances[:, 1]  # Exclude self-distance

    for distance in nearest_distances:
        if distance < cutoff:
            return True
    return False
