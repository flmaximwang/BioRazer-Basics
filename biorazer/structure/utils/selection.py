import numpy as np
from biotite import structure as bio_struct


def normalize_selection(atom_array: bio_struct.AtomArray, selection):
    """
    Normalize an optional selection into a boolean mask aligned with atom_array.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        Structure atoms used as the reference shape.
    selection : array-like or None
        If None, returns a full-True mask. Otherwise coerced to bool mask.

    Returns
    -------
    np.ndarray
        Boolean mask with the same shape as atom_array.
    """

    if selection is None:
        return np.ones(atom_array.shape, dtype=bool)

    if isinstance(selection, str) and selection in {"all", "ALL"}:
        return np.ones(atom_array.shape, dtype=bool)

    mask = np.asarray(selection, dtype=bool)
    if mask.shape != atom_array.shape:
        raise ValueError(
            "Selection mask shape mismatch: expected shape "
            f"{atom_array.shape}, got {mask.shape}."
        )

    return mask
