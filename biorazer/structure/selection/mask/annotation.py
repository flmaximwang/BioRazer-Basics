"""
Mask transformation functions for biotite AtomArray structures.

All functions in this module follow the same contract: they take an
atom array and one or more per-atom boolean masks, and return a new
boolean mask of the same shape, aligned with the input atom array.
They are used to operate on and compose existing masks.

Four families of transforms are provided:

- *composition* -- :func:`revert_mask` maps a nested mask selection
  (a mask applied to a masked subset) back into the coordinate space of
  the full atom array;
- *expansion* -- :func:`by_res_id`, :func:`by_res_name` and
  :func:`by_chain_id` extend the selection to atoms sharing a structural
  property (residue number, residue type, chain) with the
  already-selected atoms;
- *refinement* -- :func:`by_atom_name` and :func:`by_element` restrict
  the selection to atoms with a given atom name or element;
- *combination* -- :func:`invert_mask`, :func:`union_masks` and
  :func:`intersect_masks` combine masks with boolean algebra.
"""
import numpy as np
import biotite.structure as bio_struct
from biotite.structure import AtomArray


def _normalize_mask(atom_array, mask, name="mask"):
    """
    Coerce an argument into a boolean mask aligned with ``atom_array``.

    Raises a descriptive ``ValueError`` when the shape does not match,
    so that callers fail fast instead of producing silently wrong masks.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != atom_array.shape:
        raise ValueError(
            f"{name} must be a boolean mask with the same shape as "
            f"atom_array: expected {atom_array.shape}, got {mask.shape}"
        )
    return mask


def revert_mask(atom_array, mask1, mask2):
    """
    Compose two nested masks back into the space of the full atom array.

    ``mask2`` is defined on the subset ``atom_array[mask1]``. This
    function returns a mask on the full ``atom_array`` that selects
    exactly the atoms ``mask2`` picks out of that subset::

        atom_array[revert_mask(atom_array, mask1, mask2)] ==
        atom_array[mask1][mask2]

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array both masks are defined on.
    mask1 : numpy.ndarray
        Boolean mask selecting a subset of ``atom_array``.
    mask2 : numpy.ndarray
        Boolean mask defined on ``atom_array[mask1]``, selecting a
        subset of the already-selected atoms.

    Returns
    -------
    numpy.ndarray
        Boolean mask on the full ``atom_array`` such that
        ``atom_array[mask] == atom_array[mask1][mask2]``.

    Raises
    ------
    ValueError
        If ``mask1`` does not match ``atom_array``, or ``mask2`` does
        not match ``atom_array[mask1]``.
    """
    mask1 = _normalize_mask(atom_array, mask1, "mask1")
    idx1 = np.where(mask1)[0]
    mask2 = np.asarray(mask2, dtype=bool)
    if mask2.shape != idx1.shape:
        raise ValueError(
            "mask2 must be a boolean mask with the same shape as "
            f"atom_array[mask1]: expected {idx1.shape}, got {mask2.shape}"
        )
    mask = np.zeros(atom_array.shape, dtype=bool)
    mask[idx1[mask2]] = True
    return mask


def invert_mask(atom_array, mask):
    """
    Return the complement of a mask.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array.
    mask : numpy.ndarray
        Boolean mask to invert.

    Returns
    -------
    numpy.ndarray
        Boolean mask that is ``True`` exactly where ``mask`` is
        ``False``.
    """
    mask = _normalize_mask(atom_array, mask)
    return ~mask


def union_masks(atom_array, *masks):
    """
    Return the logical OR of several masks.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array all masks are aligned with.
    *masks : numpy.ndarray
        One or more boolean masks of the same shape as ``atom_array``.

    Returns
    -------
    numpy.ndarray
        Boolean mask that is ``True`` where at least one input mask is
        ``True``.
    """
    if len(masks) == 0:
        raise ValueError("union_masks requires at least one mask")
    return np.logical_or.reduce(
        [_normalize_mask(atom_array, m, f"masks[{i}]") for i, m in enumerate(masks)]
    )


def intersect_masks(atom_array, *masks):
    """
    Return the logical AND of several masks.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array all masks are aligned with.
    *masks : numpy.ndarray
        One or more boolean masks of the same shape as ``atom_array``.

    Returns
    -------
    numpy.ndarray
        Boolean mask that is ``True`` where all input masks are
        ``True``.
    """
    if len(masks) == 0:
        raise ValueError("intersect_masks requires at least one mask")
    return np.logical_and.reduce(
        [_normalize_mask(atom_array, m, f"masks[{i}]") for i, m in enumerate(masks)]
    )

def extend_by_res(atom_array, mask):
    """
    Extend a mask to all atoms of the residues it selects.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array.
    mask : numpy.ndarray
        Boolean mask whose selected atoms determine the residues to
        extend to.

    Returns
    -------
    numpy.ndarray
        Boolean mask that is ``True`` for every atom whose ``res_id``
        occurs among the residues selected by ``mask``.

    Notes
    -----
    The extension is based on ``res_id`` alone -- chain and insertion
    code are not considered. Residues with the same number in other
    chains (or a different ``ins_code``) are therefore included as well.
    """
    mask = _normalize_mask(atom_array, mask)
    res_ids_in_mask = bio_struct.get_residues(atom_array[mask])[0]
    return np.isin(atom_array.res_id, res_ids_in_mask)


def extend_by_chain(atom_array, mask):
    """
    Extend a mask to the whole chains it touches.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array.
    mask : numpy.ndarray
        Boolean mask whose selected atoms determine the chains to
        extend to.

    Returns
    -------
    numpy.ndarray
        Boolean mask that is ``True`` for every atom whose ``chain_id``
        occurs among the atoms selected by ``mask``.
    """
    mask = _normalize_mask(atom_array, mask)
    chain_ids_in_mask = np.unique(atom_array.chain_id[mask])
    return np.isin(atom_array.chain_id, chain_ids_in_mask)


def shrink_by_atom_name(atom_array, mask, names):
    """
    Restrict a mask to atoms with one of the given atom names.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array.
    mask : numpy.ndarray
        Boolean mask to restrict.
    names : str or array-like of str
        Atom name(s) to keep, e.g. ``"CA"`` or ``["N", "CA", "C", "O"]``.

    Returns
    -------
    numpy.ndarray
        ``mask & isin(atom_name, names)``.
    """
    mask = _normalize_mask(atom_array, mask)
    return mask & np.isin(atom_array.atom_name, names)


def shrink_by_element(atom_array, mask, elements):
    """
    Restrict a mask to atoms with one of the given elements.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The reference atom array.
    mask : numpy.ndarray
        Boolean mask to restrict.
    elements : str or array-like of str
        Element symbol(s) to keep, e.g. ``"N"`` or ``["N", "O", "S"]``.

    Returns
    -------
    numpy.ndarray
        ``mask & isin(element, elements)``.
    """
    mask = _normalize_mask(atom_array, mask)
    return mask & np.isin(atom_array.element, elements)

def shrink_to_bb(atom_array:AtomArray, mask: np.ndarray, bb_atom_names = ["C", "N", "O", "CA"]):
    """
    Return a mask the contains only the backbone
    """

    mask = _normalize_mask(atom_array=atom_array, mask=mask)
    return np.isin(atom_array.atom_name, bb_atom_names) & mask


def shrink_to_sc(atom_array:AtomArray, mask: np.ndarray, bb_atom_names = ["C", "N", "O", "CA"]):
    """
    Return a mask that contains only the sidechain
    """

    mask = _normalize_mask(atom_array=atom_array, mask=mask)
    return ~np.isin(atom_array.atom_name, bb_atom_names) & mask
