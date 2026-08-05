"""
Index derivation functions for biotite AtomArray structures.

The functions in this module take an atom array and return per-atom
indices -- 1D ``int`` arrays that can be used to index the atom array
directly or converted into a boolean mask via
:func:`mask_from_indices`.
"""
import numpy as np
import biotite.structure as bio_struct


def group_atoms_by_residue(atom_array):
    """
    Group the atoms of an atom array by residue.

    Returns a mapping from each residue to the indices of *all* atoms
    belonging to it::

        {("A", 1, ""): array([0, 1, 2, 3, 4]), ("A", 2, ""): array([5, 6, 7, 8])}

    A residue is identified by the tuple ``(chain_id, res_id, ins_code)``
    -- not by ``res_id`` alone, since residue numbering restarts in each
    chain and insertion codes further distinguish residues. The mapping
    entries appear in atom order of the input array (order of first
    occurrence of each residue).

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to group.

    Returns
    -------
    dict of tuple -> numpy.ndarray
        Maps each residue key ``(chain_id, res_id, ins_code)`` to a 1D
        ``int`` array with the indices of all atoms of that residue, in
        atom order of ``atom_array``.
    """
    groups = {}
    for idx, key in enumerate(
        zip(atom_array.chain_id, atom_array.res_id, atom_array.ins_code)
    ):
        groups.setdefault(key, []).append(idx)
    return {key: np.array(idxs, dtype=int) for key, idxs in groups.items()}


_BACKBONE_ATOM_NAMES = ("N", "CA", "C", "O")


def indices_of_backbone(atom_array):
    """
    Return the indices of the backbone atoms N, CA, C and O.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of all atoms whose
        ``atom_name`` is one of N, CA, C or O.
    """
    return np.where(np.isin(atom_array.atom_name, _BACKBONE_ATOM_NAMES))[0]


def indices_of_heavy_atoms(atom_array):
    """
    Return the indices of all heavy (non-hydrogen) atoms.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of all atoms whose
        ``element`` is not ``"H"``.
    """
    return np.where(atom_array.element != "H")[0]


def indices_by_atom_name(atom_array, names):
    """
    Return the indices of atoms with one of the given atom names.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.
    names : str or array-like of str
        Atom name(s) to match, e.g. ``"CA"`` or ``["N", "CA", "C", "O"]``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of matching atoms.
    """
    names = np.atleast_1d(np.asarray(names))
    return np.where(np.isin(atom_array.atom_name, names))[0]


def indices_by_element(atom_array, elements):
    """
    Return the indices of atoms with one of the given elements.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.
    elements : str or array-like of str
        Element symbol(s) to match, e.g. ``"N"`` or ``["N", "O", "S"]``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of matching atoms.
    """
    elements = np.atleast_1d(np.asarray(elements))
    return np.where(np.isin(atom_array.element, elements))[0]


def indices_by_res_name(atom_array, res_names):
    """
    Return the indices of atoms belonging to residues with one of the
    given residue names.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.
    res_names : str or array-like of str
        Residue name(s) to match, e.g. ``"ALA"`` or ``["ALA", "GLY"]``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of matching atoms.
    """
    res_names = np.atleast_1d(np.asarray(res_names))
    return np.where(np.isin(atom_array.res_name, res_names))[0]


def indices_by_residue(atom_array, residue_keys):
    """
    Return the indices of atoms belonging to one of the given residues.

    This is the flat-array counterpart of
    :func:`group_atoms_by_residue`: while the latter returns a dict
    grouping atom indices per residue, this function takes a list of
    residue keys and returns the concatenated atom indices as a single
    1D array, in atom order of ``atom_array``.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.
    residue_keys : array-like of tuple
        Residue key(s) ``(chain_id, res_id, ins_code)`` to match, as
        used by :func:`group_atoms_by_residue`, e.g. ``("A", 1, "")``
        or ``[("A", 1, ""), ("B", 2, "")]``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of all atoms of the given
        residues, in atom order of ``atom_array``. Residue keys that do
        not occur in ``atom_array`` simply contribute no atoms.
    """
    if isinstance(residue_keys, (tuple, list)) and len(residue_keys) == 0:
        return np.array([], dtype=int)
    if (
        isinstance(residue_keys, (tuple, list))
        and isinstance(residue_keys[0], (tuple, list))
    ):
        key_list = residue_keys
    else:
        key_list = [residue_keys]
    mask = np.zeros(atom_array.shape, dtype=bool)
    for key in key_list:
        chain_id, res_id, ins_code = key
        mask |= (
            (atom_array.chain_id == chain_id)
            & (atom_array.res_id == res_id)
            & (atom_array.ins_code == ins_code)
        )
    return np.where(mask)[0]


def indices_by_chain(atom_array, chain_ids):
    """
    Return the indices of atoms belonging to one of the given chains.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array to search.
    chain_ids : str or array-like of str
        Chain identifier(s) to match, e.g. ``"A"`` or ``["A", "B"]``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of matching atoms.
    """
    chain_ids = np.atleast_1d(np.asarray(chain_ids))
    return np.where(np.isin(atom_array.chain_id, chain_ids))[0]


def indices_from_mask(atom_array, mask):
    """
    Return the indices of the atoms selected by a boolean mask.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array ``mask`` is aligned with.
    mask : numpy.ndarray
        Boolean mask with the same shape as ``atom_array``.

    Returns
    -------
    numpy.ndarray
        1D ``int`` array with the indices of the ``True`` entries of
        ``mask``.

    Raises
    ------
    ValueError
        If ``mask`` does not have the same shape as ``atom_array``.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != atom_array.shape:
        raise ValueError(
            "mask must be a boolean mask with the same shape as "
            f"atom_array: expected {atom_array.shape}, got {mask.shape}"
        )
    return np.where(mask)[0]


def mask_from_indices(atom_array, indices):
    """
    Build a boolean mask that selects the given atom indices.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array the mask must be aligned with.
    indices : array-like of int
        Atom indices to mark as ``True``.

    Returns
    -------
    numpy.ndarray
        Boolean mask of shape ``atom_array.shape``, ``True`` at
        ``indices``.
    """
    mask = np.zeros(atom_array.shape, dtype=bool)
    mask[np.asarray(indices, dtype=int)] = True
    return mask
