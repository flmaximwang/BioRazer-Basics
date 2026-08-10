"""
Internal helper functions for structure operations.

The functions in this module are implementation details of the public
operations in this package (e.g. :func:`composition.replace_side_chains`)
and are not part of the public API.
"""
import numpy as np


def _ensure_common_annotations(atom_array_a, atom_array_b):
    """
    Align two atom arrays to a common set of annotation categories.

    Returns both arrays carrying the union of all annotation categories;
    categories that exist in only one array are added to the other side,
    filled with default values (zeros for numeric annotations, empty
    strings for text annotations).

    Why this is needed
    ------------------
    ``biotite.structure.concatenate()`` (and therefore ``+``) transfers
    only annotation categories that exist in *all* concatenated
    elements; a category missing from any single element is silently
    dropped. Structures loaded from different sources may carry
    different annotation sets -- e.g. a PDB-loaded array has
    ``occupancy`` and ``b_factor`` while a programmatically constructed
    array does not. Splicing atoms of such arrays without alignment
    would silently lose the ``occupancy``/``b_factor`` values of the
    whole output. This function guarantees that every downstream
    concatenation sees identical categories, so no annotation is lost.

    Parameters
    ----------
    atom_array_a, atom_array_b : biotite.structure.AtomArray
        The atom arrays to align.

    Returns
    -------
    tuple of biotite.structure.AtomArray
        The two input arrays if their annotation categories already
        match; otherwise copies of both arrays, each extended with the
        categories missing on its side, filled with default values.
    """
    cats_a = set(atom_array_a.get_annotation_categories())
    cats_b = set(atom_array_b.get_annotation_categories())
    missing_a = cats_b - cats_a
    missing_b = cats_a - cats_b
    if not missing_a and not missing_b:
        return atom_array_a, atom_array_b
    out_a = atom_array_a.copy()
    out_b = atom_array_b.copy()
    for cat in missing_a:
        dtype = atom_array_b.get_annotation(cat).dtype
        out_a.set_annotation(cat, _annotation_defaults(len(out_a), dtype))
    for cat in missing_b:
        dtype = atom_array_a.get_annotation(cat).dtype
        out_b.set_annotation(cat, _annotation_defaults(len(out_b), dtype))
    return out_a, out_b


def _annotation_defaults(length, dtype):
    """
    Create a default-filled 1D array for a missing annotation category.

    Numeric dtypes are filled with zeros, text dtypes with empty
    strings.

    Why this is needed
    ------------------
    ``biotite.structure.AtomArray.set_annotation()`` requires an actual
    array to create a new category -- there is no "undefined" state.
    When aligning arrays with different annotation sets (see
    :func:`_ensure_common_annotations`), every missing category has to
    be materialized with placeholder values. The zero representation is
    the neutral value for an annotation that carries no information, so
    a synthesized category can never be confused with real measured
    data. The placeholder also adopts the dtype the category has in the
    array where it already exists, so that the later ``concatenate()``
    does not upcast or truncate values.

    Parameters
    ----------
    length : int
        Number of atoms the returned array must cover.
    dtype : numpy.dtype
        Dtype of the annotation category in the array where it already
        exists; the placeholder adopts this dtype.

    Returns
    -------
    numpy.ndarray
        1D array of shape ``(length,)`` filled with the zero
        representation of ``dtype`` (``0`` for numerics, ``""`` for
        strings).
    """
    if np.issubdtype(dtype, np.str_):
        return np.full(length, "", dtype=dtype)
    return np.zeros(length, dtype=dtype)


def _selected_residues(groups, mask, mask_label):
    """
    Return the residues that are fully selected by a residue-level mask.

    A residue counts as selected only if *all* of its atoms are ``True``
    in ``mask``. The returned keys are ordered by first appearance,
    which defines the pairing order between the backbone and the implant
    side of a mask pair in :func:`replace_side_chains`.

    Why this is needed
    ------------------
    Besides extracting the selected residues, this function enforces the
    residue-level contract of :func:`replace_side_chains`: a partially
    selected residue indicates a mask built at atom granularity or a
    selection mistake, and grafting half a residue would produce a
    corrupted structure. Failing fast here turns such mistakes into a
    descriptive ``ValueError`` instead of silently wrong output.

    Parameters
    ----------
    groups : dict of tuple -> numpy.ndarray
        Residue-to-atom-index mapping as produced by
        :func:`biorazer.structure.selector.index.group_atoms_by_residue`.
    mask : numpy.ndarray
        1D boolean mask with one entry per atom of the source atom
        array.
    mask_label : str
        Human-readable name of the mask (e.g. ``"mask_map[0]"``),
        included in error messages to point at the offending argument.

    Returns
    -------
    list of tuple
        Residue keys ``(chain_id, res_id, ins_code)`` whose atoms are
        all selected, in order of first appearance.

    Raises
    ------
    ValueError
        If any residue is only partially selected.
    """
    selected = []
    for key, idxs in groups.items():
        values = mask[idxs]
        if values.all():
            selected.append(key)
        elif values.any():
            raise ValueError(
                f"{mask_label} is not residue-level: residue {key} is only "
                "partially selected"
            )
    return selected
