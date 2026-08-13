"""
Atom selection toolkit for biotite AtomArray structures.

This package collects functions that derive per-atom selections --
either as boolean masks or as index arrays -- from an
``biotite.structure.AtomArray`` and from existing masks:

- :mod:`~biorazer.structure.selector.mask`
    Transform an existing per-atom mask into a new mask:
    ``(atom_array, mask) -> mask``.
- :mod:`~biorazer.structure.selector.index`
    Derive per-atom indices from an atom array: ``atom_array -> indices``.
- :mod:`~biorazer.structure.selector.spatial`
    Build masks from geometry (distance cutoffs, interfaces).

All mask-returning functions return masks aligned with the input atom
array; index-returning functions return 1D ``int`` arrays.
"""
from .index import annotation

from .mask.annotation import (
    shrink_by_atom_name,
    extend_by_chain,
    shrink_by_element,
    extend_by_res,
    intersect_masks,
    invert_mask,
    revert_mask,
    union_masks,
)
from .index.annotation import (
    group_atoms_by_residue,
    indices_by_atom_name,
    indices_by_chain,
    indices_by_element,
    indices_by_res_name,
    indices_by_residue,
    indices_from_mask,
    indices_of_backbone,
    indices_of_heavy_atoms,
    mask_from_indices,
)
from .mask import annotation, spatial
from .mask.spatial import distance
