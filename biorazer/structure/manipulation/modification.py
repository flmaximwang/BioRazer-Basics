import re
import numpy as np
import biotite.structure as bio_struct
import hydride

from biorazer.database.amino_acid import AMINO_ACIDS_1TO3_UPPER
from ..selection.index.annotation import group_atoms_by_residue
from .util import (
    _ensure_common_annotations,
    _selected_residues,
)


def add_hydrogens(atom_array: bio_struct.AtomArray):
    if not hasattr(atom_array, "bonds") or not atom_array.bonds:
        bond_list = bio_struct.connect_via_residue_names(atom_array)
        atom_array.bonds = bond_list
    if not hasattr(atom_array, "charge") or not atom_array.charge:
        atom_array.set_annotation("charge", np.zeros(len(atom_array)))
    atom_array, _ = hydride.add_hydrogen(atom_array)
    atom_array.coord = hydride.relax_hydrogen(atom_array)
    return atom_array


def remove_side_chains(
    atom_array: bio_struct.AtomArray,
    mask: np.ndarray | None = None,
):
    """
    Remove side chains from residues selected by ``mask``.

    Parameters
    ----------
    atom_array
        Input atom array.
    mask
        Optional per-atom boolean mask that selects residues to be converted
        to backbone-only GLY. If ``None``, all residues are selected
        (backward-compatible behavior).

    Returns
    -------
    biotite.structure.AtomArray
        A new atom array where selected residues keep only backbone atoms
        (N/CA/C/O) and are renamed to GLY. Unselected residues are unchanged.
    """
    if mask is None:
        target_mask = np.ones(len(atom_array), dtype=bool)
    else:
        target_mask = np.asarray(mask, dtype=bool)
        if target_mask.shape != (len(atom_array),):
            raise ValueError(
                "mask must be a 1D boolean array with the same length as atom_array"
            )

    backbone_mask = np.isin(
        atom_array.get_annotation("atom_name"), ["N", "CA", "C", "O"]
    )
    keep_mask = (~target_mask) | (target_mask & backbone_mask)
    out = atom_array[keep_mask]

    if len(out) == 0:
        raise ValueError("No atoms left after side-chain removal.")

    if np.any(target_mask):
        target_chain_ids = atom_array.chain_id[target_mask]
        target_res_ids = atom_array.res_id[target_mask]
        target_residues = set(zip(target_chain_ids.tolist(), target_res_ids.tolist()))

        mutate_mask = np.array(
            [(c, r) in target_residues for c, r in zip(out.chain_id, out.res_id)],
            dtype=bool,
        )
        out.res_name[mutate_mask] = "GLY"

    return out


def replace_side_chains(
    atom_array_backbone: bio_struct.AtomArray,
    atom_array_implant: bio_struct.AtomArray,
    mask_map: list[np.ndarray],
) -> bio_struct.AtomArray:
    """
    Replace side chains of residues selected in ``atom_array_backbone`` with
    side chains taken from residues selected in ``atom_array_implant``.

    ``mask_map`` is a flat list of per-atom boolean masks grouped in pairs::

        [mask_backbone_1, mask_implant_1, mask_backbone_2, mask_implant_2, ...]

    Each pair defines a residue-to-residue relationship: the k-th residue
    selected by ``mask_backbone_i`` (in atom order) receives the side chain
    of the k-th residue selected by ``mask_implant_i``. Hence both masks of a
    pair must select the same number of residues. Masks are residue-level:
    all atoms of a residue must share the same mask value.

    For each grafted residue, the backbone atoms (N/CA/C/O, plus OXT if
    present) are kept from ``atom_array_backbone`` and renamed to the implant
    residue's ``res_name``. The remaining atoms -- the side chain -- are
    taken from the corresponding ``atom_array_implant`` residue, keeping the
    implant coordinates, and are relabeled to the backbone residue's
    ``chain_id``/``res_id``/``ins_code``. Backbone hydrogen atoms
    (H/HA/H2/H3) of grafted residues are not carried over; run
    :func:`add_hydrogens` to rebuild them. Residues not selected in the
    backbone are returned unchanged.

    The ``bonds`` annotation of the returned array is dropped, since a graft
    cannot preserve a consistent bond graph across replaced atoms. Call
    :func:`add_hydrogens` afterwards to rebuild bonds and add hydrogen atoms.

    Parameters
    ----------
    atom_array_backbone
        Atom array whose selected residues receive new side chains.
    atom_array_implant
        Atom array that provides the side chains.
    mask_map
        Flat list of per-atom boolean masks, grouped in (backbone, implant)
        pairs. Each mask must select whole residues only.

    Returns
    -------
    biotite.structure.AtomArray
        A new atom array with the grafted side chains.

    Raises
    ------
    ValueError
        If ``mask_map`` is not a non-empty even-length list, if a mask has
        the wrong length or is not residue-level, if a mask pair selects a
        different number of backbone and implant residues, or if a backbone
        residue is targeted by more than one pair.
    """
    _BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
    _SIDE_CHAIN_EXCLUDE = _BACKBONE_ATOMS | {"H", "HA", "H2", "H3"}

    if len(atom_array_backbone) == 0:
        raise ValueError("atom_array_backbone must not be empty")
    if len(atom_array_implant) == 0:
        raise ValueError("atom_array_implant must not be empty")
    if not isinstance(mask_map, (list, tuple)) or len(mask_map) == 0:
        raise ValueError("mask_map must be a non-empty list of boolean masks")
    if len(mask_map) % 2 != 0:
        raise ValueError(
            "mask_map must contain an even number of masks, "
            "grouped as (backbone_mask, implant_mask) pairs"
        )

    backbone, implant = _ensure_common_annotations(
        atom_array_backbone, atom_array_implant
    )
    backbone_groups = group_atoms_by_residue(backbone)
    implant_groups = group_atoms_by_residue(implant)

    backbone_keep_mask = np.isin(backbone.atom_name, list(_BACKBONE_ATOMS))
    implant_sc_mask = ~np.isin(implant.atom_name, list(_SIDE_CHAIN_EXCLUDE))

    # Map each grafted backbone residue to its source implant residue
    grafts = {}
    for pair_idx in range(len(mask_map) // 2):
        b_mask = np.asarray(mask_map[2 * pair_idx], dtype=bool)
        i_mask = np.asarray(mask_map[2 * pair_idx + 1], dtype=bool)
        if b_mask.shape != (len(backbone),):
            raise ValueError(
                f"mask_map[{2 * pair_idx}] must be a 1D boolean mask of length "
                f"{len(backbone)} (atom_array_backbone)"
            )
        if i_mask.shape != (len(implant),):
            raise ValueError(
                f"mask_map[{2 * pair_idx + 1}] must be a 1D boolean mask of "
                f"length {len(implant)} (atom_array_implant)"
            )

        b_selected = _selected_residues(
            backbone_groups, b_mask, f"mask_map[{2 * pair_idx}]"
        )
        i_selected = _selected_residues(
            implant_groups, i_mask, f"mask_map[{2 * pair_idx + 1}]"
        )
        if len(b_selected) != len(i_selected):
            raise ValueError(
                f"mask pair {pair_idx}: backbone mask selects "
                f"{len(b_selected)} residues but implant mask selects "
                f"{len(i_selected)}; a residue-to-residue mapping requires "
                "equal counts"
            )
        for b_key, i_key in zip(b_selected, i_selected):
            if b_key in grafts:
                raise ValueError(
                    f"backbone residue {b_key} is targeted by more than one "
                    "mask pair"
                )
            grafts[b_key] = i_key

    # Reassemble the backbone array residue by residue
    parts = []
    for b_key, b_idxs in backbone_groups.items():
        if b_key not in grafts:
            parts.append(backbone[b_idxs])
            continue

        i_key = grafts[b_key]
        i_idxs = implant_groups[i_key]
        implant_res_name = implant.res_name[i_idxs[0]]

        # Backbone atoms of the target residue, renamed to the implant type
        bb_part = backbone[b_idxs[backbone_keep_mask[b_idxs]]].copy()
        bb_part.res_name[:] = implant_res_name
        # Side-chain atoms of the implant residue, relabeled to the
        # backbone residue identity
        sc_part = implant[i_idxs[implant_sc_mask[i_idxs]]].copy()
        sc_part.chain_id[:] = b_key[0]
        sc_part.res_id[:] = b_key[1]
        sc_part.ins_code[:] = b_key[2]
        sc_part.res_name[:] = implant_res_name
        sc_part.hetero[:] = backbone.hetero[b_idxs[0]]

        parts.append(bb_part)
        parts.append(sc_part)

    out = bio_struct.concatenate(parts)
    # A grafted bond graph would be incomplete (e.g. missing CA-CB and
    # inter-residue bonds); drop it and let add_hydrogens() rebuild it.
    out.bonds = None
    return out

_MUTATION_SPEC_PATTERN = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


def mutate_without_side_chains(
    atom_array: bio_struct.AtomArray,
    mutation_spec: list[str],
) -> bio_struct.AtomArray:
    """
    Convert residues to backbone-only scaffolds of a target residue type.

    Each entry of ``mutation_spec`` names a residue and its target type in
    the form ``<source letter><res_id><target letter>``, e.g. ``"M1N"``
    (residue 1, currently MET, becomes ASN) or ``"N56E"``. For every named
    residue, all side-chain atoms are removed -- only the backbone atoms
    N/CA/C/O are kept -- and those backbone atoms are renamed to the target
    residue's standard PDB name (upper-case three-letter code, e.g. ASN).
    The result is a backbone-only scaffold (like the output of
    :func:`remove_side_chains`, but retaining the target residue name
    instead of GLY), ready for side-chain grafting via
    :func:`replace_side_chains`.

    Residues are matched by ``res_id`` only (the spec has no chain field),
    so a ``res_id`` shared by several chains selects all of them. The
    source letter is validated against the current ``res_name`` of each
    matched residue; a mismatch raises ``ValueError``, catching spec typos
    and multi-chain numbering collisions.

    Parameters
    ----------
    atom_array
        Input atom array.
    mutation_spec
        Non-empty list of mutation entries, each formatted like ``"M1N"``:
        one-letter code of the current residue, residue id, one-letter code
        of the target residue. Letters are case-insensitive.

    Returns
    -------
    biotite.structure.AtomArray
        A new atom array where the named residues keep only their backbone
        atoms (N/CA/C/O), renamed to the target residue type. All other
        residues are unchanged.

    Raises
    ------
    ValueError
        If ``mutation_spec`` is empty or contains an entry that is not
        formatted like ``"M1N"``, names an unknown residue letter, targets
        the same ``res_id`` twice, references a ``res_id`` that does not
        exist in ``atom_array``, does not match the current ``res_name`` of
        a matched residue, or targets a residue without backbone atoms.
    """
    if not isinstance(mutation_spec, (list, tuple)):
        raise TypeError("mutation_spec must be a list of strings")
    if len(mutation_spec) == 0:
        raise ValueError("mutation_spec must not be empty")

    # res_id -> (source letter, target 3-letter code)
    requested = {}
    for entry in mutation_spec:
        if not isinstance(entry, str):
            raise TypeError(
                f"mutation_spec entries must be strings, got "
                f"{type(entry).__name__}"
            )
        match = _MUTATION_SPEC_PATTERN.fullmatch(entry)
        if match is None:
            raise ValueError(
                f"invalid mutation spec {entry!r}: expected the form 'M1N' "
                "(source residue letter, residue id, target residue letter)"
            )
        src_letter, res_id_str, tgt_letter = match.groups()
        src_letter = src_letter.upper()
        tgt_letter = tgt_letter.upper()
        if src_letter not in AMINO_ACIDS_1TO3_UPPER:
            raise ValueError(
                f"unknown source residue letter {src_letter!r} in {entry!r}"
            )
        if tgt_letter not in AMINO_ACIDS_1TO3_UPPER:
            raise ValueError(
                f"unknown target residue letter {tgt_letter!r} in {entry!r}"
            )
        tgt_res_name = AMINO_ACIDS_1TO3_UPPER[tgt_letter]
        res_id = int(res_id_str)
        if res_id in requested:
            raise ValueError(
                f"residue {res_id} is targeted more than once in mutation_spec"
            )
        requested[res_id] = (src_letter, tgt_res_name)

    groups = group_atoms_by_residue(atom_array)

    # residue key (chain_id, res_id, ins_code) -> target 3-letter code
    matched = {}
    for res_id, (src_letter, tgt_res_name) in requested.items():
        keys = [key for key in groups if key[1] == res_id]
        if not keys:
            raise ValueError(
                f"mutation_spec references residue {res_id}, which does not "
                "exist in atom_array"
            )
        expected_res_name = AMINO_ACIDS_1TO3_UPPER[src_letter]
        for key in keys:
            current_res_name = atom_array.res_name[groups[key][0]]
            if current_res_name != expected_res_name:
                raise ValueError(
                    f"residue {key} is {current_res_name}, but mutation_spec "
                    f"entry expects {src_letter} ({expected_res_name})"
                )
            matched[key] = tgt_res_name

    backbone_mask = np.isin(atom_array.atom_name, ["N", "CA", "C", "O"])
    target_mask = np.zeros(len(atom_array), dtype=bool)
    for key in matched:
        idxs = groups[key]
        target_mask[idxs] = True
        if not np.any(backbone_mask[idxs]):
            raise ValueError(
                f"residue {key} has no backbone atoms (N/CA/C/O); "
                "cannot mutate it"
            )

    keep_mask = (~target_mask) | (target_mask & backbone_mask)
    out = atom_array[keep_mask]

    for key, tgt_res_name in matched.items():
        res_mask = (
            (out.chain_id == key[0])
            & (out.res_id == key[1])
            & (out.ins_code == key[2])
        )
        out.res_name[res_mask] = tgt_res_name

    return out