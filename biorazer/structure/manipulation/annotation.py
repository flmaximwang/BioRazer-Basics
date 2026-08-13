import numpy as np
import biotite.structure as bio_struct

from ..selection.index.annotation import group_atoms_by_residue


_WATER_RES_NAMES = {"HOH", "WAT", "DOD"}


def get_renumbered_res_ids(
    atom_array: bio_struct.AtomArray, res_id_start=1
) -> np.ndarray:
    """
    Get renumbered residue IDs from an AtomArray.

    This function returns a list of residue IDs that are renumbered
    based on the unique residues present in the AtomArray.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The AtomArray from which to extract residue IDs.

    Returns
    -------
    np.array
        An array of renumbered residue IDs, starting from `res_id_start`.
        The IDs are unique and ordered based on the appearance of residues
        in the AtomArray.
    """
    old_res_ids = atom_array.get_annotation("res_id")
    res_id_shift = old_res_ids[0] - res_id_start
    new_res_ids = old_res_ids - res_id_shift
    return new_res_ids

def add_bonds_to_organic(
    atom_array: bio_struct.AtomArray,
):
    """
    Auto-detect all small molecules in the AtomArray and add bonds to
    those small molecules based on bonds from the RCSB Chemical
    Component Dictionary (as exposed by ``biotite.structure.info``).

    A residue counts as a small molecule if its atoms are marked as
    hetero (HETATM). If the ``hetero`` annotation is absent (e.g. for
    programmatically built arrays), all residues are treated as
    candidates. Water residues (HOH/WAT/DOD) and monoatomic residues
    (e.g. ions) are skipped, as they cannot carry internal bonds.

    For every candidate residue, the bond definitions of its
    ``res_name`` are looked up in the RCSB Chemical Component Dictionary
    via :func:`biotite.structure.info.bonds_in_residue` -- the same
    bond data that :func:`biotite.structure.info.residue` exposes,
    without materializing atom coordinates. Residues unknown to the
    dictionary contribute no bonds. The dictionary bonds are matched to
    the structure by atom name, so bonds involving atoms not present in
    the structure are omitted. If an atom name occurs multiple times
    within a residue (e.g. altloc duplicates), all combinations are
    connected, following the convention of
    :func:`biotite.structure.connect_via_residue_names`.

    The ``bonds`` annotation of ``atom_array`` is modified in place.
    Existing bonds are preserved; where a dictionary bond already
    exists between two atoms, its bond type is updated to the
    dictionary value, since the Chemical Component Dictionary is the
    authoritative source for small-molecule bond orders.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The atom array whose small-molecule bonds are to be completed.

    Returns
    -------
    biotite.structure.AtomArray
        The input ``atom_array`` with its ``bonds`` annotation updated.

    Examples
    --------
    >>> from biotite.structure import AtomArray
    >>> array = AtomArray(4)
    >>> array.res_name = ["HEM", "HEM", "HEM", "HEM"]
    >>> array.atom_name = ["FE", "NB", "C1A", "C2A"]
    >>> array.hetero = [True, True, True, True]
    >>> array.element = ["FE", "N", "C", "C"]
    >>> array = add_bonds_to_organic(array)
    >>> print(array.bonds.get_bond_count() > 0)
    True
    """
    if len(atom_array) == 0:
        return atom_array

    if atom_array.bonds is None:
        atom_array.bonds = bio_struct.BondList(len(atom_array))

    has_hetero = "hetero" in atom_array.get_annotation_categories()
    # Cache CCD lookups, as the same ligand often occurs multiple times
    bond_dict_cache = {}
    for _, idxs in group_atoms_by_residue(atom_array).items():
        res_name = atom_array.res_name[idxs[0]]
        # Skip water and monoatomic residues (e.g. ions)
        if res_name in _WATER_RES_NAMES or len(idxs) < 2:
            continue
        # Small molecules are HETATM residues, if the flag is available
        if has_hetero and not atom_array.hetero[idxs[0]]:
            continue

        if res_name not in bond_dict_cache:
            bond_dict_cache[res_name] = bio_struct.info.bonds_in_residue(
                res_name
            )
        bond_dict = bond_dict_cache[res_name]
        if not bond_dict:
            # Not in the Chemical Component Dictionary
            continue

        # Map atom names to positions within this residue
        name_positions = {}
        for pos, name in enumerate(atom_array.atom_name[idxs]):
            name_positions.setdefault(name, []).append(pos)

        for (name1, name2), bond_type in bond_dict.items():
            for a in name_positions.get(name1, []):
                for b in name_positions.get(name2, []):
                    if a != b:
                        atom_array.bonds.add_bond(
                            idxs[a], idxs[b], int(bond_type)
                        )
    return atom_array