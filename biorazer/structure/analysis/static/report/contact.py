import numpy as np
from scipy.spatial import KDTree
from biotite import structure as bio_struct
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.utils.selection import normalize_selection
from biorazer.structure.utils.report_utils import (
    _format_atom_label,
    _normalize_interface_residues,
    _to_pymol_atom_selector,
)
from .utils import _normalize_fmt
from .. import select


def report_interface_residues(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    distance_cutoff=3.5,
    fmt="pymol",
    model_name="",
):
    """
    Report interface residues between two selections of atoms in different formats.

    Parameters
    ----------
    fmt : str
        - pymol: print PyMOL commands to visualize the interface residues
        - list: return a list of tuples (chain_id, res_id)
        - text: print interface residues in text format
    """
    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))

    interface_atom_mask_1, interface_atom_mask_2, _ = select.mask_interface_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        distance_cutoff=distance_cutoff,
    )
    interface_atom_mask = interface_atom_mask_1 | interface_atom_mask_2

    interface_residues = _normalize_interface_residues(atom_array, interface_atom_mask)

    if fmt == "list":
        return interface_residues
    elif fmt == "pymol":
        print_with_decoration("Copy the command below to PyMOL", decoration_char="#")
        print(f"select {model_name}_interface, not all")
        for chain_id, res_id in interface_residues:
            print(
                f"select {model_name}_interface, /{model_name}//{chain_id}/{res_id}/ or {model_name}_interface"
            )
        print_decoration_line(decoration_char="#")
    elif fmt == "text":
        for chain_id, res_id in interface_residues:
            print(f"Chain {chain_id}, Residue {res_id}")


def report_interface_contact_matrix(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    distance_cutoff=3.5,
    fmt="matrix",
    model_name="",
    ignore_hydrogen=True,
):
    """
    Report an atom-level interface contact matrix between two selections.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The structure containing both selections.
    selection1, selection2
        Atom selections understood by ``normalize_selection()``.
    distance_cutoff : float, optional
        Atom pairs within this cutoff are treated as interface contacts.
    fmt : str
        - matrix: return (indices_1, indices_2, contact_matrix)
        - list: return [(atom1, atom2, distance), ...]
        - text: print atom-level contacts in text format
        - pymol: print PyMOL distance commands
    model_name : str, optional
        Model name used in PyMOL output.
    ignore_hydrogen : bool, optional
        If True, only heavy atoms are considered.
    """

    fmt = _normalize_fmt(fmt, ("matrix", "list", "text", "pymol"))

    selection_mask_1 = normalize_selection(atom_array, selection1)
    selection_mask_2 = normalize_selection(atom_array, selection2)

    if ignore_hydrogen:
        heavy_mask = atom_array.element != "H"
        selection_mask_1 &= heavy_mask
        selection_mask_2 &= heavy_mask

    indices_1 = np.where(selection_mask_1)[0]
    indices_2 = np.where(selection_mask_2)[0]

    if indices_1.size == 0 or indices_2.size == 0:
        empty = np.zeros((0, 0), dtype=float)
        if fmt == "matrix":
            return np.array([], dtype=int), np.array([], dtype=int), empty
        if fmt == "list":
            return []
        return None

    tree_1 = KDTree(atom_array.coord[indices_1])
    tree_2 = KDTree(atom_array.coord[indices_2])
    neighbor_indices = tree_1.query_ball_tree(tree_2, r=distance_cutoff)

    contact_triplets = []
    used_local_1 = set()
    used_local_2 = set()

    for local_i, neighbors in enumerate(neighbor_indices):
        atom_idx_1 = int(indices_1[local_i])
        for local_j in neighbors:
            atom_idx_2 = int(indices_2[local_j])
            if atom_idx_1 == atom_idx_2:
                continue

            distance = float(
                np.linalg.norm(
                    atom_array.coord[atom_idx_1] - atom_array.coord[atom_idx_2]
                )
            )
            contact_triplets.append(
                (local_i, local_j, atom_idx_1, atom_idx_2, distance)
            )
            used_local_1.add(local_i)
            used_local_2.add(local_j)

    used_local_1 = sorted(used_local_1)
    used_local_2 = sorted(used_local_2)

    compact_indices_1 = indices_1[used_local_1]
    compact_indices_2 = indices_2[used_local_2]

    row_map = {local_i: row_i for row_i, local_i in enumerate(used_local_1)}
    col_map = {local_j: col_j for col_j, local_j in enumerate(used_local_2)}

    contact_matrix = np.full(
        (len(compact_indices_1), len(compact_indices_2)),
        np.nan,
        dtype=float,
    )

    pair_list = []
    for local_i, local_j, atom_idx_1, atom_idx_2, distance in contact_triplets:
        row_i = row_map[local_i]
        col_j = col_map[local_j]
        contact_matrix[row_i, col_j] = distance
        pair_list.append((atom_idx_1, atom_idx_2, distance))

    pair_list.sort(key=lambda item: item[2])

    if fmt == "matrix":
        return compact_indices_1, compact_indices_2, contact_matrix

    if fmt == "list":
        return [(atom_array[i], atom_array[j], dist) for i, j, dist in pair_list]

    if fmt == "text":
        print(f"Total {len(pair_list)} interface atom contacts found:")
        for idx, (i, j, distance) in enumerate(pair_list, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            print(
                f"Contact {idx}: {_format_atom_label(atom1)} -- "
                f"{_format_atom_label(atom2)}, Distance: {distance:.2f} Å"
            )
        return None

    if fmt == "pymol":
        selection_name = f"{model_name}_interface_contacts"
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {selection_name}, not all")
        for idx, (i, j, _) in enumerate(pair_list, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            atom1_selector = _to_pymol_atom_selector(model_name, atom1)
            atom2_selector = _to_pymol_atom_selector(model_name, atom2)
            print(
                f"distance {model_name}_interface_contact_{idx}, "
                f"{atom1_selector}, {atom2_selector}"
            )
            print(
                f"select {selection_name}, "
                f"byres {atom1_selector} or byres {atom2_selector} or {selection_name}"
            )
        print_decoration_line()
        return None


def report_intra_steric_clashes(
    atom_array: bio_struct.AtomArray,
    selection=None,
    cutoff=1.2,
    fmt="pymol",
    model_name="",
    ignore_hydrogen=True,
):
    """
    Report steric clashes inside one atom selection.

    Parameters
    ----------
    selection : np.ndarray, optional
        Atom mask used for intra-selection clash detection. If None, all atoms are used.
    cutoff : float, optional
        Distance below which two atoms are considered clashed.
    fmt : str, optional
        - pymol: print PyMOL commands to visualize clashing pairs
        - text: print clash pairs in text format
        - list: return list of tuples (atom1, atom2, distance)
    """

    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))

    selection_mask = normalize_selection(atom_array, selection)
    if ignore_hydrogen:
        selection_mask &= atom_array.element != "H"

    selected_indices = np.where(selection_mask)[0]
    if selected_indices.size < 2:
        return [] if fmt == "list" else None

    tree = KDTree(atom_array.coord[selected_indices])
    local_pairs = sorted(tree.query_pairs(r=cutoff))

    clashes = []
    for i, j in local_pairs:
        atom_idx_1 = int(selected_indices[i])
        atom_idx_2 = int(selected_indices[j])
        distance = float(
            np.linalg.norm(atom_array.coord[atom_idx_1] - atom_array.coord[atom_idx_2])
        )
        clashes.append((atom_idx_1, atom_idx_2, distance))

    if fmt == "list":
        return [(atom_array[i], atom_array[j], dist) for i, j, dist in clashes]

    if fmt == "text":
        print(f"Total {len(clashes)} intra steric clashes found:")
        for idx, (i, j, distance) in enumerate(clashes, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            print(
                f"Clash {idx}: {_format_atom_label(atom1)} -- {_format_atom_label(atom2)}, Distance: {distance:.2f} Å"
            )
        return None

    if fmt == "pymol":
        selection_name = f"{model_name}_intra_steric_clashes"
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {selection_name}, not all")
        for idx, (i, j, _) in enumerate(clashes, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            atom1_selector = _to_pymol_atom_selector(model_name, atom1)
            atom2_selector = _to_pymol_atom_selector(model_name, atom2)
            print(
                f"distance {model_name}_intra_clash_{idx}, {atom1_selector}, {atom2_selector}"
            )
            print(
                f"select {selection_name}, byres {atom1_selector} or byres {atom2_selector} or {selection_name}"
            )
        print_decoration_line()
        return None


def report_inter_steric_clashes(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    cutoff=1.2,
    fmt="pymol",
    model_name="",
    ignore_hydrogen=True,
):
    """
    Report steric clashes between two atom selections.

    Parameters
    ----------
    selection1 : np.ndarray
        Atom mask for the first selection.
    selection2 : np.ndarray
        Atom mask for the second selection.
    cutoff : float, optional
        Distance below which two atoms are considered clashed.
    fmt : str, optional
        - pymol: print PyMOL commands to visualize clashing pairs
        - text: print clash pairs in text format
        - list: return list of tuples (atom1, atom2, distance)
    """

    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))

    selection_mask_1 = normalize_selection(atom_array, selection1)
    selection_mask_2 = normalize_selection(atom_array, selection2)

    if ignore_hydrogen:
        heavy_mask = atom_array.element != "H"
        selection_mask_1 &= heavy_mask
        selection_mask_2 &= heavy_mask

    indices_1 = np.where(selection_mask_1)[0]
    indices_2 = np.where(selection_mask_2)[0]
    if indices_1.size == 0 or indices_2.size == 0:
        return [] if fmt == "list" else None

    tree_1 = KDTree(atom_array.coord[indices_1])
    tree_2 = KDTree(atom_array.coord[indices_2])
    neighbor_indices = tree_1.query_ball_tree(tree_2, r=cutoff)

    clashes = []
    for local_i, neighbors in enumerate(neighbor_indices):
        atom_idx_1 = int(indices_1[local_i])
        for local_j in neighbors:
            atom_idx_2 = int(indices_2[local_j])
            if atom_idx_1 == atom_idx_2:
                continue
            distance = float(
                np.linalg.norm(
                    atom_array.coord[atom_idx_1] - atom_array.coord[atom_idx_2]
                )
            )
            clashes.append((atom_idx_1, atom_idx_2, distance))

    clashes.sort(key=lambda item: item[2])

    if fmt == "list":
        return [(atom_array[i], atom_array[j], dist) for i, j, dist in clashes]

    if fmt == "text":
        print(f"Total {len(clashes)} inter steric clashes found:")
        for idx, (i, j, distance) in enumerate(clashes, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            print(
                f"Clash {idx}: {_format_atom_label(atom1)} -- {_format_atom_label(atom2)}, Distance: {distance:.2f} Å"
            )
        return None

    if fmt == "pymol":
        selection_name = f"{model_name}_inter_steric_clashes"
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {selection_name}, not all")
        for idx, (i, j, _) in enumerate(clashes, start=1):
            atom1 = atom_array[i]
            atom2 = atom_array[j]
            atom1_selector = _to_pymol_atom_selector(model_name, atom1)
            atom2_selector = _to_pymol_atom_selector(model_name, atom2)
            print(
                f"distance {model_name}_inter_clash_{idx}, {atom1_selector}, {atom2_selector}"
            )
            print(
                f"select {selection_name}, byres {atom1_selector} or byres {atom2_selector} or {selection_name}"
            )
        print_decoration_line()
        return None
