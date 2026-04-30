from __future__ import annotations

from collections import defaultdict

import numpy as np
import biotite.structure as bio_struct
import biotite.sequence as bio_seq
from scipy.spatial import KDTree

from biorazer.display import print_decoration_line, print_with_decoration
from biorazer.util.dictionary.aa_types import TYPES2AA
from ..calculate.array import sasa_array
from .util import _normalize_fmt


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}


def _iter_residue_masks(structure: bio_struct.AtomArray):
    if len(structure) == 0:
        return

    start = 0
    current_key = (
        str(structure.chain_id[0]),
        int(structure.res_id[0]),
        str(structure.res_name[0]),
    )

    for index in range(1, len(structure)):
        key = (
            str(structure.chain_id[index]),
            int(structure.res_id[index]),
            str(structure.res_name[index]),
        )
        if key != current_key:
            mask = np.zeros(len(structure), dtype=bool)
            mask[start:index] = True
            yield current_key, mask
            start = index
            current_key = key

    mask = np.zeros(len(structure), dtype=bool)
    mask[start:] = True
    yield current_key, mask


def _residue_one_letter_code(res_name: str) -> str | None:
    if len(res_name) != 3:
        return None
    try:
        return bio_seq.ProteinSequence.convert_letter_3to1(res_name)
    except KeyError:
        return None


def _is_target_sidechain_heavy_atom(atom, exclude_elements: set[str]) -> bool:
    return (
        str(atom.element) not in exclude_elements
        and str(atom.element) != "H"
        and str(atom.atom_name) not in BACKBONE_ATOM_NAMES
    )


def _collect_surface_patch_atoms(
    structure: bio_struct.AtomArray,
    target_aas: str,
    surface_atom_threshold: float,
    probe_radius: float,
    atom_filter: str | None,
    ignore_ions: bool,
    point_number: int,
    point_distr: str,
    vdw_radii: str | dict[str, float] | None,
    exclude_elements: list[str],
) -> list[dict[str, object]]:
    total_sasa_array = sasa_array(
        structure,
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
        exclude_elements=exclude_elements,
    )

    atoms = []
    target_set = set(target_aas)
    exclude_set = set(exclude_elements)
    for (chain_id, res_id, res_name), residue_mask in _iter_residue_masks(structure):
        residue_code = _residue_one_letter_code(res_name)
        if residue_code is None or residue_code not in target_set:
            continue

        residue_atoms = structure[residue_mask]
        residue_indices = np.where(residue_mask)[0]
        for local_index, atom in enumerate(residue_atoms):
            if not _is_target_sidechain_heavy_atom(atom, exclude_set):
                continue

            atom_index = int(residue_indices[local_index])
            isolated_atom_sasa = float(
                sasa_array(
                    structure[atom_index : atom_index + 1],
                    probe_radius=probe_radius,
                    atom_filter=atom_filter,
                    ignore_ions=ignore_ions,
                    point_number=point_number,
                    point_distr=point_distr,
                    vdw_radii=vdw_radii,
                    exclude_elements=exclude_elements,
                )[0]
            )
            if isolated_atom_sasa <= 0:
                continue

            atom_sasa = float(total_sasa_array[atom_index])
            atom_surface_ratio = atom_sasa / isolated_atom_sasa
            if atom_surface_ratio < surface_atom_threshold:
                continue

            atoms.append(
                {
                    "chain_id": chain_id,
                    "res_id": res_id,
                    "res_name": res_name,
                    "aa": residue_code,
                    "atom_name": str(atom.atom_name),
                    "element": str(atom.element),
                    "atom_index": atom_index,
                    "surface_ratio": atom_surface_ratio,
                    "coord": np.asarray(atom.coord, dtype=float),
                }
            )

    return atoms


def _cluster_patch_atoms(
    atoms: list[dict[str, object]],
    cluster_distance_threshold: float,
    min_patch_size: int,
) -> list[list[dict[str, object]]]:
    if not atoms:
        return []

    coords = np.asarray([atom["coord"] for atom in atoms], dtype=float)
    if len(coords) == 1:
        return [[atoms[0]]] if min_patch_size <= 1 else []

    tree = KDTree(coords)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for index in range(len(atoms)):
        adjacency[index]
    for i, j in tree.query_pairs(r=cluster_distance_threshold):
        adjacency[i].add(j)
        adjacency[j].add(i)

    clusters = []
    visited = set()
    for start in range(len(atoms)):
        if start in visited:
            continue
        stack = [start]
        component = []
        visited.add(start)
        while stack:
            current = stack.pop()
            component.append(atoms[current])
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        clusters.append(
            sorted(
                component,
                key=lambda item: (
                    item["chain_id"],
                    item["res_id"],
                    item["atom_name"],
                    item["atom_index"],
                ),
            )
        )

    clusters.sort(
        key=lambda cluster: (
            -len(cluster),
            -float(np.mean([item["surface_ratio"] for item in cluster])),
            str(cluster[0]["chain_id"]),
            int(cluster[0]["res_id"]),
            str(cluster[0]["atom_name"]),
        )
    )
    return [cluster for cluster in clusters if len(cluster) >= min_patch_size]


def _unique_residue_members(
    cluster: list[dict[str, object]],
) -> list[tuple[str, int, str]]:
    residue_members = []
    seen = set()
    for item in cluster:
        residue_key = (
            str(item["chain_id"]),
            int(item["res_id"]),
            str(item["res_name"]),
        )
        if residue_key not in seen:
            seen.add(residue_key)
            residue_members.append(residue_key)
    return residue_members


def _build_patch_summary(
    patch_clusters: list[list[dict[str, object]]],
) -> list[dict[str, object]]:
    summaries = []
    for patch_id, cluster in enumerate(patch_clusters, start=1):
        residue_members = _unique_residue_members(cluster)
        summaries.append(
            {
                "patch_id": patch_id,
                "atom_count": len(cluster),
                "residue_count": len(residue_members),
                "mean_surface_ratio": float(
                    np.mean([item["surface_ratio"] for item in cluster])
                ),
                "atom_members": [
                    (
                        str(item["chain_id"]),
                        int(item["res_id"]),
                        str(item["res_name"]),
                        str(item["atom_name"]),
                    )
                    for item in cluster
                ],
                "residue_members": residue_members,
                "members": residue_members,
            }
        )
    return summaries


def _to_pymol_atom_selector(
    model_name: str,
    chain_id: str,
    res_id: int,
    atom_name: str,
) -> str:
    if model_name:
        return f"/{model_name}//{chain_id}/{res_id}/{atom_name}"
    return f"chain {chain_id} and resi {res_id} and name {atom_name}"


def _report_patch_output(
    patch_name: str,
    patch_summary: list[dict[str, object]],
    fmt: str,
    model_name: str,
):
    if fmt == "list":
        return patch_summary

    if fmt == "text":
        print(f"Total {len(patch_summary)} {patch_name} patches found:")
        for patch in patch_summary:
            residue_members = ", ".join(
                f"{chain_id}:{res_name}{res_id}"
                for chain_id, res_id, res_name in patch["residue_members"]
            )
            print(
                f"Patch {patch['patch_id']}: atom_count={patch['atom_count']}, "
                f"residue_count={patch['residue_count']}, "
                f"mean_surface_ratio={patch['mean_surface_ratio']:.3f}, "
                f"residues=[{residue_members}]"
            )
        return None

    selection_prefix = f"{model_name}_" if model_name else ""
    print_with_decoration("Copy the command below to PyMOL")
    print(f"select {selection_prefix}{patch_name}_patches, not all")
    for patch in patch_summary:
        patch_atom_selection_name = (
            f"{selection_prefix}{patch_name}_patch_{patch['patch_id']}_atoms"
        )
        patch_residue_selection_name = (
            f"{selection_prefix}{patch_name}_patch_{patch['patch_id']}"
        )
        print(f"select {patch_atom_selection_name}, not all")
        for chain_id, res_id, _, atom_name in patch["atom_members"]:
            selector = _to_pymol_atom_selector(
                model_name,
                chain_id,
                res_id,
                atom_name,
            )
            print(
                f"select {patch_atom_selection_name}, {selector} or {patch_atom_selection_name}"
            )
        print(
            f"select {patch_residue_selection_name}, byres {patch_atom_selection_name}"
        )
        print(
            f"select {selection_prefix}{patch_name}_patches, {patch_residue_selection_name} or {selection_prefix}{patch_name}_patches"
        )
    print_decoration_line()
    return None


def report_surface_patch(
    structure: bio_struct.AtomArray,
    target_aas: str,
    patch_name: str,
    surface_atom_threshold: float = 0.25,
    cluster_distance_threshold: float = 5.0,
    min_patch_size: int = 5,
    fmt: str = "pymol",
    model_name: str = "",
    probe_radius: float = 1.4,
    atom_filter: str | None = None,
    ignore_ions: bool = True,
    point_number: int = 1000,
    point_distr: str = "Fibonacci",
    vdw_radii: str | dict[str, float] | None = "Single",
    exclude_elements: list[str] = [],
):
    """
    Report surface atom patches for a given amino-acid category.

    Parameters
    ----------
    structure : bio_struct.AtomArray
        Input structure used for SASA calculation and patch clustering.
    target_aas : str
        One-letter amino-acid codes whose sidechain heavy atoms are included in the
        patch analysis, such as hydrophobic or charged residues.
    patch_name : str
        Name used in text output and PyMOL selection names.
    surface_atom_threshold : float, default=0.25
        Minimum SASA ratio required for a sidechain heavy atom to be treated as
        surface-exposed. The ratio is computed as atom SASA in the full structure
        divided by the SASA of that isolated atom.
    cluster_distance_threshold : float, default=5.0
        Maximum sidechain heavy atom distance in Angstrom used to connect two
        candidate surface atoms into the same spatial patch.
    min_patch_size : int, default=5
        Minimum number of surface atoms required for a connected component to be
        kept as a reported patch.
    fmt : str, default="pymol"
        Output format. Supported values are "pymol", "text", and "list".
    model_name : str, default=""
        Optional model name used when composing PyMOL selectors.
    probe_radius : float, default=1.4
        Probe radius passed to the SASA calculation.
    atom_filter : str or None, optional
        Optional atom filter forwarded to the SASA calculation.
    ignore_ions : bool, default=True
        Whether ions are ignored during SASA calculation.
    point_number : int, default=1000
        Number of surface points used in the SASA estimation.
    point_distr : str, default="Fibonacci"
        Surface point distribution mode used by the SASA calculation.
    vdw_radii : dict[str, float] or None, optional
        Optional van der Waals radii overrides passed to the SASA calculation.
    exclude_elements : list[str], default=[]
        Element symbols excluded before SASA calculation, such as hydrogen.

    Returns
    -------
    list[dict[str, object]] or None
        Returns a patch summary list when ``fmt`` is "list". Otherwise prints the
        requested report and returns None.
    """
    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))
    atoms = _collect_surface_patch_atoms(
        structure=structure,
        target_aas=target_aas,
        surface_atom_threshold=surface_atom_threshold,
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
        exclude_elements=exclude_elements,
    )
    if min_patch_size < 1:
        raise ValueError(f"min_patch_size must be >= 1, got {min_patch_size}")

    patch_clusters = _cluster_patch_atoms(
        atoms,
        cluster_distance_threshold,
        min_patch_size,
    )
    patch_summary = _build_patch_summary(patch_clusters)
    return _report_patch_output(patch_name, patch_summary, fmt, model_name)


def report_hydropohbic_patch(
    structure: bio_struct.AtomArray,
    surface_atom_threshold: float = 0.25,
    cluster_distance_threshold: float = 5.0,
    min_patch_size: int = 5,
    fmt: str = "pymol",
    hydrophobic_aas=TYPES2AA["hydrophobic"],
    probe_radius: float = 1.4,
    atom_filter: str | None = None,
    ignore_ions: bool = True,
    point_number: int = 1000,
    point_distr: str = "Fibonacci",
    vdw_radii: str | dict[str, float] | None = "Single",
    exclude_elements=["H"],
    model_name: str = "",
):
    return report_surface_patch(
        structure=structure,
        target_aas=hydrophobic_aas,
        patch_name="hydrophobic",
        surface_atom_threshold=surface_atom_threshold,
        cluster_distance_threshold=cluster_distance_threshold,
        min_patch_size=min_patch_size,
        fmt=fmt,
        model_name=model_name,
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
        exclude_elements=exclude_elements,
    )


def report_charged_patch(
    structure: bio_struct.AtomArray,
    surface_atom_threshold: float = 0.25,
    cluster_distance_threshold: float = 5.0,
    min_patch_size: int = 5,
    fmt: str = "pymol",
    charged_aas=TYPES2AA["charged"],
    probe_radius: float = 1.4,
    atom_filter: str | None = None,
    ignore_ions: bool = True,
    point_number: int = 1000,
    point_distr: str = "Fibonacci",
    vdw_radii: str | dict[str, float] | None = "Single",
    exclude_elements=["H"],
    model_name: str = "",
):
    return report_surface_patch(
        structure=structure,
        target_aas=charged_aas,
        patch_name="charged",
        surface_atom_threshold=surface_atom_threshold,
        cluster_distance_threshold=cluster_distance_threshold,
        min_patch_size=min_patch_size,
        fmt=fmt,
        model_name=model_name,
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
        exclude_elements=exclude_elements,
    )
