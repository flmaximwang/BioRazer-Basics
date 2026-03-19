import numpy as np
from biotite import structure as bio_struct
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.utils.report_utils import _to_pymol_atom_selector
from .utils import _normalize_fmt
from .. import check, select


def report_hbonds(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    fmt="pymol",
    model_name="",
    cutoff_dist=2.5,
    cutoff_angle=120,
    donor_elements=("O", "N", "S"),
    acceptor_elements=("O", "N", "S"),
    periodic=False,
):
    """

    Parameters
    ----------
    format: str
        - pymol: print PyMOL commands to visualize the hydrogen bonds
        - list: return a list of tuples (donor, hydrogen, acceptor)
        - text: print hydrogen bonds in text format
    """

    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))

    model_is_hydrided = check.is_hydrided(atom_array)
    if not model_is_hydrided:
        raise Exception(
            "Model is not hydrided. Please add hydrogens before calculating hydrogen bonds."
        )

    hbonds = bio_struct.hbond(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        selection1_type="both",
        cutoff_dist=cutoff_dist,
        cutoff_angle=cutoff_angle,
        donor_elements=donor_elements,
        acceptor_elements=acceptor_elements,
        periodic=periodic,
    )

    if fmt == "pymol":
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {model_name}_hbonds, not all")
        for i, hbond in enumerate(hbonds, start=1):
            donor, h, acceptor = atom_array[hbond]
            donor_selector = _to_pymol_atom_selector(model_name, donor)
            acceptor_selector = _to_pymol_atom_selector(model_name, acceptor)
            print(
                f"distance {model_name}_hbond_{i}, {donor_selector}, {acceptor_selector}"
            )
            print(
                f"select {model_name}_hbonds, byres {donor_selector} or byres {acceptor_selector} or {model_name}_hbonds"
            )
        print(f"show sticks, byres {model_name}_hbonds")
        print_decoration_line()
    elif fmt == "text":
        print(f"Total {len(hbonds)} hydrogen bonds found:")
        for i, hbond in enumerate(hbonds, start=1):
            donor, h, acceptor = atom_array[hbond]
            print(
                f"H-bond {i}: Donor {donor.chain_id} {donor.res_name}{donor.res_id}({donor.atom_name}) -- Acceptor {acceptor.chain_id} {acceptor.res_name}{acceptor.res_id}({acceptor.atom_name}), Distance: {np.linalg.norm(donor.coord - acceptor.coord):.2f} Å"
            )
    elif fmt == "list":
        res = []
        for hbond in hbonds:
            donor, h, acceptor = atom_array[hbond]
            res.append((donor, h, acceptor))
        return res


def report_buried_unsat_hbonds(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    fmt="pymol",
    model_name="",
    **kwargs,
):
    """
    Report unsatisfied hydrogen bonds in the interface between two selections of atoms.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The biotite structure array containing the atoms.
    selection1 : np.ndarray
        The atom mask for the 1st selection.
    selection2 : np.ndarray
        The atom mask for the 2nd selection.
    format : str, optional
        The output format, by default "pymol".
    """

    fmt = _normalize_fmt(fmt, ("pymol",))

    unsat_hbond_mask_1, unsat_hbond_mask_2, _ = select.mask_buried_unsat_hbond_atoms(
        atom_array, selection1=selection1, selection2=selection2, **kwargs
    )

    unsat_hbond_mask = unsat_hbond_mask_1 | unsat_hbond_mask_2

    if fmt == "pymol":
        selection_name = f"{model_name}_unsat_hbonds"
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {selection_name}, not all")
        for atom in atom_array[unsat_hbond_mask]:
            print(
                f"select {selection_name}, /{model_name}//{atom.chain_id}/{atom.res_id}/{atom.atom_name} or {selection_name}"
            )
        print_decoration_line()
