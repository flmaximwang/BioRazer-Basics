import numpy as np
from biotite import structure as bio_struct
import hydride
from biorazer.display import print_with_decoration, print_decoration_line
from . import check, select


def report_interface_residues(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    distance_cutoff=3.5,
    fmt="pymol",
    model_name="",
):
    interface_atom_mask_1, interface_atom_mask_2, _ = select.mask_interface_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        distance_cutoff=distance_cutoff,
    )
    interface_atom_mask = interface_atom_mask_1 | interface_atom_mask_2

    interface_residues = []
    for atom in atom_array[interface_atom_mask]:
        identifier = (atom.chain_id, atom.res_id)
        if not identifier in interface_residues:
            interface_residues.append(identifier)

    if fmt == "list":
        return interface_residues
    elif fmt == "pymol":
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {model_name}_interface, not all")
        for chain_id, res_id in interface_residues:
            print(
                f"select {model_name}_interface, /{model_name}//{chain_id}/{res_id}/ or {model_name}_interface"
            )
        print_decoration_line()

    else:
        raise ValueError(f"Unsupported format {fmt}")


def report_hbonds(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    format="pymol",
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

    formats = ["pymol", "text", "list"]
    assert (
        format in formats
    ), f"Format {format} not supported. Supported formats: {formats}"

    # 需要加氢后才能计算氢键
    model_is_hydrided = check.is_hydrided(atom_array)
    if not model_is_hydrided:
        atom_array, _ = hydride.add_hydrogen(atom_array)
        atom_array.coord = hydride.relax_hydrogen(atom_array)

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

    if format == "pymol":
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select {model_name}_hbonds, not all")  # Initialize the selection
        for i, hbond in enumerate(hbonds, start=1):
            donor, h, acceptor = atom_array[hbond]
            print(
                f"distance {model_name}_hbond_{i}, /{model_name}//{donor.chain_id}/{donor.res_id}/{donor.atom_name}, /{model_name}//{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name}"
            )
            print(
                f"select {model_name}_hbonds, byres /{model_name}//{donor.chain_id}/{donor.res_id}/{donor.atom_name} or byres /{model_name}//{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name} or {model_name}_hbonds"
            )
        print(f"show sticks, byres {model_name}_hbonds")
        print_decoration_line()
    elif format == "text":
        pass
    elif format == "list":
        res = []
        for hbond in hbonds:
            donor, h, acceptor = atom_array[hbond]
            res.append((donor, h, acceptor))
        return res
    else:
        raise ValueError(f"Format {format} not supported. Supported formats: {formats}")


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
    else:
        raise ValueError(f"Unsupported format {fmt}")
