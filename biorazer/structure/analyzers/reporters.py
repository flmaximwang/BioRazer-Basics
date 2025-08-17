import numpy as np
import biotite.structure as bio_struct
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.analyzers import checkers, selectors
import hydride


def report_interface_atoms(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    distance_cutoff=3.5,
):
    pass


def report_hbonds(
    atom_array: bio_struct.AtomArray, selection1, selection2, format="pymol", **kwargs
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
    model_is_hydrided = checkers.is_hydrided(atom_array)
    if not model_is_hydrided:
        atom_array, _ = hydride.add_hydrogen(atom_array)
        atom_array.coord = hydride.relax_hydrogen(atom_array)

    hbonds = bio_struct.hbond(
        atom_array,
        selection1=selection1,
        selection2=selection2,
    )

    if format == "pymol":
        print_with_decoration("Copy the command below to PyMOL")
        print(f"select hbonds, not all")  # Initialize the selection
        for i, hbond in enumerate(hbonds, start=1):
            donor, h, acceptor = atom_array[hbond]
            if not model_is_hydrided:
                print(
                    f"distance hbond_{i}, ///{donor.chain_id}/{donor.res_id}/{donor.atom_name}, ///{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name}"
                )
            else:
                print(
                    f"distance hbond_{i}, ///{h.chain_id}/{h.res_id}/{h.atom_name}, ///{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name}"
                )
            print(
                f"select hbonds, byres ///{donor.chain_id}/{donor.res_id}/{donor.atom_name} or byres ///{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name} or hbonds"
            )
        print(f"show sticks, byres hbonds")
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


def report_unsat_hbonds(
    atom_array: bio_struct.AtomArray,
    selection1,
    selection2,
    format="pymol",
    interface_kwargs=dict(distance_cutoff=3.5),
    hbond_kwargs=dict(),
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

    if not checkers.is_hydrided(atom_array):
        atom_array, _ = hydride.add_hydrogen(atom_array)
        atom_array.coord = hydride.relax_hydrogen(atom_array)

    interface_mask_1, interface_mask_2, _ = selectors.mask_interface_atoms(
        atom_array, selection1, selection2, **interface_kwargs
    )
    interface_mask = interface_mask_1 | interface_mask_2
    interface_res_mask = selectors.by_res_id(atom_array, interface_mask)

    hbond_mask_1, hbond_mask_2, _ = selectors.mask_hbond_atoms(
        atom_array, interface_res_mask, np.ones(atom_array.shape, bool), **hbond_kwargs
    )

    NO_mask = (atom_array.element == "O") | (atom_array.element == "N")
    interface_NO_mask = interface_mask & NO_mask
    unsat_hbond_mask = interface_NO_mask & ~hbond_mask_1
