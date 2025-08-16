import numpy as np
import biotite.structure as bio_struct
from . import selectors, reporters


def sasa_no_h(atom_array: bio_struct.AtomArray, **kwargs):
    heavy_mask = atom_array.element != "H"
    atom_sasa_tmp = bio_struct.sasa(
        atom_array[heavy_mask], vdw_radii="Single", **kwargs
    )
    atom_sasa = np.zeros(atom_array.shape, dtype=float)
    atom_sasa[heavy_mask] = atom_sasa_tmp
    return atom_sasa


def unsat_hbond(
    atom_array: bio_struct.AtomArray,
    selection1: np.ndarray,
    selection2: np.ndarray,
    sasa_cutoff: float = 0.5,
    interface_kwargs: dict = dict(distance_cutoff=3.5),
    hbond_kwargs: dict = dict(),
):
    unsat_hbond_mask_1, unsat_hbond_mask_2, _ = selectors.mask_unsat_hbond_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        sasa_cutoff=sasa_cutoff,
        interface_kwargs=interface_kwargs,
        hbond_kwargs=hbond_kwargs,
    )
    unsat_hbonds_1 = int(np.sum(unsat_hbond_mask_1))
    unsat_hbonds_2 = int(np.sum(unsat_hbond_mask_2))
    return [unsat_hbonds_1 + unsat_hbonds_2, unsat_hbonds_1, unsat_hbonds_2]


def hbond(
    atom_array: bio_struct.AtomArray,
    selection1: np.ndarray,
    selection2: np.ndarray,
    **kwargs
):
    hbond_list = reporters.report_hbonds(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        format="list",
        **kwargs
    )

    return len(hbond_list)
