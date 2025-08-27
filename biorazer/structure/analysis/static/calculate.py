import numpy as np
import biotite.structure as bio_struct
from . import report, select


def sasa_no_h(
    atom_array: bio_struct.AtomArray,
    probe_radius=1.4,
    atom_filter=None,
    ignore_ions=True,
    point_number=1000,
    point_distr="Fibonacci",
    vdw_radii="Single",
):
    heavy_mask = atom_array.element != "H"
    atom_sasa_tmp = bio_struct.sasa(
        atom_array[heavy_mask],
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
    )
    atom_sasa = np.zeros(atom_array.shape, dtype=float)
    atom_sasa[heavy_mask] = atom_sasa_tmp
    return atom_sasa


def buried_unsat_hbond(
    atom_array: bio_struct.AtomArray,
    selection1: np.ndarray,
    selection2: np.ndarray,
    sasa_kwargs=dict(sasa_cutoff=0.5, probe_radius=1.4),
    interface_kwargs=dict(distance_cutoff=3.5),
    hbond_kwargs=dict(),
):
    unsat_hbond_mask_1, unsat_hbond_mask_2, _ = select.mask_buried_unsat_hbond_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        sasa_kwargs=sasa_kwargs,
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
    cutoff_dist=2.5,
    cutoff_angle=120,
    donor_elements=("O", "N", "S"),
    acceptor_elements=("O", "N", "S"),
    periodic=False,
):
    hbond_list = report.report_hbonds(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        format="list",
        cutoff_dist=cutoff_dist,
        cutoff_angle=cutoff_angle,
        donor_elements=donor_elements,
        acceptor_elements=acceptor_elements,
        periodic=periodic,
    )

    return len(hbond_list)
