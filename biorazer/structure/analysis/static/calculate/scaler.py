import numpy as np
import biotite.structure as bio_struct
import biotite.sequence as bio_seq
from .. import report, select
from .array import sasa_array, buried_unsat_hbond, hbond
from biorazer.util.dictionary.aa_types import TYPES2AA

"""
Multiple calculation can be performed on a same array to return multiple scalar values.
"""


def _normalize_sum_aas(sum_aas: str, atom_array: bio_struct.AtomArray):
    """
    Return a list of three-letter amino acid codes corresponding to the specified sum_aas.
    """
    if isinstance(sum_aas, str):
        try:
            res_str = TYPES2AA[sum_aas]
        except KeyError:
            for aa in sum_aas:
                if aa not in bio_seq.ProteinSequence.alphabet:
                    raise ValueError(f"Invalid amino acid: {aa}")
            res_str = sum_aas
    else:
        raise ValueError("sum_aas must be a string, list or None")

    return [bio_seq.ProteinSequence.convert_letter_1to3(aa) for aa in res_str]


def sasa_value(
    atom_array: bio_struct.AtomArray,
    probe_radius=1.4,
    atom_filter=None,
    ignore_ions=True,
    point_number=1000,
    point_distr="Fibonacci",
    vdw_radii="Single",
    exclude_elements: list[str] = ["H"],
    sum_aa_list: list[str] = ["all"],
):
    """
    Calculate the total solvent accessible surface area (SASA) of a structure.

    - sum_aas: If specified, only include atoms of residues with these amino acids in the sum.
        Supported amino acid types are "aliphatic", "aromatic", "polar", "charged", "hydrophobic" and "all".
        1-letter string of amino acids is also supported, e.g. "AFILMPVWY" for hydrophobic amino acids.
        3-letter string of amino acids is not supported.
    """
    my_sasa_array = sasa_array(
        atom_array,
        probe_radius=probe_radius,
        atom_filter=atom_filter,
        ignore_ions=ignore_ions,
        point_number=point_number,
        point_distr=point_distr,
        vdw_radii=vdw_radii,
        exclude_elements=exclude_elements,
    )
    result = []
    for sum_aa in sum_aa_list:
        sum_aa_normalized = _normalize_sum_aas(sum_aa, atom_array)
        my_sasa_array_sum = my_sasa_array[
            np.isin(atom_array.res_name, sum_aa_normalized)
        ]
        result.append(float(np.sum(my_sasa_array_sum)))
    return result
