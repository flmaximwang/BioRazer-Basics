"""
This module provides functions to convert core data structures across different packages in BioRazer.
"""

import biotite.sequence as bt_seq
import biotite.structure as bt_struct


def _chain_structure_to_sequence(chain_structure: bt_struct.AtomArray):
    """
    Currently this function only supports conversion to ProteinSequence.
    For NucleotideSequence support, additional logic is needed to identify nucleotide residues.
    """

    _, residue_names = bt_struct.get_residues(chain_structure)
    residue_names = [
        bt_seq.ProteinSequence.convert_letter_3to1(residue_name)
        for residue_name in residue_names
    ]
    if all(name in bt_seq.ProteinSequence.alphabet for name in residue_names):
        seq_str = "".join(residue_names)
        return bt_seq.ProteinSequence(seq_str)
    elif all(name in bt_seq.NucleotideSequence.alphabet for name in residue_names):
        seq_str = "".join(residue_names)
        return bt_seq.NucleotideSequence(seq_str)
    else:
        raise ValueError("Chain structure contains mixed or unknown residue types")


def structure2sequence(structure: bt_struct.AtomArray):
    """
    Convert a biotite.structure.AtomArray to a biotite.sequence.ProteinSequence
    or biotite.sequence.NucleotideSequence depending on the content of the structure.
    """

    if not isinstance(structure, bt_struct.AtomArray):
        raise TypeError("Input must be a biotite.structure.AtomArray")

    result = {}
    for chain_id in bt_struct.get_chains(structure):
        chain_structure = structure[structure.chain_id == chain_id]
        seq = _chain_structure_to_sequence(chain_structure)
        result[chain_id] = seq
    return result
