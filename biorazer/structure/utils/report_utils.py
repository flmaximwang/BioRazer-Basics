import numpy as np
from biotite import structure as bio_struct


def _normalize_interface_residues(
    atom_array: bio_struct.AtomArray, interface_atom_mask: np.ndarray
):
    interface_residues = []
    for atom in atom_array[interface_atom_mask]:
        identifier = (str(atom.chain_id), int(atom.res_id))
        if identifier not in interface_residues:
            interface_residues.append(identifier)
    return interface_residues


def _format_atom_label(atom):
    return f"{atom.chain_id} {atom.res_name}{atom.res_id}({atom.atom_name})"


def _to_pymol_atom_selector(model_name: str, atom):
    return f"/{model_name}//{atom.chain_id}/{atom.res_id}/{atom.atom_name}"
