import numpy as np
import biotite.structure as struc

def check_struc_type(my_struc: struc.AtomArray):
    if not isinstance(my_struc, struc.AtomArray):
        raise TypeError(f"my_struc must be an AtomArray. The given object is {type(my_struc)}")

def get_chain_atom_indices(my_struc: struc.AtomArray, chains: list):
    check_struc_type(my_struc)
    chain_ids = my_struc.get_annotation("chain_id")
    atom_indices = []
    for chain in chains:
        atom_indices.extend(
            np.where(
                (chain_ids == chain)
            )[0]
        )
    return [int(i) for i in atom_indices]

def select_chain_atoms(my_struc: struc.AtomArray, chains: list):
    return my_struc[get_chain_atom_indices(my_struc, chains)]

def get_residue_atom_indices(my_struc: struc.AtomArray, residues: list[tuple[str, int]]):
    check_struc_type(my_struc)
    chain_ids = my_struc.get_annotation("chain_id")
    res_ids = my_struc.get_annotation("res_id")
    atom_indices = []
    for chain_id, res_id in residues:
        atom_indices.extend(
            np.where(
                (chain_ids == chain_id) & (res_ids == res_id)
            )[0]
        )
    return [int(i) for i in atom_indices]

def select_residue_atoms(my_struc: struc.AtomArray, residues: list[tuple[str, int]]):
    check_struc_type(my_struc)
    return my_struc[get_residue_atom_indices(my_struc, residues)]

def get_sidechain_atom_indices(my_struc: struc.AtomArray):
    check_struc_type(my_struc)
    atom_names = my_struc.get_annotation("atom_name")
    atom_indices = np.where(
                    ~np.isin(atom_names, ["N", "CA", "C", "O"])
                )[0]
    return [int(i) for i in atom_indices]

def select_sidechain_atoms(my_struc: struc.AtomArray):
    return my_struc[get_sidechain_atom_indices(my_struc)]