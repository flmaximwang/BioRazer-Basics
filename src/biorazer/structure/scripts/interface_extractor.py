import json, os
import numpy as np
from biotite.structure.io.pdb import PDBFile
import biotite.structure.io.pdb as pdb
import pyrosetta
pyrosetta.init("-mute all")

def extract_interface_residues(pdb_file, overwrite=False):
    interface_set_json_root, pdb_ext = os.path.splitext(pdb_file)
    interface_set_json = interface_set_json_root + "_interface_annotations.json"
    if os.path.exists(interface_set_json) and not overwrite:
        return
    
    pose = pyrosetta.pose_from_file(pdb_file)
    mover = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    mover.apply(pose)
    
    res = {}
    res_index_set = mover.get_interface_set()
    res["res_indices"] = list(res_index_set)
    res["res_indices"].sort()
    res["chain_id, res_id"] = []
    model_name, pdb_ext = os.path.splitext(os.path.basename(pdb_file))
    res["pymol_selector"] = f"select {model_name}_interface_residues, "
    selector_dict = {} 
    for res_index in res["res_indices"]:
        chain_id = pose.pdb_info().chain(res_index)
        res_id = pose.pdb_info().number(res_index)
        res["chain_id, res_id"].append([chain_id, res_id])
        if not chain_id in selector_dict:
            selector_dict[chain_id] = [res_id]
        else:
            selector_dict[chain_id].append(res_id)
    select_macro_list = []
    for chain_id, res_ids in selector_dict.items():
        select_macro_list.append(f"/{model_name}//{chain_id}/{'+'.join([str(res_id) for res_id in res_ids])}")
    res["pymol_selector"] += " or ".join(select_macro_list)
    
    with open(interface_set_json, "w") as f:
        json.dump(res, f, indent=4)

def batch_extract_interface_residues(struc_dir, overwrite=False):

    for dirpath, dirnames, filenames in os.walk(struc_dir):
        for filename in filenames:
            if filename.endswith(".pdb"):
                pdb_file = os.path.join(dirpath, filename)
                extract_interface_residues(pdb_file, overwrite=overwrite)

def append_interface_residue_sidechain_atom_indices(pdb_file, interface_annotation_file):
    my_struc = pdb.get_structure(PDBFile.read(pdb_file))[0]
    interface_annotations = json.load(open(interface_annotation_file))
    interface_residues = interface_annotations["chain_id, res_id"]
    final_selector = np.zeros(my_struc.array_length(), dtype=bool)
    for chain_id, res_id in interface_residues:
        chain_id_selector = my_struc.chain_id == chain_id
        res_id_selector = my_struc.res_id == res_id
        atom_name_selector = np.ones(my_struc.array_length(), dtype=bool)
        for atom_name in ["C", "O", "N", "CA"]:
            atom_name_selector = np.logical_and(atom_name_selector, my_struc.atom_name != atom_name)
        my_selector = np.ones(my_struc.array_length(), dtype=bool)
        for selector in [chain_id_selector, res_id_selector, atom_name_selector]:
            my_selector = np.logical_and(my_selector, selector)
        final_selector = np.logical_or(final_selector, my_selector)
    sidechain_indices = np.where(final_selector)[0]
    interface_annotations["sidechain_indices"] = sidechain_indices.tolist()
    json.dump(interface_annotations, open(interface_annotation_file, "w"), indent=4)
    
def batch_append_interface_residue_sidechain_atom_indices(struc_dir):
    for dirpath, dirnames, filenames in os.walk(struc_dir):
        for filename in filenames:
            if filename.endswith(".pdb"):
                pdb_file = os.path.join(dirpath, filename)
                interface_annotation_file = os.path.splitext(pdb_file)[0] + "_interface_annotations.json"
                append_interface_residue_sidechain_atom_indices(pdb_file, interface_annotation_file)