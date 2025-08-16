# The af_utils module provides lots of useful functions to help you handle alphafold structures.

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from .pre_processors import *
from .rosetta_utils import *
from . import rosetta_utils
from biotite.structure.io import pdbx

def generate_af_file_prefix(top_dir, identifier, design):
    prediction_dir = generate_prediction_dir(top_dir, identifier)
    design = format_design_name(design)
    af_file_prefix = os.path.join(prediction_dir, design, f"fold_{design}")
    return af_file_prefix

def get_interface_residue_indices(top_dir, identifier, design, prediction_i):
    '''
    Warning: this function doesn't return atom indices, but residue indices
    This function relies on the rosetta_utils.get_interface_residues function.
    '''
    residue_identifiers = get_interface_residues(top_dir, identifier, design, prediction_i)
    meta_data = MetaFile(top_dir, identifier, design, prediction_i).get_data()
    chain_ids = meta_data["token_chain_ids"]
    res_ids = meta_data["token_res_ids"]
    residue_indices = {}
    for i, (chain_id, res_id) in enumerate(zip(chain_ids, res_ids)):
        if (chain_id, res_id) in residue_identifiers:
            if not chain_id in residue_indices:
                residue_indices[chain_id] = []
            residue_indices[chain_id].append(i)
    return residue_indices

def get_ranking_score(top_dir, identifier, design, prediction_i):
    return ConfFile(top_dir, identifier, design, prediction_i)["ranking_score"]

# def assign_pLDDT(struc: struc.AtomArray, AF_cif_file_path):
#     struc_file = pdbx.CIFFile.read(AF_cif_file_path)
#     pLDDT = struc_file.block["atom_site"]['B_iso_or_equiv'].as_array()
#     # 将 pLDDT 从 U5 转换为 float64
#     pLDDT = pLDDT.astype(np.float64)
#     struc.del_annotation("pLDDT")
#     struc.add_annotation("pLDDT", np.float64)
#     struc.set_annotation("pLDDT", pLDDT)

# def calculate_struc_average_pLDDT(struc: struc.AtomArray):
#     if not "pLDDT" in struc.get_annotation_categories():
#         raise KeyError("pLDDT not found in annotations, please assign pLDDT first")
#     return np.mean(struc.get_annotation("pLDDT"))

def read_pLDDT_from_meta_file(meta_file_path):
    with open(meta_file_path, "r") as f:
        meta = json.load(f)
    return np.array(meta["atom_plddts"])

def get_pLDDT(meta_file_path):
    with open(meta_file_path, "r") as f:
        meta = json.load(f)
    return np.array(meta["atom_plddts"])

# def get_pLDDT(design, prediction_i):
#     meta_data = read_meta_file(design, prediction_i)
#     pLDDT = np.array(meta_data["atom_plddts"])
#     return pLDDT

# def get_binder_tdomain_plddt(design, prediction_i):
#     cif_file_path = os.path.join("Predictions", design, f"fold_{design}_model_{prediction_i}.cif")
#     atoms_plddt = get_pLDDT(design, prediction_i)
#     prdct_struc = pdbx.get_structure(pdbx.CIFFile.read(cif_file_path))[0]
#     binder_tdomain_indices
#     complex_selector = np.logical_or(prdct_struc.get_annotation("domain") == "Bdr", prdct_struc.get_annotation("domain") == DOMAIN_CONVERTER[design.split("_")[2]])
#     pLDDT_vector[i] = np.mean(atoms_plddt[complex_selector])

def get_i_pTM_for_design(top_dir, identifier, design, prediction_i, chain_1_index, chain_2_index):
    conf_data = read_conf_file(top_dir, identifier, design, prediction_i)
    return conf_data['chain_pair_iptm'][chain_1_index, chain_2_index]

def get_full_pae(top_dir, identifier, design, prediction_i):
    meta_data = read_meta_file(top_dir, identifier, design, prediction_i)
    return np.array(meta_data["pae"], dtype=np.float64)

def get_i_pAE(meta_file_path):
    with open(meta_file_path, 'r') as meta_file:
        metadata = json.load(meta_file)
    identifiers = {}
    counter = 0
    chain_ids = list(set(metadata["token_chain_ids"]))
    chain_ids.sort()
    for chain_id in chain_ids:
        chain_length = metadata["token_chain_ids"].count(chain_id)
        identifiers[chain_id] = slice(counter, counter + chain_length)
        counter += chain_length
    pae = np.array(metadata["pae"])
    norm_pae = (pae - np.min(pae)) / (np.max(pae) - np.min(pae))
    i_pAE = (np.mean(norm_pae[identifiers["A"], identifiers["B"]]) + np.mean(norm_pae[identifiers["B"], identifiers["A"]])) / 2
    return i_pAE

def normalize_pae(pae: np.array):
    return (pae - np.min(pae)) / (np.max(pae) - np.min(pae))

def plot_pae(pae: np.array):
    fig, ax = plt.subplots(figsize=(20, 20))
    plt.imshow(pae, cmap='coolwarm', interpolation='nearest', vmin=0, vmax=31.75)
    plt.colorbar()
    return fig, ax

def plot_normalized_pae(pae: np.array):
    norm_pae = normalize_pae(pae)
    fig, ax = plt.subplots(figsize=(20, 20))
    plt.imshow(norm_pae, cmap='coolwarm', interpolation='nearest', vmin=0, vmax=1)
    plt.colorbar()
    return fig, ax

def plot_pae_in_meta_file(meta_file_path):
    pae = get_full_pae(meta_file_path)
    return plot_pae(pae)

def plot_normalized_pae_in_meta_file(meta_file_path):
    pae = get_full_pae(meta_file_path)
    return plot_normalized_pae(pae)

def get_interface_pae(top_dir, identifier, design, prediction_i):
    full_pae = get_full_pae(top_dir, identifier, design, prediction_i)
    interface_identifiers = get_interface_residues(top_dir, identifier, design, prediction_i)
    interface_indices = get_interface_residue_indices(top_dir, identifier, design, prediction_i)
    flattend_interface_indices = [index for sublist in interface_indices.values() for index in sublist]
    flattend_interface_indices.sort()
    interface_pae = full_pae[np.ix_(flattend_interface_indices, flattend_interface_indices)]
    return interface_pae, interface_identifiers, interface_indices

def get_cross_interface_pae(top_dir, identifier, design, prediction_i, normalize=True):
    full_pae = get_full_pae(top_dir, identifier, design, prediction_i)
    if normalize:
        full_pae = normalize_pae(full_pae)
    interface_indices = get_interface_residue_indices(top_dir, identifier, design, prediction_i)
    cross_matrices = {}
    for i in interface_indices:
        for j in interface_indices:
            if i != j:
                cross_matrices[(i, j)] = full_pae[np.ix_(interface_indices[i], interface_indices[j])]
    return cross_matrices

def get_cross_interface_pae_statistic(top_dir, identifier, design, prediction_i, statistic_func, normalize=True):
    cross_paes = get_cross_interface_pae(top_dir, identifier, design, prediction_i, normalize)
    flattend_cross_pae = []
    for key in cross_paes:
        flattend_cross_pae.extend(cross_paes[key].flatten())
    return statistic_func(flattend_cross_pae)
