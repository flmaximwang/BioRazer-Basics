import json, os
import xml.etree.ElementTree as ET
import numpy as np
import pyrosetta as pr
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector, LayerSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import TotalEnergyMetric, SasaMetric
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
import pymol
from scipy.spatial import cKDTree
from biotite.sequence import ProteinSequence
from Bio.PDB import PDBParser, Selection, Select, PDBIO

def get_interface_residues(
    pdb_file, 
    chain_analyzed, chain_interacting, 
    atom_distance_cutoff=4.0, 
    PDBParser_kwargs={'QUIET': True}
):
    '''
    使用 KDTree 算法计算两个链之间的接触残基。
    - 只要两个链中的原子之间的距离小于 atom_distance_cutoff，就认为这两个原子是接触的, 从而认为 2 个残基是接触的。
    - Supposing chain_analyzed is the chain of interest, this function returns the residues of the chain that are in contact with the chain_interacting.
    
    Returns: [(resname, resid), ...]
    '''
    
    # Load the PDB file
    parser = PDBParser(**PDBParser_kwargs)
    structure = parser.get_structure("structure", pdb_file)
    
    # Get the atoms from the chains
    ## Unfold entities to atom level
    atoms_analyzed = Selection.unfold_entities(structure[0][chain_analyzed], "A")
    atoms_interacting = Selection.unfold_entities(structure[0][chain_interacting], "A")
    # 去除氢原子
    atoms_analyzed = [atom for atom in atoms_analyzed if atom.get_id() != 'H']
    atoms_interacting = [atom for atom in atoms_interacting if atom.get_id() != 'H']
    coords_analyzed = np.array([atom.get_coord() for atom in atoms_analyzed])
    coords_interacting = np.array([atom.get_coord() for atom in atoms_interacting])
    tree_analyzed = cKDTree(coords_analyzed)
    tree_interacting = cKDTree(coords_interacting)
    
    # Get the residues at the interface
    pairs = tree_analyzed.query_ball_tree(tree_interacting, atom_distance_cutoff)
    interface_residues = set()
    for idx_analyzed, close_indices in enumerate(pairs):
        if len(close_indices) == 0:
            continue
        residue = atoms_analyzed[idx_analyzed].get_parent()
        if residue.get_resname() in ProteinSequence._dict_3to1:
            resname = ProteinSequence.convert_letter_3to1(residue.get_resname())
        else:
            resname = residue.get_resname()
        resid = residue.get_id()[1]
        interface_residues.add((resname, resid))
    interface_residues = list(interface_residues)
    interface_residues.sort(key=lambda x: x[1])
    
    return interface_residues

def score_interface_bindcraft_old(
    pdb_file, target_chain, binder_chain,
):
    '''
    按照 BindCraft 的指标对界面进行评估. 相关指标包括
    - binder_energy_score: REF2015 的标准能量函数 (kcal/mol)
    - *binder_surface_hydrophobicity: target_chain 的表面氨基酸中疏水氨基酸的比例 (归一化的), 
        通过 pyrosetta.core.select.residue_selector.LayerSelector 可以筛选出表面残基. 
        BindCraft 选择这个指标 < 0.35 的 binder 以避免设计的 binder 过于疏水
    - *interface_sc: ShapeComplementarity. 由 InterfaceAnalyzerMover 计算的形状互补性 (归一化的). 
        BindCraft 选择这个指标 > 0.55 的复合体.
    - interface_packstat: 由 InterfaceAnalyzerMover 计算的 PackStat 分数 (归一化的). 
        一般认为 PackStat 需要 > 0.65, 但是 BindCraft 没有据此进行选择
    - interface_dG: 由 InterfaceAnalyzerMover 计算的界面结合自由能 (kcal/mol)
    - interface_dSASA: 由 InterfaceAnalyzerMover 计算的界面表面积 (Å^2)
    - interface_dG_SASA_ratio: 上述 2 个指标的比值
    - interface_fraction: binder 的界面表面积占 binder 的总表面积的比例 (归一化的)
    - interface_hydrophobicity: binder 的界面氨基酸中疏水氨基酸的比例 (归一化的)
    - interface_nres: binder 的界面残基数, 显然它和 interface——dSASA 有关
    
    '''
    pose = pr.pose_from_pdb(pdb_file)
    iam = InterfaceAnalyzerMover()
    iam.set_interface(f"{target_chain}_{binder_chain}")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)
    
    # Calculating the interface hydrophobicity of the binder
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    hydrophobic_aa = 'ACFILMVWY'
    binder_interface_residues = get_interface_residues(pdb_file, binder_chain, target_chain)
    for resname, resid in binder_interface_residues:
        interface_AA[resname] += 1
    hydrophobicity_count = sum([interface_AA[aa] for aa in hydrophobic_aa])
    interface_nres = sum(interface_AA.values())
    if interface_nres == 0:
        interface_hydrophobicity = - 1
    else:
        interface_hydrophobicity = hydrophobicity_count / interface_nres * 100
    
    # Calculating buried unsat hbonds
    buns_filter_xml = ET.Element("BuriedUnsatHbonds")
    buns_filter_xml.set("report_all_heavy_atom_unsats", "true")
    buns_filter_xml.set("scorefxn", "scorefxn")
    buns_filter_xml.set("ignore_surface_res", "false")
    buns_filter_xml.set("use_ddG_style", "true")
    buns_filter_xml.set("dalphaball_sasa", "1")
    buns_filter_xml.set("probe_radius", "1.1")
    buns_filter_xml.set("burial_cutoff_apo", "0.2")
    buns_filter_xml.set("confidence", "0")
    xml_str = ET.tostring(buns_filter_xml).decode()
    buns_filter = XmlObjects.static_get_filter(xml_str)
    interface_unsat_hbonds = buns_filter.report_sm(pose)
    interface_unsat_hbonds_percentage = interface_unsat_hbonds / interface_nres * 100 if interface_nres != 0 else -1
    
    # Retrieving other statistics
    interface_score = iam.get_all_data()
    interface_sc = interface_score.sc_value # Shape Complementarity
    interface_hbonds = interface_score.interface_hbonds # Number of interface hydrogen bonds
    interface_hbond_percentage = interface_hbonds / interface_nres * 100 if interface_nres != 0 else -1 # Percentage of interface hydrogen bonds
    interface_dG = iam.get_interface_dG() # Interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # Interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # Interface packstat score
    interface_dG_SASA_ratio = interface_score.dG_dSASA_ratio * 100 # Ratio of dG / dSASA (normalized energy for interface area size)
    
    # Calculating binder energy score
    tem = TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(ChainSelector(binder_chain))
    binder_score = tem.calculate(pose)
    
    # Calculating binder SASA fraction
    bsasa = SasaMetric()
    bsasa.set_residue_selector(ChainSelector(binder_chain))
    binder_sasa = bsasa.calculate(pose)
    
    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0
    
    # Calculating surface hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains() + 1), pose.split_by_chain())}[binder_chain]
    
    layer_sel = LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(binder_pose)
    exp_apol_count = 0
    total_count = 0
    for i in range(1, len(surface_res) + 1):
        if surface_res[i]:
            res = binder_pose.residue(i)
            if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name()  == 'TYR':
                exp_apol_count += 1
            total_count += 1
    surface_hydrophobicity = exp_apol_count / total_count
    
    interface_scores = {
        'binder_score': binder_score,
        'surface_hydrophobicity': surface_hydrophobicity,
        'interface_sc': interface_sc,
        'interface_packstat': interface_packstat,
        'interface_dG': interface_dG,
        'interface_dSASA': interface_dSASA,
        'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
        'interface_fraction': interface_binder_fraction,
        'interface_hydrophobicity': interface_hydrophobicity,
        'interface_nres': interface_nres,
        'interface_hbonds': interface_hbonds,
        'interface_hbond_percentage': interface_hbond_percentage,
        'interface_delta_unsat_hbonds': interface_unsat_hbonds,
        'interface_delta_unsat_hbonds_percentage': interface_unsat_hbonds_percentage,
    }
    
    return interface_scores

def score_interface_bindcraft(
    pdb_file, target_chains='A', binder_chains='B',
):
    '''
    Args:
    - target_chains: binder 结合的链段, 可以是 1/2 个链
    - binder_chains: binder 的链段, 可以是 1/2 个链
    target_chains + binder_chains 不能超过 3 条链, 不然 InterfaceAnalyzerMover 会报错
    Return: 按照 BindCraft 的指标对界面进行评估. 相关指标包括
    - *interface_sc: ShapeComplementarity. 由 InterfaceAnalyzerMover 计算的形状互补性 (归一化的). 
        BindCraft 选择这个指标 > 0.55 的复合体.
    - interface_packstat: 由 InterfaceAnalyzerMover 计算的 PackStat 分数 (归一化的). 
        一般认为 PackStat 需要 > 0.65, 但是 BindCraft 没有据此进行选择
    - interface_dG: 由 InterfaceAnalyzerMover 计算的界面结合自由能 (kcal/mol)
    - interface_dSASA: 由 InterfaceAnalyzerMover 计算的界面表面积 (Å^2)
    - interface_dG_SASA_ratio: 上述 2 个指标的比值
    - interface_fraction: binder 的界面表面积占 binder 的总表面积的比例 (归一化的)
    - interface_hydrophobicity: binder 的界面氨基酸中疏水氨基酸的比例 (归一化的)
    - interface_nres: binder 的界面残基数, 显然它和 interface——dSASA 有关
    - interface_hbonds: binder 的界面氢键数
    - interface_hbond_percentage: binder 的界面氢键数占总界面残基数的比例. 
        理论上 1 个残基最多有 2 个氢键, 所以这个比例实际上有可能大于 1
    - *interface_delta_unsat_hbonds: binder 的界面残基中埋藏的不饱和氢键数
        BindCraft 选择这个指标 <= 3 的 binder 以避免设计的 binder 优先与水形成氢键
    - interface_delta_unsat_hbonds_percentage: binder 的界面残基中埋藏的不饱和氢键数占总界面残基数的比例
        理论上 1 个残基最多有 2 个氢键, 所以这个比例实际上有可能大于 1
    Columns: ['interface_sc', 'interface_packstat', 'interface_dG', 'interface_dSASA', 'interface_dG_SASA_ratio', 'interface_fraction', 'interface_hydrophobicity', 'interface_nres', 'interface_hbonds', 'interface_hbond_percentage', 'interface_delta_unsat_hbonds', 'interface_delta_unsat_hbonds_percentage']
    '''
    if len(target_chains) + len(binder_chains) not in [2, 3]:
        raise ValueError("target_chains and binder_chain should have a total length of 2 or 3.")
    
    pose = pr.pose_from_pdb(pdb_file)
    iam = InterfaceAnalyzerMover()
    iam.set_interface(f"{target_chains}_{binder_chains}")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True) # 这个指标对于较为准确地计算 dG 是必要的
    iam.apply(pose)
    
    # Calculating the interface hydrophobicity of the binder
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    hydrophobic_aa = 'ACFILMVWY'
    binder_interface_residues = []
    for target_chain in target_chains:
        for binder_chain in binder_chains:
            binder_interface_residues.extend(get_interface_residues(pdb_file, binder_chain, target_chain))
    for resname, resid in binder_interface_residues:
        if resname in ProteinSequence._dict_1to3:
            interface_AA[resname] += 1
    hydrophobicity_count = sum([interface_AA[aa] for aa in hydrophobic_aa])
    interface_nres = sum(interface_AA.values())
    if interface_nres == 0:
        interface_hydrophobicity = 0
    else:
        interface_hydrophobicity = hydrophobicity_count / interface_nres
    
    # Calculating buried unsat hbonds
    buns_filter_xml = ET.Element("BuriedUnsatHbonds")
    buns_filter_xml.set("report_all_heavy_atom_unsats", "true")
    buns_filter_xml.set("scorefxn", "scorefxn")
    buns_filter_xml.set("ignore_surface_res", "false")
    buns_filter_xml.set("use_ddG_style", "true") # Setting this to true makes sure that you calculate only the interface unsat hbonds
    buns_filter_xml.set("dalphaball_sasa", "1")
    buns_filter_xml.set("probe_radius", "1.1")
    buns_filter_xml.set("burial_cutoff_apo", "0.2")
    buns_filter_xml.set("confidence", "0")
    xml_str = ET.tostring(buns_filter_xml).decode()
    buns_filter = XmlObjects.static_get_filter(xml_str)
    
    tmp_pdb_path = 'tmp.pdb'
    class UnsatHbondSelect(Select):
        def __init__(self, chain_ids):
            self.chain_ids = chain_ids
            
        def accept_chain(self, chain):
            return chain.get_id() in self.chain_ids
        
        def accept_residue(self, residue):
            # 去除修饰分子
            return residue.get_resname() in list(ProteinSequence._dict_3to1) + ['HOH']
        
    interface_unsat_hbonds = 0
    for target_chain in target_chains:
        for binder_chain in binder_chains:
            # 仅仅选择 target_chain 和 binder_chain 并导出到 tmp.pdb
            my_struc = PDBParser(QUIET=True).get_structure("tmp", pdb_file)
            io = PDBIO()
            io.set_structure(my_struc)
            io.save(tmp_pdb_path, UnsatHbondSelect([target_chain, binder_chain]))
            pose = pr.pose_from_pdb(tmp_pdb_path)
            interface_unsat_hbonds += buns_filter.report_sm(pose)
            os.remove(tmp_pdb_path)
            
    interface_unsat_hbonds_percentage = interface_unsat_hbonds / interface_nres if interface_nres != 0 else 0
    
    # Retrieving other statistics
    interface_score = iam.get_all_data()
    interface_sc = interface_score.sc_value # Shape Complementarity
    interface_hbonds = interface_score.interface_hbonds # Number of interface hydrogen bonds
    interface_hbond_percentage = interface_hbonds / interface_nres if interface_nres != 0 else 0 # Percentage of interface hydrogen bonds
    interface_dG = iam.get_interface_dG() # Interface dG
    interface_dSASA = iam.get_interface_delta_sasa() # Interface dSASA (interface surface area)
    interface_packstat = iam.get_interface_packstat() # Interface packstat score
    interface_dG_SASA_ratio = interface_score.dG_dSASA_ratio # Ratio of dG / dSASA (normalized energy for interface area size)
    
    # Calculating binder SASA fraction
    bsasa = SasaMetric()
    if len(binder_chains) == 1:
        combined_selector = ChainSelector(binder_chains[0])
    elif len(binder_chains) == 2:
        combined_selector = ChainSelector(binder_chains[0]) | ChainSelector(binder_chains[1])
    else:
        raise ValueError("binder_chains should have a length of 1 or 2.")
    bsasa.set_residue_selector(combined_selector)
    binder_sasa = bsasa.calculate(pose)
    
    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa)
    else:
        interface_binder_fraction = 0
    
    interface_scores = {
        'interface_sc': interface_sc,
        'interface_packstat': interface_packstat,
        'interface_dG': interface_dG,
        'interface_dSASA': interface_dSASA,
        'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
        'interface_fraction': interface_binder_fraction,
        'interface_hydrophobicity': interface_hydrophobicity,
        'interface_nres': interface_nres,
        'interface_hbonds': interface_hbonds,
        'interface_hbond_percentage': interface_hbond_percentage,
        'interface_delta_unsat_hbonds': interface_unsat_hbonds,
        'interface_delta_unsat_hbonds_percentage': interface_unsat_hbonds_percentage,
    }
    
    return interface_scores

def binding_induced_rmsd(
    unbound_pdb, bound_pdb, 
    chain_analyzed_unbound, chain_analyzed_bound,
    export = False
):
    '''
    This function calculates the RMSD between the unbound and bound structures.
    '''
    pymol.cmd.reinitialize()
    pymol.cmd.load(unbound_pdb, 'unbound')
    pymol.cmd.load(bound_pdb, 'bound')
    rmsd = pymol.cmd.align(
        f'bound and c. {chain_analyzed_bound}', 
        f'unbound and chain {chain_analyzed_unbound}',
        cycles=0
    )[0]
    if export:
        pymol.cmd.save(f'{bound_pdb[:-4]}_aligned.pdb', 'bound')
    return rmsd

def display_interface(pdb_file, 
    chains_1, chains_2, 
    atom_distance_cutoff=4.0, 
    PDBParser_kwargs={'QUIET': True},
    verbose=True
):
    '''
    This function helps to annotate the interface between chains_1 and chains_2 in the pdb_file.
    chains_1 and chains_2 can containe multiple chains, e.g. 'AB' and 'CD'.
    '''
    pymol.cmd.reinitialize()
    pymol.cmd.load(pdb_file)
    for chain_1 in chains_1:
        for chain_2 in chains_2:
            if verbose: print(">>> Analyzing the interface between chains", chain_1, "and", chain_2)
            chain1_interface_residues = get_interface_residues(
                pdb_file, chain_1, 
                chain_2, atom_distance_cutoff, PDBParser_kwargs
            )
            chain2_interface_residues = get_interface_residues(
                pdb_file, chain_2, 
                chain_1, atom_distance_cutoff, PDBParser_kwargs
            )
            interface1_ids = list(map(lambda x: str(x[1]), chain1_interface_residues))
            interface1_ids = "+".join(interface1_ids)
            interface2_ids = list(map(lambda x: str(x[1]), chain2_interface_residues))
            interface2_ids = "+".join(interface2_ids)
            pymol.cmd.create(f"{chain_1}_{chain_2}", f"chain {chain_1} and resi {interface1_ids}")
            pymol.cmd.create(f"{chain_2}_{chain_1}", f"chain {chain_2} and resi {interface2_ids}")
            
    # Show packing
    for chain_1 in chains_1:
        for chain_2 in chains_2:
            pymol.cmd.show("sphere", f"{chain_1}_{chain_2} or {chain_2}_{chain_1}")
    pymol.cmd.save(f"{pdb_file[:-4]}_packing.pse")
    
    # Show hydrogen bonds
    pymol.cmd.hide("spheres")
    for chain_1 in chains_1:
        for chain_2 in chains_2:
            pymol.cmd.show("sticks", f"{chain_1}_{chain_2} or {chain_2}_{chain_1}")
            pymol.cmd.distance(f"{chain_1}_{chain_2}_hbonds", f"{chain_1}_{chain_2} and not hydrogens", f"{chain_2}_{chain_1} and not hydrogens", mode=2)
    pymol.cmd.save(f"{pdb_file[:-4]}_hbonds.pse")

base_cols_to_consider = [
    "pLDDT",
    "pTM",
    "i_pTM",
    "pAE",
    "i_pAE",
    "i_pLDDT",
    "ss_pLDDT",
    # "Unrelaxed_Clashes",
    # "Relaxed_Clashes",
    "Binder_Energy_Score",
    "Surface_Hydrophobicity",
    "ShapeComplementarity",
    "PackStat",
    "dG",
    "dSASA",
    "dG/dSASA",
    "Interface_SASA_%",
    "Interface_Hydrophobicity", # Some are NaN
    "n_InterfaceHbonds",
    "InterfaceHbondsPercentage",
    "n_InterfaceUnsatHbonds",
    "InterfaceUnsatHbondsPercentage",
    "Interface_Helix%",
    "Interface_BetaSheet%",
    "Interface_Loop%",
    "Binder_Helix%",
    "Binder_BetaSheet%",
    "Binder_Loop%",
    # "Hotspot_RMSD",
    # "Target_RMSD",
    "Binder_pLDDT",
    "Binder_pTM",
    "Binder_pAE",
    "Binder_RMSD",
]