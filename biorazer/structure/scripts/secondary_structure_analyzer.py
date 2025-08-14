import os
from collections import defaultdict
import numpy as np
from Bio.PDB import PDBParser, DSSP
from biorazer.structure.scripts.interface_analyzer import get_interface_residues
from sklearn.decomposition import PCA

def ss_count_to_percentage(ss_counts: dict):
    total_residues = sum(ss_counts.values())
    helix_percentage = round(ss_counts['helix'] / total_residues * 100, 2) if total_residues > 0 else 0
    sheet_percentage = round(ss_counts['sheet'] / total_residues * 100, 2) if total_residues > 0 else 0
    loop_percentage = round(ss_counts['loop'] / total_residues * 100, 2) if total_residues > 0 else 0
    return helix_percentage, sheet_percentage, loop_percentage

def get_pdb_ss(pdb_file, dssp='dssp'):
    '''
    计算 pdb_file 中的二级结构
    返回一个字典, 其中键是链 ID, 值是一个字典, 其中键是残基 ID, 值是二级结构类型
    - H: helix
    - E: sheet
    - G: 3-10 helix
    - I: pi helix
    - S: bend
    - T: turn
    - L: loop
    - C: coil
    - B: bridge
    '''
    
    model = PDBParser(QUIET=True).get_structure("structure", pdb_file)[0]
    dssp = DSSP(model, pdb_file, dssp=dssp)
    res = {}
    for (chain_id, (xx1, res_id, xx2)) in dssp.keys():
        if chain_id not in res:
            res[chain_id] = {}
        res[chain_id][res_id] = dssp[(chain_id, res_id)][2]
    return res

def extract_helices_from_dssp(ss_dict):
    '''
    根据 DSSP 字典提取螺旋信息
    返回一个字典, 其中键是链 ID, 值是一个 tuple, (start_resi, end_resi, helix_class)
    - helix_class: G (3-10 helix), H (alpha helix), I (pi helix)
    - start_resi: 螺旋的起始残基 ID
    - end_resi: 螺旋的结束残基 ID
    
    ss_dict 是一个字典, 其中键是链 ID, 值是一个字典, 其中键是残基 ID, 值是二级结构类型
    '''
    
    res = {}
    for chain_id, ss in ss_dict.items():
        res[chain_id] = []
        helix_start = None
        helix_end = None
        current_helix_type = None
        for res_id, ss_type in ss.items():
            if current_helix_type is None:
                if ss_type in 'HGI':
                    helix_start = res_id
                    current_helix_type = ss_type
                else:
                    continue
            else:
                if ss_type == current_helix_type:
                    continue
                else:
                    helix_end = res_id
                    res[chain_id].append((helix_start, helix_end, current_helix_type))
                    helix_start = None
                    helix_end = None
                    current_helix_type = None
        # 处理最后一个螺旋
        if helix_start is not None:
            helix_end = res_id
            res[chain_id].append((helix_start, helix_end, current_helix_type))
            helix_start = None
            helix_end = None
            current_helix_type = None
    return res

def extract_sheets_from_dssp(ss_dict):
    '''
    根据 DSSP 字典提取片层信息
    返回一个字典, 其中键是链 ID, 值是一个 tuple, (start_resi, end_resi, strand_num)
    - strand_num: 片层的序号
    - start_resi: 片层的起始残基 ID
    - end_resi: 片层的结束残基 ID
    
    ss_dict 是一个字典, 其中键是链 ID, 值是一个字典, 其中键是残基 ID, 值是二级结构类型
    '''
    
    res = {}
    for chain_id, ss in ss_dict.items():
        res[chain_id] = []
        sheet_start = None
        sheet_end = None
        strand_num = None
        for res_id, ss_type in ss.items():
            if ss_type == 'E':
                if sheet_start is None:
                    sheet_start = res_id
                    strand_num = 1
                sheet_end = res_id
            else:
                if sheet_start is not None:
                    res[chain_id].append((sheet_start, sheet_end, strand_num))
                    sheet_start = None
                    sheet_end = None
                    strand_num = None
        if sheet_start is not None:
            res[chain_id].append((sheet_start, sheet_end, strand_num))
    return res

def calc_interface_ss_percentage_of_af_pdb(pdb_file, chain_analyzed, chain_interacting, atom_distance_cutoff=4.0, dssp='dssp'):
    '''
    Return 
    - chain ss_percentage
        - helix_percentage
        - sheet_percentage
        - loop_percentage
    - interface_ss_percentage
        - helix_interface_percentage
        - sheet_interface_percentage
        - loop_interface_percentage
    - plddt of all interface residues (with / without secondary structures)
    - plddt of all residues with secondary structures
    '''
    
    model = PDBParser(QUIET=True).get_structure("structure", pdb_file)[0]
    dssp = DSSP(model, pdb_file, dssp=dssp)
    
    ss_counts = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface = []
    plddts_ss = []

    interface_residues = get_interface_residues(
        pdb_file, chain_analyzed, chain_interacting, 
        atom_distance_cutoff=atom_distance_cutoff
    )
    interface_res_ids = [res_id for chain_id, res_id in interface_residues]
    
    for residue in model[chain_analyzed]:
        res_id = residue.get_id()[1]
        if (chain_analyzed, res_id) in dssp:
            ss = dssp[(chain_analyzed, res_id)][2]
            if ss in 'HGI':
                ss_type = 'helix'
            elif ss in 'E':
                ss_type = 'sheet'
            else:
                ss_type = 'loop'
            ss_counts[ss_type] += 1
            
            if ss_type != 'loop':
                plddts_ss.append(np.mean([atom.bfactor for atom in residue]))
            
            if res_id in interface_res_ids:
                ss_interface_counts[ss_type] += 1
                plddts_interface.append(np.mean([atom.bfactor for atom in residue]))
            
    helix_percentage, sheet_percentage, loop_percentage = ss_count_to_percentage(ss_counts)
    helix_interface_percentage, sheet_interface_percentage, loop_interface_percentage = ss_count_to_percentage(ss_interface_counts)
    i_plddt = round(np.mean(plddts_interface), 2) if len(plddts_interface) > 0 else 0
    ss_plddt = round(np.mean(plddts_ss), 2) if len(plddts_ss) > 0 else 0
    
    return helix_percentage, sheet_percentage, loop_percentage, helix_interface_percentage, sheet_interface_percentage, loop_interface_percentage, i_plddt, ss_plddt

def generate_pdb_helix_row(
    helix_no: int, chain_id: str, 
    start_resn: str, start_resi: int, 
    end_resn: str, end_resi: int, 
    helix_class: int, 
    start_segi: str="", end_segi: str="",
    helix_id: str="", comment=""
):
    if not helix_id:
        helix_id = str(helix_no)
    helix_length = end_resi - start_resi + 1
    row_str = f"HELIX  {helix_no:>3} {helix_id:>3} {start_resn.upper()} {chain_id} {start_resi:>4}{start_segi:>1} {end_resn.upper()} {chain_id} {end_resi:>4}{end_segi:>1}{helix_class:>2}{comment:>30} {helix_length:5}"
    return row_str

def generate_pdb_sheet_row(
    strand_no: int, sheet_id: str, strand_num: int, chain_id: str,
    start_resn: str, start_resi: int,
    end_resn: str, end_resi: int,
    sense: int,
    strand_start_atom: str="", strand_start_resn: str="", strand_start_resi: int=-1,
    strand_end_atom: str="", strand_end_resn: str="", strand_end_resi: int=-1,
):
    '''
    This functionality is too complicated. Abandoned.
    '''
    if strand_start_resi == -1:
        strand_start_resi = ""
    if strand_end_resi == -1:
        strand_end_resi = ""
    row_str = f"SHEET  {strand_no:>3} {sheet_id:>3}{strand_num:2} {start_resn.upper()} {chain_id}{start_resi:>4} {end_resn.upper()} {chain_id}{end_resi:>4} {sense:>2} {strand_start_atom:>4} {strand_start_resn.upper()} {chain_id if strand_end_resi else '':>1} {strand_start_resi:>4} {strand_end_atom:>4} {strand_end_resn.upper()} {chain_id if strand_end_resi else '':>1} {strand_end_resi:>4}"
    return row_str

def generate_ss_pdb(pdb_file, dssp="dssp"):
    '''
    Input: pdb_file
    Output: pdb_file with secondary structure information
    '''
    pdb_file_out = os.path.splitext(pdb_file)[0] + "_ss.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    dssp = DSSP(structure[0], pdb_file, dssp=dssp)
    # 写入二级结构信息 (HELIX 与 SHEET 行) 以及原始结构到新的 PDB 文件 pdb_file_out 
    
