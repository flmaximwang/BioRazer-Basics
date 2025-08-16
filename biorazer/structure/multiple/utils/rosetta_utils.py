import os
from .pre_processors import *

def get_interface_annotated_pdb_path(top_dir, identifier, design, prediction_i):
    design = format_design_name(design)
    prediction_dir = generate_prediction_dir(top_dir, identifier)
    return os.path.join(prediction_dir, design, f"fold_{design}_model_{prediction_i}_rosetta_interface", f"fold_{design}_model_{prediction_i}_0001.pdb")

def get_unsat_residue_pymol_selector(top_dir, identifier, design, prediction_i):
    '''
    获取未饱和氢键的残基选择器
    '''
    interface_annotated_pdb_path = get_interface_annotated_pdb_path(top_dir, identifier, design, prediction_i)
    rosetta_pdb_file_path = interface_annotated_pdb_path
    with open(rosetta_pdb_file_path, 'r') as file:
        lines = file.readlines()
    
    flag = False
    counter = 0
    while(not flag):
        line = lines[counter]
        if "pymol-style selection for unstat hbond res" in line:
            flag = True
            selection_command = lines[counter + 1][:-1]
        else:
            counter += 1
    
    if not flag:
        raise ValueError(f"No unsatisfied hbond residues found in the pdb file {rosetta_pdb_file_path}")
    
    return selection_command

def get_interface_pymol_selector(interface_annotated_pdb_path):
    '''
    获取界面残基的选择器
    '''
    rosetta_pdb_file_path = interface_annotated_pdb_path
    with open(rosetta_pdb_file_path, 'r') as file:
        lines = file.readlines()
    
    flag = False
    counter = 0
    while(not flag):
        line = lines[counter]
        if "pymol-style selection for interface res" in line:
            flag = True
            selection_command = lines[counter + 1][:-1]
        else:
            counter += 1
    
    if not flag:
        raise ValueError(f"No interface residues found in the pdb file {rosetta_pdb_file_path}")
    
    return selection_command

def get_interface_residues(interface_annotated_pdb_path):
    '''
    :param interface_annotated_pdb_path: str
    :return: list of chain identifiers [(chain_id, res_id), ...]
    '''
    
    selection_command = get_interface_pymol_selector(interface_annotated_pdb_path)
    selection_macro = selection_command.split(', ')[1]
    macro1, macro2 = selection_macro.split(' + ')
    
    def generate_identifiers_from_macro(macro: str):
        chain = macro.split('/')[3]
        res_ids: list = macro.split('/')[4].split('+')
        while '' in res_ids:
            res_ids.remove('')
        identifiers = [(chain, int(res_id)) for res_id in res_ids]
        return identifiers
    
    identifiers_1 = generate_identifiers_from_macro(macro1)
    identifiers_2 = generate_identifiers_from_macro(macro2)
    
    return identifiers_1 + identifiers_2

def run_interface_analyzer_for(top_dir, identifier, design, verbose=True):
    '''
    使用 Rosetta 的 InterfaceAnalyzer 计算界面参数
    '''
    design = format_design_name(design)
    prediction_dir = generate_prediction_dir(top_dir, identifier)
    for i in range(5):
        name_prefix = f"{design}_model_{i}"
        input_filepath = os.path.join(prediction_dir, design, f"fold_{name_prefix}.pdb")
        out_score_path = os.path.join(prediction_dir, design, f"fold_{name_prefix}_rosetta_interface.sc")
        out_pdb_dir = os.path.join(prediction_dir, design, f"fold_{name_prefix}_rosetta_interface")
        out_pdb_path = os.path.join(out_pdb_dir, f"fold_{name_prefix}_0001.pdb")
        if not os.path.exists(out_pdb_dir):
            os.mkdir(out_pdb_dir)
        if not os.path.exists(out_score_path) or not os.path.exists(out_pdb_path):
            if os.path.exists(out_score_path):
                os.remove(out_score_path)
            if os.path.exists(out_pdb_path):
                os.remove(out_pdb_path)
            # Passing -i to source .zshrc
            # Passing -c to run the command, not accepting a command file (script)
            # By not passing -overwrite, the command will not overwrite existing files
            command_list = ["/bin/zsh", "-i", "-c", f"InterfaceAnalyzer.macosclangrelease -s {input_filepath} -out:path:pdb {out_pdb_dir} -out:file:scorefile {out_score_path} -mute core basic -overwrite"]
            print(" ".join(command_list))
            subprocess.run(command_list)
        else:
            if verbose: print(f">>> Skipping because {out_score_path} and {out_pdb_dir} already exist")

def read_interface_sc(top_dir, identifier, design, prediction_i):
    '''
    SCORE: 
    - total_score 
    - complex_normalized
    - dG_cross 
    - dG_cross/dSASAx100 
    - dG_separated 
    - dG_separated/dSASAx100 
    - dSASA_hphobic 
    - dSASA_int 
    - dSASA_polar 
    - delta_unsatHbonds 
    - hbond_E_fraction
    - hbonds_int
    - nres_all
    - nres_int   
    - packstat 
    - per_residue_energy_int               
    - sc_value side1_normalized side1_score side2_normalized side2_score description     
    '''
    design_dir = generate_prediction_dir(top_dir, identifier)
    design = format_design_name(design)
    sc_file_path = os.path.join(design_dir, design, f"fold_{design}_model_{prediction_i}_rosetta_interface.sc")
    with open(sc_file_path, 'r') as file:
        lines = file.readlines()
    
    headers = lines[1].strip().split()[1:]  # 获取标题行并去掉 "SCORE:"
    values = lines[2].strip().split()[1:]   # 获取值行并去掉 "SCORE:"
    
    score_dict = dict(zip(headers, values))
    for dict_key in score_dict:
        if not dict_key == 'description':
            score_dict[dict_key] = float(score_dict[dict_key])
    return score_dict

def get_interface_shape_complementarity(top_dir, identifier, design, prediction_i):
    return read_interface_sc(top_dir, identifier, design, prediction_i)["sc_value"]

def get_unsat_hbond_num(top_dir, identifier, design, prediction_i):
    return read_interface_sc(top_dir, identifier, design, prediction_i)["delta_unsatHbonds"]
