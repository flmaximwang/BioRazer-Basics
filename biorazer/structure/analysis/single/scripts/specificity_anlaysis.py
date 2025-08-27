import pandas as pd
import numpy as np
import collections
import os

def read_macro(macro_str):
    '''
    A macro looks like this:
    /model/segi/chain/resn`resi/name
    model: model name
    segi: segment identifier
    chain: chain identifier
    resn: residue name
    resi: residue number
    name: atom name
    '''
    import re
    macro_pattern = re.compile(r'/(?P<model>\w*)/(?P<segi>\w*)/(?P<chain>\w)/(?P<resn>\w+)`(?P<resi>\d+)/(?P<name>\w+)')
    match = macro_pattern.match(macro_str)
    if match:
        return match.groupdict()
    else:
        raise ValueError("Invalid macro string: %s" % macro_str)

def chain_priority(chain):
    return sum([ord(c) for c in chain])

def name_priority(name):
    return sum([ord(c) for c in name])

def atom_priority(atom):
    atom_macro = read_macro(atom)
    return (
        chain_priority(atom_macro['chain']),
        int(atom_macro['resi']),
        name_priority(atom_macro['name'])
    )

def analyze_specificity_with_contacting_atoms(
        contacting_atom_lists: list[str],
        output_file: str = "./Summary.csv"
    ):
    '''
    Analyzing the specificity of the contacting atoms in the contacting atom list files\n
    - contacting_atom_lists: list of file paths to the contacting atom list files\n
        - a contacting atom list file is a csv file with columns: atom1, atom2, distance (Å)\n
    '''
    basename_list = [os.path.basename(file_path) for file_path in contacting_atom_lists] 

    res_dict_temp = {}
    for file_path, basename in zip(contacting_atom_lists, basename_list):
        temp = pd.read_csv(file_path)
        for i, row in temp.iterrows():
            atom1_macro = read_macro(row['atom1'])
            atom2_macro = read_macro(row['atom2'])
            contact_tuple = (
                f"//{atom1_macro['segi']}/{atom1_macro['chain']}/{atom1_macro['resn']}`{atom1_macro['resi']}/{atom1_macro['name']}",
                f"//{atom2_macro['segi']}/{atom2_macro['chain']}/{atom2_macro['resn']}`{atom2_macro['resi']}/{atom2_macro['name']}"
            )
            if not contact_tuple in res_dict_temp:
                res_dict_temp[contact_tuple] = {}
            res_dict_temp[contact_tuple][basename] = row['distance (Å)']

    res_dict = collections.OrderedDict(sorted(res_dict_temp.items(), key=lambda x: atom_priority(x[0][0])))
    res_df = pd.DataFrame({
        "atom1": [],
        "atom2": []
    })
    for basename in basename_list:
        res_df[basename] = []
    for my_key in res_dict:
        res_df_new_line = pd.DataFrame({
            "atom1": [my_key[0]],
            "atom2": [my_key[1]],
        })
        for basename in res_dict[my_key]:
            res_df_new_line[basename] = [res_dict[my_key][basename]]
        res_df = pd.concat([res_df, pd.DataFrame(res_df_new_line)], ignore_index=True)

    # 用同一相互作用在不同模型中的出现次数作为 Specificity Index
    res_df["Specificity Index"] = np.zeros(res_df.shape[0])
    for i, row in res_df.iterrows():
        res_df.at[i, "Specificity Index"] = row[2:-1].count()
        res_df.to_csv(output_file, index=False)
    
    ordering = list(res_df.columns)
    for col in ["atom1", "atom2", "Specificity Index"]:
        ordering.remove(col)
    ordering = ["atom1", "atom2", "Specificity Index"] + ordering
    print(ordering)
    res_df = res_df[ordering]
    
    return res_df