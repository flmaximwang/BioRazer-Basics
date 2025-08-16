import os
import collections
import numpy as np
import pandas as pd
import ammolite
import biotite.structure as struc
from .utils import macro, progress

def find_interactions_between_domains(
        domain1: struc.AtomArray, domain2: struc.AtomArray, cutoff=5.0,
        interaction_name="Interactions for Domain I and Domain II",
        output_dir=".", ignore_hydrogens=True
    ):
    """
    This Python function finds interactions between atoms in two structures and saves the results to a
    CSV file.
    
    :param domain1: `domain1` is an `AtomArray` representing the first domain structure
    :type domain1: struc.AtomArray
    :param domain2: `domain2` is a parameter representing the second structure (domain) for which you
    want to find interactions with `domain1`. It is expected to be of type `struc.AtomArray`
    :type domain2: struc.AtomArray
    :param cutoff: The `cutoff` parameter in the `find_interactions_between_domains` function specifies
    the distance threshold within which atoms from `domain1` and `domain2` are considered to be
    interacting. Atoms that are closer to each other than the specified `cutoff` distance will be
    identified as interacting
    :param interaction_name: The `interaction_name` parameter is a string that specifies the name of the
    interaction being analyzed. It is used as part of the output file name and for display purposes
    during the execution of the function, defaults to Interactions for Domain I and Domain II (optional)
    :param output_dir: The `output_dir` parameter specifies the directory where the output file will be
    saved. If you do not specify a directory, the default directory will be the current working
    directory (denoted by `.`). You can provide a specific directory path where you want the results to
    be saved. For example,, defaults to . (optional)
    :param ignore_hydrogens: The `ignore_hydrogens` parameter in the `find_interactions_between_domains`
    function is a boolean flag that determines whether hydrogen atoms should be excluded from the
    analysis when finding interactions between two protein domains. If set to `True`, the function will
    filter out hydrogen atoms from both `domain1`, defaults to True (optional)
    """
    
    print(f">>> Job started for \"{interaction_name}\"...")
    
    if ignore_hydrogens:
        domain1 = domain1[domain1.get_annotation("element") != "H"]
        domain2 = domain2[domain2.get_annotation("element") != "H"]
    atom_list1 = []
    atom_list2 = []
    distance_list = []
    
    def going_through_all_atoms_to_find_interacting_molecules():
        counter = 0
        counter_max = 100000
        for atom1 in domain1:
            for atom2 in domain2:
                if abs(atom1.coord[0] - atom2.coord[0]) >= cutoff:
                    counter += 1
                    if counter < counter_max:
                        continue
                    else:
                        yield counter
                        counter = 0
                if abs(atom1.coord[1] - atom2.coord[1]) >= cutoff:
                    counter += 1
                    if counter < counter_max:
                        continue
                    else:
                        yield counter
                        counter = 0
                if abs(atom1.coord[2] - atom2.coord[2]) >= cutoff:
                    counter += 1
                    if counter < counter_max:
                        continue
                    else:
                        yield counter
                        counter = 0
                distance_between_atoms = struc.distance(atom1, atom2)
                if distance_between_atoms < cutoff:
                    atom_list1.append(atom1)
                    atom_list2.append(atom2)
                    distance_list.append((atom1, atom2, distance_between_atoms))
                counter += 1
                if counter < counter_max:
                    continue
                else:
                    yield counter
                    counter = 0
    
    print(">>> Going through all atoms to find interacting atoms...")
    pair_num = len(domain1) * len(domain2)
    my_bar = progress.MyBar(0, pair_num, "Finding...")
    counter = 0
    for i in going_through_all_atoms_to_find_interacting_molecules():
        counter += i
        my_bar.set_bar_value(counter)
    my_bar.finish()
    
    print(">>> Generating contact matrix...")
    contact_matrix = np.zeros((len(atom_list1), len(atom_list2)))
    for atom1, atom2, distance_between_atoms in distance_list:
        contact_matrix[atom_list1.index(atom1), atom_list2.index(atom2)] = distance_between_atoms
    
    print(">>> Generating contact atom dataframe...")
    contact_atom_df = pd.DataFrame({
        "Atom1": ["" for i in range(len(atom_list1) * len(atom_list2))],
        "Atom2": ["" for i in range(len(atom_list1) * len(atom_list2))],
        "Distance (Å)": np.zeros(len(atom_list1) * len(atom_list2), dtype=np.float64)
    })
    for i, atom1 in enumerate(atom_list1):
        for j, atom2 in enumerate(atom_list2):
            my_distance = contact_matrix[i, j]
            if my_distance > 0:
                contact_atom_df.loc[i * len(atom_list2) + j] = [macro.atom_to_macro(atom1), macro.atom_to_macro(atom2), contact_matrix[i, j]]
    
    # 去除所有 distance 为 0 的行
    contact_atom_df = contact_atom_df[contact_atom_df["Distance (Å)"] > 0]
    contact_atom_df.reset_index(drop=True, inplace=True)
    output_path = os.path.join(output_dir, f"{interaction_name}.csv")
    contact_atom_df.to_csv(output_path, index=False)
    
    print(f">>> Results saved to \"{output_path}\"")
    return contact_atom_df

def find_possible_hbonds(
        domain1: struc.AtomArray, domain2: struc.AtomArray, cutoff=3.5,
        interaction_name="Hydrogen bonds for Domain I and Domain II",
        output_dir="."
    ):
    
    domain1 = domain1[domain1.get_annotation("element") != "H"]
    domain1_NC_filter = list(map(lambda x, y: x or y, domain1.get_annotation("element") == "O", domain1.get_annotation("element") == "N"))
    domain1 = domain1[domain1_NC_filter]
    domain2 = domain2[domain2.get_annotation("element") != "H"]
    domain2_NC_filter = list(map(lambda x, y: x or y, domain2.get_annotation("element") == "O", domain2.get_annotation("element") == "N"))
    domain2 = domain2[domain2_NC_filter]
    
    return find_interactions_between_domains(domain1, domain2, interaction_name=interaction_name, cutoff=cutoff, output_dir=output_dir)
    

def analyze_specificity_with_contacting_atoms(
        contacting_atom_lists: list[str],
        output_file: str = "./Summary.csv"
    ):
    '''
    Analyzing the specificity of the contacting atoms in the contacting atom list files\n
    - contacting_atom_lists: list of file paths to the contacting atom list files\n
        - a contacting atom list file is a csv file with columns: Atom1, Atom2, Distance (Å)\n
    '''
    basename_list = [".".join(os.path.basename(file_path).split(".")[:-1]) for file_path in contacting_atom_lists] 

    res_dict_temp = {}
    for file_path, basename in zip(contacting_atom_lists, basename_list):
        temp = pd.read_csv(file_path)
        for i, row in temp.iterrows():
            atom1_macro = macro.read_macro(row['Atom1'])
            atom2_macro = macro.read_macro(row['Atom2'])
            contact_tuple = (
                f"//{atom1_macro['segi']}/{atom1_macro['chain']}/{atom1_macro['resn']}`{atom1_macro['resi']}/{atom1_macro['name']}",
                f"//{atom2_macro['segi']}/{atom2_macro['chain']}/{atom2_macro['resn']}`{atom2_macro['resi']}/{atom2_macro['name']}"
            )
            if not contact_tuple in res_dict_temp:
                res_dict_temp[contact_tuple] = {}
            res_dict_temp[contact_tuple][basename] = row['Distance (Å)']

    res_dict = collections.OrderedDict(sorted(res_dict_temp.items(), key=lambda x: macro.atom_priority(x[0][0])))
    res_df = pd.DataFrame({
        "Atom1": [i[0] for i in res_dict],
        "Atom2": [i[1] for i in res_dict]
    })
    for basename in basename_list:
        res_df[basename] = np.zeros(res_df.shape[0])
    for i, my_key in enumerate(res_dict):
        for basename in res_dict[my_key]:
            res_df.loc[i, basename] = res_dict[my_key][basename]

    # 用同一相互作用在不同模型中的出现次数作为 Specificity Index
    res_df["Specificity Index"] = np.zeros(res_df.shape[0])
    for i, row in res_df.iterrows():
        res_df.at[i, "Specificity Index"] = sum(row[basename_list] > 0)
        
    
    ordering = list(res_df.columns)
    for col in ["Atom1", "Atom2", "Specificity Index"]:
        ordering.remove(col)
    ordering = ["Atom1", "Atom2", "Specificity Index"] + ordering
    res_df = res_df[ordering]
    
    res_df.to_csv(output_file, index=False)
    return res_df