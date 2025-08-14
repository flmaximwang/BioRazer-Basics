# base environment is enough
import os
import re
import pandas as pd
import subprocess

def generate_design_stats_path(top_dir):
    return os.path.join(top_dir, "Designs", "Results", "final_design_stats.csv")

def generate_prediction_dir(top_dir, identifier):
    return os.path.join(top_dir, identifier)

def merge_design_stats(top_dir):
    design_stats = None
    des_res_dir = os.path.join(top_dir, "Designs", "Results")
    for my_dir in os.listdir(des_res_dir):
        if not os.path.exists(os.path.join(des_res_dir, my_dir, "final_design_stats.csv")):
            continue
        current_design_stats_path = os.path.join(des_res_dir, my_dir, "final_design_stats.csv")
        if not isinstance(design_stats, pd.DataFrame):
            design_stats = pd.read_csv(current_design_stats_path)
        else:
            design_stats = pd.concat([design_stats, pd.read_csv(current_design_stats_path)], ignore_index=True)
    design_stats.sort_values(by="Average_i_pTM", ascending=False, inplace=True, ignore_index=True)
    design_stats["Rank"] = range(1, len(design_stats) + 1)
    design_stats_path = generate_design_stats_path(top_dir)
    design_stats.to_csv(design_stats_path, index=False)
    return design_stats

def read_design_stats(top_dir):
    """
    The function reads design statistics from a CSV file if it exists, otherwise it merges the design
    statistics.
    
    :param top_dir: The `top_dir` parameter in the `read_design_stats` function is typically a string
    representing the top-level directory where the design statistics are stored or where the function
    should look for the design statistics file. This directory path is used to generate the full path to
    the design statistics file and then check if
    :return: If the `design_stats_path` does not exist, the function will return the result of the
    `merge_design_stats()` function. Otherwise, it will return a pandas DataFrame read from the CSV file
    located at the `design_stats_path`.
    """
    design_stats_path = generate_design_stats_path(top_dir)
    if not os.path.exists(design_stats_path):
        return merge_design_stats()
    return pd.read_csv(design_stats_path)

def format_design_name(design: str):
    if re.match(r".*_model[\d]+", design):
        design = "_".join(design.split("_")[:-1])
    return design.lower()

def get_design_stat(top_dir, design, stat):
    design = format_design_name(design)
    design_stats = read_design_stats(top_dir)
    return design_stats.loc[design_stats["Design"].apply(format_design_name) == design, :][stat].values[0]

def get_design_rank(top_dir, design):
    return get_design_stat(top_dir, design, "Rank")

def search_design_prediction_dir(top_dir, identifier, design):
    prediction_dir = generate_prediction_dir(top_dir, identifier)
    design = format_design_name(design)
    for my_dir in os.listdir(prediction_dir):
        dir_pattern = f"\\d*_?{design}"
        # print(dir_pattern)
        # print(my_dir)
        # print(my_dir.lower())
        if re.match(dir_pattern, my_dir.lower()):
            return os.path.join(prediction_dir, my_dir)

def remove_prediction_rank(top_dir, identifier, design: str, verbose=False):
    design = format_design_name(design)
    design_dir = search_design_prediction_dir(top_dir, identifier, design)
    if not design_dir:
        if verbose: print(">>> No prediction found for: ", design, identifier)
        return None
    
    def rename_dir_or_file(dir_or_file_path):
        # print(dir_or_file_path)
        match_res = re.match(f"(fold_)?\\d+_{design}(.*)", os.path.basename(dir_or_file_path))
        if match_res:
            prefix = match_res.group(1)
            suffix = match_res.group(2)
            old_path = dir_or_file_path
            new_path = os.path.join(os.path.dirname(dir_or_file_path), f"{prefix}{design}{suffix}")
            os.rename(old_path, new_path)
            if verbose: print(f">>> Renamed {old_path} to {new_path}")
        else:
            if verbose: print(f">>> Skipping {dir_or_file_path}")
    
    for my_dir, sub_dirs, filenames in os.walk(design_dir):
        for sub_dir in sub_dirs:
            rename_dir_or_file(os.path.join(my_dir, sub_dir))
        for filename in filenames:
            rename_dir_or_file(os.path.join(my_dir, filename))
    
    rename_dir_or_file(design_dir)

def cif_to_pdb(cif_file_path, verbose=False):
    pdb_file_path = cif_file_path.replace(".cif", ".pdb")
    command = f"obabel {cif_file_path} -O {pdb_file_path}"
    if not os.path.exists(pdb_file_path):
        if verbose: print(command)
        subprocess.run(["/bin/zsh", "-c", command], check=True)
    else:
        if verbose: print(f">>> Skipping because {pdb_file_path} already exists")
    
def convert_design_prediction_cif_to_pdb(top_dir, identifier, design, verbose=False):
    for i in range(5):
        design_prediction_dir = search_design_prediction_dir(top_dir, identifier, design)
        if not design_prediction_dir:
            if verbose: print(f">>> No prediction found for: {design}, {identifier}")
            continue
        cif_file_path = os.path.join(design_prediction_dir, f"fold_{design}_model_{i}.cif")
        if not os.path.exists(cif_file_path):
            continue
        cif_to_pdb(cif_file_path, verbose=verbose)