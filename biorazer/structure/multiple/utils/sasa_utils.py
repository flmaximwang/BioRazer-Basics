import os
import csv
import subprocess
from .af_utils import *

def calculate_pdb_sasa(pdb_file_path, verbose=True):
    sasa_file_path = pdb_file_path.replace(".pdb", ".sasa")
    command = f"freesasa {pdb_file_path} 1>{sasa_file_path}"
    if verbose: print(command)
    subprocess.run(["/bin/zsh", "-i", "-c", command])

def freesasa_raw_to_csv(raw_file_path, csv_file_path):
    with open(raw_file_path, 'r') as raw_file:
        lines = raw_file.readlines()
    
    results = {}
    for line in lines:
        if line.startswith("Total"):
            results["Total"] = float(line.split(":")[1].strip())
        elif line.startswith("Apolar"):
            results["Apolar"] = float(line.split(":")[1].strip())
        elif line.startswith("Polar"):
            results["Polar"] = float(line.split(":")[1].strip())
        elif line.startswith("CHAIN"):
            chain_info = line.split(":")
            chain_name = chain_info[0].strip()
            chain_value = float(chain_info[1].strip())
            results[chain_name] = chain_value
    
    with open(csv_file_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["", "Value"])
        for key, value in results.items():
            writer.writerow([key, value])

def run_sasa_for_prediction(top_dir, identifier, design, verbose=False):
    af_file_prefix = generate_af_file_prefix(top_dir, identifier, design)
    for i in range(5):
        pdb_file_path = f"{af_file_prefix}_model_{i}.pdb"
        if not os.path.exists(pdb_file_path):
            if verbose: print(f">>> No prediction found for {design}.")
            continue
        raw_file_path = pdb_file_path.replace(".pdb", ".sasa")
        csv_file_path = pdb_file_path.replace(".pdb", "_sasa.csv")
        if not os.path.exists(csv_file_path):
            calculate_pdb_sasa(pdb_file_path, verbose)
            freesasa_raw_to_csv(raw_file_path, csv_file_path)
        else:
            if verbose: print(f">>> Skipping because {csv_file_path} already exists.")

def get_hydrophobicity(top_dir, identifier, design, prediction_i):
    af_file_prefix = generate_af_file_prefix(top_dir, identifier, design)
    sasa_file_path = f"{af_file_prefix}_model_{prediction_i}_sasa.csv"
    with open(sasa_file_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            if row[0] == "Apolar":
                hydrophobic_sasa = float(row[1])
            elif row[0] == "Total":
                total_sasa = float(row[1])
    return hydrophobic_sasa / total_sasa