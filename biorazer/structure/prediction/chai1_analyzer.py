import os, re, shutil
import numpy as np
import pandas as pd
from pathlib import Path
from biorazer.structure.scripts.struc_file_converter import convert_cif_to_pdb
from .prediction_analyzer import SinglePrediction, FormattedSinglePrediction
import pickle

class SingleChai1(SinglePrediction):

    def format_output(self, target_dir):
        """
        Format the output of the Chai1 prediction.
        This function should return a FormattedChai1 object.
        """
        target_marker_dir = super().format_output(target_dir)

        for i, cif_model_src in enumerate(self.my_dir.glob(f"seed-*/pred.model_idx_*.cif")):
            pattern = f"{self.my_dir}/seed-(\d+)/pred\.model_idx_(\d+)\.cif"
            seed, sample = re.match(pattern, str(cif_model_src)).groups()
            target_sample_dir = target_marker_dir / f"seed-{seed}_sample-{sample}"
            if not target_sample_dir.exists():
                target_sample_dir.mkdir(parents=True)

            cif_model_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_model.cif")
            shutil.copyfile(cif_model_src, cif_model_target)
            pdb_model_target = cif_model_target.with_suffix(".pdb")
            convert_cif_to_pdb(cif_model_src, pdb_model_target)

            scores_npz = cif_model_src.parent / f"scores.model_idx_{sample}.npz"
            scores_data = np.load(scores_npz, allow_pickle=False)
            meta_data = {}
            for key in scores_data:
                meta_data[key] = scores_data[key][0]
            meta_data_pickle_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_data.pickle")
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return FormattedChai1(target_marker_dir)

class FormattedChai1(FormattedSinglePrediction):

    pass

def merge_chai1_trunks(pred_dir, output_dir=None, rank_by="ptm", seed=0):
    """
    Merge the trunks of the Chai1 protein structure.

    Args:
        pred_dir (str): Directory containing the prediction files.
        rank (str): The ranking method to use for merging.

    Returns:
        None
    """
    
    # Ensure the prediction directory exists
    if not isinstance(pred_dir, Path):
        pred_dir = Path(pred_dir)
    if not pred_dir.exists():
        raise FileNotFoundError(f"Directory {pred_dir} does not exist.")
    pred_dir: Path
    
    if output_dir is None:
        output_dir = pred_dir.parent / f"{pred_dir.name}_merged"
    elif not isinstance(output_dir, Path):
        output_dir = Path(output_dir)
    if not output_dir.exists():
        os.makedirs(output_dir)
    output_dir: Path
    
    trunk_indices = []
    for trunk_dir in pred_dir.iterdir():
        trunk_index = re.match(r"trunk_(\d+)", trunk_dir.name).group(1)
        trunk_indices.append(trunk_index)
    trunk_indices = sorted(set(trunk_indices))
    trunk_indices = np.array(trunk_indices, dtype=int)
    trunk_seeds = trunk_indices + seed
    
    results = pd.DataFrame({
        "trunk_index": [],
        "trunk_seed": [],
        "model_index": [],
        "score": [],
    })

    for trunk_index, seed in zip(trunk_indices, trunk_seeds):
        trunk_dir = pred_dir / f"trunk_{trunk_index}"
        
        score_files = sorted(list(trunk_dir.glob(f"scores.model_idx_*.npz")), key=lambda x: int(x.name.split("_")[-1].split(".")[0]))
        scores = []
        for score_file in score_files:
            scores.append(np.load(score_file)[rank_by][0])

        results = pd.concat([results, pd.DataFrame({
            "trunk_index": [trunk_index] * len(scores),
            "trunk_seed": [seed + trunk_index] * len(scores),
            "model_index": list(range(len(scores))),
            "score": scores,
        })], ignore_index=True)
    
    # Sort the results by score
    results = results.sort_values(
        by=["score", "trunk_seed"],
        ascending=[False, True],
        ignore_index=True
    )
    
    # Copy files to the output directory
    for i, row in results.iterrows():
        trunk_index = int(row["trunk_index"])
        trunk_seed = int(row["trunk_seed"])
        model_index = int(row["model_index"])
        trunk_dir = pred_dir / f"trunk_{trunk_index}"
        model_file = trunk_dir / f"pred.model_idx_{model_index}.cif"
        score_file = trunk_dir / f"scores.model_idx_{model_index}.npz"
        shutil.copy(
            model_file,
            output_dir / f"pred.rank_{i+1:0>3d}_seed_{trunk_seed}_trunk_{trunk_index}_model_idx_{model_index}.cif"
        )
        shutil.copy(
            score_file,
            output_dir / f"scores.rank_{i+1:0>3d}_seed_{trunk_seed}_trunk_{trunk_index}_model_idx_{model_index}.npz"
        )
        