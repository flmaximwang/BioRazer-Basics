import pandas as pd
import numpy as np
from pathlib import Path
import re, json, shutil, pickle
from .prediction_analyzer import SinglePrediction, FormattedSinglePrediction
from biorazer.structure.scripts.struc_file_converter import convert_cif_to_pdb

class SingleBoltz1(SinglePrediction):

    def format_output(self, target_dir):
        """
        Format the output of the Boltz1 prediction.
        This function should return a FormattedBoltz1 object.
        """
        target_marker_dir = super().format_output(target_dir)

        for i, cif_model_src in enumerate(self.my_dir.glob(f"seed-*/predictions/{self.marker}/{self.marker}_model_*.cif")):

            pattern = f"{self.my_dir}/seed-(\d+)/predictions/{self.marker}/{self.marker}_model_(\d+)\.cif"
            seed, sample = re.match(pattern, str(cif_model_src)).groups()
            target_sample_dir = target_marker_dir / f"seed-{seed}_sample-{sample}"
            if not target_sample_dir.exists():
                target_sample_dir.mkdir(parents=True)

            cif_model_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_model.cif")
            shutil.copyfile(cif_model_src, cif_model_target)
            pdb_model_target = cif_model_target.with_suffix(".pdb")
            convert_cif_to_pdb(cif_model_src, pdb_model_target)

            confidence_json = cif_model_src.parent / f"confidence_{cif_model_src.stem}.json"
            confidence_data = json.load(open(confidence_json, "r"))
            plddt_npz = cif_model_src.parent / f"plddt_{cif_model_src.stem}.npz"
            plddt_data = np.load(plddt_npz, allow_pickle=False)
            meta_data = {}
            for key in confidence_data:
                meta_data[key] = confidence_data[key]
            for key in plddt_data:
                meta_data[key] = plddt_data[key]
            meta_data_pickle_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_data.pickle")
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return FormattedBoltz1(target_marker_dir)

class FormattedBoltz1(FormattedSinglePrediction):

    pass