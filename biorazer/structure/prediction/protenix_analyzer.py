import os, re, shutil, json
import numpy as np
from pathlib import Path
from biorazer.structure.scripts.struc_file_converter import convert_cif_to_pdb
from .prediction_analyzer import SinglePrediction, FormattedSinglePrediction
import pickle

class SingleProtenix(SinglePrediction):

    def format_output(self, target_dir):
        """
        Format the output of the Chai1 prediction.
        This function should return a FormattedChai1 object.
        """
        target_marker_dir = super().format_output(target_dir)

        for i, cif_model_src in enumerate(self.my_dir.glob(f"seed_*/predictions/{self.marker}_seed_*_sample_*.cif")):
            pattern = f"{self.marker}_seed_(\d+)_sample_(\d+).cif"
            seed, sample = re.match(pattern, cif_model_src.name).groups()
            target_sample_dir = target_marker_dir / f"seed-{seed}_sample-{sample}"
            if not target_sample_dir.exists():
                target_sample_dir.mkdir(parents=True)

            cif_model_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_model.cif")
            shutil.copyfile(cif_model_src, cif_model_target)
            pdb_model_target = cif_model_target.with_suffix(".pdb")
            convert_cif_to_pdb(cif_model_src, pdb_model_target)

            summary_confidences_json = cif_model_src.parent / f"{self.marker}_seed_{seed}_summary_confidence_sample_{sample}.json"
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))
            meta_data = {}
            for key in summary_confidences_data:
                meta_data[key] = summary_confidences_data[key]
            meta_data_pickle_target = target_sample_dir / (self.my_dir.stem + "_" + f"seed-{seed}" + "_" + f"sample-{sample}" + "_data.pickle")
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)

        return FormattedProtenix(target_marker_dir)

class FormattedProtenix(FormattedSinglePrediction):

    pass