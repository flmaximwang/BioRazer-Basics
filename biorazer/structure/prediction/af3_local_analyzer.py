import pandas as pd
from pathlib import Path
import re, json, shutil, pickle
from .prediction_analyzer import SinglePrediction, FormattedSinglePrediction
from biorazer.structure.scripts.struc_file_converter import convert_cif_to_pdb

class SingleAF3Local(SinglePrediction):
    """
    This class is used to analyze the AF3 local prediction results.
    """
    
    def format_output(self, target_dir):

        target_marker_dir = super().format_output(target_dir)
        
        for i, sample in enumerate(self.my_dir.glob("seed-*_sample-*")):

            if not sample.is_dir():
                continue
            target_sample_dir = target_marker_dir / sample.stem
            if not target_sample_dir.exists():
                target_sample_dir.mkdir(parents=True)
            
            cif_model_src = sample / (self.my_dir.stem + "_" + sample.stem + "_model.cif")
            cif_model_target = target_sample_dir / (self.my_dir.stem + "_" + sample.stem + "_model.cif")
            shutil.copyfile(cif_model_src, cif_model_target)

            pdb_model_target = cif_model_target.with_suffix(".pdb")
            convert_cif_to_pdb(cif_model_src, pdb_model_target)

            confidences_json = sample / (self.my_dir.stem + "_" + sample.stem + "_confidences.json")
            confidence_data = json.load(open(confidences_json, "r"))
            summary_confidences_json = sample / (self.my_dir.stem + "_" + sample.stem + "_summary_confidences.json")
            summary_confidences_data = json.load(open(summary_confidences_json, "r"))

            meta_data = {}
            for key in summary_confidences_data:
                meta_data[key] = summary_confidences_data[key]
            for key in confidence_data:
                meta_data[key] = confidence_data[key]
            meta_data_pickle_target = target_sample_dir / (self.my_dir.stem + "_" + sample.stem + "_data.pickle")
            with open(meta_data_pickle_target, "wb") as f:
                pickle.dump(meta_data, f)
        
        return FormattedAF3Local(target_marker_dir)

class FormattedAF3Local(FormattedSinglePrediction):

    pass


def extract_af3_confidence_data(af3_dir, metric_key_list: list[str], metric_id_list: list[tuple], metric_label_list: list[str]):
    af3_dir = Path(af3_dir)
    marker = af3_dir.name
    predicted_samples = sorted(list(af3_dir.glob("seed-*_sample-*")))
    
    columns = ["Marker"]
    for metric_label in metric_label_list:
        for i in range(len(predicted_samples)):
            columns.append(f"{i+1}_{metric_label}")
    res = pd.DataFrame(columns=columns)
    res.loc[0, "Marker"] = marker

    for i, predicted_sample in enumerate(predicted_samples):
        if not predicted_sample.is_dir():
            continue
        seed, sample = re.match(r"seed-(\d+)_sample-(\d+)", predicted_sample.stem).groups()
        seed = int(seed)
        sample = int(sample)
        summary_confidences_file = predicted_sample / (af3_dir.stem + "_" + predicted_sample.stem + "_summary_confidences.json")
        with open(summary_confidences_file, "r") as f:
            summary_confidences_data = json.load(f)
        for metric_key, metric_id, metric_label in zip(metric_key_list, metric_id_list, metric_label_list):
            metric_data = summary_confidences_data[metric_key]
            if not isinstance(metric_data, list):
                data = metric_data
            else:
                data = None
                for j in metric_id:
                    if data is None:
                        data = metric_data[j]
                    else:
                        data = data[j]
            res.loc[0, f"{i+1}_{metric_label}"] = data
    return res