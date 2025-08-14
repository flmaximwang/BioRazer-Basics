import os, json, re, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, Selection
import pymol
from pathlib import Path

def rename_af3_pre_dir(pre_dir, new_name):
    dir_of_pre_dir = os.path.dirname(pre_dir)
    old_name = os.path.basename(pre_dir)
    for filename in os.listdir(pre_dir):
        if old_name in filename:
            os.rename(
                os.path.join(pre_dir, filename),
                os.path.join(pre_dir, filename.replace(old_name, new_name))
            )
    os.rename(pre_dir, os.path.join(dir_of_pre_dir, new_name))

def collect_af3_model_0(af3_pred_batch_dir, target_dir=None):
    '''
    收集 AF3 预测的第一个模型, 也就是 rank_001 的模型,
    '''
    af3_pred_batch_dir = Path(af3_pred_batch_dir)
    if target_dir is None:
        target_dir = af3_pred_batch_dir.parent / f"{af3_pred_batch_dir.name}_model_0"
    else:
        target_dir = Path(target_dir)
    if not target_dir.exists():
        target_dir.mkdir(parents=True, exist_ok=True)
    for pred_dir in af3_pred_batch_dir.iterdir():
        if not pred_dir.is_dir():
            continue
        for model_file in pred_dir.iterdir():
            if re.match(r"fold_(.*)_model_rank_001_seed_\d+\.cif", model_file.name):
                shutil.copyfile(
                    model_file,
                    target_dir / model_file.name
                )
                break

def collect_best_models(af3_res: pd.DataFrame, dir_of_preds: str, marker_col="Marker", struc_file_fmt=".cif", target_dir=None, map_col_suffix="map"):
    if target_dir is None:
        target_dir = dir_of_preds + "_best"
    if not os.path.exists(target_dir):
        os.makedirs(target_dir, exist_ok=True)
    if not f"1_{map_col_suffix}" in af3_res.columns:
        af3_res[f"1_{map_col_suffix}"] = 1
    for rank, row_i in enumerate(af3_res.index, start=1):
        marker = af3_res.loc[row_i, marker_col]
        best_model = af3_res.loc[row_i, f"1_{map_col_suffix}"]
        model_file_pattern = f"fold_{marker}_model_rank_{int(best_model):0>3}_seed_\d+\{struc_file_fmt}"
        # print(model_file_pattern)
        for filename in os.listdir(os.path.join(dir_of_preds, marker)):
            if re.match(model_file_pattern, filename):
                model_file = os.path.join(dir_of_preds, marker, filename)
                break
        shutil.copyfile(
            model_file,
            os.path.join(target_dir, f"rank_{rank:0>3}_{marker}_best_model{struc_file_fmt}")
        )

def merge_af3_duplicates(dir_of_preds, pred_prefix, target_dir = None, dup_num: int=1):
    '''
    使用多个 Seed 对同一个 protein 进行多次预测, 你可以使用这个函数将这些预测合并在一起, 输出到一个新的目录中.
    这些预测的目录名应当以 pred_prefix 开头, 例如 "pred_1", "pred_2", ...
    '''
    if target_dir is None:
        target_dir = dir_of_preds
    os.makedirs(os.path.join(target_dir, pred_prefix), exist_ok=True)
    af3_res = []
    for dirname in os.listdir(dir_of_preds):
        matched = re.match(f"{pred_prefix}_(\d+)", dirname)
        if matched and int(matched.group(1)) in range(dup_num):
            with open(os.path.join(dir_of_preds, dirname, f"fold_{dirname}_job_request.json")) as f:
                job_meta = json.load(f)
            seed = job_meta[0]["modelSeeds"][0]
            for i in range(5):
                with open(os.path.join(dir_of_preds, dirname, f"fold_{dirname}_summary_confidences_{i}.json")) as f:
                    confidence_data  = json.load(f)
                model_ranking_score = confidence_data["ranking_score"]
                af3_res.append([dirname, i, seed, model_ranking_score])
            old_name = f"fold_{dirname}_job_request.json"
            assert os.path.exists(os.path.join(dir_of_preds, dirname, old_name))
            new_name = f"fold_{pred_prefix}_seed_{seed:0>10}_job_request.json"
            shutil.copyfile(
                os.path.join(dir_of_preds, dirname, old_name),
                os.path.join(target_dir, pred_prefix, new_name)
            )
    
    af3_res.sort(key=lambda x: x[-1], reverse=True)
    for i, row in enumerate(af3_res, start=1):
        dirname, model_i, seed, ranking_score = row
        # print(ranking_score)
        for insertion, ext_name in zip(["full_data", "summary_confidences", "model"], ["json", "json", "cif"]):
            old_name = f"fold_{dirname}_{insertion}_{model_i}.{ext_name}"
            assert os.path.exists(os.path.join(dir_of_preds, dirname, old_name))
            new_name = f"fold_{pred_prefix}_{insertion}_rank_{i:0>3}_seed_{seed:0>10}.{ext_name}"
            shutil.copyfile(
                os.path.join(dir_of_preds, dirname, old_name),
                os.path.join(target_dir, pred_prefix, new_name)
            )

def batch_merge_af3_duplicates(dir_of_preds, target_dir = None, dup_num: int=1):
    if target_dir is None:
        target_dir = dir_of_preds + "_merged"
    os.makedirs(target_dir, exist_ok=True)
    pred_prefixes = set()
    for dirname in os.listdir(dir_of_preds):
        matched = re.match(r"(.*)_(\d+)", dirname)
        if matched:
            pred_prefixes.add(matched.group(1))
    for pred_prefix in pred_prefixes:
        # print(pred_prefix)
        merge_af3_duplicates(dir_of_preds, pred_prefix, target_dir, dup_num=dup_num)

def get_files_with_pattern(af3_dir, insertion_str: str):
    '''
    获取指定目录下, 所有符合 pattern 的文件 (一个预测)
    '''

    af3_dir = Path(af3_dir)
    marker = af3_dir.name
    files_with_pattern = []

    counter = 1
    while True:
        rank_exist = False
        for filename in af3_dir.iterdir():
            if re.match(f"fold_{marker}_{insertion_str}_rank_{counter:0>3}_seed_\d+.*", filename.name):
                files_with_pattern.append(filename)
                rank_exist = True
                break
        if not rank_exist:
            break
        counter += 1
    
    return files_with_pattern

def get_files_of_summary_confidences(af3_dir):
    return get_files_with_pattern(af3_dir, "summary_confidences")

def get_files_of_full_data(af3_dir):
    return get_files_with_pattern(af3_dir, "full_data")

def get_files_of_model(af3_dir):
    return get_files_with_pattern(af3_dir, "model")

def extract_af3_plddt_mean(
    af3_dir, 
    chain_id: list[str],
    metric_label_list: list[str]
):
    af3_dir = Path(af3_dir)
    marker = af3_dir.name
    columns = ["Marker"]
    full_data_files = []

    counter = 1
    while True:
        full_data_pattern = f"fold_{marker}_full_data_rank_{counter:0>3}_seed_\d+\.json"
        rank_exist = False
        for filename in af3_dir.iterdir():
            if re.match(full_data_pattern, filename.name):
                full_data_files.append(filename)
                rank_exist = True
                break
        if not rank_exist:
            break
        counter += 1
    
    for metric_label in metric_label_list:
        for i in range(1, len(full_data_files) + 1):
            columns.append(f"{i}_{metric_label}")
    
    res = pd.DataFrame(columns=columns)
    res.loc[0, "Marker"] = marker
    for i, confidence_data_file in enumerate(full_data_files):
        with open(confidence_data_file) as f:
            full_data = json.load(f)
        for metric_key, subset_idx, metric_label in zip(metric_key_list, subset_idx_list, metric_label_list):
            metric_data = np.array(full_data[metric_key])
            data = np.mean(metric_data[subset_idx])
            res.loc[0, f"{i+1}_{metric_label}"] = data
    
    return res

def extract_af3_confidence_data(af3_dir, metric_key_list: list[str], metric_id_list: list[tuple], metric_label_list: list[str]):
    
    af3_dir = Path(af3_dir)
    confidence_data_files = get_files_of_summary_confidences(af3_dir)
    marker = af3_dir.name
    columns = ["Marker"]

    for metric_label in metric_label_list:
        for i in range(1, len(confidence_data_files) + 1):
            columns.append(f"{i}_{metric_label}")
    
    res = pd.DataFrame(columns=columns)
    res.loc[0, "Marker"] = marker
    for i, confidence_data_file in enumerate(confidence_data_files):
        with open(confidence_data_file) as f:
            confidence_data = json.load(f)
        for metric_key, metric_id, metric_label in zip(metric_key_list, metric_id_list, metric_label_list):
            metric_data = confidence_data[metric_key]
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

def extract_af3_confidence_data_default(af3_dir):

    return extract_af3_confidence_data(
        af3_dir,
        metric_key_list=["fraction_disordered", "has_clash", "iptm", "num_recycles", "ptm", "ranking_score"],
        metric_id_list=[(), (), (), (), (), ()]
    )

def batch_extract_af3_meta(dir_of_af3_dir, extract_func, extract_func_kwargs={}):
    counter = 0
    dir_of_af3_dir = Path(dir_of_af3_dir)
    for af3_dir in dir_of_af3_dir.iterdir():
        if not af3_dir.is_dir():
            continue
        if counter == 0:
            res = extract_func(af3_dir, **extract_func_kwargs)
        else:
            res = pd.concat([res, extract_func(af3_dir, **extract_func_kwargs)], ignore_index=True)
        counter += 1
    
    return res

AF3_CONFIDENCE_EXAMPLE = {
 "chain_iptm": [
  0.75,
  0.75
 ],
 "chain_pair_iptm": [
  [
   0.81,
   0.75
  ],
  [
   0.75,
   0.82
  ]
 ],
 "chain_pair_pae_min": [
  [
   0.76,
   2.13
  ],
  [
   2.07,
   0.76
  ]
 ],
 "chain_ptm": [
  0.81,
  0.82
 ],
 "fraction_disordered": 0.0,
 "has_clash": 0.0,
 "iptm": 0.75,
 "num_recycles": 10.0,
 "ptm": 0.81,
 "ranking_score": 0.76
}

AF3_FULL_DATA_EXAMPLE = {
    "atom_chain_ids": ["A", "A", "B", "B"],
    "atom_plddts": [57.24,59.68,61.59,58.29],
    "contact_probs": [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ],
    "token_chain_ids": ["A", "B"],
    "token_res_ids": [1, 2],
    "pae": [
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ],
}