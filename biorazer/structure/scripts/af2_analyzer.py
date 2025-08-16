import json, os, re
import pandas as pd
import numpy as np

def extract_af2_meta(pred_dir):
    marker = os.path.basename(pred_dir)
    counter = 1
    all_files = os.listdir(pred_dir)
    full_marker_list = []
    while True:
        found_model = False
        for f in all_files:
            f: str
            if f.startswith(f"{marker}_scores_rank_{counter:03d}"):
                found_model = True
                full_marker_list.append(f.split("_scores_rank_"))
                break
        if not found_model:
            break
        counter += 1
    model_num = counter - 1
    
    base_cols = ["pLDDT", "pAE", "pTM", "i_pTM", "ranking_score"]
    sort_col = ["Marker"]
    for key in base_cols:
        sort_col.append(f"Average_{key}")
        for i in range(1, model_num + 1):
            sort_col.append(f"{i}_{key}")
    res = pd.DataFrame(columns=sort_col)
    res.loc[0, "Marker"] = marker
    
    for i in range(1, model_num + 1):
        f = f"{full_marker_list[i-1][0]}_scores_rank_{full_marker_list[i-1][1]}"
        with open(os.path.join(pred_dir, f), "r") as file:
            try:
                meta_data = json.load(file)
            except:
                print(f"There's some problem with {f}")
            full_plddt = meta_data['plddt']
            plddt = np.mean(full_plddt)
            full_pae = np.array(meta_data['pae']) / 31.0
            pae = np.mean(full_pae)
            ptm = meta_data['ptm']
            if 'iptm' in meta_data:
                # 当预测的链数为 1 时, meta_data 中不会有 iptm
                iptm = meta_data['iptm']
            else:
                iptm = 0
        res.loc[0, f"{i}_pLDDT"] = plddt
        res.loc[0, f"{i}_pAE"] = pae
        res.loc[0, f"{i}_pTM"] = ptm
        res.loc[0, f"{i}_i_pTM"] = iptm
        if iptm > 0:
            res.loc[0, f"{i}_ranking_score"] = 0.8 * iptm + 0.2 * ptm
        else:
            res.loc[0, f"{i}_ranking_score"] = ptm
    
    for key in base_cols:
        res[f"Average_{key}"] = res[[f"{i}_{key}" for i in range(1, model_num+1)]].mean(axis=1)
    return res

def batch_extract_af2_meta(dir_of_af3_dir, extract_func = extract_af2_meta, exclude=[".DS_Store"], extract_func_kwargs={}):
    counter = 0
    for af2_dir in os.listdir(dir_of_af3_dir):
        if not os.path.isdir(os.path.join(dir_of_af3_dir, af2_dir)):
            continue
        if af2_dir in exclude:
            continue
        if counter == 0:
            res = extract_func(os.path.join(dir_of_af3_dir, af2_dir), **extract_func_kwargs)
        else:
            new_line = extract_func(os.path.join(dir_of_af3_dir, af2_dir), **extract_func_kwargs)
            res.loc[counter] = new_line.loc[0]
        counter += 1
    
    return res


def collect_AF2_results(res_dir: str):
    job_files = {}
    for marker in os.listdir(res_dir):
        if os.path.isdir(os.path.join(res_dir, marker)):
            continue
        matched = re.match(f"(.*)\.a3m", marker)
        if matched:
            job_name = matched.group(1)
            job_files[job_name] = []
    print(">>> Found {} jobs".format(len(job_files)))
    for filename in os.listdir(res_dir):
        for job_name in job_files:
            flag = False
            if re.match(f"{job_name}_coverage\.png", filename):
                flag = True
            elif re.match(f"{job_name}_pae\.png", filename):
                flag = True
            elif re.match(f"{job_name}_plddt\.png", filename):
                flag = True
            elif re.match(f"{job_name}_predicted_aligned_error_.*\.json", filename):
                flag = True
            elif re.match(f"{job_name}_scores_.*\.json", filename):
                flag = True
            elif re.match(f"{job_name}_unrelaxed_.*\.pdb", filename):
                flag = True
            elif re.match(f"{job_name}_relaxed_.*\.pdb", filename):
                flag = True
            elif re.match(f"{job_name}\.a3m", filename):
                flag = True
            elif re.match(f"{job_name}\.done\.txt", filename):
                flag = True
            else:
                pass
            if flag:
                job_files[job_name].append(filename)
                break
    for job_name in job_files:
        os.makedirs(os.path.join(res_dir, job_name))
        for filename in job_files[job_name]:
            os.rename(
                os.path.join(res_dir, filename),
                os.path.join(res_dir, job_name, filename)
            )