import pandas as pd

def merge_multiple_af3_res(
    af3_res: pd.DataFrame, dup_num: int,
    marker_prefix: str
):
    af3_dups: list[pd.DataFrame] = []
    af3_res = af3_res.loc[:, ~af3_res.columns.str.startswith('Average')]
    for key in ["Length", "Seed", "MPNN", "AF3_Dup"]:
        assert key in af3_res.columns
    for i in range(dup_num):
        af3_dups.append(af3_res[af3_res["AF3_Dup"] == i].copy())
    cols_to_rename = list(af3_res)
    for col in ["Length", "Seed", "MPNN", "AF3_Dup", "Marker"]:
        cols_to_rename.remove(col)
    for i in range(dup_num):
        af3_dups[i].rename(
            columns={
                f"{col}": f"{int(col[0])+5*i}_{col[2:]}" for col in cols_to_rename
            },
            inplace=True
        )
        af3_dups[i].drop(['Marker', 'AF3_Dup'], axis=1, inplace=True)
    f_res = af3_dups[0]
    for i in range(1, dup_num):
        f_res = f_res.merge(af3_dups[i], on=["Length", "Seed", "MPNN"])
    f_cols = ["Marker"]
    f_res['Marker'] = marker_prefix + "_" + f_res['Length'].astype(str) + "_" + f_res['Seed'].astype(str) + "_" + f_res['MPNN'].astype(str)
    base_cols = set()
    for col in cols_to_rename:
        base_cols.add(col[2:])
    base_cols = list(base_cols)
    for base_col in base_cols:
        f_res[f"Average_{base_col}"] = f_res[[f"{i}_{base_col}" for i in range(5*dup_num)]].mean(axis=1)
        f_cols.append(f"Average_{base_col}")
    for base_col in base_cols:
        for i in range(5*dup_num):
            f_cols.append(f"{i}_{base_col}")
    f_cols.extend(["Length", "Seed", "MPNN"])
    f_res = f_res[f_cols]
    f_res.sort_values(by = ["Length", "Seed", "MPNN"], inplace=True)
    
    return f_res
