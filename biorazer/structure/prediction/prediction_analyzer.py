from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, re, pickle, pymol


class SinglePrediction:

    def __init__(self, my_dir):

        self.my_dir = Path(my_dir)
        self.marker = self.my_dir.name

    def format_output(self, target_dir):
        """
        Implement this function to format the output of your prediction.
        This function should return a FormattedSinglePrediction object.
        At the same time, it should create a directory named after the marker in the target_dir,
        with the following structure:
        target_dir/
            marker/
                seed-1_sample-1/
                    my_marker_seed-1_sample-1_model.cif
                    my_marker_seed-1_sample-1_model.pdb
                    my_marker_seed-1_sample-1_data.pickle
                seed-1_sample-2/
                    ...
                seed-2_sample-1/
                    ...
                ...
        """

        target_dir = Path(target_dir)
        target_marker_dir = target_dir / self.marker
        if not target_marker_dir.exists():
            target_marker_dir.mkdir(parents=True)

        return target_marker_dir


class FormattedSinglePrediction(SinglePrediction):

    def __init__(self, my_dir):
        super().__init__(my_dir)
        samples = list(self.my_dir.glob("seed-*_sample-*"))
        self.samples = pd.DataFrame(
            columns=["seed", "sample", "cif_model", "pdb_model", "meta_data"]
        )
        for i, predicted_sample in enumerate(samples):
            seed, sample = re.match(
                r"seed-(\d+)_sample-(\d+)", predicted_sample.stem
            ).groups()
            self.samples.loc[i, "seed"] = int(seed)
            self.samples.loc[i, "sample"] = int(sample)
            self.samples.loc[i, "cif_model"] = predicted_sample / (
                self.my_dir.stem + "_" + predicted_sample.stem + "_model.cif"
            )
            self.samples.loc[i, "pdb_model"] = predicted_sample / (
                self.my_dir.stem + "_" + predicted_sample.stem + "_model.pdb"
            )
            self.samples.loc[i, "meta_data"] = predicted_sample / (
                self.my_dir.stem + "_" + predicted_sample.stem + "_data.pickle"
            )
        self.samples.sort_values(by=["seed", "sample"], inplace=True, ignore_index=True)

    def get_meta_data_keys(self):
        """
        Get the keys of the meta data.
        """
        tmp_meta_data = pickle.load(open(self.samples.loc[0, "meta_data"], "rb"))
        return list(tmp_meta_data.keys())

    def get_meta_data(self, seed, sample, meta_data_key):
        """
        Get the meta data of a specific sample.
        seed: seed of the sample
        sample: sample number
        meta_data_key: key of the meta data
        """
        selector = (self.samples["seed"] == seed) & (self.samples["sample"] == sample)
        selected_sample = self.samples[selector]
        meta_data = pickle.load(open(selected_sample.iloc[0]["meta_data"], "rb"))
        return meta_data[meta_data_key]

    def add_simple_metric(self, metric_name, metric_idx=None, metric_label=None):
        """
        Add a simple metric to the samples dataframe.
        - metric_name: name of the metric
        - metric_idx: index of the metric. for chain_pairs, like (0, 1). for pLDDT, like None
        - metric_label: label of the metric
        """
        if metric_label is None:
            metric_label = metric_name
        for i in self.samples.index:
            meta_data = pickle.load(open(self.samples.loc[i, "meta_data"], "rb"))
            if not metric_name in meta_data:
                raise ValueError(
                    f"{metric_name} not found in {self.samples.loc[i, 'meta_data']}"
                )

            if metric_idx is None:
                self.samples.loc[i, metric_label] = meta_data[metric_name]
            else:
                tmp_data = meta_data[metric_name]
                for idx in metric_idx:
                    tmp_data = tmp_data[idx]
                self.samples.loc[i, metric_label] = tmp_data

    def flatten_samples(self):

        results = pd.DataFrame()
        results.loc[0, "Marker"] = self.marker
        properties = [i for i in self.samples.columns if not i in ["seed", "sample"]]

        for property in properties:
            for i in self.samples.index:
                seed = self.samples.loc[i, "seed"]
                sample = self.samples.loc[i, "sample"]
                results.loc[0, f"{property}_seed-{seed}_sample-{sample}"] = (
                    self.samples.loc[i, property]
                )

        return results


def rerank_models(
    af3_res: pd.DataFrame, model_num: int, rerank_by: str, ascending: bool = False
):
    """
    使用 merge_multiple_af3_res 函数合并多个 AF3 预测后, 你可以使用这个函数对所有的 model 重新进行排序.
    这个函数会在 dataframe 中添加新的列 i_map, 例如 0_map = 1, 1_map =2, 2_map = 3, ...
    它表示现在的 model 0 对应原来的 model 1, model 1 对应 model 2, ...
    - ascending: 是否升序排列 (default: False)
    """
    base_cols = set()
    for col in af3_res.columns:
        matched = re.match(r"(\d+)_(.*)", col)
        if matched:
            base_cols.add(matched.group(2))
    if not rerank_by in base_cols:
        raise ValueError(
            f"rerank_by {rerank_by} not found in columns. Please use one of {base_cols}"
        )

    for row_i in af3_res.index:
        cols_for_reranking = [f"{i+1}_{rerank_by}" for i in range(model_num)]
        new_order = np.argsort(af3_res.loc[row_i, cols_for_reranking].values)
        if not ascending:
            new_order = new_order[::-1]
        for base_col in base_cols:
            related_cols = [f"{i+1}_{base_col}" for i in range(model_num)]
            af3_res.loc[row_i, related_cols] = af3_res.loc[row_i, related_cols].values[
                new_order
            ]
        af3_res.loc[row_i, [f"{i+1}_map" for i in range(model_num)]] = new_order


def plot_rank_differences(
    data_list: list[pd.DataFrame],
    marker_col: str,
    markers: list[str],
    labels: list[str],
    ymin: int,
    ymax: int,
    ylabel="Rank",
    figsize=(10, 10),
    show_diff=False,
    colors=[],
    col_to_show_number="Average_i_pTM",
    col_to_show_number_short="ipTM",
    legend_cols=1,
):
    plt.figure(figsize=figsize)
    df_num = len(data_list)
    if len(colors) == 0:
        colors = [f"C{i}" for i in range(len(markers))]
    elif len(colors) != len(markers):
        raise ValueError(
            "The number of colors should be the same as the number of designs."
        )
    else:
        pass

    for i, marker in enumerate(markers):
        if not show_diff:
            label = (
                marker
                + f": {col_to_show_number_short} "
                + " vs. ".join(
                    [
                        f"{data.loc[data[marker_col] == marker, col_to_show_number].values[0]:.2f}"
                        for data in data_list
                    ]
                )
            )
            plt.plot(
                range(df_num),
                [
                    data.loc[data[marker_col] == marker, :].index[0] + 1
                    for data in data_list
                ],
                "o-",
                label=label,
                color=colors[i],
            )
        else:
            if df_num != 2:
                raise ValueError("Only support two dataframes for now.")
            before = (
                data_list[0].loc[data_list[0][marker_col] == marker, :].index[0] + 1
            )
            after = data_list[1].loc[data_list[1][marker_col] == marker, :].index[0] + 1
            label = (
                marker
                + f": {col_to_show_number_short} "
                + " vs. ".join(
                    [
                        f"{data.loc[data[marker_col] == marker, col_to_show_number].values[0]:.2f}"
                        for data in data_list[0:2]
                    ]
                )
            )
            if after - before < 0:
                plt.plot([0, 1], [before, after], "o-", label=label, color="red")
            else:
                plt.plot([0, 1], [before, after], "o-", label=label, color="green")

    plt.ylim(ymax, 0)
    # 隐藏 x 轴
    plt.xticks(range(df_num), labels)
    # plt.yticks(np.arange(0, ymax+1, 1))
    plt.yticks([])
    plt.ylabel(ylabel)
    plt.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0, ncol=legend_cols
    )


def plot_metric_differences(
    data_list: list[pd.DataFrame],
    marker_col: str,
    markers: list[str],
    labels: list[str],
    ymin: float,
    ymax: float,
    metric_name: str,
    ylabel="Metric Name",
    show_diff=False,
    colors=[],
    figsize=(10, 10),
    metric_name_short=None,
    step=0.025,
    legend_bbox_to_anchor=(1.05, 1),
    legend_cols=1,
):
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    df_num = len(data_list)
    if not metric_name_short:
        metric_name_short = metric_name
    if len(colors) == 0:
        colors = [f"C{i}" for i in range(len(markers))]
    elif len(colors) != len(markers):
        raise ValueError(
            "The number of colors should be the same as the number of designs."
        )
    else:
        pass
    for i, design in enumerate(markers):
        if not show_diff:
            label = (
                design
                + f": {metric_name_short} "
                + " vs. ".join(
                    [
                        f"{data.loc[data[marker_col] == design, metric_name].values[0]:.2f}"
                        for data in data_list
                    ]
                )
            )
            ax.plot(
                range(df_num),
                [
                    data.loc[data[marker_col] == design, metric_name].values[0]
                    for data in data_list
                ],
                "o-",
                label=label,
                color=colors[i],
            )
        else:
            if df_num != 2:
                raise ValueError("Only support two dataframes for now.")
            before = (
                data_list[0]
                .loc[data_list[0][marker_col] == design, metric_name]
                .values[0]
            )
            after = (
                data_list[1]
                .loc[data_list[1][marker_col] == design, metric_name]
                .values[0]
            )
            label = (
                design
                + f": {metric_name_short} "
                + " vs. ".join(
                    [
                        f"{data.loc[data[marker_col] == design, metric_name].values[0]:.2f}"
                        for data in data_list[0:2]
                    ]
                )
            )
            if after - before < 0:
                plt.plot([0, 1], [before, after], "o-", label=label, color="red")
            else:
                plt.plot([0, 1], [before, after], "o-", label=label, color="green")

    ax.set_ylim(ymin, ymax)
    # 隐藏 x 轴
    ax.set_xticks(range(df_num), labels)
    ax.set_yticks(np.arange(ymin, ymax + 0.001 * step / abs(step), step))
    ax.set_ylabel(ylabel)
    ax.legend(
        bbox_to_anchor=legend_bbox_to_anchor,
        loc="upper left",
        borderaxespad=0.0,
        ncol=legend_cols,
    )
    return fig, ax


def plot_metric_distribution(
    pred_stat: pd.DataFrame,
    metric_name: str,
    ax: plt.Axes = None,
    y_shift=0,
    y_fluctuation: float = 0.1,
    subplots_kwargs={},
    label="",
    scatter_kwargs={},
    errorbar_kwargs={"fmt": ",", "alpha": 0.5, "capsize": 3},
):
    if not ax:
        fig, ax = plt.subplots(1, 1, **subplots_kwargs)
    else:
        fig = ax.get_figure()

    metric_name_type = -1
    for column_name in pred_stat.columns:
        if re.match(rf"\d+_{metric_name}", column_name):
            metric_name_type = 0
        elif re.match(rf"{metric_name}_seed-\d+_sample-\d+", column_name):
            metric_name_type = 1
        else:
            continue
        break
    if metric_name_type == -1:
        raise ValueError(
            f"No supported metric name found in {pred_stat.columns}. Please use the format of 1_metric_name or metric_name_seed-*_sample-*"
        )
    if metric_name_type == 0:
        related_metrics = pd.DataFrame()
        counter = 0
        for column_name in pred_stat.columns:
            matched = re.match(rf"(\d)+_{metric_name}", column_name)
            related_metrics.loc[counter, "index"] = int(matched.group(1))
            related_metrics.loc[counter, "metric_name"] = column_name
            counter += 1
        related_metrics.sort_values(by=["index"], inplace=True, ignore_index=True)
        related_metric_names = related_metrics["metric_name"].values
    elif metric_name_type == 1:
        related_metrics = pd.DataFrame()
        counter = 0
        for column_name in pred_stat.columns:
            matched = re.match(rf"{metric_name}_seed-(\d+)_sample-(\d+)", column_name)
            if matched:
                related_metrics.loc[counter, "seed"] = int(matched.group(1))
                related_metrics.loc[counter, "sample"] = int(matched.group(2))
                related_metrics.loc[counter, "metric_name"] = column_name
                counter += 1
        related_metrics.sort_values(
            by=["seed", "sample"], inplace=True, ignore_index=True
        )
        related_metric_names = related_metrics["metric_name"].values
    else:
        raise ValueError(
            f"Metric name type changed for unknown reason. Please check the code."
        )

    model_num = len(related_metric_names)
    y = np.array(
        [j for j in range(pred_stat.shape[0]) for _ in range(model_num)]
    ) + np.random.normal(0, y_fluctuation, pred_stat.shape[0] * model_num)
    y += y_shift

    x = np.array(
        [
            pred_stat.loc[pred_stat.index[i], j]
            for i in range(pred_stat.shape[0])
            for j in related_metric_names
        ]
    )
    a = ax.scatter(x, y, label=label, **scatter_kwargs)

    x_mean = pred_stat[related_metric_names].mean(axis=1)
    x_err = pred_stat[related_metric_names].std(axis=1)

    y = np.array([j for j in range(pred_stat.shape[0])])
    y += y_shift
    color = a.get_facecolor()[0]
    ax.errorbar(x_mean, y, xerr=x_err, color=color, **errorbar_kwargs)

    return fig, ax


def format_metric_distribution_plot(
    pred_stat: pd.DataFrame,
    ax: plt.Axes,
    xlabel: str,
    ylabel: str,
    title: str,
    legend_loc: str = "upper left",
    marker_col: str = "Marker",
):
    ax.set_yticks(range(pred_stat.shape[0]))
    ax.set_yticklabels(pred_stat[marker_col])
    ax.set_ylim(pred_stat.shape[0], -1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc=legend_loc)
