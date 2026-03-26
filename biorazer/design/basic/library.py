from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .entry import Entry
from ._shared import _normalize_dir


@dataclass
class Library:
    """
    Properties
    ----------
    entries: list[Entry]
        list of designs in the library
    """

    entries: list[Entry] = field(default_factory=list)
    library_dir: str | None = None
    dataframe: pd.DataFrame | None = None

    @classmethod
    def from_dir(cls, library_dir: str | Path, entry_type=Entry):
        """
        Create a Library object from a directory.
        This is useful when you want to create a Library object from a directory without an
        associated samples DataFrame.
        """
        library_dir_path = _normalize_dir(library_dir)

        entries = []
        for entry_dir in library_dir_path.glob("[!.]*"):
            if entry_dir.is_dir():
                entry = entry_type.from_dir(entry_dir)
                entries.append(entry)

        library = cls(entries=entries, library_dir=str(library_dir_path))
        return library

    @classmethod
    def from_samples(cls, samples: pd.DataFrame):
        """
        Create a Library object from a samples DataFrame.
        This is useful when you want to create a Library object from a samples DataFrame
        without an associated directory.
        """
        if "entry_marker" not in samples.columns:
            raise ValueError(
                "The samples DataFrame must contain a 'entry_marker' column."
            )

        library = cls(entries=[], library_dir=None, dataframe=samples.copy())
        for marker in samples["entry_marker"].unique():
            entry_samples = samples[samples["entry_marker"] == marker].reset_index(
                drop=True
            )
            entry = cls.entry_type.from_dataframe(entry_samples, marker)
            library.entries.append(entry)
        return library

    def __getitem__(self, identifier: str | int):
        if isinstance(identifier, int):
            return self.entries[identifier]
        elif isinstance(identifier, str):
            for entry in self.entries:
                if entry.marker == identifier:
                    return entry
            raise KeyError(f"No entry found with marker {identifier}.")
        else:
            raise TypeError(
                "Identifier must be either an integer index or a string marker."
            )

    def __iter__(self):
        for entry in self.entries:
            yield entry

    def generate_dataframe(self):
        """
        This function merges samples from all entries into a single DataFrame.
        The merged DataFrame will have an additional column "marker" to indicate the source
        entry.
        """
        if not self.entries:
            raise ValueError("No entries in the library to generate dataframe from.")

        dataframe_list = []
        for entry in self.entries:
            entry.generate_dataframe()
            samples = entry.dataframe.copy()
            dataframe_list.append(samples)

        self.dataframe = pd.concat(dataframe_list, ignore_index=True)

    def flatten_dataframe(self):
        if self.dataframe is None:
            raise ValueError("Samples are not available to flatten.")
        self.dataframe = self.dataframe.groupby("entry_marker").agg(list)
        self.dataframe.reset_index(inplace=True)

    def explode_dataframe(self, columns: list[str]):
        if self.dataframe is None:
            raise ValueError("Samples are not available to explode.")

        columns_to_explode = [
            col for col in self.dataframe.columns if col != "entry_marker"
        ]
        self.dataframe = self.dataframe.explode(columns_to_explode).reset_index(
            drop=True
        )

    def add_simple_metrics(self, *simple_metric_list):
        """
        Add simple metrics to the samples DataFrame by applying the provided metric functions.

        Parameters
        ----------
        simple_metric_list: list of callable
            A list of functions that take a sample row as input and return a metric value.
        """
        for entry in self.entries:
            try:
                entry.add_simple_metrics(*simple_metric_list)
            except ValueError as e:
                raise ValueError(
                    f"Error adding simple metrics to entry {entry.marker}: {e}"
                )

    def map_entry_marker_to_y(self, row_selector=None):
        """
        Call this function before plot_metric_distribution to get a mapping dict from marker
        to y value.
        """
        if self.dataframe is None:
            raise ValueError(
                "Samples are not available. Please call merge_samples() first."
            )
        if row_selector is None:
            row_selector = self.dataframe.index
        df = self.dataframe.loc[row_selector, :].reset_index(drop=True)
        _, uniques = pd.factorize(df["entry_marker"])
        mapping_dict = {marker: float(i) + 1.0 for i, marker in enumerate(uniques)}
        return mapping_dict

    def plot_metric_distribution(
        self,
        ax: plt.Axes,
        metric_name: str,
        mapping_dict: dict,
        y_shift=0.0,
        y_fluctuation: float = 0.1,
        scatter_kwargs=dict(s=5),
        errorbar_kwargs=dict(fmt=",", alpha=0.5, capsize=3),
    ):
        """
        Plot the distribution of a metric for each entry in the library.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes on which to plot the metric distribution.
        metric_name : str
            The name of the metric to plot.
        mapping_dict : dict
            A mapping from entry markers to y-axis values.
        y_shift : float, optional
            A value to shift the y-axis positions, by default 0.0
        y_fluctuation : float, optional
            The standard deviation of the random noise added to the y-axis positions, by default 0.1
        scatter_kwargs : dict, optional
            Additional keyword arguments for the scatter plot, by default dict(s=5)
        errorbar_kwargs : dict, optional
            Additional keyword arguments for the error bars, by default dict(fmt=",", alpha=0.5, capsize=3)

        """
        if self.dataframe is None:
            raise ValueError(
                "Samples are not available. Please call merge_samples() first."
            )

        row_selector = np.isin(
            self.dataframe["entry_marker"], list(mapping_dict.keys())
        )
        df = self.dataframe.loc[row_selector, :].reset_index(drop=True)
        y = df["entry_marker"].map(mapping_dict).astype(float)
        y += np.random.normal(0, y_fluctuation, len(df)) + y_shift
        x = pd.to_numeric(df[metric_name], errors="coerce")
        scatter = ax.scatter(x, y, label=metric_name, **scatter_kwargs)

        df_summary = df.groupby(by="entry_marker").agg({metric_name: ["mean", "std"]})
        y = np.asarray([mapping_dict[i] for i in df_summary.index], dtype=float)
        x_mean = df_summary[metric_name]["mean"].values
        x_std = df_summary[metric_name]["std"].values
        y += y_shift
        color = scatter.get_facecolor()[0]
        ax.errorbar(x_mean, y, xerr=x_std, color=color, **errorbar_kwargs)

    def annotate_designs(self, ax: plt.Axes, mapping_dict: dict):
        """
        Annotate the y ticks with design markers.

        Parameters
        ----------
        mapping_dict: dict
            The mapping dict from marker to y value, which can be obtained by calling
            self.map_entry_marker_to_y(). This is useful when you want plot a subset of designs.
        """
        yticks, yticklabels = [], []
        for key, value in mapping_dict.items():
            yticks.append(value)
            yticklabels.append(key)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)

    @staticmethod
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
                after = (
                    data_list[1].loc[data_list[1][marker_col] == marker, :].index[0] + 1
                )
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
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
            ncol=legend_cols,
        )

    @staticmethod
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

    @staticmethod
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
