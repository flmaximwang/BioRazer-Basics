from dataclasses import dataclass, field
from pathlib import Path
from abc import abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class Entry:
    """
    An entry is one design.
    An entry may contain multiple samples (e.g., data from different seeds, different experimental repeats, etc.).
    An entry is associated with a directory, within which related scores, models, and metadata are stored.
    An entry is identified by a marker, which is the name of the directory.

    Properties
    ----------
    marker: str
        the name of the entry, usually the name of the directory
    parent_dir: str
        the parent directory of the entry directory. Necessary when you want read metadata files from the directory.
    samples: pd.DataFrame
        a DataFrame containing the samples of the entry.
    """

    marker: str
    parent_dir: str = None
    samples: pd.DataFrame = None

    @classmethod
    def from_dir(cls, entry_dir, formatted=True):
        """
        Create an Entry object from a directory.
        This is useful when you want to create an Entry object from a directory without an associated samples DataFrame.

        Parameters
        ----------
        my_dir: str or Path
            the path to the directory associated with the entry.

        Returns
        -------
        entry: Entry
            an Entry object with the given directory.
        """
        entry_dir_path = Path(entry_dir)
        entry_parent_dir_path = entry_dir_path.parent
        marker = entry_dir_path.name
        entry: Entry = cls(
            marker=marker, parent_dir=str(entry_parent_dir_path), samples=None
        )
        if formatted:
            entry.collect_samples()
            entry.get_basic_info()
        return entry

    @classmethod
    def from_samples(cls, samples: pd.DataFrame, marker: str):
        """
        Create an Entry object from a marker and a samples DataFrame.
        This is useful when you want to create an Entry object from a samples DataFrame without an associated directory.

        Parameters
        ----------
        marker: str
            the name of the entry
        samples: pd.DataFrame
            a DataFrame containing the samples of the entry.

        Returns
        -------
        entry: Entry
            an Entry object with the given marker and samples.
        """
        entry: Entry = cls(marker=marker, parent_dir=None, samples=samples)
        entry.get_basic_info()
        return entry

    def to_parent_dir(self, parent_dir, new_marker=None):
        """
        Implement this function to save the entry to a directory in a formatter manner.
        The directory structure should be like this:
        ```
        parent_dir/
            marker/
                ... # your output files
        ```
        """
        if new_marker is None:
            target_dir_path = Path(parent_dir) / self.marker
        else:
            target_dir_path = Path(parent_dir) / new_marker
        if not target_dir_path.exists():
            target_dir_path.mkdir(parents=True)

        return target_dir_path

    @abstractmethod
    def is_formatted(self) -> bool:
        """
        Implement this function to check if the entry is already formatted.
        Return True if the entry is formatted, False otherwise.
        """

    @abstractmethod
    def collect_samples(self, overwrite=False):
        """
        Implement this function to collect samples from your entries.
        The collected samples should be stored in self.samples as a pandas DataFrame.
        """

    @abstractmethod
    def get_basic_info(self):
        """
        Implement this function to get basic information about your entry and
        stored them in self as needed. This function may rely on self.samples, so call it
        after collect_samples.
        """

    def get_marker(self) -> str:
        return self.marker

    def get_dir_path(self):
        return Path(self.parent_dir) / self.marker

    @abstractmethod
    def get_sample(self, *args):
        """
        Implement this function to get a specific sample from your entry.
        The returned sample should be a pandas Series.
        """

    def get_samples(self) -> pd.DataFrame:
        """
        Return
        ------
        samples: pd.DataFrame
            A DataFrame containing the samples of the entry.
        """
        return self.samples

    def set_samples(self, samples: pd.DataFrame):
        self.samples = samples

    def set_new_marker(self, new_marker):
        """
        Calling this function before format_output will change the marker name used in the output directory.
        """
        self.new_marker = new_marker

    def flatten_samples(self):
        """
        Normally, self.samples has multiple rows, each row is a sample from a specific seed and a specific experimental repeat.
        This function should flatten the samples into a single row
        """
        self.set_samples(self.get_samples().groupby("marker").agg(list))
        self.get_samples().reset_index(inplace=True)

    def explode_samples(self, columns):
        """
        Explode the specified columns in self.samples.
        This is useful when you want to convert list-like entries in the specified columns into separate rows.

        Parameters
        ----------
        columns: list of str
            The columns to be exploded.
        """
        self.set_samples(self.get_samples().explode(columns).reset_index(drop=True))

    @abstractmethod
    def extract_simple_metrics(self, *args):
        """
        You may not want to extract all metrics from the samples at the beginning of your analysis.
        Call this function to add extra simple metrics to self.samples.
        """


@dataclass
class Library:
    """

    Properties
    ----------
    entries: list[Entry]
        list of designs in the library
    """

    entries: list[Entry] = field(default_factory=list)
    library_dir: str = None
    samples: pd.DataFrame = None
    entry_type = Entry

    @classmethod
    def from_dir(cls, library_dir, formatted=True):
        """
        Create a Library object from a directory.
        This is useful when you want to create a Library object from a directory without an associated samples DataFrame.

        Parameters
        ----------
        library_dir: str or Path
            the path to the directory associated with the library.

        Returns
        -------
        library: Library
            a Library object with the given directory.
        """
        library_dir_path = Path(library_dir)
        if not library_dir_path.exists():
            raise ValueError(f"The directory {library_dir} does not exist.")
        if not library_dir_path.is_dir():
            raise ValueError(f"The path {library_dir} is not a directory.")
        entries = []
        for entry_dir in library_dir_path.glob("[!.]*"):
            if entry_dir.is_dir():
                entry = cls.entry_type.from_dir(entry_dir, formatted=formatted)
                entries.append(entry)
        library: Library = cls(
            entries=entries, library_dir=str(library_dir_path), samples=None
        )
        if formatted:
            library.merge_samples()
        return library

    @classmethod
    def from_samples(cls, samples: pd.DataFrame):
        """
        Create a Library object from a samples DataFrame.
        This is useful when you want to create a Library object from a samples DataFrame without an associated directory.

        Parameters
        ----------
        samples: pd.DataFrame
            a DataFrame containing the samples of the library.

        Returns
        -------
        library: Library
            a Library object with the given samples.
        """
        if not "marker" in samples.columns:
            raise ValueError("The samples DataFrame must contain a 'marker' column.")
        library: Library = cls(entries=[], library_dir=None, samples=samples)
        for marker in samples["marker"].unique():
            entry_samples = samples[samples["marker"] == marker].reset_index(drop=True)
            entry = cls.entry_type.from_samples(entry_samples, marker)
            library.entries.append(entry)
        return library

    def __getitem__(self, idx):
        return self.entries[idx]

    def __iter__(self):
        for entry in self.entries:
            yield entry

    def merge_samples(self):
        """
        This function merges samples from all entries into a single DataFrame.
        The merged DataFrame will have an additional column "marker" to indicate the source entry.
        """

        samples_list = []
        markers = []
        for entry in self.entries:
            samples = entry.samples
            samples_list.append(samples)
            markers.extend([entry.marker] * len(samples))
        all_samples = samples_list[0].copy()
        for i in range(1, len(samples_list)):
            all_samples = pd.concat([all_samples, samples_list[i]], ignore_index=True)
        all_sample_cols = all_samples.columns.to_list()
        all_samples["marker"] = markers
        all_samples = all_samples.loc[:, ["marker"] + all_sample_cols]
        self.samples = all_samples
        return all_samples

    def get_samples(self) -> pd.DataFrame:
        """
        After calling collect_samples, this function becomes available to return the collected samples.
        """
        return self.samples

    def set_samples(self, samples: pd.DataFrame):
        self.samples = samples

    def map_marker_to_y(self, row_selector=None):
        """
        Call this function before plot_metric_distribution to get a mapping dict from marker to y value.

        Parameters
        ----------
        - row_selector: a boolean array to select rows from self.samples. If None, all rows are used.

        Returns
        -------
        - The returned mapping dict is like {"marker1": 1.0, "marker2": 2.0, ...}
        """
        if row_selector is None:
            row_selector = self.get_samples().index
        df = self.get_samples().loc[row_selector, :].reset_index(drop=True)
        codes, uniques = pd.factorize(df["marker"])
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
        row_selector = np.isin(self.get_samples()["marker"], list(mapping_dict.keys()))
        df = self.get_samples().loc[row_selector, :].reset_index(drop=True)
        y = df["marker"].map(mapping_dict)
        y += np.random.normal(0, y_fluctuation, len(df)) + y_shift
        x = df[metric_name]
        scatter = ax.scatter(x, y, label=metric_name, **scatter_kwargs)

        df_summary = df.groupby(by="marker").agg({metric_name: ["mean", "std"]})
        y = np.asarray([mapping_dict[i] for i in df_summary.index])
        x_mean = df_summary[metric_name]["mean"].values
        x_std = df_summary[metric_name]["std"].values
        y += y_shift
        color = scatter.get_facecolor()[0]
        errorbar = ax.errorbar(x_mean, y, xerr=x_std, color=color, **errorbar_kwargs)

    def annotate_designs(self, ax: plt.Axes, mapping_dict: dict):
        """
        Annotate the y ticks with design markers.

        Parameters
        ----------
        mapping_dict: dict
            The mapping dict from marker to y value, which can be obtained by calling self.map_marker_to_y().
            This is useful when you want plot a subset of designs.
        """
        yticks, yticklabels = [], []
        for key, value in mapping_dict.items():
            yticks.append(value)
            yticklabels.append(key)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
