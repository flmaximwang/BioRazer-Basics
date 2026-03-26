from __future__ import annotations

import pickle
from abc import abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd
from .single_test import SingleTest
from ._shared import _normalize_dir


@dataclass
class EntryProperty:
    marker: str
    dir_path: Path | None = None
    tests: list[SingleTest] | None = None
    test_type: SingleTest = SingleTest
    dataframe: pd.DataFrame | None = None

    @property
    def parent_dir(self):
        if self.dir_path is None:
            return None
        return self.dir_path.parent

    def load_pickle(self):
        with open(self.pickle, "rb") as f:
            self._pickle_cache = pickle.load(f)

    def unload_pickle(self):
        self._pickle_cache = None


@dataclass
class EntryFormatter(EntryProperty):

    @classmethod
    def _guess_marker_from_dir(cls, dir_path: str | Path) -> str:
        """
        Implement this function to guess the marker name from the directory path.
        By default, it returns the name of the directory. You can override this function
        if your directory structure is more complicated and the marker name cannot be directly obtained from the directory name.
        """
        dir_path = Path(dir_path)
        return dir_path.stem

    @classmethod
    def _prepare_format_paths(
        cls, src_dir: str | Path, target_dir: str | Path, keep_marker=True
    ):
        src_dir = Path(src_dir)
        target_dir = Path(target_dir)
        old_marker = cls._guess_marker_from_dir(src_dir)
        if keep_marker:
            new_marker = old_marker
            target_dir = target_dir / new_marker
        else:
            new_marker = target_dir.stem
        if not src_dir.exists():
            raise FileNotFoundError(f"Source directory {src_dir} does not exist.")
        if not target_dir.exists():
            target_dir.mkdir(parents=True, exist_ok=True)
        return src_dir, target_dir, old_marker, new_marker

    @classmethod
    @abstractmethod
    def format(cls, src_dir: str | Path, target_dir: str | Path, keep_marker=True):
        """
        程序输出的目录结构可能不符合我们分析的要求。实现这个函数来将程序输出的目录结构格式化成我们分析需要的目录结构。
        keep_marker:
            True, 以 src_dir 的目录名, 并调用 _guess_marker_from_dir 来猜测 marker 名, 作为 target_dir 的子目录;
            False, 以 target_dir.stem 作为 marker 名, 直接将 src_dir 的内容复制到 target_dir 中。
        """
        src_dir, target_dir, old_marker, new_marker = cls._prepare_format_paths(
            src_dir, target_dir, keep_marker
        )
        raise NotImplementedError(
            "Please implement the format function for your entry."
        )

    @abstractmethod
    def is_formatted(self) -> bool:
        """
        Implement this function to check if the entry is already formatted.
        Return True if the entry is formatted, False otherwise.
        """

    @abstractmethod
    def _collect_samples(self):
        """
        Implement this function to collect samples from your entries.
        The collected samples should be stored in self.samples as a pandas DataFrame.
        """

    @abstractmethod
    def get_basic_info(self):
        """
        Implement this function to get basic information about your entry and
        store them in self as needed. This function may rely on self.samples, so call it
        after collect_samples.
        """


class EntryPropertyFormatted(EntryFormatter):

    def generate_dataframe(self):
        if self.tests is None:
            raise ValueError("Tests are not provided for this entry.")
        data = [test.to_dict() for test in self.tests]
        self.dataframe = pd.DataFrame(data)


@dataclass
class EntryIO(EntryPropertyFormatted):

    @classmethod
    def from_dir(cls, entry_dir: str | Path):
        """
        Create an Entry object from a directory.
        This is useful when you want to create an Entry object from a directory without an
        associated samples DataFrame.
        """
        entry_dir = _normalize_dir(entry_dir)
        marker = entry_dir.stem
        single_test_dirs = sorted(list(entry_dir.glob(f"{marker}*")))
        tests = []
        for test_dir in single_test_dirs:
            test = cls.test_type.from_dir(test_dir, entry_marker=marker)
            tests.append(test)
        entry = cls(marker=marker, dir_path=entry_dir, tests=tests)

        return entry

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame, marker: str):
        """
        Create an Entry object from a marker and a samples DataFrame.
        This is useful when you want to create an Entry object from a samples DataFrame
        without an associated directory.
        """
        entry = cls(marker=marker, dir_path=None, dataframe=dataframe.copy())
        entry.get_basic_info()
        return entry

    def to_dir(self, dir_path: str | Path, add_marker=True):
        """
        Implement this function to save the entry to a directory in a formatted manner.
        if add_marker is True, the directory structure should be like this:
        ```
        parent_dir/
            marker/
                ... # your output files
        ```
        If add_marker is False, the directory structure should be like this:
        ```
        parent_dir/
                ... # your output files
        ```
        """
        if add_marker:
            target_dir_path = Path(dir_path) / self.marker
        else:
            target_dir_path = Path(dir_path)
        target_dir_path.mkdir(parents=True, exist_ok=True)
        return target_dir_path


@dataclass
class Entry(EntryIO):
    """
    An entry is one design.
    An entry may contain multiple samples (e.g., data from different seeds, different
    experimental repeats, etc.). An entry is associated with a directory, within which
    related scores, models, and metadata are stored. An entry is identified by a marker,
    which is the name of the directory.

    Properties
    ----------
    marker: str
        the name of the entry, usually the name of the directory
    parent_dir: str
        the parent directory of the entry directory. Necessary when you want read metadata
        files from the directory.
    samples: pd.DataFrame
        a DataFrame containing the samples of the entry.
        Each sample (row) corresponds to a specific condition (e.g., different seeds,
        different experimental repeats).
        Each column (metric/metadata) contains the value of a specific metric/metadata for
        all samples.
    """

    def flatten_dataframe(self):
        """
        Normally, self.samples has multiple rows, each row is a sample from a specific seed
        and a specific experimental repeat. This function flattens the samples into a
        single row and merge data from each cell into a list. After flattening, self.samples will have only one row,
        and each cell will contain a list of values from the original samples.

        Call self.explode_samples to explode the specified columns after flattening if you want to restore the original samples.
        """
        if self.dataframe is None:
            raise ValueError("Samples are not available to flatten.")
        self.dataframe = self.dataframe.groupby("entry_marker").agg(list)
        self.dataframe.reset_index(inplace=True)

    def explode_dataframe(self):
        """
        Explode the specified columns in self.dataframe.

        Parameters
        ----------
        columns: list of str
            The columns to be exploded.
        """
        if self.dataframe is None:
            raise ValueError("Samples are not available to explode.")

        columns_to_explode = [
            col for col in self.dataframe.columns if col != "entry_marker"
        ]
        self.dataframe = self.dataframe.explode(columns_to_explode).reset_index(
            drop=True
        )

    @property
    def simple_metrics(self):
        test_type_instance: SingleTest = self.test_type(
            entry_marker=self.marker, marker="temp"
        )
        return test_type_instance.simple_metrics

    def add_simple_metrics(self, *metric_key_list):
        """
        You may not want to extract all metrics from the samples at the beginning of your
        analysis. Call this function to add extra simple metrics to self.samples based on
        data in your entry directory.
        """

        for single_test in self.tests:
            for metric_key in metric_key_list:
                single_test.add_simple_metric(metric_key)
