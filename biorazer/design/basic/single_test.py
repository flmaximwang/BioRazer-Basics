import pickle
from dataclasses import dataclass, field
from pathlib import Path


def _normalize_paths(test_dir: str | Path, metadata_files):
    if test_dir is not None and not isinstance(test_dir, Path):
        test_dir = Path(test_dir)

    normalized_metadata_files = metadata_files.copy()
    for key, value in normalized_metadata_files.items():
        if value is not None and not isinstance(value, Path):
            normalized_metadata_files[key] = Path(value)

    return test_dir, normalized_metadata_files


def _validate_paths(test_dir: Path, metadata_files: dict[str, Path]):
    if test_dir is not None and not test_dir.exists():
        raise FileNotFoundError(f"Test directory {test_dir} does not exist.")

    for path in metadata_files.values():
        if path is not None and not path.exists():
            raise FileNotFoundError(f"Metadata file {path} does not exist.")


@dataclass
class SingleTest:
    """
    This class is used to analyze the results of a single test (e.g., a single design, a single seed, etc.).
    Directory Structure:
    - marker/
        - marker_data.pickle
        - marker_*.* (other files related to the test, e.g., confidence json, cif files, etc.)
    """

    entry_marker: str
    marker: str
    test_dir: Path | None = None  # Path to the directory containing the test results
    metadata_files: dict = field(
        default_factory=lambda: {
            "pickle": None,
        }
    )  # Store any additional metadata related to the test
    metadata_file_suffixes: dict = field(
        default_factory=lambda: {
            "pickle": "_data.pkl",
        }
    )  # Expected suffixes for the metadata files, used for auto-discovery
    pickle_data: dict | None = None  # Cache for the loaded pickle
    # Map to tell how to extract simple metrics from the pickle cache.
    # The key is the metric name, and the value is a list of keys to access the metric data in the pickle cache.
    simple_metrics: dict = field(default_factory=dict)
    metrics: dict = field(
        default_factory=dict
    )  # Store the extracted metrics for the test

    def __post_init__(self):
        self.marker = str(self.marker)
        self.test_dir, self.metadata_files = _normalize_paths(
            self.test_dir, self.metadata_files
        )
        _validate_paths(self.test_dir, self.metadata_files)

    @property
    def pickle(self):
        return self.metadata_files.get("pickle", None)

    @pickle.setter
    def pickle(self, value):
        self.metadata_files["pickle"] = None if value is None else Path(value)

    @classmethod
    def from_dir(cls, test_dir: str | Path, entry_marker: str):
        """
        Create a SingleTest object from a directory.
        This is useful when you want to create a SingleTest object from a directory without an
        associated samples DataFrame.
        """
        test_dir_path = Path(test_dir)
        merged_marker = test_dir_path.stem
        if not merged_marker.startswith(entry_marker):
            raise ValueError(
                f"The directory name {merged_marker} does not seem to start with the entry marker {entry_marker}. Please check the directory name and the entry marker."
            )
        test_marker = merged_marker[len(entry_marker) :]
        res_obj = cls(
            marker=test_marker, test_dir=test_dir_path, entry_marker=entry_marker
        )
        metadata_files = res_obj.metadata_files
        metadata_file_suffixes = res_obj.metadata_file_suffixes
        for key, suffix in metadata_file_suffixes.items():
            expected_file = test_dir_path / f"{merged_marker}{suffix}"
            if expected_file.exists():
                metadata_files[key] = expected_file
            else:
                raise FileNotFoundError(
                    f"No metadata file with suffix {suffix} is found in the directory {test_dir_path} for the test {test_marker}. Please check the directory and the expected metadata files. If you want to use a different metadata file structure, please specify the metadata files explicitly when creating the SingleTest object or implement a custom from_dir method."
                )
        return res_obj

    def load_pickle(self):
        with open(self.pickle, "rb") as f:
            self.pickle_data = pickle.load(f)

    def unload_pickle(self):
        self.pickle_data = None

    def add_simple_metric(self, metric_key):
        if self.pickle_data is None:
            self.load_pickle()
        if metric_key not in self.simple_metrics:
            raise ValueError(f"Metric {metric_key} is not found in the pickle cache.")
        for i in self.simple_metrics[metric_key]:
            try:
                temp = self.pickle_data[i]
                self.metrics[metric_key] = temp
            except KeyError:
                self.metrics[metric_key] = None

    def to_dict(self):
        result = {"entry_marker": self.entry_marker, "marker": self.marker}
        result.update({k: str(v) for k, v in self.metadata_files.items()})
        result.update(self.metrics)
        return result

    def rename_files(self, new_marker: str | None = None):
        if new_marker is not None:
            self.marker = str(new_marker)

        old_test_dir = self.test_dir
        new_marker = self.marker
        old_marker = None if old_test_dir is None else old_test_dir.stem
        new_test_dir = (
            None if old_test_dir is None else old_test_dir.with_name(new_marker)
        )

        if (
            old_test_dir is not None
            and old_test_dir.exists()
            and new_test_dir != old_test_dir
            and new_test_dir.exists()
        ):
            raise FileExistsError(
                f"Target test directory {new_test_dir} already exists."
            )

        rename_plan = []
        for key, path in self.metadata_files.items():
            if path is None:
                continue

            source_path = path
            if old_test_dir is not None:
                try:
                    relative_path = path.relative_to(old_test_dir)
                except ValueError:
                    relative_path = None
                else:
                    if new_test_dir is not None:
                        source_path = new_test_dir / relative_path

            target_name = path.name
            if old_marker is not None and path.name.startswith(old_marker):
                target_name = f"{new_marker}{path.name[len(old_marker):]}"
            target_path = source_path.with_name(target_name)

            if source_path != target_path and target_path.exists():
                raise FileExistsError(
                    f"Target metadata file {target_path} already exists."
                )

            rename_plan.append((key, source_path, target_path))

        if (
            old_test_dir is not None
            and old_test_dir.exists()
            and new_test_dir != old_test_dir
        ):
            old_test_dir.rename(new_test_dir)

        updated_metadata_files = self.metadata_files.copy()
        for key, source_path, target_path in rename_plan:
            if source_path.exists() and source_path != target_path:
                source_path.rename(target_path)
            updated_metadata_files[key] = target_path

        self.test_dir = new_test_dir
        self.metadata_files = updated_metadata_files
