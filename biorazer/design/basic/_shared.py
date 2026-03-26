from pathlib import Path


def _normalize_dir(dir_path: str | Path) -> Path:
    dir_path = Path(dir_path)
    if not dir_path.is_dir():
        raise ValueError(
            f"The directory {dir_path} does not exist or is not a directory."
        )
    if not dir_path.exists():
        raise FileNotFoundError(f"The directory {dir_path} does not exist.")
    return dir_path
