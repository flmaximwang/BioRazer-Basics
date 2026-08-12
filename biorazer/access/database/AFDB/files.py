"""AlphaFold DB 预测结构文件下载。"""
from pathlib import Path

from ....logger import initialize_logger
from .info import AFDBEntry
from .query import uniprot_to_entries

# AlphaFold DB 提供的文件类型 (与 AFDBEntry.file_types 一致)
SUPPORTED_FMT = ["bcif", "cif", "pdb", "pdbImage", "plddtDoc", "paeDoc"]


def fetch(
    uniprot_id: str,
    fmt: str = "cif",
    download_dir: str | Path = ".",
    overwrite=False,
    logger=None,
):
    """下载 AlphaFold DB 中该 UniProt 登录号对应的预测结构文件。

    Parameters
    ----------
    uniprot_id : str
        UniProt 登录号, 例如 "P0DP23"。
    fmt : str
        文件类型, 见 SUPPORTED_FMT。
    download_dir : str | Path
        下载目录。
    overwrite : bool
        已存在同名文件时是否覆盖。
    logger : logging.Logger
        日志器, 缺省自动初始化。

    Returns
    -------
    Path
        下载得到的文件路径。
    """
    if not logger:
        logger = initialize_logger(__name__)
    if fmt not in SUPPORTED_FMT:
        raise ValueError(
            f"Unsupported format: {fmt} (supported: {', '.join(SUPPORTED_FMT)})"
        )
    entries = uniprot_to_entries(uniprot_id)
    if not isinstance(entries, list) or not entries or not isinstance(entries[0], AFDBEntry):
        raise ValueError(f"No AFDB entry found for {uniprot_id}")
    entry = entries[0]
    file_path = Path(download_dir) / f"{entry.id}.{fmt}"
    if file_path.exists() and not overwrite:
        logger.warning(f"{file_path} already exists, skipping")
        return file_path
    Path(download_dir).mkdir(parents=True, exist_ok=True)
    entry.download(fmt, folder_dir=str(download_dir))
    logger.info(f"{entry.id}.{fmt} downloaded to {download_dir}")
    return file_path
