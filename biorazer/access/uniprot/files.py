import requests
from pathlib import Path
from ..logger import initialize_logger

REST_API = "https://rest.uniprot.org"


def fetch(
    uniprot_id: str,
    fmt: str,
    download_dir: str | Path = ".",
    overwrite=False,
    logger=None,
):
    if not logger:
        logger = initialize_logger(__name__)
    if not overwrite:
        file_path = Path(download_dir) / f"{uniprot_id}.{fmt}"
        if file_path.exists():
            return file_path
    if fmt in ["fasta", "fa"]:
        return _fetch_fasta(
            uniprot_id, download_dir, overwrite=overwrite, logger=logger
        )
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def _fetch_fasta(uniprot_id: str, download_dir: str | Path, overwrite, logger):
    logger.info(f"Fetching {uniprot_id} in FASTA format from UniProt")
    url = f"{REST_API}/uniprotkb/{uniprot_id}.fasta"
    headers = {"accept": "text/x-fasta"}
    r = requests.get(url, headers=headers)
    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)
    file_path = download_dir / f"{uniprot_id}.fasta"
    with open(file_path, "w") as f:
        f.write(r.text)
    if logger:
        logger.info(f"{uniprot_id}.fasta downloaded to {download_dir}")
    return file_path
