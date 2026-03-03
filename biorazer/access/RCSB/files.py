import logging, requests, re
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET
from html import unescape


@dataclass
class PDBStructure:
    pdb_code: str = None

    def download(
        self,
        fmt="pdb",
        folder_dir=".",
        overwrite=False,
        logger: logging.Logger = None,
    ):
        if not logger:
            logger = logging.getLogger(__name__)
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s - %(levelname)s - %(message)s",
            )
        logger.info(f"Downloading {self.pdb_code}.{fmt} to {folder_dir}")
        if Path(f"{folder_dir}/{self.pdb_code}.{fmt}").exists() and not overwrite:
            logger.warning(
                f"{self.pdb_code}.{fmt} already exists in {folder_dir}, skipping"
            )
        else:
            r = requests.get(f"https://files.rcsb.org/download/{self.pdb_code}.{fmt}")
            if r.status_code != 200:
                logger.warning(f"\n{r.text}")
            if not Path(folder_dir).exists():
                Path(folder_dir).mkdir(parents=True)
            with open(Path(f"{folder_dir}") / f"{self.pdb_code}.{fmt}", "wb") as f:
                f.write(r.content)
            if logger:
                logger.info(f"{self.pdb_code}.{fmt} downloaded to {folder_dir}")
