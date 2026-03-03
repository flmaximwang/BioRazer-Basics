import requests
from pathlib import Path
from dataclasses import dataclass


@dataclass
class AFDBEntry:

    data: dict = None

    @property
    def id(self):
        return self.data["entryId"]

    @property
    def file_types(self):
        return ["bcif", "cif", "pdb", "pdbImage", "plddtDoc", "paeDoc"]

    @property
    def file_Urls(self):
        result = {}
        for file_type in self.file_types:
            key = f"{file_type}Url"
            result[key] = self.data.get(key, None)
        return result

    def download(self, file_type: str, folder_dir="."):

        r = requests.get(self.file_Urls[f"{file_type}Url"])
        if r.status_code != 200:
            raise ValueError(f"File type {file_type} not found for entry {self.id}")
        with open(Path(folder_dir) / f"{self.id}.{file_type}", "wb") as f:
            f.write(r.content)
