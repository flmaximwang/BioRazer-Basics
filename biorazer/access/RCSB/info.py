import requests, json
from abc import abstractmethod
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET
from .files import PDBStructure


@dataclass
class PDBObject:
    """
    Abstract base class for RCSB PDB objects
    """

    id: str = None

    def dump(self, data_dir=".", indent=2):
        if not Path(data_dir).exists():
            Path(data_dir).mkdir(parents=True, exist_ok=True)
        json.dump(self.data, open(f"{data_dir}/{self.id}.json", "w"), indent=indent)

    @classmethod
    def load(cls, json_file: str | Path):
        data = json.load(open(json_file).read())
        entry = cls(id=data["rcsb_id"])
        entry.data = data
        return entry

    @abstractmethod
    def _fetch(self):
        """Fetch data from RCSB API and store it in self.data"""

    def fetch(self, retry: int = 3):
        while retry > 0:
            try:
                self._fetch()
                if "status" in self.data:
                    raise Exception(f"[{self.data['status']}] {self.data['message']}")
                return
            except requests.exceptions.ConnectionError as e:
                print(f"ConnectionError: {e}, retrying... ({retry} retries left)")
                retry -= 1
        raise ConnectionError(
            f"Failed to fetch data from RCSB API after {retry} retries"
        )


@dataclass
class PDBEntry(PDBObject):
    """
    Entry 是 RCSB PDB 中最高级的条目, 一个 Entry 包括多个 Entity, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    """

    id: str = None

    def _fetch(self):
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{self.id.upper()}")
        self.data = json.loads(r.text)

    def download(self, fmt="pdb", folder_dir=".", overwrite=False, logger=None):
        structure = PDBStructure(pdb_code=self.id)
        structure.download(
            fmt=fmt, folder_dir=folder_dir, overwrite=overwrite, logger=logger
        )

    @property
    def info(self):
        return self.data["rcsb_entry_info"]

    @property
    def title(self):
        return self.data["struct"]["title"]

    @property
    def entity_counts(self):
        return {
            "polymer_entity_count": (
                self.info["polymer_entity_count"]
                if "polymer_entity_count" in self.info
                else 0
            ),
            "nonpolymer_entity_count": (
                self.info["nonpolymer_entity_count"]
                if "nonpolymer_entity_count" in self.info
                else 0
            ),
            "branched_entity_count": (
                self.info["branched_entity_count"]
                if "branched_entity_count" in self.info
                else 0
            ),
        }

    @property
    def entity_ids(self):
        if "polymer_entity_ids" in self.data["rcsb_entry_container_identifiers"]:
            polymer_entity_ids: list = self.data["rcsb_entry_container_identifiers"][
                "polymer_entity_ids"
            ]
        else:
            polymer_entity_ids = []
        if "non_polymer_entity_ids" in self.data["rcsb_entry_container_identifiers"]:
            non_polymer_entity_ids: list = self.data[
                "rcsb_entry_container_identifiers"
            ]["non_polymer_entity_ids"]
        else:
            non_polymer_entity_ids = []
        if "branched_entity_ids" in self.data["rcsb_entry_container_identifiers"]:
            branched_polymer_entity_ids: list = self.data[
                "rcsb_entry_container_identifiers"
            ]["branched_entity_ids"]
        else:
            branched_polymer_entity_ids = []
        all_entity_ids = []
        for entity_id in polymer_entity_ids:
            all_entity_ids.append((PDBPolymerEntity, f"{self.id}_{entity_id}"))
        for entity_id in non_polymer_entity_ids:
            all_entity_ids.append((PDBNonpolymerEntity, f"{self.id}_{entity_id}"))
        for entity_id in branched_polymer_entity_ids:
            all_entity_ids.append((PDBBranchedEntity, f"{self.id}_{entity_id}"))
        all_entity_ids.sort(key=lambda x: int(x[1].split("_")[1]))
        return all_entity_ids

    @property
    def citation(self):
        return self.data["rcsb_primary_citation"]


@dataclass
class PDBEntity(PDBObject):
    """
    Entity 是 RCSB PDB 中第二级的条目, 一个 Entity 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    输入的 entity_id 形如 "6LU7_1"
    """

    @property
    def entry_id(self):
        return self.id.split("_")[0]

    @property
    def entity_index(self):
        return self.id.split("_")[1]

    @property
    @abstractmethod
    def asym_ids(self):
        """Return the list of asym_ids associated with this entity"""

    @property
    @abstractmethod
    def auth_asym_ids(self):
        """Return the list of auth_asym_ids associated with this entity"""

    @property
    @abstractmethod
    def database_accessions(self):
        """Return the list of database accessions associated with this entity"""


@dataclass
class PDBPolymerEntity(PDBEntity):
    """
    PDBPolymerEntity 是 RCSB PDB 中的聚合物实体, 例如蛋白质, 核酸等
    输入的 entity_id 形如 "6LU7_1"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{self.entry_id}/{self.entity_index}"
        )
        self.data = json.loads(r.text)

    @property
    def sequence(self):
        return self.data["entity_poly"]["pdbx_seq_one_letter_code_can"]

    @property
    def asym_ids(self):
        return self.data["rcsb_polymer_entity_container_identifiers"]["asym_ids"]

    @property
    def auth_asym_ids(self):
        return self.data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]

    @property
    def database_accessions(self):
        return self.data["rcsb_polymer_entity_container_identifiers"][
            "reference_sequence_identifiers"
        ]

    @property
    def src(self):
        return self.data["entity_src_nat"]


@dataclass
class PDBNonpolymerEntity(PDBEntity):
    """
    PDBNonpolymerEntity 是 RCSB PDB 中的非聚合物实体, 例如小分子配体, 离子等
    输入的 entity_id 形如 "6LU7_4"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{self.entry_id}/{self.entity_index}"
        )
        self.data = json.loads(r.text)

    @property
    def asym_ids(self):
        return self.data["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"]

    @property
    def auth_asym_ids(self):
        return self.data["rcsb_nonpolymer_entity_container_identifiers"][
            "auth_asym_ids"
        ]

    @property
    def database_accessions(self):
        return [
            self.data["rcsb_nonpolymer_entity_container_identifiers"]["chem_ref_def_id"]
        ]


@dataclass
class PDBBranchedEntity(PDBEntity):
    """
    PDBBranchedEntity 是 RCSB PDB 中的支链实体, 例如多糖等
    输入的 entity_id 形如 "6LU7_3"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/branched_entity/{self.entry_id}/{self.entity_index}"
        )
        self.data = json.loads(r.text)

    @property
    def asym_ids(self):
        return self.data["rcsb_branched_entity_container_identifiers"]["asym_ids"]

    @property
    def auth_asym_ids(self):
        return self.data["rcsb_branched_entity_container_identifiers"]["auth_asym_ids"]

    @property
    def database_accessions(self):
        return self.data["rcsb_branched_entity_container_identifiers"][
            "chem_comp_monomers"
        ]


@dataclass
class PDBChemicalComponent(PDBObject):
    """
    Chemical Component 是 RCSB PDB 中的化学成分, 例如氨基酸, 核苷酸, 小分子配体等
    输入的 chem_comp_id 形如 "ATP"
    """

    def _fetch(self):
        r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{self.id}")
        self.data = json.loads(r.text)


@dataclass
class PDBInstance(PDBObject):
    """
    Entity Instance 是 RCSB PDB 中第三级的条目, 一个 Entity Instance 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity_instance, nonpolymer_entity_instance, branched_entity_instance,
    输入的 entity_instance_id 形如 "6LU7_A"
    """

    @property
    def entry_id(self):
        return self.id.split("_")[0]

    @property
    def entity_instance_index(self):
        return self.id.split("_")[1]

    @abstractmethod
    def _fetch(self):
        """Fetch data from RCSB API and store it in self.data"""


@dataclass
class PDBPolymerInstance(PDBInstance):
    """
    PDBPolymerEntityInstance 是 RCSB PDB 中的聚合物实体实例, 例如蛋白质链, 核酸链等
    输入的 entity_instance_id 形如 "6LU7_A"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{self.entry_id}/{self.entity_instance_index}"
        )
        self.data = json.loads(r.text)


@dataclass
class PDBNonPolymerInstance(PDBInstance):
    """
    PDBNonpolymerEntityInstance 是 RCSB PDB 中的非聚合物实体实例, 例如小分子配体实例, 离子实例等
    输入的 entity_instance_id 形如 "6LU7_B"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{self.entry_id}/{self.entity_instance_index}"
        )
        self.data = json.loads(r.text)


@dataclass
class PDBBranchedInstance(PDBInstance):
    """
    PDBBranchedEntityInstance 是 RCSB PDB 中的支链实体实例, 例如多糖链等
    输入的 entity_instance_id 形如 "6LU7_C"
    """

    def _fetch(self):
        r = requests.get(
            f"https://data.rcsb.org/rest/v1/core/branched_entity_instance/{self.entry_id}/{self.entity_instance_index}"
        )
        self.data = json.loads(r.text)


def PDBPolymer(id: str):
    if len(id.split("_")) == 1:
        return PDBEntry(id)
    pdb_code, extra_id = id.split("_")
    if extra_id.isdigit():
        return PDBPolymerEntity(id)
    elif extra_id.isalpha():
        return PDBPolymerInstance(id)
    else:
        raise ValueError(f"Invalid PDBPolymer id: {id}")
