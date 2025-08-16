import os, json, requests
from rcsbsearchapi import AttributeQuery

def query_uniprot(uniprot_id, return_type="entry"):
    '''
    This function queries all RCSB codes that contain the target protein denoted by the Uniprot accession ID
    Args:
    - return_type: str, one of ["entry", "polymer_entity", "nonpolymer_entity", "branched_entity"]
    '''
    assert return_type in ["entry", "polymer_entity", "nonpolymer_entity", "branched_entity"]
    q_uniprot = AttributeQuery(
        attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        operator="exact_match",
        value= uniprot_id
    )
    res = list(q_uniprot(return_type))
    return res

def get_entry(entry_id: str, dump=False, data_dir = "."):
    '''
    Entry 是 RCSB PDB 中最高级的条目, 一个 Entry 包括多个 Entity, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    '''
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{entry_id.upper()}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entry_id}.json", "w"), indent=2)
    return data

def get_polymer_entity(entity_id: str, dump=False, data_dir="."):
    '''
    Entity 是 RCSB PDB 中第二级的条目, 一个 Entity 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    输入的 entity_id 形如 "6LU7_1"
    '''
    entry_id, entity_index = entity_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry_id}/{entity_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_id}.json", "w"), indent=2)
    return data

def get_nonpolymer_entity(entity_id: str, dump=False, data_dir="."):
    '''
    Entity 是 RCSB PDB 中第二级的条目, 一个 Entity 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    输入的 entity_id 形如 "6LU7_1"
    '''
    entry_id, entity_index = entity_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{entry_id}/{entity_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_id}.json", "w"), indent=2)
    return data

def get_branched_entity(entity_id: str, dump=False, data_dir="."):
    '''
    Entity 是 RCSB PDB 中第二级的条目, 一个 Entity 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity, nonpolymer_entity, branched_entity,
    输入的 entity_id 形如 "6LU7_1"
    '''
    entry_id, entity_index = entity_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/branched_entity/{entry_id}/{entity_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_id}.json", "w"), indent=2)
    return data

def get_chemical_component(chem_comp_id: str, dump=False, data_dir="."):
    chem_comp_id = chem_comp_id.upper()
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{chem_comp_id}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{chem_comp_id}.json", "w"), indent=2)
    return data

def get_polymer_entity_instance(entity_instance_id: str, dump=False, data_dir="."):
    '''
    Entity Instance 是 RCSB PDB 中第三级的条目, 一个 Entity Instance 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity_instance, nonpolymer_entity_instance, branched_entity_instance,
    输入的 entity_instance_id 形如 "6LU7_A"
    '''
    entry_id, entity_instance_index = entity_instance_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{entry_id}/{entity_instance_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_instance_id}.json", "w"), indent=2)
    return data

def get_nonpolymer_entity_instance(entity_instance_id: str, dump=False, data_dir="."):
    '''
    Entity Instance 是 RCSB PDB 中第三级的条目, 一个 Entity Instance 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity_instance, nonpolymer_entity_instance, branched_entity_instance,
    输入的 entity_instance_id 形如 "6LU7_A"
    '''
    entry_id, entity_instance_index = entity_instance_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{entry_id}/{entity_instance_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_instance_id}.json", "w"), indent=2)
    return data

def get_branched_entity_instance(entity_instance_id: str, dump=False, data_dir="."):
    '''
    Entity Instance 是 RCSB PDB 中第三级的条目, 一个 Entity Instance 包括多个 Entity Instance, 并且有多类, 包括 polymer_entity_instance, nonpolymer_entity_instance, branched_entity_instance,
    输入的 entity_instance_id 形如 "6LU7_A"
    '''
    entry_id, entity_instance_index = entity_instance_id.split("_")
    r = requests.get(f"https://data.rcsb.org/rest/v1/core/branched_entity_instance/{entry_id}/{entity_instance_index}")
    data = json.loads(r.text)
    if dump: 
        json.dump(data, open(f"{data_dir}/{entity_instance_id}.json", "w"), indent=2)
    return data

def download_struc(pdb_code, format, data_dir, verbose=True):
    if os.path.exists(f"{data_dir}/{pdb_code}.{format}"):
        if verbose:
            print(f">>> {pdb_code}.{format} already exists in {data_dir}")
        return
    if verbose:
        print(f">>> Downloading {pdb_code}.{format}")
    r = requests.get(f"https://files.rcsb.org/download/{pdb_code}.{format}")
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
    with open(os.path.join(f"{data_dir}", f"{pdb_code}.{format}"), "wb") as f:
        f.write(r.content)
    if verbose:
        print(f">>> {pdb_code}.{format} downloaded to {data_dir}")

def download_struc_entity_info(pdb_code, data_dir):
    '''
    Download entity information for a PDB entry
    '''
    if os.path.exists(f"{data_dir}/{pdb_code}.json"):
        print(f">>> {pdb_code}.json already exists in {data_dir}")
        return
    print(f">>> Downloading information of {pdb_code}...")
    res = {}
    entry = get_entry(pdb_code)
    for entity_id in entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        print(f"Processing {pdb_code}_{entity_id}")
        entity = get_polymer_entity(f"{pdb_code}_{entity_id}")
        res[entity_id] = {
            "chains": entity["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"],
            "sequence": entity["entity_poly"]["pdbx_seq_one_letter_code_can"],
            "annotations": []
        }
        try:
            for annotation in entity["rcsb_polymer_entity_annotation"]:
                res[entity_id]["annotations"].append(annotation["name"])
        except KeyError:
            pass
    json.dump(res, open(f"{data_dir}/{pdb_code}.json", "w"), indent=2)
    print(f">>> {pdb_code}.json downloaded to {data_dir}")

def breakdown_entry_to_entity_accession_ids(pdb_code: str, verbose=True):
    '''
    This function focuses on an entry. All entity instances in the entry are considered.
    A dictionary is returned, with keys being entity IDs and values being the corresponding (usually) UniProt accession IDs.
    Returns: dict[str<chain_id>, str<uniprot_id>]
    - asym_ids: segment ids, 用于 RCSB 中检索 Instances
    - auth_asym_ids: chain ids, 不可用于 RCSB 中检索 Instances
    '''
    
    entry_data = get_entry(pdb_code)
    if "polymer_entity_ids" in entry_data["rcsb_entry_container_identifiers"]:
        polymer_entity_ids: list = entry_data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    else:
        polymer_entity_ids = []
    if "non_polymer_entity_ids" in entry_data["rcsb_entry_container_identifiers"]:
        non_polymer_entity_ids: list = entry_data["rcsb_entry_container_identifiers"]["non_polymer_entity_ids"]
    else:
        non_polymer_entity_ids = []
    if "branched_entity_ids" in entry_data["rcsb_entry_container_identifiers"]:
        branched_polymer_entity_ids: list = entry_data["rcsb_entry_container_identifiers"]["branched_entity_ids"]
    else:
        branched_polymer_entity_ids = []
    all_entity_ids = []
    for entity_id in polymer_entity_ids:
        all_entity_ids.append((0, f"{pdb_code}_{entity_id}"))
    for entity_id in non_polymer_entity_ids:
        all_entity_ids.append((1, f"{pdb_code}_{entity_id}"))
    for entity_id in branched_polymer_entity_ids:
        all_entity_ids.append((2, f"{pdb_code}_{entity_id}"))
    all_entity_ids.sort(key=lambda x: int(x[1].split("_")[1]))
    res = {}
    for entity_type, entity_id in all_entity_ids:
        if verbose: print(f">>> Analyzing {entity_id}")
        res_entry = {
                "asym_ids": [],
                "auth_asym_ids": [],
                "database_accessions": [],
            }
        if entity_type == 0:
            entity_data = get_polymer_entity(entity_id)
            res_entry["asym_ids"] = entity_data["rcsb_polymer_entity_container_identifiers"]["asym_ids"]
            res_entry["auth_asym_ids"] = entity_data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
            if "reference_sequence_identifiers" in entity_data["rcsb_polymer_entity_container_identifiers"]:
                for database in entity_data["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"]:
                    res_entry["database_accessions"].append(database["database_name"] + ":" + database["database_accession"])
            res_entry["database_accessions"].append(entity_data['rcsb_polymer_entity']['pdbx_description'])
        elif entity_type == 1:
            entity_data = get_nonpolymer_entity(entity_id)
            res_entry["asym_ids"] = entity_data["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"]
            res_entry["auth_asym_ids"] = entity_data["rcsb_nonpolymer_entity_container_identifiers"]["auth_asym_ids"]
            res_entry["database_accessions"].append(entity_data["rcsb_nonpolymer_entity_container_identifiers"]["chem_ref_def_id"])
        elif entity_type == 2:
            entity_data = get_branched_entity(entity_id)
            res_entry["asym_ids"] = entity_data["rcsb_branched_entity_container_identifiers"]["asym_ids"]
            res_entry["auth_asym_ids"] = entity_data["rcsb_branched_entity_container_identifiers"]["auth_asym_ids"]
            res_entry["database_accessions"] = entity_data['rcsb_branched_entity_container_identifiers']['chem_comp_monomers']
        else:
            raise ValueError("Invalid entity type")
        res[entity_id] = res_entry
    return res