from biotite.structure.io import pdb, pdbx
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
import os, json

def read_struc(struc_file):
    # 已经确认了里面所有的结构都只有 2 条链, A 是 target, B 是 binder
    if struc_file.endswith(".pdb"):
        my_atom_array = pdb.get_structure(pdb.PDBFile.read(struc_file))[0]
        res = {}
        for i in bio_struc.get_chains(my_atom_array):
            seq_i = "".join(
                list(
                    map(
                        lambda x: bio_seq.ProteinSequence.convert_letter_3to1(x),
                        bio_struc.get_residues(my_atom_array[my_atom_array.chain_id == i])[1]
                    )
                )
            )
            res[i] = seq_i
    elif struc_file.endswith(".cif"):
        my_atom_array = pdbx.get_structure(pdbx.CIFFile.read(struc_file))[0]
        res = {}
        for i in bio_struc.get_chains(my_atom_array):
            seq_i = "".join(
                list(
                    map(
                        lambda x: bio_seq.ProteinSequence.convert_letter_3to1(x),
                        bio_struc.get_residues(my_atom_array[my_atom_array.chain_id == i])[1]
                    )
                )
            )
            res[i] = seq_i
    return res

def generate_protein_chain(sequence):
    return {
        "proteinChain": {
            "sequence": sequence,
            "count": 1
        }
    }

def generate_entry(name, sequences, modelSeeds=[]):
    return {
        "name": name,
        "modelSeeds": modelSeeds,
        "sequences": sequences,
    }

def generate_af3_json(struc_marker: str, sort_key=lambda x: x["name"], recognize_format_in_dir='pdb'):
    '''
    根据 struc_marker 生成提交 AF3 预测的 JSON 文件. struc_marker 可以是多种类型
    - 如果是一个目录, 则会读取目录下所有的 recognize_format 的文件, 每个 pdb 文件对应一个 entry
    - 如果是一个文件, 则会读取这个文件, 生成一个 entry
    regonize_format 是用来识别文件格式的, 默认是 pdb, 只有在 struc_marker 是目录时才会用到
    '''
    res = []
    if not os.path.exists(struc_marker):
        raise FileNotFoundError(f"File or directory {struc_marker} not found")
    if os.path.isdir(struc_marker):
        for filename in os.listdir(struc_marker):
            if filename.endswith(".pdb"):
                entry_sequence_dict = read_struc(os.path.join(struc_marker, filename))
                entry_sequence_list = [
                    generate_protein_chain(entry_sequence_dict["A"]),
                    generate_protein_chain(entry_sequence_dict["B"])
                ]
                entry = generate_entry(".".join(filename.split(".")[:-1]), entry_sequence_list)
                res.append(entry)
        res.sort(key=lambda x: x["name"])
        return res
    else:
        _, file_ext = os.path.splitext(struc_marker)
        if file_ext in ['.pdb', '.cif']:
            entry_sequence_dict = read_struc(struc_marker)
            entry_sequence_list = []
            for my_key in entry_sequence_dict:
                entry_sequence_list.append(generate_protein_chain(entry_sequence_dict[my_key]))
            entry = generate_entry(os.path.basename(struc_marker).split(".")[0], entry_sequence_list)
            return [entry]
        else:
            raise ValueError(f"Unrecognized file format {file_ext}")

def remove_predicted_AF3_jobs(af3_jobs: list, res_dir: str):
    '''
    这个函数用于删除已经预测过的 job, 以避免重复提交
    '''
    jobs_to_remove = []
    for i, job_request in enumerate(af3_jobs):
        job_name = job_request['name']
        job_dir = os.path.join(res_dir, f"{job_name}")
        if os.path.exists(job_dir):
            jobs_to_remove.append(i)
    for i in reversed(jobs_to_remove):
        af3_jobs.pop(i)

def dump_AF3_jobs(af3_jobs: list, prefix: str, job_num_per_file: int = 30):
    for i in range(0, len(af3_jobs), job_num_per_file):
        with open(f"{prefix}_{i//job_num_per_file}.json", "w") as f:
            json.dump(af3_jobs[i:i+job_num_per_file], f, indent=4)

def split_AF3_jobs(af3_jobs: list, job_num_per_file: int = 20):
    res = []
    for i in range(0, len(af3_jobs), job_num_per_file):
        res.append(af3_jobs[i:i+job_num_per_file])

def split_AF3_job_file(job_file: str, job_num_per_file: int = 20):
    '''
    这个函数用于将 AF3 的 job 文件分割成多个小文件, 每个文件包含 job_num_per_file 个 job
    通过这个函数, 你可以将众多任务进行拆分, 从而使用多个账号并行提交任务
    '''
    with open(job_file) as f:
        jobs = json.load(f)
    for i in range(0, len(jobs), job_num_per_file):
        with open(f"{job_file[:-5]}_{i//job_num_per_file}.json", "w") as f:
            json.dump(jobs[i:i+job_num_per_file], f, indent=4)