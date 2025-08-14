import os
import numpy as np
import biotite.structure as struc
from biotite.structure.io import pdbx
from biotite.application.dssp import DsspApp
from biotite.application.localapp import LocalApp
from .pre_processors import *
from .af_utils import *
# 用 dssp 为每个预测结构决定二级结构, 并导出到 mmcif 文件

class DsspAppMac(DsspApp):
    
    def run(self):
        in_file = pdbx.cif.CIFFile()
        pdbx.convert.set_structure(in_file, self._array)
        in_file.write(self._in_file)
        self._in_file.flush()
        self.set_arguments([self._in_file.name, self._out_file.name])
        LocalApp.run(self)
    
    @staticmethod
    def annotate_sse(atom_array, bin_path="mkdssp"):
        app = DsspAppMac(atom_array, bin_path)
        app.start()
        app.join()
        return app.get_sse()

def export_sse_cif(struc: struc.AtomArray):
    sse = struc.get_annotation("sse")
    chain_id_array = struc.get_annotation("chain_id")
    res_id_array = struc.get_annotation("res_id")
    sse_cif_file = pdbx.CIFFile()
    sse_cif_block = pdbx.CIFBlock()
    sse_cif_file.blocks.append(sse_cif_block)
    sse_cif_block["data_"] = "SSE"
    sse_cif_block["loop_"] = ["_struct_conf.conf_type_id", "_struct_conf.beg_label_asym_id", "_struct_conf.beg_label_seq_id", "_struct_conf.end_label_asym_id", "_struct_conf.end_label_seq_id", "_struct_conf.pdbx_PDB_helix_class"]
    for i, (chain_id, res_id, sse_type) in enumerate(zip(chain_id_array, res_id_array, sse)):
        if sse_type == "H":
            sse_type = "HELX_P"
        elif sse_type == "E":
            sse_type = "STRN_P"
        elif sse_type == "C":
            sse_type = "TURN_P"
        else:
            sse_type = "NA"
        sse_cif_block.add_row([sse_type, chain_id, res_id, chain_id, res_id, sse_type])
    return sse_cif_file

def assign_sse_to_struc(struc: struc.AtomArray, sse_cif_path=None):
    if not sse_cif_path:
        sse_per_res = DsspAppMac.annotate_sse(struc)
    else:
        sse_cif_file = pdbx.CIFFile.read(sse_cif_path)
        sse_per_res = sse_cif_file.block["dssp_struct_summary"]["secondary_structure"].as_array()
        # 将其中的 . 转换为 C
        sse_per_res = list(map(lambda x: "C" if x == "." else x, sse_per_res))
    struc.del_annotation("sse")
    struc.add_annotation("sse", 'U1')
    chain_id_array = struc.get_annotation("chain_id")
    res_id_array = struc.get_annotation("res_id")
    res_ss_index = -1
    identifier_initial = ('', '')
    for i, (chain_id, res_id) in enumerate(zip(chain_id_array, res_id_array)):
        if (chain_id, res_id) != identifier_initial:
            identifier_initial = (chain_id, res_id)
            res_ss_index += 1
        struc.sse[i] = sse_per_res[res_ss_index]
        
def test_equality(struc: struc.AtomArray, sse_cif_path=None):
    '''
    Test if the assigned sse from the exported mmCIF file is equal to the one obtained from DSSP
    '''
    assign_sse_to_struc(struc, sse_cif_path)
    sse_1 = struc.get_annotation("sse").copy()
    assign_sse_to_struc(struc)
    sse_2 = struc.get_annotation("sse").copy()
    print(np.equal(sse_1, sse_2).all())

def run_dssp_for_design(top_dir, identifier, design, verbose=False):
    af_file_prefix = generate_af_file_prefix(top_dir, identifier, design)
    design = format_design_name(design)
    for i in range(5):
        input_filepath = f"{af_file_prefix}_model_{i}.cif"
        annotated_cif_path = f"{af_file_prefix}_model_{i}_dssp.mmcif"
        if not os.path.exists(input_filepath):
            if verbose: print(f">>> No prediction found for {design}.")
            continue
        if not os.path.exists(annotated_cif_path):
            command = f"mkdssp --output-format mmcif {input_filepath} {annotated_cif_path}"
            if verbose: print(command)
            os.system(command)
        else:
            if verbose: print(f">>> Skipping because {annotated_cif_path} already exists.")

def read_res_wise_sse_for_prediction(top_dir, identifier, design, prediction_i):
    af_file_prefix = generate_af_file_prefix(top_dir, identifier, design)
    dssp_mmcif_path = f"{af_file_prefix}_model_{prediction_i}_dssp.mmcif"
    sse_cif_file = pdbx.CIFFile.read(dssp_mmcif_path)
    sse_per_res = sse_cif_file.block["dssp_struct_summary"]["secondary_structure"].as_array()
    # 将其中的 . 转换为 C
    sse_per_res = list(map(lambda x: "C" if x == "." else x, sse_per_res))
    return sse_per_res

def get_sse_ratio(top_dir, identifier, design, prediction_i, sse_list):
    sse_per_res = read_res_wise_sse_for_prediction(top_dir, identifier, design, prediction_i)
    res = {}
    possible_sse = ["G", "H", "I", "T", "E", "B", "S", "C"]
    for sse in sse_list:
        if not sse in possible_sse:
            raise ValueError(f"SSE should be in {possible_sse}, but got {sse}")
        res[sse] = sse_per_res.count(sse) / len(sse_per_res)
    return res

def get_helix_ratio(top_dir, identifier, design, prediction_i):
    res = get_sse_ratio(top_dir, identifier, design, prediction_i, ["G", "H", "I"])
    return sum(res.values()), res