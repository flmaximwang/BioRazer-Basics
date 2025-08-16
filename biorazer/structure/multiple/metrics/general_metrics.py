import json

from .. import *
from ..utils import af_utils, rosetta_utils, select_utils

class PartMetric(Metric):
    
    '''
    PartMetric is a metric that is generated from part of a structure
    '''
    
    def __init__(
        self,
        project,
        identifiers,
        metric_name="Alphafold metric",
        metric_files={
            "meta_data": "0/./fold_{sn}_full_data_{i}.json",
            "cif": "0/./fold_{sn}_model_{i}.cif",
            "conf_data": "0/./fold_{sn}_model_{i}.conf",
            "part_indices": "0/./fold_{sn}_part_indices_{i}.json"
        }
    ):
        super().__init__(
            project=project,
            identifiers=identifiers,
            metric_name=metric_name,
            metric_files=metric_files
        )
        
    def _generate_part_indices(self, sequence_name, sequence_index):
        '''
        Implement this method to define specific atom indices for which part to choose for pLDDT calculaion
        '''
        
    def generate_part_indices(self):
        my_project = self.get_project()
        sequence_names = my_project.get_sequence_names()
        for sequence_name in sequence_names:
            for i in range(5):
                self._generate_part_indices(sequence_name, i)

class PartAtomMetric(PartMetric):
    
    '''
    PartResidueMetrics are metrics that are generated from residues of a structure
    '''
    
    def __init__(
        self,
        project: Project,
        identifiers,
        metric_name = "Residue metric for specific atoms",
        metric_files={
            "meta_data": "0/./fold_{sn}_full_data_{i}.json",
            "cif": "0/./fold_{sn}_model_{i}.cif",
            "conf_data": "0/./fold_{sn}_model_{i}.conf",
            "part_indices": "0/./fold_{sn}_part_indices_{i}.json"
        }
    ):
        super().__init__(
            project=project,
            identifiers=identifiers,
            metric_name=metric_name,
            metric_files=metric_files
        )
    
    def get_part_atom_indices(self):
        part_residue_indices_file = self.get_metric_file_path("part_atom_indices")
        with open(part_residue_indices_file) as f:
            part_residue_indices_file = json.load(f)
        return part_residue_indices_file["part_atom_indices"]
        
class PartResidueMetric(PartMetric):
    
    '''
    PartAtomMetrics are metrics that are generated from atoms of a structure
    '''
    
    def __init__(
        self,
        project: Project,
        identifiers,
        metric_name = "Atom metric for specific atoms",
        metric_files={
            "meta_data": "0/./fold_{sn}_full_data_{i}.json",
            "cif": "0/./fold_{sn}_model_{i}.cif",
            "conf_data": "0/./fold_{sn}_model_{i}.conf",
            "part_residue_indices": "0/./fold_{sn}_part_residue_indices_{i}.json"
        }
    ):
        super().__init__(
            project=project,
            identifiers=identifiers,
            metric_name=metric_name,
            metric_files=metric_files
        )
    
    def get_part_residue_indices(self, sequence_name, sequence_index):
        part_residue_indices_file = self.get_metric_file_path("part_residue_indices", sequence_name, sequence_index)
        with open(part_residue_indices_file) as f:
            part_residue_indices_file = json.load(f)
        return part_residue_indices_file["part_residue_indices"]