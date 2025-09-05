import json
import numpy as np

class MetaFile:
    
    def __init__(self, meta_file_path):
        self.meta_file_path = meta_file_path
        with open(self.meta_file_path) as f:
            self.meta_data = json.load(f)
        for key in self.meta_data:
            self.meta_data[key] = np.array(self.meta_data[key])
            
    def get_data(self):
        '''
        Alphafold meta files contain 
        - atom_chain_ids (length of all your atoms)
        - atom_plddts (length of all your atoms)
        - pae (length of your residues, 2d matrix)
        - token_chain_ids (length of your protein)
        - tok (en_res_ids (length of your protein)
        '''
        return self.meta_data
    
    def __getitem__(self, key):
        return self.get_data()[key]

class ConfFile: 
    
    def __init__(self, conf_file_path):
        self.conf_file_path = conf_file_path
        with open(self.conf_file_path) as f:
            self.conf_data = json.load(f)
        self.conf_data['chain_pair_iptm'] = np.array(self.conf_data['chain_pair_iptm'])
        self.conf_data['chain_pair_pae_min'] = np.array(self.conf_data['chain_pair_pae_min'])
        
    def get_data(self):
        '''
        Alphaffold conf files contain
        - 'chain_iptm': length of your chains, 2d matrix, interchain pTM
        - 'chain_pair_iptm',
        - 'chain_pair_pae_min',
        - 'chain_ptm',
        - 'fraction_disordered',
        - 'has_clash',
        - 'iptm',
        - 'num_recycles',
        - 'ptm',
        - 'ranking_score']
        - 
        '''
        return self.conf_data

    def __getitem__(self, key):
        return self.get_data()[key]