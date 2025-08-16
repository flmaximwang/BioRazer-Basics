import requests
import sys
import time
import os
import shutil
import urllib

from collections import OrderedDict

import pandas as pd

from selenium import webdriver
from selenium.webdriver.remote.webdriver import WebDriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import *

from . import rosetta_utils
from . import af_utils
from . import pre_processors
from . import summary

class CASTpJob:
    
    def __init__(self, 
            pdb_path,
            probe_radius=1.4,
            seconds_to_wait_after_submitting: int=1,
            seconds_for_navigating_time_out: float=10,
            seconds_to_wait_after_navigating: float=0.5,
            seconds_to_wait_for_downloading: int=10,
            mode = "a",
        ):
        if not os.path.exists(pdb_path):
            raise FileNotFoundError(f"File {pdb_path} not found")
        if not pdb_path.endswith(".pdb"):
            raise ValueError("PDB file must end with .pdb")
        self.pdb_path = pdb_path
        self.res_zip_path = pdb_path[:-4] + ".zip"
        self.res_dir = pdb_path[:-4] + "_castp"
        self.res_name = os.path.basename(self.res_dir)
        if probe_radius < 0 or probe_radius > 10:
            raise ValueError("Probe radius must be between 0 and 10")
        self.probe_radius=str(probe_radius)
        if not isinstance(seconds_to_wait_after_submitting, int):
            raise TypeError("seconds_to_wait_after_submitting must be an integer")
        self.seconds_to_wait_after_submitting = seconds_to_wait_after_submitting
        self.seconds_for_navigating_time_out = seconds_for_navigating_time_out
        self.seconds_to_wait_after_navigating = seconds_to_wait_after_navigating
        if not isinstance(seconds_to_wait_for_downloading, int):
            raise TypeError("seconds_to_wait_for_downloading must be an integer")
        self.seconds_to_wait_for_downloading = seconds_to_wait_for_downloading
        self.mode = mode
        self.job_id = None
        self.driver: WebDriver = None
        self.job_path = None
        
    def test_connectivity(self):
        try:
            requests.get('http://sts.bioe.uic.edu/castp/index.html')
        except requests.exceptions.ConnectionError:
            raise requests.exceptions.ConnectionError("Error: No internet connection")
        
    def submit(self):
        file = open(self.pdb_path, "rb")
        self.probe_radius = str(self.probe_radius)
        url = 'http://sts.bioe.uic.edu/castp/submit_calc.php'
        headers = {
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:102.0) Gecko/20100101 Firefox/102.0',
            'Accept': '*/*',
            'Accept-Language': 'en-US,en;q=0.5',
            'Accept-Encoding': 'gzip, deflate',
            'X-Requested-With': 'XMLHttpRequest',
            'Origin': 'http://sts.bioe.uic.edu',
            'DNT': '1',
            'Connection': 'keep-alive',
            'Referer': 'http://sts.bioe.uic.edu/castp/calculation.html'
        }
        files = {
            'file': (self.pdb_path, file, 'application/x-aportisdoc'),
        }
        data = {
            'probe': self.probe_radius,
            'email': 'null',
        }

        self.test_connectivity()
        print(f">>> Submitting {self.pdb_path} to CASTp server...")
        response = requests.post(url, headers=headers, data=data, files=files)
        if response.status_code != 200:
            raise Exception(f"Error: {response.status_code}")
        
        self.job_id = str(response.text)
        print(f">>> Job {self.job_id} successfully submitted to CASTp")
        file.close()
        
        for i in range(self.seconds_to_wait_after_submitting):
            print(f">>> Waiting {self.seconds_to_wait_after_submitting - i} seconds for the job to be processed on the server...", end="\r")
            time.sleep(1)
        print("")
        self.generate_job_path()
    
    def generate_job_path(self):
        self.job_path =  os.path.join(os.path.expanduser("~"), "Downloads", f"{self.job_id}.zip")

    @staticmethod
    def page_is_loaded(driver: WebDriver):
        document_completed = driver.execute_script("return document.readyState") == "complete"
        if not document_completed: return False
        try:
            driver.find_element(By.XPATH, '//*[@id="web_header"]/div[1]/img')
            return True
        except NoSuchElementException:
            return False

    def navigate_to_job(self):
        print(f">>> Navigating to job {self.job_id}...")
        self.test_connectivity()
        while True:
            self.driver.get(f'http://sts.bioe.uic.edu/castp/index.html?{self.job_id}')
            driver_wait = WebDriverWait(self.driver, self.seconds_for_navigating_time_out)
            driver_wait.until(self.page_is_loaded)
            time.sleep(self.seconds_to_wait_after_navigating)
            try:
                toast_container = self.driver.find_element(By.XPATH, '//*[@id="toast-container"]')
                job_finished = True
                for toast in toast_container.find_elements(By.XPATH, 'div'):
                    if 'Cannot find' in toast.text:
                        job_finished = False
                if not job_finished:
                    print(f">>> Job {self.job_id} in progress, refreshing...")
                    time.sleep(1)
                    continue
                else: 
                    break
            except NoSuchElementException:
                break
        print(f">>> Job {self.job_id} finished...")

    def job_is_downloaded(self):
        return os.path.exists(self.job_path)
    
    def download_res(self):
        print(f">>> Downloading results for job {self.job_id}")
        js = 'document.getElementById("downloadbtn").click();'
        self.driver.execute_script(js)
        for i in range(self.seconds_to_wait_for_downloading):
            if self.job_is_downloaded():
                print(f">>> Results for job {self.job_id} successfully downloaded")
                break
            time.sleep(1)

    def unzip_job(self):
        os.system(f"unzip {self.job_path} -d {self.res_dir} &>/dev/null")
        os.remove(self.job_path)
        print(f">>> Results for job {self.job_id} unzipped")
        for filename in os.listdir(self.res_dir):
            new_filename = filename.replace(self.job_id, self.res_name)
            os.rename(os.path.join(self.res_dir, filename), os.path.join(self.res_dir, new_filename))
        
    def download_res_for_job_to_path(self):
        if os.path.exists(self.job_path):
            print(f">>> Results for job {self.job_id} already downloaded, skipping...")
            return
        print(f">>> Downloading results for job {self.job_id}...")
        self.driver = webdriver.Safari()
        self.driver.set_window_position(0, 0)
        self.driver.set_window_size(1080, 1415)
        self.navigate_to_job()
        self.download_res()
        self.driver.close()
        self.unzip_job()

    def run(self):
        if os.path.exists(self.res_dir) and self.mode == "a":
            print(f">>> Results for {self.pdb_path} already calculated, skipping...")
            return
        self.submit()
        self.download_res_for_job_to_path()

def get_pocket_composition(castp_res_dir):
    '''
    Parse .poc file from CASTP web server.  Lines are PDB format with pocket number at end:

    ATOM    165  NH2 ARG A  45      -9.395 -35.779   0.675  1.00 32.68   2  POC
    ATOM    194  O   VAL B  27     -12.428  -9.493   5.152  1.00 13.76   1  POC
    ATOM    196  CG1 VAL B  27     -11.788  -7.554   2.677  1.00 14.96   1  POC
    This script is modified from https://rbvi.github.io/chimerax-recipes/castp/castp.html
    Return: {pocket_id: [{"chain_id": chain_id, "res_num": res_num, "atom_name": atom_name}, ...]}
    '''

    poc_file_path = os.path.join(castp_res_dir, f"{os.path.basename(castp_res_dir)}.poc")
    with open(poc_file_path, 'r') as f:
        lines = f.readlines()

    pocket_composition: OrderedDict[int: list[dict]] = {}
    for line in lines:
        if line.startswith('ATOM  '):
            atom_name, chain_id, res_num = line[12:16].strip(), line[21:22], int(line[22:26])
            atom_dict = {'chain_id': chain_id, 'res_num': res_num, 'atom_name': atom_name}
            pocket_num = int(line[66:71])
            if pocket_num in pocket_composition:
                pocket_composition[pocket_num].append(atom_dict)
            else:
                pocket_composition[pocket_num] = [atom_dict]

    pocket_composition = OrderedDict(sorted(pocket_composition.items(), key=lambda x: x[0]))
    return pocket_composition

def get_interface_pocket_ids(castp_res_dir, interface_annotated_pdb_path):
    pocket_composition = get_pocket_composition(castp_res_dir)
    interface_residues = rosetta_utils.get_interface_residues(interface_annotated_pdb_path)
    interface_pockets = []
    for pocket_id, atom_identifiers in pocket_composition.items():
        # 只有当组成 pocket 的所有残基都在界面上 (位于 interface_residues 中) 时，才认为这个 pocket 是界面 pocket
        is_interface_pocket = True
        for atom_identifier in atom_identifiers:
            if (atom_identifier['chain_id'], atom_identifier['res_num']) not in interface_residues:
                is_interface_pocket = False
                break
        if is_interface_pocket:
            interface_pockets.append(pocket_id)
    return interface_pockets

def get_interface_pocket_composition(castp_res_dir):
    pocket_composition = get_pocket_composition(castp_res_dir)
    interface_pockets = []
    for pocket_id, atom_identifiers in pocket_composition.items():
        chains = set([atom['chain_id'] for atom in atom_identifiers])
        if len(chains) > 1:
            interface_pockets.append(pocket_id)
    interface_pockets.sort()
    interface_pocket_composition = OrderedDict()
    for pocket_id in interface_pockets:
        interface_pocket_composition[pocket_id] = pocket_composition[pocket_id]
        
    return interface_pocket_composition

def get_pocket_info(castp_res_dir):
    '''
    Parse .pocInfo file from CASTP web server.  This file looks like
    POC:	Molecule	ID	N_mth	Area_sa	Area_ms	Vol_sa	Vol_ms	Lenth	cnr
    POC:	j_673c86babd4d6	1	3	3866.369	3991.202	12090.182	17541.664	2344.235	885
    POC:	j_673c86babd4d6	2	1	205.816	292.521	164.714	514.249	181.495	86
    POC:	j_673c86babd4d6	3	1	75.203	159.931	41.722	188.781	68.046	34
    POC:	j_673c86babd4d6	4	1	4.373	5.662	37.361	44.936	7.750	2
    POC:	j_673c86babd4d6	5	2	65.288	135.212	32.516	149.889	61.813	37
    POC:	j_673c86babd4d6	6	3	98.594	258.104	26.776	230.779	93.897	41
    '''
    pocket_info_file_path = os.path.join(castp_res_dir, f"{os.path.basename(castp_res_dir)}.pocInfo")
    pocket_info = pd.read_csv(pocket_info_file_path, sep='\t')
    return pocket_info

def get_interface_pocket_info(castp_res_dir, interface_annotated_pdb_path):
    pocket_info = get_pocket_info(castp_res_dir)
    interface_pocket_ids = get_interface_pocket_ids(castp_res_dir, interface_annotated_pdb_path)
    interface_pocket_info = pocket_info[pocket_info["ID"].isin(interface_pocket_ids)]
    return interface_pocket_info
    
def run_castp_for_prediction(top_dir, identifier, design, probe_radius=1.4, mode="a"):
    '''
    Run CASTp for a specific prediction. The results will be saved in the prediction directory.
    '''
    design = pre_processors.format_design_name(design)
    af_file_prefix = af_utils.generate_af_file_prefix(top_dir, identifier, design)
    for i in range(5):
        pdb_path = f"{af_file_prefix}_model_{i}.pdb"
        castp_job = CASTpJob(pdb_path, probe_radius=probe_radius, mode=mode)
        castp_job.run()



class SolventAccessibleVolumeForInterfaceEmptySpaceMetric(summary.Metric):
    
    def __init__(
        self,
        parent_prediction_collection_summary: summary.PredictionCollectionSummary,
        related_prediction_cateogory: str,
        metric_name = "Solvent Accessible Volume For Interface Empty Space"
    ):
        assert isinstance(parent_prediction_collection_summary, summary.PredictionCollectionSummary)
        assert isinstance(related_prediction_cateogory, str)
        super().__init__(
            parent_prediction_collection_summary=parent_prediction_collection_summary,
            related_prediction_categories=[related_prediction_cateogory],
            metric_name=metric_name,
            metric_files={
                "castp_res_dir": "0/./_model_{i}_castp",
                "mouth_atoms": "0/./_model_{i}_castp/_model_{i}_castp.mouth",
                "mouth_info": "0/./_model_{i}_castp/_model_{i}_castp.mouthInfo",
                "pdb": "0/./_model_{i}_castp/_model_{i}_castp.mouthInfo",
                "pocket_atoms": "0/./_model_{i}_castp/_model_{i}_castp.poc",
                "pocket_info": "0/./_model_{i}_castp/_model_{i}_castp.pocInfo",
                "spheres": "0/./_model_{i}_castp/_model_{i}_castp.spheres.json",
                "interface_annotated_pdb": "0/./_model_{i}_rosetta_interface/_model_{i}_0001.pdb"
            }
        )
        
    def _summarize_metric_vector_for_prediction_i(self, prediction_name, i):
        castp_res_dir = self.get_metric_file_path("castp_res_dir", prediction_name, i)
        interface_annotated_pdb_path = self.get_metric_file_path("interface_annotated_pdb", prediction_name, i)
        interface_pocket_info = get_interface_pocket_info(castp_res_dir, interface_annotated_pdb_path)
        return interface_pocket_info["Vol_sa"].sum()

class dSASAAveragedSolventAccessibleVolumeForInterfaceEmptySpaceMetric(summary.Metric):

    def __init__(
        self,
        parent_prediction_collection_summary: summary.PredictionCollectionSummary,
        related_prediction_cateogory: str,
        metric_name = "dSASA Averaged Solvent Accessible Volume For Interface Empty Space"
    ):
        assert isinstance(parent_prediction_collection_summary, summary.PredictionCollectionSummary)
        assert isinstance(related_prediction_cateogory, str)
        super().__init__(
            parent_prediction_collection_summary=parent_prediction_collection_summary,
            related_prediction_categories=[related_prediction_cateogory],
            metric_name=metric_name,
            metric_files={
                "castp_res_dir": "0/./_model_{i}_castp",
                "mouth_atoms": "0/./_model_{i}_castp/_model_{i}_castp.mouth",
                "mouth_info": "0/./_model_{i}_castp/_model_{i}_castp.mouthInfo",
                "pdb": "0/./_model_{i}_castp/_model_{i}_castp.mouthInfo",
                "pocket_atoms": "0/./_model_{i}_castp/_model_{i}_castp.poc",
                "pocket_info": "0/./_model_{i}_castp/_model_{i}_castp.pocInfo",
                "spheres": "0/./_model_{i}_castp/_model_{i}_castp.spheres.json",
                "rosetta_interface_score": "0/./_model_{i}_rosetta_interface.sc",
            }
        )
    
    def _summarize_metric_vector_for_prediction_i(self, prediction_name, i):
        my_df = self.get_parent_prediction_summary()
        top_dir = my_df.get_top_dir()
        castp_res_dir = self.get_metric_file_path("castp_res_dir", prediction_name, i)
        interface_pocket_info = get_interface_pocket_info(castp_res_dir)
        interface_vol_sa =  interface_pocket_info["Vol_sa"].sum()
        interface_dsasa = rosetta_utils.read_interface_sc(top_dir, "Complex", prediction_name, i)["dSASA_int"]
        return interface_vol_sa / interface_dsasa