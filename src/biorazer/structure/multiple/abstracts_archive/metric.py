import os
import operator
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from biotite.structure.io import pdbx
# from .pre_processors import *
from ..utils import af_utils
from ..abstracts import StructureStack

class Metric:
    '''
    A metric is a class that helps you to calculate and summarize a metric for StructureStack.
    
    To give a specific metric, you need to inherit this class and implement the following methods:
    - generate_metric_files(self)
    - _summarize_metric_vector_for_prediction_i(self, prediction_name, i)
    
    To use a metric, your first create it by attach it to
    - your project
    - the identifiers of structures related to the metric. For example, 
    inteface related metrics are related to "complex" structures, by not to "monomer" structures.
    
    :param metric_files: A dictionary that contains the file paths of the metric files. it looks like
    {key: "{0}/fold_{sn}_model_{i}_castp"}.
    - {0} is replaced with the directory path of the 1st structure stack you give in generate_metric_files,
    - {sn} is replaced with the folder name of the structure stack,
    - {i} is replaced with the indices.
    '''
    
    def __init__(
            self,
            metric_name: str,
            metric_files: dict[str, str],
            indices = [0, 1, 2, 3, 4]
        ):
        self.metric_name = metric_name
        self.metric_files = metric_files
        self.indices = indices
        
    def generate_metric_files(
            self,
            structure_stacks: list[StructureStack],
        ):
        pass
        
    def update_metric_to_summary(self):
        my_df = self.get_project()
        for subscript in ["mean", "std", "0", "1", "2", "3", "4"]:
            actual_col = f"{self.metric_name}, {subscript}"
            if not actual_col in my_df.columns:
                my_df[actual_col] = np.nan
    
    def remove_metric_from_summary(self):
        my_df = self.get_project()
        for subscript in ["mean", "std", "0", "1", "2", "3", "4"]:
            actual_col = f"{self.metric_name}, {subscript}"
            if actual_col in my_df.columns:
                my_df.drop(actual_col, axis=1, inplace=True)
        
    def get_project(self):
        return self.project
    
    def get_metric_index(self):
        return self.metric_index
    
    def get_metric_name(self):
        return self.metric_name
    
    def get_metric_categories(self):
        return list(self.metric_files.keys())
    
    def check_metric_file_path_pattern(self, metric_file_path_pattern):
        if not isinstance(metric_file_path_pattern, str):
            raise ValueError("metric_file_path_pattern should be a string")
        metric_file_path_parts = metric_file_path_pattern.split("/")
        if metric_file_path_parts[0] not in self.identifiers:
            raise ValueError(f"identifier should be one of {self.identifiers}")
        if metric_file_path_parts[1] != ".":
            raise ValueError("The second part of the metric file path should be '.'")
        if not "{sn}" in metric_file_path_pattern:
            raise ValueError("The metric file path should contain '{sn}'")
        if not "{i}" in metric_file_path_pattern:
            raise ValueError("The metric file path should contain '{i}'")
    
    def parse_metric_file_path_pattern(self,
            metric_file_path_pattern: str,
            struc_stacks: list[StructureStack], 
            sequence_index,
        ):
        """
        Parses the metric file path pattern and replaces placeholders with actual values.
        Args:
            metric_file_path_pattern (str): The pattern of the metric file path containing placeholders.
            sequence_name (str): The name of the sequence to replace the placeholder.
            sequence_index (int): The index of the sequence to replace the placeholder.
        Returns:
            str: The fully constructed metric file path with all placeholders replaced by actual values.
        """
        
        metric_file_path_parts = metric_file_path_pattern.split("/")
        struc_stack_names = [struc_stack.get_name() for struc_stack in struc_stacks]
        structure_batch_dir = self.get_project().get_top_dir()
        metric_file_path_parts[0] = self.identifiers[int(metric_file_path_parts[0])]
        metric_file_path_parts[1] = sequence_name 
        for i in range(2, len(metric_file_path_parts)):
            a = metric_file_path_parts[i]
            b = a.replace("{sn}", sequence_name)
            c = b.replace("{i}", str(sequence_index))
            metric_file_path_parts[i] = c
        res_path_list = [structure_batch_dir]
        res_path_list.extend(metric_file_path_parts)
        return os.path.join(*res_path_list)
    
    def get_metric_file_path(self, metric_file_category, sequence_name, sequence_index):
        
        metric_file_path_pattern = self.metric_files[metric_file_category]
        return self.parse_metric_file_path_pattern(metric_file_path_pattern, sequence_name, sequence_index)
    
    def change_metric_file_path_pattern(self, old_metric_name, new_metric_file_path_pattern):
        for sequence_name in self.get_project().get_sequence_names():
            for i in range(5):
                old_metric_file_path = self.get_metric_file_path(old_metric_name, sequence_name, i)
                new_metric_file_path = self.parse_metric_file_path_pattern(new_metric_file_path_pattern, sequence_name, i)
                os.rename(old_metric_file_path, new_metric_file_path)
        self.metric_files[old_metric_name] = new_metric_file_path_pattern
    
    def open_metric_file(self, metric_file_category, prediction_name, prediction_i, mode="r"):
        metric_file_path = self.get_metric_file_path(metric_file_category, prediction_name, prediction_i)
        return open(metric_file_path, mode)
    
    def check_summarization_mode(self, mode):
        if not mode in ["r", "w", "a"]:
            raise ValueError("mode should be 'r' for read, 'w' for write or 'a' for append")
    
    def _msg_for_summary(self, sequence_name, mode):
        """
        This function prints messages related to analyzing a design and calculating or skipping a metric
        based on the mode provided.
        
        :param design: The `design` parameter represents the specific design that is being analyzed. It
        could refer to a particular model, layout, blueprint, or any other design element that is being
        evaluated
        :param metric_index: Metric_index is the index of the metric being analyzed or calculated. It is
        used to access the metric name from the self.design_metrics list for displaying purposes in the function
        :param mode: The `mode` parameter is used to determine the action to be taken in the function. It
        can have two possible values: 0 for calculation, 1 for skipping
        """
        print(f">>> Analyzing design {sequence_name}")
        if mode == 0:
            print(f">>> Calculating metric: {self.metric_name}")
        else:
            print(f">>> Skipping calculated metric: {self.metric_name}")
            
    def _set_metric_vector_for_prediction(self, prediction_name, vector):
        my_df = self.get_project()
        for i, vector_i in enumerate(vector):
            my_df.loc[my_df.get_sequence_names() == prediction_name, f"{self.metric_name}, {i}"] = vector_i
        my_df.loc[my_df.get_sequence_names() == prediction_name, f"{self.metric_name}, mean"] = np.mean(vector)
        my_df.loc[my_df.get_sequence_names() == prediction_name, f"{self.metric_name}, std"] = np.std(vector)
        
    def get_metric_vector_for_design(self, design):
        my_df = self.get_project()
        return np.array(my_df.loc[my_df.get_sequence_names() == design, [f"{self.metric_name}, {i}" for i in range(5)]], dtype=np.float64)
    
    def get_metric_mean_for_design(self, design):
        my_df = self.get_project()
        return my_df.loc[my_df.get_sequence_names() == design, f"{self.metric_name}, mean"].values[0]
    
    def get_metric_std_for_design(self, design):
        my_df = self.get_project()
        return my_df.loc[my_df.get_sequence_names() == design, f"{self.metric_name}, std"].values[0]
    
    def _summarize_metric_vector_for_i(self, sequence_name, sequence_index):
        '''
        Implement this method to summarize the i-th model in a prediction.
        Just return the metric value of the i-th model
        '''
        self._msg_for_summary(sequence_name=sequence_name, mode=0)
    
    def summarize_metric_vector_for_prediction(self, prediction_name):
        vector_length = 5
        metric_vector = np.zeros(vector_length)
        for i in range(5):
            metric_vector[i] = self._summarize_metric_vector_for_i(prediction_name, i)
        self._set_metric_vector_for_prediction(prediction_name, metric_vector)
    
    def summarize_based_on_mode(self, mode="a"):
        self.check_summarization_mode(mode)
        self.update_metric_to_summary()
        my_df = self.get_project()
        prediction_names = my_df.get_sequence_names()
        if mode == "w":
            for prediction_name in prediction_names:
                self._msg_for_summary(prediction_name, 0)
                self.summarize_metric_vector_for_prediction(prediction_name)
        elif mode == "a":
            for prediction_name in prediction_names:
                row_selector = prediction_names == prediction_name
                col_selector = [f"{self.metric_name}, {i}" for i in ["mean", "std", 0, 1, 2, 3, 4]]
                values_accessed = my_df.loc[row_selector, col_selector].values
                if np.isnan(values_accessed).any():
                    self._msg_for_summary(prediction_name, 0)
                    self.summarize_metric_vector_for_prediction(prediction_name)
                else:
                    self._msg_for_summary(prediction_name, 1)
        else:
            raise NotImplementedError(f"Read mode '{mode}' not implemented")
    
    def get_grouped_bar_x(self, group_i, group_len, index_i, barwidth=1):
        index_i = float(index_i)
        return (group_i * (group_len + 0.5) + index_i) * barwidth
    
    
    def generate_metric_files(self):
        '''
        Run apps to generate metric files
        '''
        pass