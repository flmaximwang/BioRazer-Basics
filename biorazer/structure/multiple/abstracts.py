import os, re
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class StructureStack:
    
    def __init__(self,
        struc_dir: str,
    ):
        self.struc_dir = struc_dir
        self.metrics = {}
    
    def __setitem__(self, key, value):
        self.metrics[key] = value
        
    def __getitem__(self, key):
        return self.metrics[key]
    
    def __repr__(self):
        return repr(self.struc_dir)
    
    def get_name(self):
        return os.path.split(self.struc_dir)[-1]

class Metric:
    '''
    A metric is a class that helps you to calculate and summarize a metric for StructureStack.
    
    :param metric_name: The name of the metric.
    
    :param metric_files: A dictionary that contains the file paths of the metric files. it looks like
    {key: "{0}/fold_{name}_model_{i}.cif"}.
    - {0} is replaced with the directory path of the 1st structure stack you give in generate_metric_files,
    - {name} is replaced with the folder name of the structure stack,
    - {i} is replaced with the indices.
    
    :param indices: The indices of multiple structures. Usually used for AF predicted structures
    '''
    
    def __init__(
            self,
            metric_name: str,
            input_file_patterns: dict[str, str],
            output_file_patterns: dict[str, str],
            struc_indices = [0, 1, 2, 3, 4]
        ):
        self.metric_name = metric_name
        self.input_file_patterns = input_file_patterns
        self.output_file_patterns = output_file_patterns
        self.struc_indices = struc_indices
        
    def generate_metric_files(
            self,
            structure_stacks: list[StructureStack],
        ):
        pass
    
    def pattern_to_path(
        self,
        path_pattern: str,
        struc_stacks: list[StructureStack], 
        struc_index: int,
    ):
        """
        
        Args:
        - file_path_pattern: The pattern of the file path containing placeholders.
        - struc_stacks: The list of structure stacks to replace the placeholders.
        - struc_index: The index of the structure stack to replace the placeholder.
        
        Returns:
        - The fully constructed file path with all placeholders replaced by actual values.
        """
        
        # Check file path pattern
        if not isinstance(path_pattern, str):
            raise ValueError("file_path_pattern should be a string")
        path_parts = path_pattern.split("/")
        
        struc_stack_index = int(re.match(r"\{(\d+)\}", path_parts[0])[1])
        struc_stack_names = [struc_stack.get_name() for struc_stack in struc_stacks]
        struc_stack_dir = struc_stacks[struc_stack_index].struc_dir
        struc_stack_name = struc_stack_names[struc_stack_index]
        
        path_parts[0] = struc_stack_dir
        for i in range(1, len(path_parts)):
            a = path_parts[i]
            b = a.replace("{name}", struc_stack_name)
            c = b.replace("{i}", str(struc_index))
            path_parts[i] = c
        return os.path.join(*path_parts)
    
    def get_input_file_path(
        self,
        input_file_key: str,
        struc_stacks: list[StructureStack], 
        struc_index: int,
    ):
        
        path_pattern = self.input_file_patterns[input_file_key]
        return self.pattern_to_path(path_pattern, struc_stacks, struc_index)
    
    def get_output_file_path(
        self,
        output_file_key: str,
        struc_stacks: list[StructureStack], 
        struc_index: int,
    ):
        
        path_pattern = self.output_file_patterns[output_file_key]
        return self.pattern_to_path(path_pattern, struc_stacks, struc_index)
    
    def calculate(self, struc_stacks: list[StructureStack]):
        pass
    
    def open_metric_file(self, metric_file_category, prediction_name, prediction_i, mode="r"):
        metric_file_path = self.get_input_files(metric_file_category, prediction_name, prediction_i)
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

class Entry:
    '''
    An entry is a collection of structure stacks and metrics.
    '''
    
    def __init__(self, struc_stack_dict: dict[str, StructureStack] = {}):
        self.struc_stack_dict = struc_stack_dict
        self.metrics: dict[str, ] = {}
        
    def __repr__(self):
        return repr(self.struc_stack_dict)
    
    def add_struc_stack(self, struc_type, struc_dir):
        self.struc_stack_dict[struc_type] = StructureStack(struc_dir)
        
    def remove_struc_stack(self, struc_type):
        self.struc_stack_dict.pop(struc_type)
    
    def calculate_metric(self, metric: Metric, struc_types: list[str], overwrite=False):
        res = metric.calculate([self.struc_stack_dict[struc_type] for struc_type in struc_types], overwrite=False)
        self.metrics[metric.metric_name] = {"method": metric, "values": res}

class Project:
    '''
    A project collects a series of protein sequences and their structures, and the metrics calculated from sequences and structures.
    
    Data structure:
    - top_dir
        -  summary_filename (.csv file)
        -  Binder (identifier name)
            -  ncam_d1_binder_l108_s51378_mpnn6
        - Complex
            - ncam_d1_binder_l108_s51378_mpnn6
    
    You have a bundle of protein sequences stored in top_dir/summary_filename (a .csv file)
    - Names of your sequences are stored in summary_file[prediction_name_col] (for example, the 'Design' column)
    - Structures of your sequences are obtained with different methods, including AlphaFold, Rosetta, experiments, etc.
    They are classified by identifiers (for example, 'Binder', 'Complex') and stored in top_dir/<identifier>/<sequence_name>
    Usually multiple structures are generated for each sequence, and they are stored in the same folder for easy access
    - Structures are further analyzed with different apps to generate files stored in each prediction.
    
    To 
    '''
    
    def __init__(
        self,
        project_dir,
        summary_filename="Default_Project",
        struc_types=["Default"],
        exclude_dirs = [".DS_Store"]
    ):
        self.project_dir = project_dir
        self.summary_filename = summary_filename
        summary_file_path = self.get_summary_file_path()
        if not os.path.exists(summary_file_path):
            self.data = pd.DataFrame({})
        else:
            self.data = pd.read_csv(summary_file_path)
        
        self.entries: OrderedDict[str, Entry] = {}
        for struc_type in struc_types:
            for struc_stack_name in os.listdir(os.path.join(project_dir, struc_type)):
                struc_stack_dir = os.path.join(project_dir, struc_type, struc_stack_name)
                if struc_stack_name in exclude_dirs:
                    exclude_dirs.remove(struc_stack_name)
                    continue
                if not struc_stack_name in self.entries:
                    self.entries[struc_stack_name] = Entry(struc_stack_dict={})
                self.entries[struc_stack_name].add_struc_stack(struc_type, struc_stack_dir)
        
        self.metrics = {}
        
    def sort_entries(self, key):
        self.entries = OrderedDict(sorted(self.entries.items(), key=key))
    
    def get_project_dir(self):
        return self.project_dir
    
    def set_summary_filename(self, summary_filename):
        self.summary_filename = summary_filename
    
    def get_summary_file_path(self):
        return os.path.join(self.project_dir, f"{self.summary_filename}.csv")
    
    def export(self):
        self.data.to_csv(self.get_summary_file_path())
    
    def calculate_metric(self, metric: Metric, struc_types: list[str], overwrite=False):
        for entry in self.entries.values():
            entry.calculate_metric(metric, struc_types, overwrite=False)
    
    def get_prediction_dir_identifier(self, identifier):
        return os.path.join(self.project_dir, identifier)
    
    def get_sequence_names(self):
        return self[self.sequence_name_col]
    
    def add_extra_cols(self, extra_cols, values=np.nan, overwrite=False):
        for col in extra_cols:
            if (not col in self.columns) or overwrite:
                self[col] = values
            else:
                print(f"Column {col} already exists, skipping. set overwrite=True to overwrite")
        self.extra_cols.extend(extra_cols)
    
    def get_metric_names(self):
        all_cols = list(self.columns)
        cols_to_ignore = [self.sequence_name_col, self.sequence_col] + self.extra_cols
        all_cols.remove(self.sequence_name_col)
        all_cols.remove(self.sequence_col)
        metric_names = []
        for metric_name in all_cols:
            metric_names.append(",".join(metric_name.split(",")[:-1]))
        metric_names = list(set(metric_names))
        metric_names.sort()
        return metric_names
        
    def check_summarization_mode(self, mode):
        if not mode in ["r", "w", "a"]:
            raise ValueError("mode should be 'r' for read, 'w' for write or 'a' for append")
    
    def get_metric_vector_for_design(self, metric_name, prediction_name):
        row_selector = self.get_sequence_names() == prediction_name
        col_selector = [f"{metric_name}, {i}" for i in range(5)]
        return np.array(self.loc[row_selector, col_selector], dtype=np.float64)
    
    def get_metric_mean_for_prediction(self, metric_name, prediction_name):
        row_selector = self.get_sequence_names() == prediction_name
        col_selector = f"{metric_name}, mean"
        return self.loc[row_selector, col_selector].values[0]
    
    def get_metric_std_for_prediction(self, metric_name, prediction_name):
        row_selector = self.get_sequence_names() == prediction_name
        col_selector = f"{metric_name}, std"
        return self.loc[row_selector, col_selector].values[0]
    
    def get_grouped_bar_x(self, group_i, group_len, index_i, barwidth=1):
        index_i = float(index_i)
        return (group_i * (group_len + 0.5) + index_i) * barwidth

    def plot_metrics(self, metric_names, figsize=None, ylims: list[tuple]=[(0, 1)], thresholds=None):
        
        metric_num = len(metric_names)
        entry_names = self.entries.keys()
        prediction_length = len(entry_names)
        
        def calculate_figsize(design_length, metric_indices_length):
            return (design_length * (metric_indices_length + 1) * 0.2, 5)
        if not figsize:
            figsize = calculate_figsize(prediction_length, metric_num)
            
        fig, ax = plt.subplots(figsize=figsize)
        bar_width = 1
        
        # 绘制多指标的柱状图
        if len(ylims) == 1:
            ylims = ylims * metric_num
        
        for i, metric_name in enumerate(metric_names):
            x_positions = [self.get_grouped_bar_x(group_i, metric_num, i, bar_width) for group_i in range(prediction_length)]
            y_positions = []
            
            for prediction_name in entry_names:
                y_positions.append(self.get_metric_mean_for_prediction(metric_name, prediction_name))
            y_pos_std = []
            for prediction_name in entry_names:
                y_pos_std.append(self.get_metric_std_for_prediction(metric_name, prediction_name))
            ax.bar(x_positions, y_positions, yerr=y_pos_std, label=metric_name, capsize=3, color=f"C{i}", width=bar_width, alpha=0.5)
            x_pos_for_vector = []
            for x in x_positions:
                for vector_i in range(5):
                    x_pos_for_vector.append(x + (vector_i - 2) * bar_width / 10)
            y_vector_positions = []
            for prediction_name in entry_names:
                y_vector_positions.extend(self.get_metric_vector_for_design(metric_name, prediction_name))
            ax.scatter(x_pos_for_vector, y_vector_positions, color=f"C{i}", marker=".")
        ax.set_ylim(*ylim)
        ax.set_xlabel(self.sequence_name_col)
        if metric_num == 1:
            ax.set_ylabel(metric_name)
        else:
            ax.set_ylabel("Metric Value")
        xticks = [self.get_grouped_bar_x(group_i, metric_num, metric_num / 2 - 0.5, bar_width) for group_i in range(prediction_length)]
        ax.set_xticks(xticks)
        ax.set_xticklabels(entry_names, rotation=90)
        # 将图例放在图的右上方, 但是不在内部
        ax.legend(loc="lower right", bbox_to_anchor=(1, 1), ncol=1)
        
        if thresholds:
            plot_threshhold(ax, thresholds)
        
        return fig, ax
    
    def plot_raw_data_correlation(self, metric1, metric2, xlim=(0, 1), ylim=(0, 1)):
        fig, ax = plt.subplots()
        xdata = []
        ydata = []
        for i in range(5):
            actual_metric1 = f"{metric1}, {i}"
            actual_metric2 = f"{metric2}, {i}"
            xdata.extend(self[actual_metric1])
            ydata.extend(self[actual_metric2])
        plt.scatter(xdata, ydata)
        ax.set_xlabel(metric1)
        ax.set_ylabel(metric2)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        return fig, ax
    
    def filter(self, filters: dict, stat="mean", summary_filename = "summary_filtered.csv"):
        '''
        Filters look like:
        {
            metric_name: (">", 0.5),
        }
        - The key is the name of the metric, which you can get from get_metric_names().
        - The value is a tuple of two elements:
            - The first element is the operator. It can be one of the following:
                - ">", "<", ">=", "<=", "=="
            - The second element is the threshold value.
        :param stat str: The statistic to be used for filtering. It can be one of the following:
            - "mean": The mean of the metric values.
            - "std": The median of the metric values.
            - 0~4: The i-th part of the metric vector.
        '''
        
        def get_operator_func(op_str):
            ops = {
                ">": operator.gt,
                "<": operator.lt,
                ">=": operator.ge,
                "<=": operator.le,
                "==": operator.eq
            }
            return ops[op_str]
        
        row_selectors = np.ones(len(self.get_sequence_names()), dtype=bool)
        for metric_name, (operator_str, threshold) in filters.items():
            op_func = get_operator_func(operator_str)
            metric_values = self[f'{metric_name}, {stat}']
            row_selectors = np.logical_and(row_selectors, op_func(metric_values, threshold))
        filtered_res = self.loc[row_selectors, :]
        return self.__class__(
            self.project_dir,
            filtered_res,
            summary_filename=summary_filename,
            sequence_name_col=self.sequence_name_col
        )
    
    def remove_metric(self, metric_name):
        self.drop([col for col in self.columns if metric_name in col], axis=1, inplace=True)
        
    def rename_metric(self, old_metric_name, new_metric_name):
        for col in self.columns:
            if col.startswith(old_metric_name):
                self.rename(columns={col: col.replace(old_metric_name, new_metric_name)}, inplace=True)
                
    def plot_sequence_metrics(self, sequence_name, metrics, figsize=None, ticksize=10, ylim=(0, 1), **kwargs):
        metric_num = len(metrics)
        
        fig_kw = {"subplot_kw": {"projection": "polar"}}
        if figsize: fig_kw["figsize"] = figsize
        fig, ax = plt.subplots(**fig_kw)
        
        angles = [n / float(metric_num) * 2 * math.pi for n in range(metric_num)]
        angles += angles[:1] # 让 polar plot 绕回来, 形成一个闭环
        plt.xticks(angles[:-1], metrics, color='black', size=ticksize)
        plt.ylim(*ylim)
        for i in range(5):
            values = list(self.loc[self.get_sequence_names() == sequence_name, [f"{metric}, {i}" for metric in metrics]].values[0])
            values += values[:1]
            ax.plot(angles, values, linewidth=2, linestyle='solid', label=f"Model {i}")
            # ax.fill(angles, values, alpha=0.1)
        ax.legend(loc="center left", bbox_to_anchor=(1.5, 0.5), ncol=1)
        ax.set_title(sequence_name)
        return fig, ax
    
    def plot_sequences_mean_metrics(self, sequence_names, metrics, figsize=None, ticksize=10, ylim=(0, 1), **kwargs):
        metric_num = len(metrics)
        
        fig_kw = {"subplot_kw": {"projection": "polar"}}
        if figsize: fig_kw["figsize"] = figsize
        fig, ax = plt.subplots(**fig_kw)
        
        angles = [n / float(metric_num) * 2 * math.pi for n in range(metric_num)]
        angles += angles[:1]
        plt.xticks(angles[:-1], metrics, color='black', size=ticksize)
        plt.ylim(*ylim)
        for sequence_name in sequence_names:
            values = list(self.loc[self.get_sequence_names() == sequence_name, [f"{metric}, mean" for metric in metrics]].values[0])
            values += values[:1]
            ax.plot(angles, values, linewidth=2, linestyle='solid', label=sequence_name)
            # ax.fill(angles, values, alpha=0.1)
        ax.legend(loc="center left", bbox_to_anchor=(1.5, 0.5), ncol=1)
        ax.set_title(sequence_name)
        return fig, ax
    
    def rescale_metric(self, metric_name, rescale_min, rescale_max):
        for col in self.columns:
            if col.startswith(metric_name) and not col.startswith(f"{metric_name}, {rescale_min}~{rescale_max}"):
                col_parts = col.split(", ")
                new_col_parts = col_parts[:-1] + [f"{rescale_min}~{rescale_max}"] + [col_parts[-1]]
                new_col = ", ".join(new_col_parts)
                self[new_col] = (self[col] - rescale_min) / (rescale_max - rescale_min)