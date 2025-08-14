import os
import pandas as pd

def batch_extract_prediction(batch_dir, extract_func, exclude=[".DS_Store"]):
    """
    Extracts predictions from a directory of predictions.
    :param prediction_dir: Directory containing predictions.
    :param extract_func: Function to extract predictions from a single file.
    :param ignore_dirs: List of directories to ignore.
    :return: List of extracted predictions.
    """
    counter = 0
    for pred_dir in os.listdir(batch_dir):
        if pred_dir in exclude:
            continue
        if counter == 0:
            res = extract_func(os.path.join(batch_dir, pred_dir))
        else:
            res.loc[counter, :] = extract_func(os.path.join(batch_dir, pred_dir)).loc[0, :]
        counter += 1
    
    return res