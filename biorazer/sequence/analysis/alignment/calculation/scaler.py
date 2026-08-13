from numbers import Number
import numpy as np
from scipy.stats import entropy

def calculate_entropy(array_like: np.ndarray | list, base=2):
    """
    Calculate the entropy of a 1-D array
    """
    if isinstance(array_like, list):
        for i in array_like:
            if not isinstance(i, Number):
                raise TypeError("Input list should contain only numeric elements")
        array = np.array(array_like)
    elif isinstance(array_like, np.ndarray):
        if array_like.shape.__len__() != 1:
            raise ValueError("Input array should be 1-dimensional")
        array = array_like
    else:
        raise TypeError(f"Unsupported input array {type(array_like)}")

    return entropy(array, base=2)