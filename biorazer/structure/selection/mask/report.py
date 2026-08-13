import numpy as np
from biotite.structure import AtomArray
from .annotation import extend_by_res

def report_mask_by_res(atom_array: AtomArray, mask: np.ndarray):

    mask_by_res = extend_by_res(atom_array=atom_array, mask=mask)
    
