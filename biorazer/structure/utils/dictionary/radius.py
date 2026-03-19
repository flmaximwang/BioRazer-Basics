import numpy as np

vdw_dict = {
    "H": np.float32(0),
    "C": np.float32(1.70),
    "N": np.float32(1.60),
    "O": np.float32(1.50),
    "P": np.float32(1.80),
    "S": np.float32(1.80),
}


def vdw_radii(element):
    return vdw_dict[element]


vdw_radii = np.vectorize(vdw_radii)
