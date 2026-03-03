import warnings
from pathlib import Path
import numpy as np
from biotite.structure.io import pdb, pdbx
from biotite.structure import AtomArray
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biorazer.io import Converter
from biorazer.sequence.io import SEQ2FASTA

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    module="biotite.structure.io.pdbx",
    message="Attribute .* not found within .* category.*",
)


class CIF2STRUCT(Converter):

    def read(self, **kwargs) -> AtomArray:
        cif_file = pdbx.CIFFile.read(self.input_file)
        resn = cif_file.block["chem_comp_atom"]["comp_id"].data.array
        name = cif_file.block["chem_comp_atom"]["atom_id"].data.array
        element = cif_file.block["chem_comp_atom"]["type_symbol"].data.array
        charge = cif_file.block["chem_comp_atom"]["charge"].data.array
        x = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_x_ideal"].data.array
        y = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_y_ideal"].data.array
        z = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_z_ideal"].data.array
        structure = bio_struc.AtomArray(len(name))
        structure.res_name = resn
        structure.atom_name = name
        structure.element = element
        structure.charge = charge
        structure.coord = np.column_stack((x, y, z))
        return structure
