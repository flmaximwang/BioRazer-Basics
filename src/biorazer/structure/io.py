from pathlib import Path
from biotite.structure.io import pdb, pdbx
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biorazer.io import Converter


def sequence_from_file(struc_file, protein_chains=[]):
    """
    This function reads a structure file and returns a dictionary with the chain IDs as keys and the sequences as values.
    """
    struc_file = Path(struc_file)
    if struc_file.suffix not in [".pdb", ".cif"]:
        raise ValueError("File must be in PDB or CIF format.")
    if not struc_file.exists():
        raise FileNotFoundError(f"File {struc_file} does not exist.")
    if struc_file.suffix == ".pdb":
        my_atom_array = pdb.get_structure(pdb.PDBFile.read(struc_file))[0]
    elif struc_file.suffix == ".cif":
        my_atom_array = pdbx.get_structure(pdbx.CIFFile.read(struc_file))[0]
    res = {}
    if len(protein_chains) == 0:
        protein_chains = list(bio_struc.get_chains(my_atom_array))
    for i in protein_chains:
        seq_i = "".join(
            list(
                map(
                    lambda x: bio_seq.ProteinSequence.convert_letter_3to1(x),
                    bio_struc.get_residues(my_atom_array[my_atom_array.chain_id == i])[
                        1
                    ],
                )
            )
        )
        res[i] = seq_i
    return res


class PDB2STRUCT(Converter):

    def read(self) -> bio_struc.AtomArray:
        return pdb.get_structure(pdb.PDBFile.read(self.input_file))[0]


class CIF2STRUCT(Converter):

    def read(self):
        return pdbx.get_structure(pdbx.CIFFile.read(self.input_file))[0]


class STRUCT2CIF(Converter):

    def write(self, tmp):
        output_file_obj = pdbx.CIFFile()
        pdbx.set_structure(output_file_obj, tmp)
        output_file_obj.write(self.output_file)


class STRUCT2PDB(Converter):

    def write(self, tmp):
        output_file_obj = pdb.PDBFile()
        pdb.set_structure(output_file_obj, tmp)
        output_file_obj.write(self.output_file)


class PDB2CIF(PDB2STRUCT, STRUCT2CIF):

    pass


class CIF2PDB(CIF2STRUCT, STRUCT2PDB):

    pass
