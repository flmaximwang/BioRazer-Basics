from pathlib import Path
from biotite.structure.io import pdb, pdbx
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biorazer.io import Converter
from biorazer.sequence.io import SEQ2FASTA


class PDB2STRUCT(Converter):
    def read(self, **kwargs) -> bio_struc.AtomArray:
        return pdb.get_structure(pdb.PDBFile.read(self.input_file), **kwargs)[0]


class CIF2STRUCT(Converter):
    def read(self, **kwargs) -> bio_struc.AtomArray:
        return pdbx.get_structure(pdbx.CIFFile.read(self.input_file), **kwargs)[0]


class STRUCT2CIF(Converter):
    def write(self, tmp, **kwargs):
        output_file_obj = pdbx.CIFFile()
        pdbx.set_structure(output_file_obj, tmp, **kwargs)
        output_file_obj.write(self.output_file)


class STRUCT2PDB(Converter):
    def write(self, tmp, **kwargs):
        output_file_obj = pdb.PDBFile()
        pdb.set_structure(output_file_obj, tmp, **kwargs)
        output_file_obj.write(self.output_file)


class PDB2CIF(PDB2STRUCT, STRUCT2CIF):
    pass


class CIF2PDB(CIF2STRUCT, STRUCT2PDB):
    pass


class CIF2CIF(CIF2STRUCT, STRUCT2CIF):
    pass


class PDB2PDB(PDB2STRUCT, STRUCT2PDB):
    pass


class PDB2SEQ(Converter):
    """
    Converts a PDB file to a sequence dictionary.
    """

    def read(self, **kwargs) -> dict:
        structure = PDB2STRUCT(self.input_file, self.output_file).read(**kwargs)
        chain_ids = bio_struc.get_chains(structure)
        res = {}
        for chain_id in chain_ids:
            chain_structure = structure[structure.chain_id == chain_id]
            res_ids, res_names = bio_struc.get_residues(chain_structure)
            one_char_res_names = list(
                map(lambda x: bio_seq.ProteinSequence.convert_letter_3to1(x), res_names)
            )
            sequence = "".join(one_char_res_names)
            res[chain_id] = sequence
        return res


class PDB2FASTA(PDB2SEQ, SEQ2FASTA):
    pass
