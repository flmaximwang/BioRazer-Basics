from pathlib import Path
from biotite.structure.io import pdb, pdbx
from biotite.structure import AtomArray
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biorazer.io import Converter
from biorazer.sequence.io import SEQ2FASTA


class PDB2STRUCT(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdb.get_structure(pdb.PDBFile.read(self.input_file), **kwargs)[0]


class CIF2STRUCT(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdbx.get_structure(pdbx.CIFFile.read(self.input_file), **kwargs)[0]


class STRUCT2CIF(Converter):
    def write(self, tmp, **kwargs):
        output_file_obj = pdbx.CIFFile()
        pdbx.set_structure(output_file_obj, tmp, **kwargs)
        output_file_obj.write(self.output_file)


class STRUCT2PDB(Converter):
    def write(self, tmp: AtomArray, **kwargs):
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
            res_ids = list(res_ids)
            res_names = list(res_names)
            one_char_res_names = []
            for i in range(min(res_ids), max(res_ids) + 1):
                try:
                    idx = res_ids.index(i)
                except ValueError:
                    # Missing residue
                    one_char_res_names.append("X")
                    continue
                res_name = res_names[idx]
                if len(res_name) != 3:
                    # Nucleotides
                    break
                try:
                    one_char_res_names.append(
                        bio_seq.ProteinSequence.convert_letter_3to1(res_name)
                    )
                except KeyError:
                    # Non-standard amino acid or ligand
                    one_char_res_names.append("X")

            sequence = "".join(one_char_res_names)
            res[chain_id] = sequence

        return res


class PDB2FASTA(PDB2SEQ, SEQ2FASTA):
    pass
