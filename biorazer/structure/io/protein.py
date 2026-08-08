from pathlib import Path
import numpy as np
from biotite.structure.io import pdb, pdbx
from biotite.structure import AtomArray
import biotite.structure as bio_struc
import biotite.sequence as bio_seq
from biotite.structure.io.pdb.hybrid36 import encode_hybrid36
from biorazer.io import Converter
from biorazer.sequence.io import SEQ2FASTA


class PDB_STRUCT(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdb.get_structure(pdb.PDBFile.read(self.input_io), **kwargs)[0]


class CIF_STRUCT(Converter):
    def read(self, **kwargs) -> AtomArray:
        return pdbx.get_structure(pdbx.CIFFile.read(self.input_io), **kwargs)[0]


class STRUCT_CIF(Converter):
    def write(self, tmp, **kwargs):
        output_file_obj = pdbx.CIFFile()
        pdbx.set_structure(output_file_obj, tmp, **kwargs)
        output_file_obj.write(self.output_io)


# Same limit as biotite's PDB writer (_PDB_MAX_RESIDUES), so that the
# residue sequence numbers in LINK/SSBOND match the ATOM records.
_PDB_MAX_RESIDUES = 9999


def _format_pdb_atom_name(atom_name: str, element: str) -> str:
    """
    Format a biotite atom name into the 4-column PDB atom-name field
    (names of atoms with a one-letter element are shifted right by one).
    """
    name = atom_name.strip()
    if len(element) == 1 and len(name) < 4:
        name = " " + name
    return name.ljust(4)


def _format_pdb_res_id(res_id: int, hybrid36: bool) -> str:
    """Format a residue sequence number exactly like the ATOM records."""
    if hybrid36:
        return encode_hybrid36(res_id, 4)
    if res_id > 0:
        res_id = ((res_id - 1) % _PDB_MAX_RESIDUES) + 1
    return f"{res_id:>4}"


def _format_link_record(array: AtomArray, protein_i: int, ligand_j: int, hybrid36: bool) -> str:
    """Build a LINK record for a covalent bond between a protein atom and a ligand atom."""
    dist = np.linalg.norm(array.coord[protein_i] - array.coord[ligand_j])
    return (
        "LINK  "
        + " " * 6
        + _format_pdb_atom_name(array.atom_name[protein_i], array.element[protein_i])
        + " "
        + f"{array.res_name[protein_i]:>3}"
        + " "
        + f"{array.chain_id[protein_i]:<1}"
        + _format_pdb_res_id(array.res_id[protein_i], hybrid36)
        + f"{array.ins_code[protein_i]:<1}"
        + " " * 15
        + _format_pdb_atom_name(array.atom_name[ligand_j], array.element[ligand_j])
        + " "
        + f"{array.res_name[ligand_j]:>3}"
        + " "
        + f"{array.chain_id[ligand_j]:<1}"
        + _format_pdb_res_id(array.res_id[ligand_j], hybrid36)
        + f"{array.ins_code[ligand_j]:<1}"
        + "  "
        + f"{'1555':>6}"
        + " "
        + f"{'1555':>6}"
        + " "
        + f"{dist:>5.2f}"
        + "  "
    )


def _format_ssbond_record(
    array: AtomArray, cys_i: int, cys_j: int, serial: int, hybrid36: bool
) -> str:
    """Build an SSBOND record for a CYS SG - CYS SG disulfide bond."""
    dist = np.linalg.norm(array.coord[cys_i] - array.coord[cys_j])
    return (
        "SSBOND"
        + " "
        + f"{serial:>3}"
        + " "
        + "CYS"
        + " "
        + f"{array.chain_id[cys_i]:<1}"
        + " "
        + _format_pdb_res_id(array.res_id[cys_i], hybrid36)
        + f"{array.ins_code[cys_i]:<1}"
        + "   "
        + "CYS"
        + " "
        + f"{array.chain_id[cys_j]:<1}"
        + " "
        + _format_pdb_res_id(array.res_id[cys_j], hybrid36)
        + f"{array.ins_code[cys_j]:<1}"
        + " " * 23
        + f"{'1555':>6}"
        + " "
        + f"{'1555':>6}"
        + " "
        + f"{dist:>5.2f}"
        + "  "
    )


def _format_link_records(array: AtomArray, hybrid36: bool) -> list[str]:
    """Find covalent bonds between a canonical amino acid and a small molecule."""
    protein_mask = bio_struc.filter_canonical_amino_acids(array)
    ligand_mask = ~protein_mask & ~bio_struc.filter_solvent(array)
    records = []
    for atom_i, atom_j, _ in array.bonds.as_array():
        if protein_mask[atom_i] and ligand_mask[atom_j]:
            records.append(_format_link_record(array, atom_i, atom_j, hybrid36))
        elif protein_mask[atom_j] and ligand_mask[atom_i]:
            records.append(_format_link_record(array, atom_j, atom_i, hybrid36))
    return records


def _format_ssbond_records(array: AtomArray, hybrid36: bool) -> list[str]:
    """Find disulfide bonds (SG of two CYS residues) and format them as SSBOND."""
    cys_sg = (array.res_name == "CYS") & (np.char.strip(array.atom_name) == "SG")
    records = []
    for atom_i, atom_j, _ in array.bonds.as_array():
        if atom_i != atom_j and cys_sg[atom_i] and cys_sg[atom_j]:
            records.append(
                _format_ssbond_record(array, atom_i, atom_j, len(records) + 1, hybrid36)
            )
    return records


class STRUCT_PDB(Converter):
    def write(self, tmp: AtomArray, hybrid36: bool = False):
        """
        Write a PDB file. Intermolecular covalent bonds between a protein
        and a small molecule are additionally written as LINK records, and
        disulfide bonds as SSBOND records, appended at the end of the file.
        """
        output_file_obj = pdb.PDBFile()
        pdb.set_structure(output_file_obj, tmp, hybrid36=hybrid36)
        if tmp.bonds is not None and tmp.coord.ndim == 2:
            output_file_obj.lines.extend(
                _format_link_records(tmp, hybrid36)
                + _format_ssbond_records(tmp, hybrid36)
            )
        output_file_obj.write(self.output_io)


class PDB_SEQ(Converter):
    """
    Converts a PDB file to a sequence dictionary.
    """

    def read(self, **kwargs) -> dict:
        structure = PDB_STRUCT(self.input_io, self.output_io).read(**kwargs)
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
