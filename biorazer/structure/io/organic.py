import io
import warnings
from pathlib import Path
import numpy as np

from rdkit import Chem
from rdkit.Chem import Mol, SDWriter

from biotite.structure.io import pdb, pdbx
from biotite.structure.io.mol import SDFile, set_structure
import biotite.structure as bio_struc
from biotite.structure import AtomArray, BondList
import biotite.sequence as bio_seq

from biorazer.io import Converter
from biorazer.sequence.io import SEQ2FASTA

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    module="biotite.structure.io.pdbx",
    message="Attribute .* not found within .* category.*",
)

def _sdf_block(tmp: AtomArray) -> str:
    """
    Serialize an AtomArray as a single V2000 SDF entry.

    The serialization is delegated to biotite's
    :class:`~biotite.structure.io.mol.SDFile`, which maps biotite bond
    types -- including the aromatic types -- to V2000 bond types and
    writes formal charges both as inline V2000 charge codes and as
    ``M CHG`` records. If the array has no bonds annotation, bonds are
    inferred from interatomic distances with biotite's
    :func:`connect_via_distances` (single bonds).
    """
    if tmp.bonds is None:
        tmp = tmp.copy()
        tmp.bonds = bio_struc.connect_via_distances(
            tmp, default_bond_type=bio_struc.BondType.SINGLE
        )
    sdf = SDFile()
    set_structure(sdf, tmp, record_name="")
    text_io = io.StringIO()
    sdf.write(text_io)
    return text_io.getvalue()


class CIF2STRUCT(Converter):
    """
    Read small molecule cif files to biotite AtomArray
    """

    def read(self, **kwargs) -> AtomArray:
        cif_file = pdbx.CIFFile.read(self.input_io)
        resn = cif_file.block["chem_comp_atom"]["comp_id"].data.array
        name = cif_file.block["chem_comp_atom"]["atom_id"].data.array
        element = cif_file.block["chem_comp_atom"]["type_symbol"].data.array
        charge = cif_file.block["chem_comp_atom"]["charge"].data.array
        x = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_x_ideal"].data.array
        y = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_y_ideal"].data.array
        z = cif_file.block["chem_comp_atom"]["pdbx_model_Cartn_z_ideal"].data.array
        structure = bio_struc.AtomArray(len(name))
        structure.add_annotation("charge", dtype=int)
        structure.res_name = resn
        structure.atom_name = name
        structure.element = element
        structure.charge = np.asarray(charge, dtype=int)
        structure.coord = np.column_stack((x, y, z))
        return structure


class SDF2MOL(Converter):

    def read(self, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True, **kwargs) -> Mol:
        if isinstance(self.input_io, io.StringIO):
            # ForwardSDMolSupplier requires a binary file object
            supplier = Chem.ForwardSDMolSupplier(
                io.BytesIO(self.input_io.getvalue().encode("utf-8")),
                sanitize=sanitize,
                removeHs=removeHs,
                strictParsing=strictParsing,
            )
            return next(iter(supplier), None)
        mol = Chem.MolFromMolFile(str(self.input_io), sanitize=sanitize, removeHs=removeHs, strictParsing=strictParsing)
        return mol


class STRUCT2SDF(Converter):
    """
    Write a biotite AtomArray to an SDF file or an ``io.StringIO``.

    The atom array is serialized as a single V2000 SDF entry without
    RDKit: element, coordinates, formal charges and bonds are taken
    directly from the array. If the array carries no bonds annotation
    (e.g. freshly read small-molecule CIFs), single bonds are inferred
    from interatomic distances. For chemically correct bond orders, run
    the array through
    :func:`biorazer.structure.operation.annotation.add_bonds_to_organic`
    first.
    """

    def write(self, tmp: AtomArray) -> str | None:
        """
        Parameters
        ----------
        tmp : AtomArray
            The atom array to write.

        Returns
        -------
        str or None
            The SDF text if `self.output_io` is an ``io.StringIO``,
            otherwise None.
        """
        text = _sdf_block(tmp)
        if isinstance(self.output_io, io.StringIO):
            self.output_io.write(text)
            return self.output_io.getvalue()
        Path(self.output_io).write_text(text)
        return None

class STRUCT_MOL2(Converter):

    def write(self, tmp: AtomArray, molecule_name="UNK") -> None:
        bonds: BondList = tmp.bonds
        n_atoms = len(tmp)
        n_bonds = len(bonds.as_array()) if bonds is not None else 0
        
        with self._text_io(self.output_io, "w") as f:
            # @<TRIPOS>MOLECULE block
            f.write("@<TRIPOS>MOLECULE\n")
            f.write(f"{molecule_name}\n")
            f.write(f"{n_atoms} {n_bonds} 0 0 0\n")
            f.write("SMALL\n")
            f.write("USER_CHARGES\n\n")

            f.write("@<TRIPOS>ATOM\n")
            for i in range(n_atoms):
                atom = tmp[i]
                x, y, z = atom.coord
                atom_name = (
                    atom.atom_name if hasattr(atom, "atom_name") else atom.element
                )
                atom_type = (
                    atom.element
                )  # Map to SYBYL atom types if needed (e.g., C.3, H)
                res_name = (
                    atom.res_name if hasattr(atom, "res_name") else "RES"
                )
                res_id = atom.res_id if hasattr(atom, "res_id") else 1
                charge = atom.charge if hasattr(atom, "charge") else 0.0

                f.write(
                    f"{i+1:7d} {atom_name:<8s} {x:10.4f} {y:10.4f} {z:10.4f}"
                    f" {atom_type:<5s} {res_id:5d} {res_name:<8s} {charge:10.4f}\n"
                )

            # @<TRIPOS>BOND block (if bonds exist)
            if n_bonds > 0:
                f.write("\n@<TRIPOS>BOND\n")

                _BOND_ORDER_MAP = {
                    1: "1",
                    2: "2",
                    3: "3",
                    5: "ar",
                    6: "ar"
                }

                for i, bond in enumerate(bonds.as_array()):
                    # bond is (atom1_index, atom2_index, bond_order)

                    f.write(
                        f"{i+1:6d} {bond[0]+1:6d} {bond[1]+1:6d} {_BOND_ORDER_MAP[int(bond[2])]:>5}\n"
                    )


class MOL2SDF(SDF2MOL):
    """
    Write an RDKit Mol to an SDF file or an ``io.StringIO``.

    The read side is inherited from :class:`SDF2MOL`.
    """

    def write(self, tmp: Mol) -> str | None:
        """
        Parameters
        ----------
        tmp : Mol
            The RDKit molecule to write.

        Returns
        -------
        str or None
            The SDF text if `self.output_io` is an ``io.StringIO``,
            otherwise None.
        """
        target = self.output_io if isinstance(self.output_io, io.StringIO) else str(self.output_io)
        writer = SDWriter(target)
        try:
            writer.write(tmp)
        finally:
            writer.close()
        if isinstance(self.output_io, io.StringIO):
            return self.output_io.getvalue()
        return None


class CONF2SDF(SDF2MOL):
    """
    Write conformer(s) of an RDKit Mol to an SDF file or an
    ``io.StringIO``, one SDF entry per conformer.

    The read side is inherited from :class:`SDF2MOL`.
    """

    def write(
        self,
        tmp: Mol,
        conf_id: int | list[int] | None = None,
    ) -> str | None:
        """
        Parameters
        ----------
        tmp : Mol
            The RDKit molecule whose conformers are written.
        conf_id : int, list of int or None
            The conformer id(s) to write. If None, all conformers of
            `tmp` are written.

        Returns
        -------
        str or None
            The SDF text if `self.output_io` is an ``io.StringIO``,
            otherwise None.
        """
        target = self.output_io if isinstance(self.output_io, io.StringIO) else str(self.output_io)
        conformer_ids = [conf.GetId() for conf in tmp.GetConformers()]
        if conf_id is None:
            conf_ids = conformer_ids
            if not conf_ids:
                raise ValueError("Mol has no conformers, nothing to write")
        elif isinstance(conf_id, int):
            if conf_id not in conformer_ids:
                raise ValueError(f"Mol has no conformer with id {conf_id}")
            conf_ids = [conf_id]
        else:
            conf_ids = list(conf_id)

        writer = SDWriter(target)
        try:
            for cid in conf_ids:
                writer.write(tmp, confId=cid)
        finally:
            writer.close()
        if isinstance(self.output_io, io.StringIO):
            return self.output_io.getvalue()
        return None