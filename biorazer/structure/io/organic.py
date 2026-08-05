import io
import warnings
from pathlib import Path
import numpy as np

from rdkit import Chem
from rdkit.Chem import Mol, SDWriter

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

#: V2000 formal-charge codes (molfile spec) for charges in [-3, 3]
_SDF_CHARGE = {3: 1, 2: 2, 1: 3, 0: 0, -1: 5, -2: 6, -3: 7}

#: biotite BondType -> V2000 bond type (1 single, 2 double, 3 triple, 4 aromatic)
_SDF_BOND = {1: 1, 2: 2, 3: 3, 4: 4}


def _sdf_block(tmp: AtomArray) -> str:
    """
    Serialize an AtomArray as a single V2000 SDF entry.

    The entry contains all atoms of the array (element, coordinates and
    formal charge) and its bonds. If the array has no bonds annotation,
    bonds are inferred from interatomic distances with biotite's
    :func:`connect_via_distances` (single bonds).

    Charges are written both as inline V2000 charge codes (range
    [-3, 3]) and as ``M CHG`` records, which also cover charges outside
    that range. Aromatic bonds (biotite ``BondType.AROMATIC``) are
    written as V2000 bond type 4.
    """
    n_atoms = len(tmp)
    if tmp.bonds is None:
        bonds = bio_struc.connect_via_distances(tmp)
    else:
        bonds = tmp.bonds
    n_bonds = bonds.get_bond_count()
    has_charge = "charge" in tmp.get_annotation_categories()

    lines = [
        "",  # molecule name
        "  biorazer          3D",
        "",  # comment line
    ]
    lines.append(f"{n_atoms:3d}{n_bonds:3d}  0  0  0  0  0  0  0  0999 V2000")

    charged_atoms = []  # (1-based atom index, charge) for the M CHG records
    for i in range(n_atoms):
        x, y, z = tmp.coord[i]
        q = int(tmp.charge[i]) if has_charge else 0
        if q != 0:
            charged_atoms.append((i + 1, q))
        charge_code = _SDF_CHARGE.get(q, 0)
        element = tmp.element[i] or "*"
        lines.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} {element:<3s}  0{charge_code:3d}"
            "  0  0  0  0  0  0  0  0  0  0"
        )

    for a, b, order in bonds.as_array():
        lines.append(
            f"{a + 1:3d}{b + 1:3d}{_SDF_BOND.get(int(order), 1):3d}  0  0  0  0"
        )

    # Each M CHG line holds up to 8 charged atoms (V2000 property block)
    for start in range(0, len(charged_atoms), 8):
        batch = charged_atoms[start : start + 8]
        lines.append(
            f"M  CHG{len(batch):3d}"
            + "".join(f" {atom_i:3d} {c:3d}" for atom_i, c in batch)
        )

    lines.append("M  END")
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


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

    def read(self) -> AtomArray:
        raise NotImplementedError(
            "STRUCT2SDF is a write-only converter, "
            "use CIF2STRUCT or SDF2MOL for reading"
        )

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
