import numpy as np
from biotite import structure as bio_struct


def _normalize_fmt(fmt, supported_formats):
    """
    Normalize and validate report output format.

    Parameters
    ----------
    fmt : str
        Raw format input from caller.
    supported_formats : tuple[str, ...]
        Supported canonical format names.
    """

    if not isinstance(fmt, str):
        raise TypeError(f"Format must be str, got {type(fmt).__name__}.")

    normalized_fmt = fmt.strip().lower()
    fmt_alias = {
        "txt": "text",
        "lst": "list",
        "table": "pandas",
    }
    normalized_fmt = fmt_alias.get(normalized_fmt, normalized_fmt)

    if normalized_fmt not in supported_formats:
        supported = ", ".join(supported_formats)
        raise ValueError(f"Format {fmt} not supported. Supported formats: {supported}")

    return normalized_fmt


def _normalize_interface_residues(
    atom_array: bio_struct.AtomArray, interface_atom_mask: np.ndarray
):
    interface_residues = []
    for atom in atom_array[interface_atom_mask]:
        identifier = (str(atom.chain_id), int(atom.res_id))
        if identifier not in interface_residues:
            interface_residues.append(identifier)
    return interface_residues


def _format_atom_label(atom):
    return f"{atom.chain_id} {atom.res_name}{atom.res_id}({atom.atom_name})"


def _to_pymol_atom_selector(model_name: str, atom):
    return f"/{model_name}//{atom.chain_id}/{atom.res_id}/{atom.atom_name}"
