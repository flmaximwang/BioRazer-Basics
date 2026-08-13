"""
Report functions for per-atom boolean masks.

The functions in this module turn a per-atom mask into a
residue-level report, following the same ``fmt`` conventions
(``pymol`` / ``text`` / ``list``) as the static analysis report
modules.
"""
from biotite import structure as bio_struct
import numpy as np
from biotite.structure import AtomArray
from biorazer.display import print_with_decoration, print_decoration_line
from biorazer.structure.util.report import (
    _normalize_fmt,
)
from .annotation import extend_by_res


def report_mask_by_res(
    atom_array: AtomArray,
    selection: np.ndarray,
    fmt="pymol",
    pymol_model_name="",
    pymol_selection_name=None,
):
    """
    Report the residues covered by an atom mask in different formats.

    The mask is first extended to whole residues via
    :func:`extend_by_res`, so any atom selected by ``mask`` marks its
    entire residue as covered. Note that the extension is based on
    ``res_id`` alone -- residues sharing the same number in another
    chain (or with a different ``ins_code``) are marked as well.

    Parameters
    ----------
    atom_array : biotite.structure.AtomArray
        The structure the mask is defined on.
    mask : numpy.ndarray
        Boolean atom mask aligned with ``atom_array``.
    fmt : str
        - pymol: print PyMOL commands to select the masked residues
        - list: return a list of tuples (chain_id, res_id)
        - text: print the masked residues in text format
    pymol_model_name : str, optional
        PyMOL object (model) 名称, 用于构造 /model//chain/res/ 选择器。
    pymol_selection_prefix : str, optional
        PyMOL selection 的命名前缀, 生成的 selection 名为 {prefix}_mask。
        用于在同一 model 上区分多个 mask; 缺省回退为 pymol_model_name。

    Returns
    -------
    list of tuple or None
        A list of ``(chain_id, res_id)`` tuples when ``fmt == "list"``,
        otherwise ``None`` (the report is printed).
    """
    fmt = _normalize_fmt(fmt, ("pymol", "text", "list"))

    mask_by_res = extend_by_res(atom_array=atom_array, mask=selection)

    residues = []
    for atom in atom_array[selection]:
        identifier = (str(atom.chain_id), int(atom.res_id))
        if identifier not in residues:
            residues.append(identifier)

    if fmt == "list":
        return residues
    elif fmt == "pymol":
        selection_name = f"{pymol_selection_name or pymol_model_name}"
        print_with_decoration("Copy the command below to PyMOL", decoration_char="#")
        print(f"select {selection_name}, not all")
        for chain_id, res_id in residues:
            print(
                f"select {selection_name}, /{pymol_model_name}//{chain_id}/{res_id}/ or {selection_name}"
            )
        print_decoration_line(decoration_char="#")
    elif fmt == "text":
        for chain_id, res_id in residues:
            print(f"Chain {chain_id}, Residue {res_id}")
