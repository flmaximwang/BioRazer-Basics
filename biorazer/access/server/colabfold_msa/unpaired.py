"""colabfold_msa.unpaired - 非配对 (unpaired) MSA 的提交与结果整理。

提交各链到 ticket/msa 端点 (服务器各链独立搜索), 返回后把各链段
按链位置补 '-' 合并为整复合物 a3m (colabfold pad_sequences 的展开)。
"""

import os
import sys
from typing import List, Tuple

from .api import _submit_and_download, _unpaired_mode
from .io import _extract_and_clean, _read_a3m_lines, _split_a3m_sections


def _merge_unpaired_a3m(lines: List[str], chain_lengths: List[int]) -> str:
    """把多链 unpaired 结果合并为整复合物 a3m (colabfold pad_sequences)。

    各链段保留自己的 header，序列行按链位置前后补 '-' 到复合物全长，
    便于作为多链 a3m 直接使用。
    """
    sections = _split_a3m_sections(lines)
    blanks = ["-" * L for L in chain_lengths]
    out: List[str] = []
    for n, sec in enumerate(sections):
        for line in sec:
            if line.startswith(">"):
                out.append(line)
            else:
                out.append("".join(blanks[:n] + [line] + blanks[n + 1:]))
    return "\n".join(out) + "\n"


def run_unpaired(chains: List[str], tmp: str, job_dir: str,
                 use_env: bool, use_filter: bool,
                 host: str, ua: str,
                 tar_name: str = "out.tar.gz") -> Tuple[List[str], List[int]]:
    """提交非配对任务 (ticket/msa) 并整理各数据库 a3m 到 job_dir。

    与 paired 同跑 (paired+unpaired) 时传 tar_name="out_unpaired.tar.gz"
    避免与配对结果在同一个 tmp 目录里同名冲突。

    Returns
    -------
    (a3m_files, Ms)  a3m_files 为写出的文件列表 (可为空), Ms 为提交编号
    """
    tar_gz, Ms = _submit_and_download(
        chains, tmp, "ticket/msa", _unpaired_mode(use_env, use_filter),
        host, ua, tar_name=tar_name)
    _extract_and_clean(tar_gz, tmp)
    db_names = ["uniref.a3m"]
    if use_env:
        db_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")
    chain_lengths = [len(c) for c in chains]
    a3m_files: List[str] = []
    for db in db_names:
        src = os.path.join(tmp, db)
        if not os.path.isfile(src):
            continue
        lines = _read_a3m_lines(src)
        if len(chains) > 1:
            text = _merge_unpaired_a3m(lines, chain_lengths)
        else:
            text = "\n".join(lines) + "\n"
        dst = os.path.join(job_dir, db)
        with open(dst, "w") as f:
            f.write(text)
        a3m_files.append(dst)
    return a3m_files, Ms