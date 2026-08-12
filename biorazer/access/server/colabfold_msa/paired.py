"""colabfold_msa.paired - 多链配对 (paired) MSA 的提交与结果合并。

ColabFold 约定: 一条记录内以 ':' 分隔各链; 提交时拆成独立记录
(>101, >102, ...) 到 ticket/pair 端点, 服务器端配对; 返回后把各链段
按行拼接为整复合物一份 a3m (colabfold pair_sequences 的展开)。
"""

import os
import sys
from typing import List, Tuple

from .api import _paired_mode, _submit_and_download
from .io import _extract_and_clean, _read_a3m_lines, _split_a3m_sections


def _merge_paired_a3m(lines: List[str], n_chains: int) -> str:
    """把多链配对结果 (每条链一段) 合并为整复合物 a3m。

    与 colabfold 的 pair_sequences 一致: 服务器保证各段行数一致且行
    对应 (同一配对的链序列; 缺失方以 DUMMY/- 填充)，逐行拼接——
    header 从第二段起 '>' 换成 '\\t'，序列行直接拼接。
    单链 (仅一段) 时原样返回。
    """
    sections = _split_a3m_sections(lines)
    if len(sections) <= 1:
        return "\n".join(lines) + "\n"
    if len(sections) != n_chains:
        print(f"  [!] 服务器返回 {len(sections)} 段, 期望 {n_chains} 段, 按实际段数合并",
              file=sys.stderr)
    n_rows = min(len(s) for s in sections)
    merged: List[str] = []
    for i in range(n_rows):
        parts = []
        for j, sec in enumerate(sections):
            line = sec[i]
            if line.startswith(">") and j > 0:
                line = "\t" + line[1:]
            parts.append(line)
        merged.append("".join(parts))
    return "\n".join(merged) + "\n"


def run_paired(chains: List[str], tmp: str, job_dir: str,
               use_env: bool, pair_strategy: str,
               host: str, ua: str) -> Tuple[List[str], List[int]]:
    """提交配对任务 (ticket/pair) 并整理 pair.a3m 到 job_dir。

    服务器返回的 tar 解压到 tmp 后原地清理。多链 (len(chains) > 1)
    时把各链段合并为整复合物 a3m，单链时原样写出。

    Returns
    -------
    (a3m_files, Ms)  a3m_files 为写出的文件列表 (可为空), Ms 为提交编号
    """
    tar_gz, Ms = _submit_and_download(
        chains, tmp, "ticket/pair", _paired_mode(use_env, pair_strategy),
        host, ua, tar_name="out_paired.tar.gz")
    _extract_and_clean(tar_gz, tmp)
    src = os.path.join(tmp, "pair.a3m")
    if not os.path.isfile(src):
        return [], Ms
    text = _merge_paired_a3m(_read_a3m_lines(src), len(chains))
    dst = os.path.join(job_dir, "pair.a3m")
    with open(dst, "w") as f:
        f.write(text)
    return [dst], Ms