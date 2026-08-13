"""colabfold_msa.paired - 多链配对 (paired) MSA 的提交与结果合并。

ColabFold 约定: 一条记录内以 ':' 分隔各链; 提交时拆成独立记录
(>101, >102, ...) 到 ticket/pair 端点, 服务器端配对; 返回后把各链段
按行拼接为整复合物一份 a3m (colabfold pair_sequences 的展开) 写入
paired/ 子目录, 同时保留每链原始段 paired_0.a3m, paired_1.a3m, ...。

模式未请求 (active=False) 或服务器无结果时, 仅写 query 序列的 fallback
文件, 保证下游路径恒定。
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


def _write_query_only(chains: List[str], mode_dir: str, M0: int = 101
                      ) -> Tuple[List[str], List[int]]:
    """模式未请求/无结果时, 仅写 query 序列的 fallback 文件。

    单链 → paired.a3m; 多链 → paired.a3m (合并形式) + paired_0.a3m, ...
    """
    files: List[str] = []
    # paired.a3m: 合并形式的 query-only (多链时 header tab 连接)
    synthetic: List[str] = []
    for i, seq in enumerate(chains):
        synthetic.append(f">{M0 + i}")
        synthetic.append(seq)
    p = os.path.join(mode_dir, "paired.a3m")
    with open(p, "w") as f:
        f.write(_merge_paired_a3m(synthetic, len(chains)))
    files.append(p)
    if len(chains) > 1:
        for i, seq in enumerate(chains):
            p = os.path.join(mode_dir, f"paired_{i}.a3m")
            with open(p, "w") as f:
                f.write(f">{M0 + i}\n{seq}\n")
            files.append(p)
    return files, list(range(M0, M0 + len(chains)))


def run_paired(chains: List[str], tmp: str, job_dir: str,
               use_env: bool, pair_strategy: str,
               host: str, ua: str,
               active: bool = True) -> Tuple[List[str], List[int], bool]:
    """提交配对任务 (ticket/pair) 并整理 a3m 到 job_dir/paired/。

    服务器返回一个 pair.a3m (tar 内原始文件名) 含所有链段 (>101, >102,
    ...), 逐行拼接为整复合物 paired.a3m; 多链时每链原始段另拆分为
    paired_0.a3m, paired_1.a3m, ... (单链无需序号后缀)。

    Returns
    -------
    (a3m_files, Ms, had_hits)  a3m_files 为写出的文件列表, Ms 为提交编号,
    had_hits 为服务器是否返回真实命中 (决定是否出图)
    """
    mode_dir = os.path.join(job_dir, "paired")
    os.makedirs(mode_dir, exist_ok=True)

    if not active:
        files, Ms = _write_query_only(chains, mode_dir)
        return files, Ms, False

    tar_gz, Ms = _submit_and_download(
        chains, tmp, "ticket/pair", _paired_mode(use_env, pair_strategy),
        host, ua, tar_name="out_paired.tar.gz")
    _extract_and_clean(tar_gz, tmp)
    src = os.path.join(tmp, "pair.a3m")  # 服务器 tar 内原始文件名
    if not os.path.isfile(src):
        files, Ms = _write_query_only(chains, mode_dir)
        return files, Ms, False

    lines = _read_a3m_lines(src)
    text = _merge_paired_a3m(lines, len(chains))
    dst = os.path.join(mode_dir, "paired.a3m")
    with open(dst, "w") as f:
        f.write(text)
    files = [dst]

    sections = _split_a3m_sections(lines)
    had_hits = len(lines) > len(chains)
    if len(chains) > 1:
        for i in range(len(chains)):
            if i < len(sections):
                p = os.path.join(mode_dir, f"paired_{i}.a3m")
                with open(p, "w") as f:
                    f.write("\n".join(sections[i]) + "\n")
                files.append(p)
    return files, Ms, had_hits
