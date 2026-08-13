"""colabfold_msa.unpaired - 非配对 (unpaired) MSA 的提交与结果整理。

提交各链到 ticket/msa 端点 (服务器各链独立搜索, 所有链的结果在同一个
uniref.a3m 里按段拼接), 返回后写入 unpaired/ 子目录:
- 服务器返回的原始文件原样保留: uniref.a3m, bfd.*.a3m (use_env 时);
- 按段拆分的每链文件: unpaired_0.a3m, unpaired_1.a3m, ...
  (各库对应链段拼接, 不 padding; 单链时即 unpaired.a3m, 无序号后缀)。

不产 padded 的整复合物合并 a3m — unpaired 按链独立搜索, 合并无意义。

模式未请求 (active=False) 或服务器无结果时, 仅写 query 序列的 fallback
文件, 保证下游路径恒定。
"""

import os
import shutil
from typing import List, Tuple

from .api import _submit_and_download, _unpaired_mode
from .io import _extract_and_clean, _read_a3m_lines, _split_a3m_sections


def _write_query_only(chains: List[str], mode_dir: str, M0: int = 101
                      ) -> Tuple[List[str], List[int]]:
    """模式未请求/无结果时, 仅写 query 序列的 fallback 文件。

    单链 → unpaired.a3m; 多链 → unpaired_0.a3m, unpaired_1.a3m, ...
    """
    files: List[str] = []
    if len(chains) == 1:
        p = os.path.join(mode_dir, "unpaired.a3m")
        with open(p, "w") as f:
            f.write(f">{M0}\n{chains[0]}\n")
        files.append(p)
    else:
        for i, seq in enumerate(chains):
            p = os.path.join(mode_dir, f"unpaired_{i}.a3m")
            with open(p, "w") as f:
                f.write(f">{M0 + i}\n{seq}\n")
            files.append(p)
    return files, list(range(M0, M0 + len(chains)))


def run_unpaired(chains: List[str], tmp: str, job_dir: str,
                 use_env: bool, use_filter: bool,
                 host: str, ua: str,
                 tar_name: str = "out.tar.gz",
                 active: bool = True) -> Tuple[List[str], List[int], bool]:
    """提交非配对任务 (ticket/msa) 并整理每链 a3m 到 job_dir/unpaired/。

    与 paired 同跑 (paired+unpaired) 时传 tar_name="out_unpaired.tar.gz"
    避免与配对结果在同一个 tmp 目录里同名冲突。

    服务器返回一个 uniref.a3m 含所有链段 (>101, >102, ...), 按段拆分
    为每链文件 (各数据库的对应段拼接, 不 padding)。

    Returns
    -------
    (a3m_files, Ms, had_hits)  a3m_files 为写出的文件列表, Ms 为提交编号,
    had_hits 为服务器是否返回真实命中 (决定是否出图)
    """
    mode_dir = os.path.join(job_dir, "unpaired")
    os.makedirs(mode_dir, exist_ok=True)

    if not active:
        files, Ms = _write_query_only(chains, mode_dir)
        return files, Ms, False

    tar_gz, Ms = _submit_and_download(
        chains, tmp, "ticket/msa", _unpaired_mode(use_env, use_filter),
        host, ua, tar_name=tar_name)
    _extract_and_clean(tar_gz, tmp)

    db_names = ["uniref.a3m"]
    if use_env:
        db_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")

    # 保留服务器返回的原始文件 (uniref.a3m / bfd.*.a3m, 原样复制,
    # 不做 padding/合并), 同时按段拆分出每链文件
    raw_paths: List[str] = []
    sections_by_db: List[List[List[str]]] = []
    had_hits = False
    for db in db_names:
        src = os.path.join(tmp, db)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(mode_dir, db)
        shutil.copy2(src, dst)
        raw_paths.append(dst)
        lines = _read_a3m_lines(src)
        sections = _split_a3m_sections(lines)
        sections_by_db.append(sections)
        # 每个库至少每链 1 行 query; 多于链数说明有命中
        if len(lines) > len(chains):
            had_hits = True

    if not sections_by_db:
        # 服务器 tar 无期望的数据库文件 → query-only fallback
        files, Ms = _write_query_only(chains, mode_dir)
        return files, Ms, False

    files: List[str] = list(raw_paths)
    if len(chains) == 1:
        # 单链: 无需序号后缀
        dst = os.path.join(mode_dir, "unpaired.a3m")
        with open(dst, "w") as f:
            for sections in sections_by_db:
                if sections:
                    f.write("\n".join(sections[0]) + "\n")
        files.append(dst)
    else:
        # 每链一段: 各库该链段拼接 (不 padding)
        for i in range(len(chains)):
            parts: List[str] = []
            for sections in sections_by_db:
                if i < len(sections):
                    parts.append("\n".join(sections[i]) + "\n")
            if parts:
                p = os.path.join(mode_dir, f"unpaired_{i}.a3m")
                with open(p, "w") as f:
                    f.write("".join(parts))
                files.append(p)
    return files, Ms, had_hits
