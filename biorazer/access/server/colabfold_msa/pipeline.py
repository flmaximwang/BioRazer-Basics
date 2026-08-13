"""colabfold_msa.pipeline - MSA 任务流水线编排与结果整理。

单链:
  from biorazer.access.server.colabfold_msa.pipeline import run_search
  result = run_search([("my_protein", "MTSENLYFQG...")], "msa_out/")
  # result.per_seq["my_protein"].files  -> 各数据库 a3m 路径
  # result.per_seq["my_protein"].plots  -> [logo.png]

多链 complex (一条记录, ':' 分隔各链, ColabFold 约定; 用于 paired):
  result = run_search([("complex", "SEQ1...:SEQ2...")], "msa_out/",
                       pair_mode="paired")
  # 提交时每条链拆成独立记录 (>101, >102, ...), 服务器端配对,
  # 结果合并为整复合物的一份 a3m: msa_out/complex/paired/paired.a3m

多条记录 = 多个独立任务 (各自提交, 互不配对):
  result = run_search([("prot_A", "SEQ1..."), ("prot_B", "SEQ2...")], "msa_out/")
  # 输出 msa_out/prot_A/, msa_out/prot_B/

输出 (按模式分目录, 模式未请求/无结果时该目录仍建, 仅含 query 序列;
a3m 均为每链拆分文件, 不产合并形态):
  msa_out/
    ├── seq_name/
    │   ├── unpaired/            (unpaired 模式结果)
    │   │   ├── uniref.a3m       (服务器原始返回, 原样保留)
    │   │   ├── bfd.*.a3m        (环境数据库原始返回, 启用时)
    │   │   ├── unpaired_0.a3m, unpaired_1.a3m, ...  (每链拆分; 单链时即 unpaired.a3m)
    │   │   ├── logo_0.png / coverage_0.png, ...     (每链各一张; 单链时 logo.png)
    │   │   ├── pdb70.m8         (服务器原始返回, 原样保留)
    │   │   ├── pdb70_0.m8, pdb70_1.m8, ...  (模板 hit 按链拆分; 单链时无)
    │   │   └── templates/       (模板 CIF, 按链命名 {pdb}_{chain}.cif)
    │   ├── paired/              (paired 模式结果)
    │   │   ├── paired.a3m       (整复合物配对 MSA, 各链段逐行拼接)
    │   │   ├── paired_0.a3m, paired_1.a3m, ...  (每链原始段; 单链时不拆分)
    │   │   ├── logo.png / coverage.png   (整复合物, 便于对比各链 MSA 差异)
    │   │   ├── logo_0.png / coverage_0.png, ...  (每链各一张; 单链时无 _N)
    │   └── msa.sh
    └── .template_cache/         (临时缓存, 自动删除)
"""

import os
import shutil
import sys
import tempfile
from typing import Dict, List, NamedTuple, Optional, Tuple

import matplotlib.pyplot as plt
from biotite.sequence.graphics import plot_sequence_logo
from biotite.sequence.profile import SequenceProfile

from biorazer.sequence.analysis.alignment.plot import plot_msa_coverage
from biorazer.sequence.io import A3m_Alignment

from .http import DEFAULT_HOST, DEFAULT_UA
from .paired import run_paired
from .template import _handle_templates
from .unpaired import run_unpaired


# ── 返回类型 ──────────────────────────────────────────
class SeqResult(NamedTuple):
    """每条序列的搜索结果"""
    files: List[str]      # [uniref.a3m, bfd.*.a3m, merged.a3m]
    plots: List[str]      # [logo.png, ...]
    report: str           # report.txt 路径 (空字符串表示无)

class SearchResult(NamedTuple):
    """完整的搜索结果"""
    per_seq: Dict[str, SeqResult]  # seq_name → SeqResult
    merged: str                    # 顶层 merged.a3m 路径 (空字符串表示无)
    templates: Optional[dict]      # pdb70.m8 映射


# ── 出图 ──────────────────────────────────────────────
def _generate_logo(a3m_path: str, out_path: str,
                   cols_per_row: int = 50,
                   tick_every: int = 10,
                   row_height: float = 2.5,
                   row_width: float = 18) -> str:
    """从 a3m 文件生成多行 sequence logo 图片。

    每行显示 cols_per_row 个位置，每 tick_every 个位置一个刻度。

    Returns
    -------
    成功时返回 logo 图片路径，失败时返回空字符串。
    """
    try:
        alignment = A3m_Alignment(a3m_path).read()
        profile = SequenceProfile.from_alignment(alignment)
        # Add gap column so gap-only positions display as '-' in the logo
        import numpy as np
        from biotite.sequence.alphabet import LetterAlphabet
        from biotite.sequence.graphics import get_color_scheme
        orig_alph = list(profile.alphabet)
        alph_with_gap = LetterAlphabet(orig_alph + ['-'])
        n_aa = len(orig_alph)
        sym_with_gap = np.zeros((profile.symbols.shape[0], n_aa + 1), dtype=int)
        sym_with_gap[:, :n_aa] = profile.symbols
        sym_with_gap[:, n_aa] = profile.gaps
        gap_color = '#aaaaaa'
        colors = list(get_color_scheme('flower', profile.alphabet)) + [gap_color]
        profile = SequenceProfile(sym_with_gap, np.zeros_like(profile.gaps), alph_with_gap)
        total_pos = profile.symbols.shape[0]
        nrows = max(1, (total_pos + cols_per_row - 1) // cols_per_row)

        fig, axes = plt.subplots(nrows, 1, figsize=(row_width, nrows * row_height),
                                 squeeze=False)
        for i in range(nrows):
            ax = axes[i, 0]
            start = i * cols_per_row
            end = min(start + cols_per_row, total_pos)
            chunk = profile[start:end]
            n_actual = end - start
            if n_actual < cols_per_row:
                pad_width = cols_per_row - n_actual
                pad_sym = np.zeros((pad_width, chunk.symbols.shape[1]))
                padded_sym = np.vstack([chunk.symbols, pad_sym])
                pad_gaps = np.zeros(pad_width, dtype=chunk.gaps.dtype)
                padded_gaps = np.concatenate([chunk.gaps, pad_gaps])
                chunk = SequenceProfile(padded_sym, padded_gaps, chunk.alphabet)
            plot_sequence_logo(ax, chunk, scheme=colors)

            # x-axis: tick every tick_every positions, label = global position
            # plot_sequence_logo sets x from 0..(cols_per_row-1)
            ticks_all = list(range(tick_every, cols_per_row + 1, tick_every))
            ticks, labels = [], []
            for t in ticks_all:
                pos = start + t
                if pos <= total_pos:
                    ticks.append(t)
                    labels.append(str(pos))
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels, fontsize=6)

            if i < nrows - 1:
                ax.set_xlabel("")
            # preserve y info; the logo itself shows bit-conservation height

        plt.tight_layout()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    except Exception as e:
        print(f"  [!] logo 生成失败: {e}", file=sys.stderr)
        return ""


def _generate_coverage_plot(a3m_path: str, out_path: str) -> str:
    """从 a3m 文件生成 coverage 热图。

    Returns
    -------
    成功时返回图片路径，失败时返回空字符串。
    """
    try:
        alignment = A3m_Alignment(a3m_path).read()
        plot_msa_coverage(alignment)
        fig = plt.gcf()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    except Exception as e:
        print(f"  [!] coverage 图生成失败: {e}", file=sys.stderr)
        return ""


# ── 单任务 (一条记录) ────────────────────────────────
def _run_one_job(name: str, seq: str, out_dir: str,
                 pair_mode: str, use_env: bool, use_filter: bool,
                 host: str, ua: str, pair_strategy: str,
                 local_rcsb_database: Optional[str],
                 fetch_templates: bool = True) -> SeqResult:
    """提交一条序列 (记录) 并整理该任务的结果到 out_dir/<name>/。

    多链 complex (':' 分隔) 时把每条链拆成独立记录提交 (ColabFold
    服务器在记录之间做配对)，返回后把各链段合并为整复合物 a3m。
    """
    chains = seq.split(":") if ":" in seq else [seq]
    job_dir = os.path.join(out_dir, name)
    os.makedirs(job_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(dir=out_dir, prefix=".tmp_")
    try:
        a3m_files: List[str] = []
        Ms: List[int] = []
        paired_had_hits = False
        unpaired_had_hits = False

        # ── paired 部分 ──────────────────────────────
        if pair_mode in ("paired", "paired+unpaired"):
            files, Ms, paired_had_hits = run_paired(
                chains, tmp, job_dir, use_env, pair_strategy, host, ua)
            a3m_files.extend(files)
        else:
            # 模式未请求: 仍建 paired/ 子目录, 仅写 query 序列
            files, _Ms, _ = run_paired(chains, tmp, job_dir, use_env,
                                       pair_strategy, host, ua, active=False)
            a3m_files.extend(files)

        # ── unpaired 部分 ────────────────────────────
        if pair_mode in ("unpaired", "paired+unpaired"):
            tar_name = "out_unpaired.tar.gz" if pair_mode == "paired+unpaired" else "out.tar.gz"
            files, Ms, unpaired_had_hits = run_unpaired(
                chains, tmp, job_dir, use_env, use_filter, host, ua,
                tar_name=tar_name)
            a3m_files.extend(files)
        else:
            # 模式未请求: 仍建 unpaired/ 子目录, 仅写 query 序列
            files, _Ms, _ = run_unpaired(chains, tmp, job_dir, use_env,
                                         use_filter, host, ua, active=False)
            a3m_files.extend(files)

        # ── 模板 (仅 msa ticket 附带 pdb70.m8) ────────
        _handle_templates(tmp, out_dir, Ms, [(name, seq)], host, ua,
                          local_rcsb_database, fetch_templates)

        # ── msa.sh ────────────────────────────────────
        msa_sh = os.path.join(tmp, "msa.sh")
        if os.path.isfile(msa_sh):
            shutil.copy2(msa_sh, os.path.join(job_dir, "msa.sh"))

        # ── logo + coverage (拆分后的每链文件 + 整复合物 paired.a3m) ──
        result = SeqResult(files=list(a3m_files), plots=[], report="")
        for f in a3m_files:
            mode_dir = os.path.dirname(f)
            base = os.path.basename(f)
            # 出图对象: 每链拆分文件 (unpaired_N/paired_N, 单链无后缀)
            # + 整复合物 paired.a3m (便于对比各链 MSA 差异)
            if not (base.startswith("unpaired_") or base.startswith("paired_")
                    or base == "unpaired.a3m" or base == "paired.a3m"):
                continue
            # 无真实命中 (仅 query 序列) 不出图
            if base.startswith("unpaired_") or base == "unpaired.a3m":
                if not unpaired_had_hits:
                    continue
            else:
                if not paired_had_hits:
                    continue
            logo_path = os.path.join(mode_dir, "logo.png")
            cov_path = os.path.join(mode_dir, "coverage.png")
            # 每链拆分文件: 图名带链序号, 如 logo_0.png
            if base.startswith("unpaired_") or base.startswith("paired_"):
                idx = base.rsplit("_", 1)[1].split(".")[0]
                logo_path = os.path.join(mode_dir, f"logo_{idx}.png")
                cov_path = os.path.join(mode_dir, f"coverage_{idx}.png")
            if _generate_logo(f, logo_path):
                result.plots.append(logo_path)
            if _generate_coverage_plot(f, cov_path):
                result.plots.append(cov_path)
        return result
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_search(
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    pair_mode: str = "unpaired",
    use_env: bool = True,
    use_filter: bool = True,
    host: str = DEFAULT_HOST,
    ua: str = DEFAULT_UA,
    pair_strategy: str = "greedy",
    local_rcsb_database: Optional[str] = None,
    fetch_templates: bool = True,
) -> SearchResult:
    """
    调用 MMseqs2 API，返回结构化搜索结果。

    Parameters
    ----------
    named_seqs : [(name, seq), ...]
        蛋白序列列表（每项为 (序列名, 序列)）。每条记录 = 一个独立任务
        (各自提交):
        - 单链: 序列本身;
        - 多链 complex: 序列以 ':' 分隔各链 (ColabFold 约定, 如
          "AAA:BBB")。提交时每条链拆成一个记录 (>101, >102, ...),
          服务器端配对, 结果合并为该 complex 的一份 a3m。
        注意: 多条记录之间互不配对 (各自独立搜索)。
    out_dir : str
        输出目录
    pair_mode : str
        "unpaired" / "paired" / "paired+unpaired"
    use_env : bool
        是否使用环境数据库，默认 True
    use_filter : bool
        是否过滤 MSA，默认 True
    host : str
        MMseqs2 服务器地址
    ua : str
        User-Agent
    pair_strategy : str
        "greedy"（快）或 "complete"（全）
    local_rcsb_database : str, optional
        本地 RCSB 镜像根目录 (divided PDB 布局，如 /mnt/data/public/RCSB)。
        给出时模板结构直接从本地读取，不从 RCSB 下载。
    fetch_templates : bool
        是否获取模板结构 (RCSB 下载或本地镜像读取)，默认 True。
        False 时只保留服务器返回的 pdb70.m8 (原样 + 按链拆分),
        不下载/不读取模板 CIF, 不产出 templates/ 目录。
    """
    os.makedirs(out_dir, exist_ok=True)

    per_seq_result: Dict[str, SeqResult] = {}
    for name, seq in named_seqs:
        print(f"[→] 任务 {name}: {seq.count(':') + 1} 条链", file=sys.stderr)
        per_seq_result[name] = _run_one_job(
            name, seq, out_dir,
            pair_mode=pair_mode, use_env=use_env, use_filter=use_filter,
            host=host, ua=ua, pair_strategy=pair_strategy,
            local_rcsb_database=local_rcsb_database,
            fetch_templates=fetch_templates,
        )
    return SearchResult(per_seq=per_seq_result, merged="", templates=None)