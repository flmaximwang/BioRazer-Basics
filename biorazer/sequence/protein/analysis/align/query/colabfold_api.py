"""
colabfold_api - 调用 ColabFold MMseqs2 公共 API 生成蛋白 MSA (A3M)

单链:
  from biorazer.sequence.protein.analysis.align.query import run_search
  result = run_search([("my_protein", "MTSENLYFQG...")], "msa_out/")
  # result.per_seq["my_protein"].files  -> 各数据库 a3m 路径
  # result.per_seq["my_protein"].plots  -> [logo.png]

多链 (unpaired, 各链独立搜索):
  result = run_search([("chain_A", "SEQ1..."), ("chain_B", "SEQ2...")], "msa_out/")

多链 (paired, 用于 OpenDDE/AF3 多聚体):
  result = run_search([("chain_A", "SEQ1..."), ("chain_B", "SEQ2...")],
                       "msa_out/", pair_mode="paired")

输出:
  msa_out/
    ├── chain_A/
    │   ├── uniref.a3m       (unpaired 时)
    │   ├── bfd.*.a3m        (环境数据库，启用时)
    │   ├── merged.a3m       (该链所有数据库合并)
    │   ├── logo.png         (sequence logo)
    │   ├── coverage.png     (coverage 热图)
    │   ├── pdb70.m8         (该链模板搜索结果)
    │   └── templates/       (该链的模板 CIF, 已按链拆分)
    ├── chain_B/
    │   └── ...
    └── .template_cache/     (临时缓存，自动删除)
"""

import argparse
import json
import os
import random
import shutil
import sys
import tarfile
import time
from typing import Dict, List, NamedTuple, Optional, Tuple
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

import matplotlib.pyplot as plt
from biotite.sequence.graphics import plot_sequence_logo
from biotite.sequence.profile import SequenceProfile

from ..io import A3M2ALIGN
from ..plot import plot_msa_coverage

# ── 默认值 ──────────────────────────────────────────
DEFAULT_HOST = "https://api.colabfold.com"
DEFAULT_UA = "colabfold_msa/2.0 (standalone; https://github.com/sokrypton/ColabFold)"


# ── HTTP 调用（纯 stdlib）────────────────────────────
def _post(url: str, data: dict, ua: str) -> dict:
    body = "&".join(f"{k}={v}" for k, v in data.items()).encode()
    req = Request(url, data=body, method="POST",
                  headers={"User-Agent": ua,
                           "Content-Type": "application/x-www-form-urlencoded"})
    with urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _get(url: str, ua: str) -> dict:
    req = Request(url, headers={"User-Agent": ua})
    with urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _download(url: str, out_path: str, ua: str) -> None:
    req = Request(url, headers={"User-Agent": ua})
    with urlopen(req, timeout=300) as resp:
        with open(out_path, "wb") as f:
            while True:
                chunk = resp.read(8192)
                if not chunk:
                    break
                f.write(chunk)


# ── API 调用 ──────────────────────────────────────────
def _submit(seqs: List[str], host: str, endpoint: str, mode: str,
            ua: str, N: int = 101) -> Tuple[dict, List[int]]:
    """提交 MSA 任务，带指数退避重试。返回 (响应, Ms) 其中 Ms = [N, N+1, ...]"""
    query = "\n".join(f">{N + i}\n{seq}" for i, seq in enumerate(seqs))
    Ms = list(range(N, N + len(seqs)))
    wait = 2
    for attempt in range(10):
        try:
            return _post(f"{host}/{endpoint}",
                         {"q": query, "mode": mode}, ua), Ms
        except (URLError, HTTPError, OSError) as e:
            print(f"  [!] 提交错误: {e} (尝试 {attempt + 1}/10)", file=sys.stderr)
            time.sleep(wait)
            wait = min(wait * 2, 60)
    raise RuntimeError("提交失败，已达最大重试次数")


def _poll(ticket_id: str, host: str, ua: str) -> dict:
    """查询任务状态，带指数退避重试"""
    wait = 2
    for attempt in range(10):
        try:
            return _get(f"{host}/ticket/{ticket_id}", ua)
        except (URLError, HTTPError, OSError) as e:
            print(f"  [!] 状态查询错误: {e} (尝试 {attempt + 1}/10)", file=sys.stderr)
            time.sleep(wait)
            wait = min(wait * 2, 60)
    raise RuntimeError("状态查询失败，已达最大重试次数")


def _download_ticket(ticket_id: str, out_path: str, host: str, ua: str) -> None:
    """下载结果 tar.gz"""
    wait = 2
    for attempt in range(10):
        try:
            return _download(f"{host}/result/download/{ticket_id}",
                             out_path, ua)
        except (URLError, HTTPError, OSError) as e:
            print(f"  [!] 下载错误: {e} (尝试 {attempt + 1}/10)", file=sys.stderr)
            time.sleep(wait)
            wait = min(wait * 2, 60)
    raise RuntimeError("下载失败，已达最大重试次数")


def _submit_and_download(seqs: List[str], out_dir: str,
                         endpoint: str, mode: str,
                         host: str, ua: str,
                         tar_name: str = "out.tar.gz") -> Tuple[str, List[int]]:
    """提交 → 轮询 → 下载，完整流水线。返回 (tar_gz 路径, Ms)"""
    tar_gz = os.path.join(out_dir, tar_name)

    print(f"  [→] 提交 {len(seqs)} 条序列 [mode={mode}]...", file=sys.stderr)
    N = 101
    timeout_count = 0
    while True:
        out, Ms = _submit(seqs, host, endpoint, mode, ua, N)

        # 速率限制 / 未知 → 等待后重提
        while out.get("status") in ("UNKNOWN", "RATELIMIT"):
            s = 5 + random.randint(0, 5)
            print(f"  [!] {out['status']}，等待 {s}s 后重提...", file=sys.stderr)
            time.sleep(s)
            out, Ms = _submit(seqs, host, endpoint, mode, ua, N)

        if out.get("status") == "ERROR":
            raise RuntimeError("MMseqs2 API 返回 ERROR")
        if out.get("status") == "MAINTENANCE":
            raise RuntimeError("MMseqs2 API 维护中，稍后再试")

        tid = out["id"]
        stat = out
        prev_status = ""
        while stat.get("status") in ("UNKNOWN", "RUNNING", "PENDING"):
            status = stat.get("status", "UNKNOWN")
            if status != prev_status:
                print(f"  [~] MSA: {status}", file=sys.stderr)
                prev_status = status
            time.sleep(5 + random.randint(0, 5))
            stat = _poll(tid, host, ua)

        if stat.get("status") == "COMPLETE":
            print(f"  [✓] MSA: COMPLETE", file=sys.stderr)
            _download_ticket(tid, tar_gz, host, ua)
            return tar_gz, Ms
        elif stat.get("status") == "ERROR":
            timeout_count += 1
            if timeout_count >= 3:
                raise RuntimeError("MMseqs2 任务 3 次超时，放弃重试")
            print(f"  [!] 任务出错/超时 (尝试 {timeout_count}/3)，递增种子重提...",
                  file=sys.stderr)
            N += 1
            continue


# ── 模板下载与链拆分 ────────────────────────────────
def _extract_chain_from_cif(cif_path: str, chain: str, out_path: str) -> bool:
    """从 mmCIF 文件中提取指定链/Entity 的原子记录，写入新 CIF。

    chain 参数可以是:
    - auth_asym_id (如 "C", "a") — 匹配 auth_asym_id 列
    - entity_id (数字, 如 "1") — 提取该 entity 的所有 chain

    自动处理大小写不匹配。
    """
    try:
        with open(cif_path) as f:
            lines = f.readlines()

        # === 第一遍：扫描结构 ===
        # 1) 解析 entity → strand 映射
        entity_strands: dict = {}  # entity_id → [auth_asym_id, ...]
        current_eid = None
        for line in lines:
            if line.startswith("_entity_poly.entity_id"):
                val = line.split(None, 1)
                if len(val) >= 2:
                    current_eid = val[1].strip().rstrip(";")
            elif line.startswith("_entity_poly.pdbx_strand_id") and current_eid is not None:
                val = line.split(None, 1)
                if len(val) >= 2:
                    strand_str = val[1].strip().rstrip(";")
                    entity_strands[current_eid] = [s.strip() for s in strand_str.split(",")]

        # 2) 解析 atom_site 列索引
        label_col = auth_col = -1
        in_atom_loop = False
        cols: List[str] = []
        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; cols = []; continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    col_name = line.split()[0]
                    cols.append(col_name)
                    if col_name == "_atom_site.label_asym_id": label_col = len(cols) - 1
                    if col_name == "_atom_site.auth_asym_id":  auth_col = len(cols) - 1
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False

        # 3) 收集所有可用的 auth_asym_id
        auth_set: set = set()
        in_atom_loop = False; idx = 0
        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; idx = 0; continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    idx += 1
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False
                elif auth_col >= 0 and line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) > auth_col:
                        auth_set.add(parts[auth_col])

        # === 匹配目标 ===
        # chain 参数可能是 auth_asym_id (字母) 或 entity_id (数字)
        target_chains: list = []

        # 先尝试作为 entity_id
        if chain in entity_strands:
            target_chains = entity_strands[chain]
        else:
            # 尝试作为 auth_asym_id (精确匹配)
            if chain in auth_set:
                target_chains = [chain]
            else:
                # 尝试忽略大小写匹配 auth_asym_id
                for ac in auth_set:
                    if ac.lower() == chain.lower():
                        target_chains.append(ac)
                # 最后尝试 label_asym_id 忽略大小写
                if not target_chains and label_col >= 0:
                    label_set: set = set()
                    in_atom_loop = False; idx = 0
                    for line in lines:
                        if line.startswith("loop_"):
                            in_atom_loop = True; idx = 0; continue
                        if in_atom_loop:
                            if line.startswith("_atom_site."): idx += 1
                            elif line.strip()=="" or line.startswith("#"): in_atom_loop=False
                            elif line.startswith(("ATOM","HETATM")):
                                parts=line.split()
                                if len(parts)>label_col: label_set.add(parts[label_col])
                    for lc in label_set:
                        if lc.lower() == chain.lower() and lc not in target_chains:
                            target_chains.append(lc)

        if not target_chains:
            print(f"  [!] CIF {cif_path} 中未找到 chain/entity '{chain}' "
                  f"(auth: {sorted(auth_set)})", file=sys.stderr)
            return False

        # === 第二遍：筛选原子行 ===
        in_atom_loop = False; idx = 0
        auth_col2 = label_col2 = -1
        cols2: List[str] = []
        out: List[str] = []

        def _should_keep(parts: list) -> bool:
            for tc in target_chains:
                if auth_col2 >= 0 and len(parts) > auth_col2 and parts[auth_col2] == tc:
                    return True
                if label_col2 >= 0 and len(parts) > label_col2 and parts[label_col2] == tc:
                    return True
            return False

        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; idx = 0; cols2 = []; auth_col2 = label_col2 = -1
                out.append(line); continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    col_name = line.split()[0]; cols2.append(col_name)
                    if col_name == "_atom_site.label_asym_id": label_col2 = len(cols2) - 1
                    if col_name == "_atom_site.auth_asym_id":  auth_col2 = len(cols2) - 1
                    out.append(line)
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False; out.append(line)
                else:
                    parts = line.split()
                    if _should_keep(parts):
                        out.append(line)
            else:
                out.append(line)

        atom_count = sum(1 for l in out if l.startswith(("ATOM", "HETATM")))
        if atom_count == 0:
            print(f"  [!] CIF {cif_path}: 匹配 chain(s) {target_chains} 后无原子记录",
                  file=sys.stderr)
            return False

        with open(out_path, "w") as f:
            f.writelines(out)
        return True
    except Exception as e:
        print(f"  [!] CIF 链过滤失败 {cif_path} [{chain}]: {e}", file=sys.stderr)
        return False


def _download_and_split_templates(
    chain_templates: Dict[int, List[Tuple[str, Optional[str]]]],
    Ms: List[int],
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    host: str, ua: str,
) -> None:
    """从 RCSB PDB 下载 CIF → 按链拆分到各序列的 templates/ 目录。

    使用 https://files.rcsb.org/download/{pdb_id}.cif 下载。
    自动通过 http_proxy/https_proxy 环境变量使用代理。
    """
    # 收集所有不重复的 PDB ID
    all_pdbs: set = set()
    for entries in chain_templates.values():
        for pdb_id, _ in entries:
            all_pdbs.add(pdb_id.lower())
    if not all_pdbs:
        return

    pdb_list = sorted(all_pdbs)
    total = len(pdb_list)
    print(f"  [→] 下载模板: {total} 个 PDB (RCSB)", file=sys.stderr)

    cache = os.path.join(out_dir, ".template_cache")
    os.makedirs(cache, exist_ok=True)

    # 逐个从 RCSB 下载 CIF
    downloaded = 0
    failed = 0
    rcsb_url = "https://files.rcsb.org/download"
    for pdb_id in pdb_list:
        dst = os.path.join(cache, f"{pdb_id}.cif")
        if os.path.isfile(dst):
            downloaded += 1
            continue
        url = f"{rcsb_url}/{pdb_id}.cif"
        ok = False
        for attempt in range(3):
            try:
                req = Request(url, headers={"User-Agent": ua})
                with urlopen(req, timeout=120) as resp:
                    with open(dst, "wb") as f:
                        while True:
                            chunk = resp.read(8192)
                            if not chunk:
                                break
                            f.write(chunk)
                # 验证是有效 CIF（含 data_ 行）
                with open(dst) as f:
                    first = f.read(200)
                if "data_" in first:
                    ok = True
                    break
                else:
                    os.remove(dst)
                    raise ValueError("不是有效 CIF 文件")
            except Exception as e:
                if attempt < 2:
                    print(f"    [!] {pdb_id} 下载失败 (尝试 {attempt + 1}/3): {e}",
                          file=sys.stderr)
                    time.sleep(3)
        if ok:
            downloaded += 1
        else:
            failed += 1
        # 进度
        done = downloaded + failed
        if done % max(1, total // 10) == 0 or done == total:
            print(f"    [{done}/{total}] 已下载 {downloaded}，失败 {failed}",
                  file=sys.stderr)

    # 为每条序列从 CIF 中提取对应链
    chain_downloaded = 0
    for i, (name, _) in enumerate(named_seqs):
        M = Ms[i] if i < len(Ms) else -1
        if M not in chain_templates:
            continue
        seq_dir = os.path.join(out_dir, name)
        tpl_dir = os.path.join(seq_dir, "templates")
        os.makedirs(tpl_dir, exist_ok=True)
        chain_count = 0
        for pdb_id, chain in chain_templates[M]:
            cif_name = f"{pdb_id}.cif"
            src = os.path.join(cache, cif_name)
            if not os.path.isfile(src):
                continue
            if chain:
                dst = os.path.join(tpl_dir, f"{pdb_id}_{chain}.cif")
                if _extract_chain_from_cif(src, chain, dst):
                    chain_count += 1
            else:
                dst = os.path.join(tpl_dir, cif_name)
                # copy2 会保留元数据，但 CIF 内容更重要
                with open(src) as f_in:
                    content = f_in.read()
                with open(dst, "w") as f_out:
                    f_out.write(content)
                chain_count += 1
        chain_downloaded += chain_count
        if chain_count > 0:
            print(f"    [✓] {name}: {chain_count} 个模板 -> {tpl_dir}", file=sys.stderr)

    # 清理缓存
    shutil.rmtree(cache, ignore_errors=True)
    if chain_downloaded > 0:
        print(f"  [✓] 模板: {chain_downloaded} 个 CIF (按链拆分)", file=sys.stderr)


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


# ── 辅助函数 ──────────────────────────────────────────
def _extract_and_clean(tar_gz: str, out_dir: str) -> None:
    """解压 tar.gz 到 out_dir，解压后删除 tar.gz"""
    if not os.path.isfile(tar_gz):
        raise FileNotFoundError(f"缓存 tar.gz 不存在: {tar_gz}")
    with tarfile.open(tar_gz) as tar:
        tar.extractall(out_dir)
    os.remove(tar_gz)


def _split_a3m_by_chain(a3m_path: str, Ms: List[int],
                         output_dirs: List[str],
                         db_label: str = "uniref") -> Dict[int, str]:
    """将 a3m 文件按链 header 编号拆分成链独立文件。

    Parameters
    ----------
    a3m_path : str      未拆分的 a3m 文件路径
    Ms : List[int]      各链对应的 header 编号（如 [101, 102, 103]）
    output_dirs : List[str]  每条序列对应的输出目录（与 Ms 等长）
    db_label : str      数据库名称，用于文件命名

    Returns
    -------
    {M: chain_a3m_path}  映射
    """
    if not os.path.isfile(a3m_path):
        return {}

    # 逐行读，按 M 分组
    lines_by_M: Dict[int, List[str]] = {}
    current_M: Optional[int] = None
    with open(a3m_path) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                try:
                    current_M = int(line[1:].rstrip().split()[0])
                except ValueError:
                    # 非数字 header（如 >seq_name）— 丢弃
                    continue
                if current_M not in lines_by_M:
                    lines_by_M[current_M] = []
                lines_by_M[current_M].append(line)
            elif current_M is not None:
                # 跳过空字节
                clean_line = line.replace("\x00", "")
                lines_by_M[current_M].append(clean_line)

    # 构建 M → output_dir 映射
    m_to_dir: Dict[int, str] = {}
    for idx, M in enumerate(Ms):
        if idx < len(output_dirs):
            m_to_dir[M] = output_dirs[idx]

    # 写入链独立文件
    written = {}
    for M, chain_dir in m_to_dir.items():
        if M not in lines_by_M:
            continue
        out_path = os.path.join(chain_dir, f"{db_label}.a3m")
        with open(out_path, "w") as f:
            f.writelines(lines_by_M[M])
        written[M] = out_path

    return written


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
        alignment = A3M2ALIGN(a3m_path).read()
        profile = SequenceProfile.from_alignment(alignment)
        total_pos = profile.symbols.shape[0]
        nrows = max(1, (total_pos + cols_per_row - 1) // cols_per_row)

        fig, axes = plt.subplots(nrows, 1, figsize=(row_width, nrows * row_height),
                                 squeeze=False)
        for i in range(nrows):
            ax = axes[i, 0]
            start = i * cols_per_row
            end = min(start + cols_per_row, total_pos)
            chunk = profile[start:end]
            plot_sequence_logo(ax, chunk, scheme="flower")

            # x-axis: tick every tick_every positions, label = global position
            # plot_sequence_logo sets x from 0..(end-start-1)
            ticks = list(range(0, end - start, tick_every))
            labels = [str(start + t + 1) for t in ticks]
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
        alignment = A3M2ALIGN(a3m_path).read()
        plot_msa_coverage(alignment)
        fig = plt.gcf()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    except Exception as e:
        print(f"  [!] coverage 图生成失败: {e}", file=sys.stderr)
        return ""


def run_search(
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    pair_mode: str = "unpaired",
    use_env: bool = True,
    use_filter: bool = True,
    host: str = DEFAULT_HOST,
    ua: str = DEFAULT_UA,
    pair_strategy: str = "greedy",
) -> SearchResult:
    """
    调用 MMseqs2 API，返回结构化搜索结果。

    Parameters
    ----------
    named_seqs : [(name, seq), ...]
        蛋白序列列表（每项为 (序列名, 序列)）
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
    """
    os.makedirs(out_dir, exist_ok=True)

    seqs_only = [s for _, s in named_seqs]

    # ── 确定 endpoint / mode ──────────────────────────
    def _unpaired_mode() -> str:
        if use_filter:
            return "env" if use_env else "all"
        return "env-nofilter" if use_env else "nofilter"

    def _paired_mode() -> str:
        m = "paircomplete" if pair_strategy == "complete" else "pairgreedy"
        return m + "-env" if use_env else m

    # ── 收集结果 ──────────────────────────────────────
    per_seq_result: Dict[str, SeqResult] = {}
    templates_map: Optional[dict] = None

    if pair_mode == "paired":
        tar_gz, Ms = _submit_and_download(seqs_only, out_dir, "ticket/pair",
                                          _paired_mode(), host, ua,
                                          tar_name="out_paired.tar.gz")
        _extract_and_clean(tar_gz, out_dir)
        db_names = ["pair.a3m"]
    elif pair_mode == "paired+unpaired":
        # unpaired
        tar_gz_u, Ms_u = _submit_and_download(seqs_only, out_dir, "ticket/msa",
                                              _unpaired_mode(), host, ua,
                                              tar_name="out_unpaired.tar.gz")
        _extract_and_clean(tar_gz_u, out_dir)
        # paired
        tar_gz_p, Ms_p = _submit_and_download(seqs_only, out_dir, "ticket/pair",
                                              _paired_mode(), host, ua,
                                              tar_name="out_paired.tar.gz")
        _extract_and_clean(tar_gz_p, out_dir)
        # Ms 用于拆分，unpaired 和 paired 的编号一致（都从 N=101 开始）
        Ms = Ms_u
        db_names = ["uniref.a3m", "pair.a3m"]
        if use_env:
            db_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")
    else:  # unpaired
        tar_gz, Ms = _submit_and_download(seqs_only, out_dir, "ticket/msa",
                                          _unpaired_mode(), host, ua)
        _extract_and_clean(tar_gz, out_dir)
        db_names = ["uniref.a3m"]
        if use_env:
            db_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")

    # ── 模板解析 ──────────────────────────────────────
    pdb70_m8 = os.path.join(out_dir, "pdb70.m8")
    chain_templates: Dict[int, List[Tuple[str, Optional[str]]]] = {}  # M → [(pdb_id, chain)]
    chain_m8_lines: Dict[int, List[str]] = {}                         # M → [raw_line]
    if os.path.isfile(pdb70_m8):
        with open(pdb70_m8) as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    M = int(parts[0])
                    target = parts[1]
                    # pdb70 数据库目标 ID 格式: "1abc" 或 "1abc_A"
                    chain: Optional[str] = None
                    pdb_id = target
                    if "_" in target:
                        pdb_id, chain = target.split("_", 1)
                    chain_templates.setdefault(M, []).append((pdb_id, chain))
                    chain_m8_lines.setdefault(M, []).append(line)
        os.remove(pdb70_m8)  # 删除未拆分的 .m8
        print(f"  [→] 模板搜索: {sum(len(v) for v in chain_templates.values())} 个 hit",
              file=sys.stderr)

    # ── 创建每条序列的子目录 ──────────────────────────
    seq_dirs: List[str] = []
    for i, (name, _) in enumerate(named_seqs):
        seq_dir = os.path.join(out_dir, name)
        os.makedirs(seq_dir, exist_ok=True)
        seq_dirs.append(seq_dir)
        per_seq_result[name] = SeqResult(files=[], plots=[], report="")

    # ── 复制 msa.sh 到每条序列的子目录 ────────────────
    msa_sh = os.path.join(out_dir, "msa.sh")
    if os.path.isfile(msa_sh):
        for seq_dir in seq_dirs:
            shutil.copy2(msa_sh, os.path.join(seq_dir, "msa.sh"))
        os.remove(msa_sh)

    # ── 按数据库拆分成链独立文件 ──────────────────────
    for db_name in db_names:
        db_path = os.path.join(out_dir, db_name)
        if not os.path.isfile(db_path):
            continue
        written = _split_a3m_by_chain(db_path, Ms, seq_dirs, db_label=db_name.replace(".a3m", ""))
        for i, (name, _) in enumerate(named_seqs):
            M = Ms[i] if i < len(Ms) else -1
            if M in written:
                result = per_seq_result[name]
                result.files.append(written[M])
        # 清理原始未拆分文件
        os.remove(db_path)

    # ── 拆分 .m8 到每条序列的子目录 ──────────────────
    for i, (name, _) in enumerate(named_seqs):
        M = Ms[i] if i < len(Ms) else -1
        if M not in chain_m8_lines:
            continue
        seq_dir = os.path.join(out_dir, name)
        m8_path = os.path.join(seq_dir, "pdb70.m8")
        with open(m8_path, "w") as f:
            for raw_line in chain_m8_lines[M]:
                f.write(raw_line + "\n")

    # ── 下载模板 → 按链拆分 ─────────────────────────
    _download_and_split_templates(chain_templates, Ms, named_seqs, out_dir, host, ua)

    # ── 生成 per-seq merged.a3m + 绘图 ──────────────
    for i, (name, _) in enumerate(named_seqs):
        result = per_seq_result[name]
        if not result.files:
            continue
        seq_dir = os.path.join(out_dir, name)
        merged_per_seq = os.path.join(seq_dir, "merged.a3m")
        parts = []
        for f in result.files:
            with open(f) as fh:
                parts.append(fh.read().rstrip() + "\n")
        with open(merged_per_seq, "w") as f:
            f.write("".join(parts))
        result.files.append(merged_per_seq)

        # logo
        logo_path = os.path.join(seq_dir, "logo.png")
        _generate_logo(merged_per_seq, logo_path)
        if os.path.isfile(logo_path):
            result.plots.append(logo_path)

        # coverage
        cov_path = os.path.join(seq_dir, "coverage.png")
        _generate_coverage_plot(merged_per_seq, cov_path)
        if os.path.isfile(cov_path):
            result.plots.append(cov_path)

    return SearchResult(per_seq=per_seq_result,
                        merged="",
                        templates=chain_templates if chain_templates else None)


# ── FASTA / 验证 ─────────────────────────────────────
def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """解析 FASTA 文本，返回 [(名称, 序列), ...]

    裸字符串输入 → [("default", 序列)]
    """
    seqs: List[Tuple[str, str]] = []
    cur_name = "default"
    cur_seq: List[str] = []
    for line in text.strip().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if cur_seq:
                seqs.append((cur_name, "".join(cur_seq)))
                cur_seq = []
            # 取 header 第一个空白前的部分作为名称
            cur_name = s[1:].split()[0] if s[1:].strip() else "default"
        else:
            cur_seq.append(s.upper())
    if cur_seq:
        seqs.append((cur_name, "".join(cur_seq)))
    return seqs


def validate(seqs: List[Tuple[str, str]]) -> None:
    """验证序列只含合法氨基酸字符 (接受 [(name, seq)] 格式)"""
    aas = set("ACDEFGHIKLMNPQRSTVWY")
    for i, (name, s) in enumerate(seqs, 1):
        bad = set(s.upper()) - aas
        if bad:
            print(f"错误: 序列 '{name}' 含非法字符: {bad}", file=sys.stderr)
            sys.exit(1)


def merge_a3m(file_list: List[str]) -> str:
    """合并多个 A3M 文件为一个字符串"""
    parts = []
    for f in file_list:
        with open(f) as fh:
            parts.append(fh.read().rstrip() + "\n")
    return "".join(parts)


# ── CLI 子命令注册 ──────────────────────────────────
def register_subcommand(sub) -> argparse.ArgumentParser:
    """在 argparse subparsers 上注册 colabfold-msa 子命令"""
    p = sub.add_parser(
        "colabfold-msa",
        help="调用 ColabFold MMseqs2 API 生成 MSA (A3M)",
        description="调用 ColabFold MMseqs2 公共 API 生成蛋白 MSA (A3M)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True,
                    help="输入: FASTA 文件路径 或 直接传序列字符串 (单链或 > 开头 inline FASTA)")
    p.add_argument("-o", "--output", required=True, help="输出目录")
    p.add_argument("--pair-mode", default="unpaired",
                    choices=["unpaired", "paired", "paired+unpaired"],
                    help="MSA 配对模式 (默认 unpaired)")
    p.add_argument("--no-env", action="store_true", help="不使用环境数据库")
    p.add_argument("--no-filter", action="store_true", help="不做 MSA 过滤")
    p.add_argument("--pair-strategy", default="greedy",
                    choices=["greedy", "complete"], help="配对策略 (默认 greedy)")
    p.add_argument("--host", default=DEFAULT_HOST, help="API 服务器地址 (默认 %(default)s)")
    p.add_argument("--debug", action="store_true", help="打印调试信息")
    p.set_defaults(func=_run_colabfold)
    return p


def _run_colabfold(args) -> None:
    """执行 colabfold-msa 子命令"""
    raw = args.input
    if os.path.isdir(raw):
        entries = os.listdir(raw)
        print(f"错误: '{raw}' 是目录，请指定目录内的序列文件:", file=sys.stderr)
        for e in sorted(entries):
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)
    elif os.path.isfile(raw):
        with open(raw) as f:
            named_seqs = parse_fasta(f.read())
    elif raw.lstrip().startswith(">"):
        named_seqs = parse_fasta(raw)
    else:
        named_seqs = [("default", raw.upper())]
    validate(named_seqs)
    print(f"[→] 读取 {len(named_seqs)} 条序列", file=sys.stderr)

    result = run_search(
        named_seqs=named_seqs,
        out_dir=args.output,
        pair_mode=args.pair_mode,
        use_env=not args.no_env,
        use_filter=not args.no_filter,
        host=args.host,
        pair_strategy=args.pair_strategy,
    )

    print(f"\n[✓] 完成!", file=sys.stderr)
    for name, seq_result in result.per_seq.items():
        print(f"  序列: {name}", file=sys.stderr)
        for f in seq_result.files:
            print(f"    A3M: {f}", file=sys.stderr)
        for p in seq_result.plots:
            print(f"    图:  {p}", file=sys.stderr)
    if result.templates:
        total = sum(len(v) for v in result.templates.values())
        print(f"  模板: {total} 个 PDB", file=sys.stderr)
