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
    │   └── logo.png         (sequence logo)
    ├── chain_B/
    │   └── ...
    ├── merged.a3m           (所有链 × 所有数据库合并)
    └── templates/
        └── pdb70.m8          (模板搜索结果)
"""

import argparse
import json
import os
import random
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
        while stat.get("status") in ("UNKNOWN", "RUNNING", "PENDING"):
            time.sleep(5 + random.randint(0, 5))
            stat = _poll(tid, host, ua)

        if stat.get("status") == "COMPLETE":
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


def _write_templates(templates_map: dict, out_dir: str,
                     host: str, ua: str) -> str:
    """
    从模板服务器下载 PDB 模板，返回模板根目录路径。
    按 API 返回的每组模板提取独立的 tar.gz。
    """
    tpl_dir = os.path.join(out_dir, "templates")
    os.makedirs(tpl_dir, exist_ok=True)
    symlink_paths = {}

    for M, pdb_ids in templates_map.items():
        # 取前 20 个模板（ColabFold 默认限制）
        ids = pdb_ids[:20]
        if not ids:
            continue
        cache_key = "_".join(ids)
        cache_file = os.path.join(tpl_dir, f"templates_{M}.done")
        sub_dir = os.path.join(tpl_dir, f"templates_{M}")
        if os.path.isfile(cache_file):
            symlink_paths[M] = sub_dir
            continue

        tids = ",".join(ids)
        print(f"  [→] 下载模板 [组 {M}]: {len(ids)} 个 PDB", file=sys.stderr)
        wait = 2
        for attempt in range(5):
            try:
                req = Request(f"{host}/template/{tids}",
                              headers={"User-Agent": ua})
                with urlopen(req, timeout=120) as resp:
                    os.makedirs(sub_dir, exist_ok=True)
                    with tarfile.open(fileobj=resp, mode="r|gz") as tar:
                        tar.extractall(path=sub_dir)
                # 创建 A3M 兼容的 symlink（AF3 pipeline 需要）
                a3m_ffindex = os.path.join(sub_dir, "pdb70_a3m.ffindex")
                cs_ffindex = os.path.join(sub_dir, "pdb70_cs219.ffindex")
                cs_ffdata = os.path.join(sub_dir, "pdb70_cs219.ffdata")
                if os.path.isfile(a3m_ffindex) and not os.path.isfile(cs_ffindex):
                    os.symlink("pdb70_a3m.ffindex", cs_ffindex)
                if not os.path.isfile(cs_ffdata):
                    with open(cs_ffdata, "w") as f:
                        f.write("")
                # 标记完成
                with open(cache_file, "w") as f:
                    f.write("done")
                symlink_paths[M] = sub_dir
                print(f"  [✓] 模板组 {M}: {sub_dir}", file=sys.stderr)
                break
            except (URLError, HTTPError, OSError) as e:
                print(f"  [!] 模板下载错误: {e} (尝试 {attempt + 1}/5)",
                      file=sys.stderr)
                time.sleep(wait)
                wait = min(wait * 2, 30)
    return tpl_dir


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


def _generate_logo(a3m_path: str, out_path: str, logo_width_per_aa: float = 0.35,
                   fig_height: float = 3.0) -> str:
    """从 a3m 文件生成 sequence logo 图片。

    使用 biotite 的 SequenceProfile 和 plot_sequence_logo。

    Returns
    -------
    成功时返回 logo 图片路径，失败时返回空字符串。
    """
    try:
        alignment = A3M2ALIGN(a3m_path).read()
        profile = SequenceProfile.from_alignment(alignment)
        seq_len = profile.symbols.shape[0]
        fig_width = max(seq_len * logo_width_per_aa, 4)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        plot_sequence_logo(ax, profile, scheme="flower")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    except Exception as e:
        print(f"  [!] logo 生成失败: {e}", file=sys.stderr)
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
    if os.path.isfile(pdb70_m8):
        templates_map = {}
        with open(pdb70_m8) as f:
            for line in f:
                parts = line.rstrip().split()
                if len(parts) >= 4:
                    M = int(parts[0])
                    pdb = parts[1]
                    templates_map.setdefault(M, []).append(pdb)

    # ── 创建每条序列的子目录 ──────────────────────────
    seq_dirs: List[str] = []
    for i, (name, _) in enumerate(named_seqs):
        seq_dir = os.path.join(out_dir, name)
        os.makedirs(seq_dir, exist_ok=True)
        seq_dirs.append(seq_dir)
        per_seq_result[name] = SeqResult(files=[], plots=[], report="")

    # ── 按数据库拆分成链独立文件 ──────────────────────
    all_merged_parts: List[str] = []
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
                # 同时收集到顶层合并
                with open(written[M]) as f:
                    content = f.read().rstrip()
                    header = f">{name}|{db_name.replace('.a3m','')}\n"
                    all_merged_parts.append(header + content)
        # 清理原始未拆分文件
        os.remove(db_path)

    # ── 生成 per-seq merged.a3m ──────────────────────
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

        # ── 生成 logo ──────────────────────────────────
        logo_path = os.path.join(seq_dir, "logo.png")
        _generate_logo(merged_per_seq, logo_path)
        if os.path.isfile(logo_path):
            result.plots.append(logo_path)

    # ── 生成顶层 merged.a3m ──────────────────────────
    merged_top = ""
    if all_merged_parts:
        merged_top = os.path.join(out_dir, "merged.a3m")
        with open(merged_top, "w") as f:
            f.write("\n".join(all_merged_parts) + "\n")

    return SearchResult(per_seq=per_seq_result,
                        merged=merged_top,
                        templates=templates_map)


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
        if seq_result.report:
            print(f"    报告: {seq_result.report}", file=sys.stderr)
    if result.merged:
        print(f"  合并: {result.merged}", file=sys.stderr)
    if result.templates:
        total = sum(len(v) for v in result.templates.values())
        print(f"  模板: {total} 个 PDB", file=sys.stderr)
