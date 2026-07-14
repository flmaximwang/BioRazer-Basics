"""
colabfold_api - 调用 ColabFold MMseqs2 公共 API 生成蛋白 MSA (A3M)

单链:
  from biorazer.sequence.protein.analysis.align.query import run_search
  files, tmpl_map = run_search(["MTSENLYFQG..."], "msa_out/")

多链 (unpaired, 各链独立搜索):
  files, tmpl_map = run_search(["CHAIN1", "CHAIN2"], "msa_out/")

多链 (paired, 用于 OpenDDE/AF3 多聚体):
  files, tmpl_map = run_search(["CHAIN1", "CHAIN2"], "msa_out/", pair_mode="paired")

多链 (paired+unpaired, 两者同时生成):
  files, tmpl_map = run_search(["CHAIN1", "CHAIN2"], "msa_out/", pair_mode="paired+unpaired")

输出:
  msa_out/
    ├── uniref.a3m            (单链 或 unpaired)
    ├── bfd.mgnify30.*.a3m    (环境数据库，启用时)
    ├── pair.a3m              (多链 paired 时生成)
    └── pdb70.m8              (模板搜索结果)
"""

import argparse
import json
import os
import random
import sys
import tarfile
import time
from typing import List, Optional, Tuple
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError

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
            ua: str, N: int = 101) -> dict:
    """提交 MSA 任务，带指数退避重试"""
    query = "\n".join(f">{N + i}\n{seq}" for i, seq in enumerate(seqs))
    wait = 2
    for attempt in range(10):
        try:
            return _post(f"{host}/{endpoint}",
                         {"q": query, "mode": mode}, ua)
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
                         tar_name: str = "out.tar.gz") -> str:
    """提交 → 轮询 → 下载，完整流水线。返回 tar_gz 路径。"""
    tar_gz = os.path.join(out_dir, tar_name)

    if os.path.isfile(tar_gz):
        print(f"  [✓] 缓存: {tar_gz}", file=sys.stderr)
        return tar_gz

    print(f"  [→] 提交 {len(seqs)} 条序列 [mode={mode}]...", file=sys.stderr)
    N = 101
    timeout_count = 0
    while True:
        out = _submit(seqs, host, endpoint, mode, ua, N)

        # 速率限制 / 未知 → 等待后重提
        while out.get("status") in ("UNKNOWN", "RATELIMIT"):
            s = 5 + random.randint(0, 5)
            print(f"  [!] {out['status']}，等待 {s}s 后重提...", file=sys.stderr)
            time.sleep(s)
            out = _submit(seqs, host, endpoint, mode, ua, N)

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
            return tar_gz
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


def run_search(seqs: List[str], out_dir: str,
               pair_mode: str = "unpaired",
               use_env: bool = True, use_filter: bool = True,
               host: str = DEFAULT_HOST, ua: str = DEFAULT_UA,
               pair_strategy: str = "greedy") -> Tuple[List[str], Optional[dict]]:
    """
    调用 MMseqs2 API，返回 (A3M 文件路径列表, 模板映射)。
    templates_map 格式: {序列索引: [pdb_id, ...]}

    Parameters
    ----------
    seqs : List[str]
        蛋白序列列表（多链时为各链序列）
    out_dir : str
        输出目录
    pair_mode : str
        配对模式: "unpaired"（单链独立搜）、"paired"（仅配对 MSA）、
        "paired+unpaired"（两者同时生成）
    use_env : bool
        是否使用环境数据库 (BFD/MGnify)，默认 True
    use_filter : bool
        是否做 MSA 过滤，默认 True
    host : str
        MMseqs2 服务器地址
    ua : str
        User-Agent 标识
    pair_strategy : str
        配对策略: "greedy"（快，默认）或 "complete"（全）
    """
    os.makedirs(out_dir, exist_ok=True)

    # ── 确定 endpoint / mode ──────────────────────────
    def _unpaired_mode() -> str:
        if use_filter:
            return "env" if use_env else "all"
        return "env-nofilter" if use_env else "nofilter"

    def _paired_mode() -> str:
        m = "paircomplete" if pair_strategy == "complete" else "pairgreedy"
        return m + "-env" if use_env else m

    def _extract_tar(tar_gz: str) -> None:
        """解压 tar.gz 到 out_dir（如尚未提取）"""
        if not os.path.isfile(tar_gz):
            raise FileNotFoundError(f"缓存 tar.gz 不存在: {tar_gz}")
        with tarfile.open(tar_gz) as tar:
            tar.extractall(out_dir)

    def _parse_templates() -> Optional[dict]:
        """从 pdb70.m8 解析模板映射"""
        pdb70_m8 = os.path.join(out_dir, "pdb70.m8")
        if not os.path.isfile(pdb70_m8):
            return None
        tmpl = {}
        with open(pdb70_m8) as f:
            for line in f:
                parts = line.rstrip().split()
                if len(parts) >= 4:
                    M = int(parts[0])
                    pdb = parts[1]
                    tmpl.setdefault(M, []).append(pdb)
        return tmpl

    if pair_mode == "paired":
        _submit_and_download(seqs, out_dir, "ticket/pair",
                             _paired_mode(), host, ua,
                             tar_name="out_paired.tar.gz")
        _extract_tar(os.path.join(out_dir, "out_paired.tar.gz"))
        need_names = ["pair.a3m"]

    elif pair_mode == "paired+unpaired":
        # unpaired
        _submit_and_download(seqs, out_dir, "ticket/msa",
                             _unpaired_mode(), host, ua,
                             tar_name="out_unpaired.tar.gz")
        _extract_tar(os.path.join(out_dir, "out_unpaired.tar.gz"))
        # paired
        _submit_and_download(seqs, out_dir, "ticket/pair",
                             _paired_mode(), host, ua,
                             tar_name="out_paired.tar.gz")
        _extract_tar(os.path.join(out_dir, "out_paired.tar.gz"))
        need_names = ["uniref.a3m", "pair.a3m"]
        if use_env:
            need_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")

    else:  # unpaired
        _submit_and_download(seqs, out_dir, "ticket/msa",
                             _unpaired_mode(), host, ua)
        _extract_tar(os.path.join(out_dir, "out.tar.gz"))
        need_names = ["uniref.a3m"]
        if use_env:
            need_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")

    # ── 收集结果 ────────────────────────────────────────
    found = [f for f in (os.path.join(out_dir, n) for n in need_names)
             if os.path.isfile(f)]
    templates_map = _parse_templates()
    return found, templates_map


# ── FASTA / 验证 ─────────────────────────────────────
def parse_fasta(text: str) -> List[str]:
    """解析 FASTA 文本，返回序列列表"""
    seqs = []
    cur = []
    for line in text.strip().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(s.upper())
    if cur:
        seqs.append("".join(cur))
    return seqs


def validate(seqs: List[str]) -> None:
    """验证序列只含合法氨基酸字符"""
    aas = set("ACDEFGHIKLMNPQRSTVWY")
    for i, s in enumerate(seqs, 1):
        bad = set(s.upper()) - aas
        if bad:
            print(f"错误: 序列 {i} 含非法字符: {bad}", file=sys.stderr)
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
    # 检测输入类型: 文件 / 目录 / inline FASTA / 裸序列
    raw = args.input
    if os.path.isdir(raw):
        entries = os.listdir(raw)
        print(f"错误: '{raw}' 是目录，请指定目录内的序列文件:", file=sys.stderr)
        for e in sorted(entries):
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)
    elif os.path.isfile(raw):
        with open(raw) as f:
            seqs = parse_fasta(f.read())
    elif raw.lstrip().startswith(">"):
        seqs = parse_fasta(raw)
    else:
        seqs = [raw.upper()]
    validate(seqs)
    print(f"[→] 读取 {len(seqs)} 条序列", file=sys.stderr)

    # 执行搜索
    files, tmpl_map = run_search(
        seqs=seqs,
        out_dir=args.output,
        pair_mode=args.pair_mode,
        use_env=not args.no_env,
        use_filter=not args.no_filter,
        host=args.host,
        pair_strategy=args.pair_strategy,
    )

    print(f"\n[✓] 完成!", file=sys.stderr)
    print(f"    A3M 文件:", file=sys.stderr)
    for f in files:
        print(f"      {f}", file=sys.stderr)
    if tmpl_map:
        total = sum(len(v) for v in tmpl_map.values())
        print(f"    模板: {total} 个 PDB", file=sys.stderr)
