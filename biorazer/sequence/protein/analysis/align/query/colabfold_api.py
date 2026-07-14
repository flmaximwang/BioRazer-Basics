"""
colabfold_api - 调用 ColabFold MMseqs2 公共 API 生成蛋白 MSA (A3M)

单链:
  from biorazer.sequence.protein.analysis.align.query import run_search
  files, tmpl_map = run_search(["MTSENLYFQG..."], "msa_out/", paired=False)

多链 (配对, 用于 OpenDDE/AF3 多聚体):
  files, tmpl_map = run_search(["CHAIN1", "CHAIN2"], "msa_out/", paired=True)

输出:
  msa_out/
    ├── uniref.a3m            (单链 或 unpaired)
    ├── bfd.mgnify30.*.a3m    (环境数据库，启用时)
    ├── pair.a3m              (多链 paired 时生成)
    └── pdb70.m8              (模板搜索结果)
"""

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


def run_search(seqs: List[str], out_dir: str, paired: bool,
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
    paired : bool
        是否做配对 MSA（多链多聚体）
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

    # 确定 endpoint 和 mode
    if paired:
        endpoint = "ticket/pair"
        if pair_strategy == "complete":
            mode = "paircomplete"
        else:
            mode = "pairgreedy"
        if use_env:
            mode += "-env"
    else:
        endpoint = "ticket/msa"
        if use_filter:
            mode = "env" if use_env else "all"
        else:
            mode = "env-nofilter" if use_env else "nofilter"

    tar_gz = os.path.join(out_dir, "out.tar.gz")

    if os.path.isfile(tar_gz):
        print(f"  [✓] 缓存: {tar_gz}", file=sys.stderr)
    else:
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
            job_time = 0
            while stat.get("status") in ("UNKNOWN", "RUNNING", "PENDING"):
                time.sleep(5 + random.randint(0, 5))
                stat = _poll(tid, host, ua)
                job_time += 5 + random.randint(0, 5)

            if stat.get("status") == "COMPLETE":
                _download_ticket(tid, tar_gz, host, ua)
                break
            elif stat.get("status") == "ERROR":
                # 超时未完成 → 递增种子重试
                timeout_count += 1
                if timeout_count >= 3:
                    raise RuntimeError("MMseqs2 任务 3 次超时，放弃重试")
                print(f"  [!] 任务出错/超时 (尝试 {timeout_count}/3)，递增种子重提...",
                      file=sys.stderr)
                N += 1
                continue

    # 解压 A3M 文件
    if paired:
        need_names = ["pair.a3m"]
    else:
        need_names = ["uniref.a3m"]
        if use_env:
            need_names.append("bfd.mgnify30.metaeuk30.smag30.a3m")

    need = [os.path.join(out_dir, n) for n in need_names]
    missing = [f for f in need if not os.path.isfile(f)]
    if missing:
        if not os.path.isfile(tar_gz):
            raise FileNotFoundError(f"缓存 tar.gz 不存在: {tar_gz}")
        with tarfile.open(tar_gz) as tar:
            tar.extractall(out_dir)

    # 检查是否生成 pdb70.m8（模板搜索结果，解压后也存在 tar 中）
    pdb70_m8 = os.path.join(out_dir, "pdb70.m8")
    templates_map = None
    if os.path.isfile(pdb70_m8):
        templates_map = {}
        with open(pdb70_m8) as f:
            for line in f:
                parts = line.rstrip().split()
                if len(parts) >= 4:
                    M = int(parts[0])
                    pdb = parts[1]
                    if M not in templates_map:
                        templates_map[M] = []
                    templates_map[M].append(pdb)

    found = [f for f in need if os.path.isfile(f)]
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
