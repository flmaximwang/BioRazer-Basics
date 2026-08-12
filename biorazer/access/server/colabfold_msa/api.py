"""colabfold_msa.api - MMseqs2 API 任务生命周期: 提交 → 轮询 → 下载。

mode 字符串构造 (_unpaired_mode / _paired_mode) 也在此: 它们是 MMseqs2
API 的协议参数。HTTP 传输细节在 http.py, 重试策略在本模块。
"""

import os
import random
import sys
import time
from typing import List, Tuple
from urllib.error import HTTPError, URLError

from .http import _download, _get, _post


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


# ── endpoint / mode ──────────────────────────────────
def _unpaired_mode(use_env: bool, use_filter: bool) -> str:
    if use_filter:
        return "env" if use_env else "all"
    return "env-nofilter" if use_env else "nofilter"


def _paired_mode(use_env: bool, pair_strategy: str) -> str:
    m = "paircomplete" if pair_strategy == "complete" else "pairgreedy"
    return m + "-env" if use_env else m