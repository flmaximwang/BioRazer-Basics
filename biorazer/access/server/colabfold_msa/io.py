"""colabfold_msa.io - FASTA / A3M 的文件读入、解析与合并。

public: parse_fasta / validate / merge_a3m (CLI 与外部调用的纯函数 API);
private: tar 解压、a3m 行读取、按链分段/拆分等流水线内部工具。
"""

import os
import re
import sys
import tarfile
from typing import List, Optional, Tuple


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
    """验证序列只含合法氨基酸字符 (接受 [(name, seq)] 格式)

    多链序列以 ':' 分隔各链 (ColabFold 多链约定，如 "AAA:BBB")，逐链验证。
    """
    aas = set("ACDEFGHIKLMNPQRSTVWY")
    for i, (name, s) in enumerate(seqs, 1):
        for chain in s.upper().split(":"):
            bad = set(chain) - aas
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


# ── tar.gz 解压 ──────────────────────────────────────
def _extract_and_clean(tar_gz: str, out_dir: str) -> None:
    """解压 tar.gz 到 out_dir，解压后删除 tar.gz"""
    if not os.path.isfile(tar_gz):
        raise FileNotFoundError(f"缓存 tar.gz 不存在: {tar_gz}")
    with tarfile.open(tar_gz) as tar:
        tar.extractall(out_dir)
    os.remove(tar_gz)


# ── A3M 读取 / 分段 / 多链拆分 ───────────────────────
def _read_a3m_lines(path: str) -> List[str]:
    """读取 a3m，去除 \\x00 前缀/后缀与空行。

    服务器返回的 a3m 中多链段头 (如 >102) 前带有 \\x00，行尾也可能有
    \\x00 (mmseqs2 的 null 终止符)，统一清理后再解析。
    """
    lines: List[str] = []
    with open(path, "rb") as f:
        for raw in f:
            l = raw.decode("utf-8", "ignore").rstrip("\n").strip("\x00").lstrip("\x00")
            if l.strip():
                lines.append(l)
    return lines


def _split_a3m_sections(lines: List[str]) -> List[List[str]]:
    """按 query header (纯数字 >N) 把 a3m 分成多链段。

    每次提交的链对应一个段 (提交顺序 = 段顺序)，段内为该链的
    query + 各行 (配对/搜索命中，每行 header + 序列)。
    """
    sections: List[List[str]] = []
    cur: Optional[List[str]] = None
    for l in lines:
        if re.fullmatch(r">\d+", l):
            cur = []
            sections.append(cur)
        if cur is not None:
            cur.append(l)
    return sections