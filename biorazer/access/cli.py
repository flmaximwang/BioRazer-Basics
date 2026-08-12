"""biorazer access fetch 子命令: 从 RCSB / UniProt / AlphaFold DB 下载条目文件。

用法: biorazer fetch <id> --source RCSB --fmt pdb -o 3e2c.pdb
-o 缺省下载到当前目录; 传已存在目录或以 / 结尾视为目录, 否则视为目标文件名。
"""

import argparse
import shutil
import sys
import tempfile
from pathlib import Path

import requests

from biorazer.access.database.AFDB.files import fetch as afdb_fetch
from biorazer.access.database.RCSB.files import fetch as rcsb_fetch
from biorazer.access.database.uniprot.files import fetch as uniprot_fetch
from biorazer.logger import initialize_logger

# 数据来源 -> (fetch 函数, 缺省格式)
SOURCES = {
    "RCSB": (rcsb_fetch, "pdb"),
    "UNIPROT": (uniprot_fetch, "fasta"),
    "AFDB": (afdb_fetch, "cif"),
}


def _normalize_source(source: str) -> str:
    """来源参数大小写不敏感。"""
    return source.upper()


def _resolve_target(identifier: str, fmt: str, output: str | None) -> Path:
    """-o 语义: 省略 -> ./{id}.{fmt}; 目录(已存在或以 / 结尾) -> 目录内 {id}.{fmt}; 否则视为完整文件路径。"""
    if output is None:
        return Path(f"{identifier}.{fmt}")
    raw = str(output)
    if raw.endswith(("/", "\\")) or Path(output).is_dir():
        return Path(output) / f"{identifier}.{fmt}"
    return Path(output)


def _download_to(source_fetch, identifier, fmt, target: Path, overwrite: bool, logger) -> Path:
    """经源 fetch 下载到临时目录, 再移动到目标路径。

    各源 fetch 均写入 {id}.{fmt} 风格文件名, 临时目录中转使 -o 自定义文件名
    无需改动 access 层接口。
    """
    if target.exists() and not overwrite:
        logger.warning(f"{target} already exists, skipping")
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="biorazer_fetch_") as tmp:
        source_fetch(identifier, fmt=fmt, download_dir=tmp, overwrite=True, logger=logger)
        produced = list(Path(tmp).iterdir())
        if len(produced) != 1:
            raise RuntimeError(
                f"fetch 在临时目录产生了 {len(produced)} 个文件, 预期 1 个: {produced}"
            )
        shutil.move(str(produced[0]), target)
    return target


def _run_fetch(args) -> None:
    """执行 fetch 子命令。"""
    fetch, default_fmt = SOURCES[args.source]
    fmt = args.fmt or default_fmt
    target = _resolve_target(args.id, fmt, args.output)
    logger = initialize_logger(__name__)
    try:
        _download_to(fetch, args.id, fmt, target, args.overwrite, logger)
    except ValueError as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
    except requests.exceptions.RequestException as e:
        print(f"网络错误: {e}", file=sys.stderr)
        sys.exit(1)
    print(target)


def _add_fetch_parser(sub) -> argparse.ArgumentParser:
    """在 argparse subparsers 上注册 fetch 子命令。"""
    p = sub.add_parser(
        "fetch",
        help="从 RCSB / UniProt / AlphaFold DB 下载条目文件",
        description="从 access 支持的各数据库下载条目文件到本地",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "id",
        help="数据库条目 ID (如 PDB 代码 3e2c、UniProt 登录号 P0DP23)",
    )
    p.add_argument(
        "-s", "--source",
        required=True,
        type=_normalize_source,
        choices=["RCSB", "UNIPROT", "AFDB"],
        help="数据来源 (缺省格式: RCSB=pdb, UNIPROT=fasta, AFDB=cif)",
    )
    p.add_argument(
        "--fmt",
        default=None,
        help="下载格式, 缺省按来源取默认 (RCSB=pdb / UNIPROT=fasta / AFDB=cif)",
    )
    p.add_argument(
        "-o", "--output",
        default=None,
        help="输出路径: 省略 -> 当前目录 ./{id}.{fmt}; 目录(已存在或以 / 结尾) -> 下载到目录内; 否则视为目标文件名",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="已存在同名文件时强制重新下载",
    )
    p.set_defaults(func=_run_fetch)
    return p


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册 access fetch 子命令。"""
    _add_fetch_parser(sub)
