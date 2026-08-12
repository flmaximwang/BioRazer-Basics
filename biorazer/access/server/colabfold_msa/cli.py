"""colabfold_msa.cli - colabfold-msa 子命令的 argparse 定义与执行。

上层 biorazer 入口 (biorazer/cli.py) 通过 register_subcommand(sub) 挂载,
parser/runner 的实现按仓库惯例以 _add_*_parser / _run_* 命名。
"""

import argparse
import os
import sys

from .http import DEFAULT_HOST
from .io import parse_fasta, validate
from .pipeline import run_search


def _add_colabfold_parser(sub) -> argparse.ArgumentParser:
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
    p.add_argument("--local-rcsb-database", default=None, metavar="DIR",
                    help="本地 RCSB 镜像根目录 (divided PDB 布局，如 /mnt/data/public/RCSB)。"
                         "给出时模板结构直接从本地读取，不从 RCSB 下载 (默认 %(default)s)")
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
        local_rcsb_database=args.local_rcsb_database,
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


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册 colabfold-msa 子命令。"""
    _add_colabfold_parser(sub)