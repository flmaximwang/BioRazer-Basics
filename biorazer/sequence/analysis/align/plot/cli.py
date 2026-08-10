"""biorazer plot 子命令: plot-seqlogo / plot-msa / plot-msa-coverage"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment

from biorazer.sequence.analysis.align.plot import (
    plot_msa,
    plot_msa_coverage,
    plot_seqlogo,
)
from biorazer.sequence.io import FASTA_PROFILE, Fasta_StrDict


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册 plot 相关子命令"""
    _add_seqlogo_parser(sub)
    _add_msa_parser(sub)
    _add_msa_coverage_parser(sub)


def _parse_int_list(spec, what):
    """'1,2,3' -> list[int]; None -> None; 非法输入直接报错退出"""
    if spec is None:
        return None
    try:
        return [int(t) for t in spec.split(",")]
    except ValueError:
        sys.exit(f"error: {what} 需为逗号分隔的整数 (如 1,2,3), got {spec!r}")


def _parse_figsize(spec):
    """'WxH' -> (W, H) inches; None -> None"""
    if spec is None:
        return None
    try:
        w, h = spec.lower().split("x")
        return (float(w), float(h))
    except ValueError:
        sys.exit(f"error: --figsize 需为 WxH 格式 (如 10x5), got {spec!r}")


def _gapped_to_alignment(seqs):
    """等长(可含 '-')序列列表 -> biotite Alignment,供 plot-msa-coverage 使用"""
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        sys.exit("error: 输入序列必须等长 (已比对, 可含 '-' 空位)")
    ungapped = [s.replace("-", "") for s in seqs]
    trace = np.full((length, len(seqs)), -1, dtype=int)
    for i, s in enumerate(seqs):
        pos = 0
        for k, ch in enumerate(s):
            if ch != "-":
                trace[k, i] = pos
                pos += 1
    return Alignment([ProteinSequence(s) for s in ungapped], trace)


def _add_seqlogo_parser(sub):
    p = sub.add_parser(
        "plot-seqlogo",
        help="从等长 FASTA 生成 sequence logo",
        description="从等长蛋白序列 FASTA 生成 sequence logo"
                    " (info/freq 两种模式, 多链以 ':' 分隔)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, metavar="FASTA",
                   help="输入 FASTA (等长序列, 多链以 ':' 分隔)")
    p.add_argument("-o", "--output", metavar="PNG",
                   help="输出 PNG (默认 <input>.logo.png)")
    p.add_argument("--title", default="", help="图标题")
    p.add_argument("--per-line", type=int, default=50, metavar="N",
                   help="每个子图最多 N 个残基, 超长链自动分段 (默认 50)")
    p.add_argument("--width", type=float, default=None,
                   help="画布总宽度 (inch, 与 --height 同时指定)")
    p.add_argument("--height", type=float, default=None,
                   help="画布总高度 (inch, 与 --width 同时指定)")
    p.add_argument("--dpi", type=int, default=300, help="输出分辨率 (默认 300)")
    p.add_argument("--chains", metavar="LIST", default=None,
                   help="只画指定链: 0-based 序号, 逗号分隔 (默认全部)")
    p.add_argument("--ignore-seqs", metavar="LIST", default=None,
                   help="忽略 FASTA 中指定序号的序列: 0-based, 逗号分隔"
                        " (默认全部保留)")
    p.add_argument("--mode", choices=("info", "freq"), default="info",
                   help="堆高含义: info = 信息量 bits; freq = 残基频率,"
                        " 每列堆高恒为 1 (默认 info)")
    p.add_argument("--res-id-shift", metavar="N[,N...]", default=None,
                   help="各链第一个残基的真实编号, 逗号分隔逐链, 单值用于全部链"
                        " (默认 1)")
    p.add_argument("--mark-res-ids", metavar="N[,N...]", default=None,
                   help="标记这些残基 (shift 后真实编号, 逗号分隔;"
                        " 编号跨链重叠时各链均标记)")
    p.set_defaults(func=_run_seqlogo)
    return p


def _run_seqlogo(args) -> None:
    """执行 plot-seqlogo 子命令"""
    ignore = _parse_int_list(args.ignore_seqs, "--ignore-seqs")
    chains = _parse_int_list(args.chains, "--chains")
    try:
        profiles = FASTA_PROFILE(input_io=str(args.input)).read(
            ignore_seqs=ignore, chains=chains
        )
    except (OSError, ValueError) as e:
        sys.exit(f"error: 读取 {args.input} 失败: {e}")
    if not profiles:
        sys.exit("error: 没有可用的序列/链 (检查 --ignore-seqs / --chains)")

    shifts = _parse_int_list(args.res_id_shift, "--res-id-shift")
    if shifts is not None and len(shifts) == 1:
        shifts = shifts[0]  # 单值 -> 全部链 (plot_seqlogo 语义)
    figsize = None
    if args.width is not None or args.height is not None:
        if args.width is None or args.height is None:
            sys.exit("error: --width 与 --height 需同时指定 (或都不指定)")
        figsize = (args.width, args.height)

    try:
        fig, axes = plot_seqlogo(
            profiles,
            mode=args.mode,
            per_line=args.per_line,
            res_id_shifts=shifts,
            mark_res_ids=_parse_int_list(args.mark_res_ids, "--mark-res-ids"),
            title=args.title,
            figsize=figsize,
            dpi=args.dpi,
        )
    except ValueError as e:
        sys.exit(f"error: {e}")
    out = args.output or (str(Path(args.input)) + ".logo.png")
    fig.savefig(out, dpi=args.dpi)
    plt.close(fig)
    print(f"logo 已保存: {out}")


def _add_msa_parser(sub):
    p = sub.add_parser(
        "plot-msa",
        help="多序列比对图",
        description="对输入的未比对序列做多序列比对 (biotite align_multiple) 并绘图",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, metavar="FASTA",
                   help="输入 FASTA (未比对序列)")
    p.add_argument("-o", "--output", required=True, metavar="PNG",
                   help="输出 PNG")
    p.add_argument("--plot-method", choices=("similarity", "type", "default"),
                   default="similarity", help="比对图样式 (默认 similarity)")
    p.add_argument("--symbols-per-line", type=int, default=50, metavar="N",
                   help="每行符号数 (默认 50)")
    p.add_argument("--labels", metavar="LIST", default=None,
                   help="序列标签, 逗号分隔 (默认 1,2,3...)")
    p.add_argument("--figsize", metavar="WxH", default=None,
                   help="画布尺寸 (inch, 如 10x5)")
    p.add_argument("--dpi", type=int, default=300, help="输出分辨率 (默认 300)")
    p.set_defaults(func=_run_msa)
    return p


def _run_msa(args) -> None:
    """执行 plot-msa 子命令"""
    try:
        seq_dict = Fasta_StrDict(input_io=str(args.input)).read()
    except (OSError, ValueError) as e:
        sys.exit(f"error: 读取 {args.input} 失败: {e}")
    sequences = list(seq_dict.values())
    if not sequences:
        sys.exit(f"error: {args.input} 中没有序列")
    labels = None
    if args.labels is not None:
        labels = [t.strip() for t in args.labels.split(",")]
    fig, ax = plot_msa(
        sequences=sequences,
        labels=labels,
        symbols_per_line=args.symbols_per_line,
        plot_method=args.plot_method,
        figsize=_parse_figsize(args.figsize) or (10, 5),
    )
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)
    print(f"比对图已保存: {args.output}")


def _add_msa_coverage_parser(sub):
    p = sub.add_parser(
        "plot-msa-coverage",
        help="MSA 覆盖度图",
        description="可视化多序列比对 (MSA) 的覆盖度 (与查询序列的相似度热图)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input", required=True, metavar="FASTA",
                   help="输入 FASTA (已比对等长序列, 可含 '-' 空位)")
    p.add_argument("-o", "--output", required=True, metavar="PNG",
                   help="输出 PNG")
    p.add_argument("--part-lengths", metavar="N,N,...", default=None,
                   help="各链长度 (如 120,34)")
    p.add_argument("--no-sort", action="store_true",
                   help="不按与查询序列的相似度排序")
    p.add_argument("--figsize", metavar="WxH", default=None,
                   help="画布尺寸 (inch, 如 8x5)")
    p.add_argument("--dpi", type=int, default=300, help="输出分辨率 (默认 300)")
    p.set_defaults(func=_run_msa_coverage)
    return p


def _run_msa_coverage(args) -> None:
    """执行 plot-msa-coverage 子命令"""
    try:
        seq_dict = Fasta_StrDict(input_io=str(args.input)).read()
    except (OSError, ValueError) as e:
        sys.exit(f"error: 读取 {args.input} 失败: {e}")
    seqs = list(seq_dict.values())
    if not seqs:
        sys.exit(f"error: {args.input} 中没有序列")
    msa = _gapped_to_alignment(seqs)
    part_lengths = _parse_int_list(args.part_lengths, "--part-lengths") or []
    plot_msa_coverage(
        msa,
        part_lengths=part_lengths,
        sort_lines=not args.no_sort,
        figsize=_parse_figsize(args.figsize) or (8, 5),
        dpi=args.dpi,
    )
    plt.savefig(args.output, dpi=args.dpi)
    plt.close(plt.gcf())
    print(f"覆盖度图已保存: {args.output}")
