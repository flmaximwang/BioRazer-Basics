import argparse
import sys

import numpy as np
import matplotlib.pyplot as plt
from biotite.sequence import ProteinSequence

from biorazer.sequence.io import Fasta_StrDict
from ..util import Alignment


def plot_msa_coverage(
    msa: Alignment,
    part_lengths: list[int] = [],
    sort_lines=True,
    figsize=(8, 5),
    dpi=100,
):
    """
    Visualize the coverage of a multiple sequence alignment (MSA).

    Parameters
    ----------
    msa: Alignment
        Biotite Alignment object containing the MSA data
    part_lengths: list of int
        List of lengths for each segment in the MSA.
        Sometimes your query is merged by multiple chains,
        and this list indicates the lengths of each chain.
    sort_lines: bool
        Whether to sort the sequences based on their similarity to the query sequence
    dpi: int
        Dots per inch for the output figure resolution
    """
    msa_seqs = list(map(str, msa.get_gapped_sequences()))
    if len(part_lengths) == 0:
        part_lengths = [len(msa_seqs[0])]
    Ls = part_lengths

    # 获取查询序列
    seq = msa_seqs.pop(0)
    seq = np.array(list(seq))
    msa_seqs = np.array([list(seq) for seq in msa_seqs])

    # 计算链的累积长度
    Ln = np.cumsum([0] + Ls)
    N = msa_seqs.shape[0]

    # 处理MSA数据
    gap = msa_seqs != "-"  # 检测非空位（21代表空位）
    qid = msa_seqs == seq  # 与查询序列相同的位点

    # 计算链级别的覆盖情况
    gapid = np.stack([gap[:, Ln[i] : Ln[i + 1]].max(-1) for i in range(len(Ls))], -1)

    # 组织绘图数据
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]

        # 计算序列相似度
        seqid = np.stack(
            [qid_[:, Ln[i] : Ln[i + 1]].mean(-1) for i in range(len(Ls))], -1
        ).sum(-1) / (g.sum(-1) + 1e-8)

        # 处理非空位数据
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan

        # 排序处理
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1, None]

        Nn.append(len(lines_))
        lines.append(lines_)

    # 合并数据
    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)

    # 绘制图像
    plt.figure(figsize=figsize, dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(
        lines,
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
        extent=(0, lines.shape[1], 0, lines.shape[0]),
    )

    # 添加链分隔线
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    # 添加序列分组线
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    # 绘制覆盖率曲线
    plt.plot((np.isnan(lines) == False).sum(0), color="black")
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")


# --- plot-msa-coverage 子命令 ---

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
