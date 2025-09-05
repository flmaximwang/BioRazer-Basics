import numpy as np
import matplotlib.pyplot as plt
from .alignment import Alignment


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
