from biotite.sequence import ProteinSequence, SequenceProfile
from biotite.sequence.align import Alignment, align_multiple, SubstitutionMatrix
from biotite.sequence.graphics import (
    plot_alignment_type_based,
    plot_alignment_similarity_based,
    plot_alignment,
    plot_sequence_logo,
)
from biotite.sequence.graphics.colorschemes import get_color_scheme
from biotite.visualize import plot_scaled_text
import numpy as np
import matplotlib.pyplot as plt
from .util import Alignment


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


def plot_msa(
    sequences: list[str | ProteinSequence] = None,
    msa: Alignment | None = None,
    labels=None,
    symbols_per_line=50,
    symbol_spacing=10,
    number_functions=None,
    msa_params: dict = {},
    plot_method: str = "similarity",
    plot_params: dict = {},
    figsize=(10, 5),
):
    """
    Generate a multiple-sequence alignment image

    Parameters
    ----------
    msa_params: dict
        See https://www.biotite-python.org/latest/apidoc/biotite.sequence.align.align_multiple.html
    number_functions: list[callable]
        Example: [lambda i: i+0, lambda i: i+2] # Keep original and add 2 to the position numbers
    plot_params: dict
        See https://www.biotite-python.org/latest/apidoc/biotite.sequence.graphics.plot_alignment_similarity_based.html#biotite.sequence.graphics.plot_alignment_similarity_based
        See https://www.biotite-python.org/latest/apidoc/biotite.sequence.graphics.plot_alignment_type_based.html#biotite.sequence.graphics.plot_alignment_type_based
        See https://www.biotite-python.org/latest/apidoc/biotite.sequence.graphics.plot_alignment.html#biotite.sequence.graphics.plot_alignment
    """

    method_map = {
        "similarity": plot_alignment_similarity_based,
        "type": plot_alignment_type_based,
        "default": plot_alignment,
    }

    # Check Input
    if msa is None and sequences is None:
        raise ValueError("Either 'sequences' or 'msa' must be provided.")
    if msa is not None and sequences is not None:
        raise ValueError("Only one of 'sequences' or 'msa' should be provided.")
    if plot_method not in method_map:
        raise ValueError(
            f"Invalid plot_method '{plot_method}'. Valid options are: {list(method_map.keys())}"
        )

    # Format Input
    if labels is None:
        labels = [f"{i+1}" for i in range(len(sequences))]
        plot_params["labels"] = labels
    plot_params["labels"] = labels
    plot_params["symbols_per_line"] = symbols_per_line
    plot_params["symbol_spacing"] = symbol_spacing
    plot_params["show_numbers"] = True
    plot_params["number_functions"] = number_functions

    formatted_sequences = []
    for sequence in sequences:
        if isinstance(sequence, str):
            formatted_sequences.append(ProteinSequence(sequence))
        elif isinstance(sequence, ProteinSequence):
            formatted_sequences.append(sequence)
        else:
            raise ValueError(
                "Each item in 'sequences' must be either a string or a ProteinSequence."
            )
    print(formatted_sequences)

    # Get MSA
    if msa is not None:
        return msa
    else:

        msa, order, tree, distance_matrix = align_multiple(
            formatted_sequences,
            SubstitutionMatrix.std_protein_matrix(),
            **(msa_params or {}),
        )

    plot_func = method_map[plot_method]
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    plot_func(ax, msa, **plot_params)

    return fig, ax

def _plot_sequence_logo_freq(ax, profile):
    """频率 logo: 字母高度 = 该残基频率,每列堆高恒为 1.0。

    biotite 的 plot_sequence_logo 把堆高硬编码为信息量(bits),
    频率模式复用其低层原语 plot_scaled_text + flower 配色,
    字母高度直接取频率,无小样本偏差,但堆高不再反映保守性。
    """
    colors = get_color_scheme("flower", profile.alphabet)
    counts = profile.symbols.astype(float)
    freq = counts / counts.sum(axis=1, keepdims=True)
    order = np.argsort(freq, axis=1)  # 矮在下、高在上(与 info 模式一致)
    for i in range(len(profile)):
        bottom = 0.0
        for j in order[i]:
            h = freq[i, j]
            if h > 0:
                plot_scaled_text(
                    ax,
                    profile.alphabet.decode(int(j)),
                    i + 0.5,
                    bottom,
                    width=1,
                    height=h,
                    color=colors[j],
                )
                bottom += h
    ax.set_xlim(0.5, len(profile.symbols) + 0.5)
    ax.set_ylim(0, 1.0)


def plot_seqlogo(
    profiles: (
        SequenceProfile
        | Alignment
        | dict[int, SequenceProfile]
        | list[SequenceProfile]
    ),
    mode: str = "info",
    per_line: int = 50,
    res_id_shifts: int | list[int] | None = None,
    mark_res_ids: int | list[int] | None = None,
    title: str = "",
    figsize: tuple[float, float] | None = None,
    dpi: int | None = None,
):
    """
    Draw a sequence logo for one or more chains.

    Each chain becomes a column of the figure; chains longer than
    ``per_line`` positions are split into row segments of at most
    ``per_line`` positions. The layout follows the reference script
    ``seqlogo.py`` of the gammaPFD 6_PocketDesign pipeline.

    Parameters
    ----------
    profiles : SequenceProfile, Alignment, dict[int, SequenceProfile] or list[SequenceProfile]
        The input data. A single :class:`SequenceProfile` (or
        :class:`Alignment`, converted via
        ``SequenceProfile.from_alignment``) is drawn as one chain. A
        dict maps chain indices to profiles (e.g. the return value of
        ``FASTA_PROFILE.read``); chains are drawn in dict order. A list
        is drawn as chains 0, 1, ... in list order.
    mode : {'info', 'freq'}, optional
        ``'info'``: letter heights encode information content in bits
        (Schneider & Stephens 1990, same as biotite's
        ``plot_sequence_logo``): stack height = log2(|alphabet|) -
        positional entropy, so a fully conserved column reaches about
        4.58 bits for the 24-symbol protein alphabet. ``'freq'``:
        letter height = residue frequency and every column sums to
        1.0, reflecting composition rather than conservation. Neither
        mode applies a small-sample correction.
    per_line : int, optional
        Maximum number of positions per row segment (default 50).
    res_id_shifts : int or list of int, optional
        Real residue number of the first position of each chain. A
        single int applies to all chains; a list must contain one
        entry per chain (default: all chains start at 1).
    mark_res_ids : int or list of int, optional
        Residue numbers (after ``res_id_shifts``) to mark with a red
        triangle above the stack. Ids falling into several chains'
        ranges are marked in each of them (default: no marks).
    title : str, optional
        Figure title (default: no title).
    figsize : tuple of float, optional
        Figure size in inches. Defaults to ``per_line * 0.24`` inches
        per column and ~1.07 inches per row segment.
    dpi : int, optional
        Figure resolution; used when the figure is saved without an
        explicit dpi (default: matplotlib default).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure.
    axes : numpy.ndarray of matplotlib.axes.Axes
        Axes grid of shape (n_row_segments, n_chains); cells without
        data are turned off.

    Raises
    ------
    ValueError
        If ``profiles`` is empty or contains an empty chain, ``mode``
        is not ``'info'`` or ``'freq'``, ``per_line`` is not positive,
        or ``res_id_shifts``/``mark_res_ids`` contain non-int values or
        the wrong number of shifts.

    Examples
    --------
    >>> profiles = FASTA_PROFILE(input_io="design.fa").read(chains=[0, 1])
    >>> fig, axes = plot_seqlogo(
    ...     profiles, res_id_shifts=[95, 34], mark_res_ids=[95, 96, 106]
    ... )
    >>> fig.savefig("design_logo.png", dpi=300)
    """
    if mode not in ("info", "freq"):
        raise ValueError(
            f"Invalid mode '{mode}'. Valid options are: 'info', 'freq'"
        )
    if per_line <= 0:
        raise ValueError(f"per_line must be positive, got {per_line}")

    # 归一化输入: 统一的 (chain_id, profile) 有序列表
    if isinstance(profiles, Alignment):
        profiles = SequenceProfile.from_alignment(profiles)
    if isinstance(profiles, SequenceProfile):
        chain_items = [(0, profiles)]
    elif isinstance(profiles, dict):
        chain_items = list(profiles.items())
    else:
        chain_items = list(enumerate(profiles))
    if not chain_items:
        raise ValueError("profiles must contain at least one chain")

    n_chains = len(chain_items)
    chain_lens = [len(profile) for _, profile in chain_items]
    if any(length == 0 for length in chain_lens):
        raise ValueError("profiles contain an empty chain")

    # 各链首个残基的真实编号(单值 -> 全部链)
    if res_id_shifts is None:
        shifts = [1] * n_chains
    elif isinstance(res_id_shifts, int):
        shifts = [res_id_shifts] * n_chains
    else:
        shifts = list(res_id_shifts)
        if len(shifts) != n_chains:
            raise ValueError(
                f"res_id_shifts must have {n_chains} entries (one per "
                f"chain), got {len(shifts)}"
            )
        if any(not isinstance(s, int) for s in shifts):
            raise ValueError(f"res_id_shifts entries must be ints, got {shifts}")

    # 要标记的残基编号(shift 后真实编号)
    if mark_res_ids is None:
        mark_set = set()
    else:
        if isinstance(mark_res_ids, int):
            mark_res_ids = [mark_res_ids]
        mark_set = set(mark_res_ids)
        if any(not isinstance(m, int) for m in mark_set):
            raise ValueError(
                f"mark_res_ids entries must be ints, got {sorted(mark_set)}"
            )

    # 分段: 每行最多 per_line 个残基,列 = 链,行 = 段
    segments = []  # (chain_index, start0, end0)
    n_segs_per_chain = [0] * n_chains
    for ci, length in enumerate(chain_lens):
        for s in range(0, length, per_line):
            segments.append((ci, s, min(s + per_line, length)))
            n_segs_per_chain[ci] += 1
    n_rows = max(n_segs_per_chain)
    n_cols = n_chains

    if figsize is None:
        # 每格 ~0.24 in 宽,行高 = 3.2 的 1/3
        figsize = (
            n_cols * max(6.0, per_line * 0.24),
            n_rows * (3.2 / 3.0),
        )

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        dpi=dpi,
        constrained_layout=True,
        squeeze=False,
    )

    for ci, (chain_id, profile) in enumerate(chain_items):
        col_segs = [(s, e) for (cc, s, e) in segments if cc == ci]
        for r, (s, e) in enumerate(col_segs):
            ax = axes[r, ci]
            sub = profile[s:e]  # SequenceProfile 切片,与整链字母宽度一致
            if mode == "freq":
                _plot_sequence_logo_freq(ax, sub)
            else:
                plot_sequence_logo(ax, sub)  # 默认配色 flower
            # 不足 per_line 的末段也按满宽画,字母宽度与其它段一致
            ax.set_xlim(0.5, per_line + 0.5)
            ax.set_title(
                f"Chain {chain_id} ({shifts[ci] + s}-{shifts[ci] + e - 1})",
                fontsize=9,
            )
            # 字母 (i) 左下角在 i+0.5、宽 1 -> 中心 = i+1;编号 = shift + 段偏移 + i
            ax.set_xticks([i + 1.0 for i in range(e - s)])
            ax.set_xticklabels(
                [str(shifts[ci] + s + i) for i in range(e - s)], fontsize=7
            )
            if mode == "freq":
                ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
            else:
                max_entropy = np.log2(len(profile.alphabet))
                ax.set_yticks(range(0, int(np.ceil(max_entropy)) + 1))
            ax.tick_params(axis="y", labelsize=7)
            ax.spines[["top", "right"]].set_visible(False)
            # 标记指定残基(shift 后编号): 字母上方红色倒三角
            if mark_set:
                if mode == "freq":
                    top_y = 1.0 + 0.08
                    tri_y = top_y - 0.05
                else:
                    top_y = np.log2(len(profile.alphabet)) + 0.4
                    tri_y = top_y - 0.25
                ax.set_ylim(0, top_y)
                for i in range(e - s):
                    if shifts[ci] + s + i in mark_set:
                        ax.plot(
                            i + 1.0,
                            tri_y,
                            marker="v",
                            color="#d62728",
                            markersize=7,
                            zorder=5,
                            clip_on=False,
                        )
            if r == 0:
                ax.set_ylabel(
                    "Frequency" if mode == "freq" else "Information (bits)"
                )
            if r == n_rows - 1:
                ax.set_xlabel("Position")

    # 关闭不足最高列的行数对应的空白子图
    for ci in range(n_chains):
        for r in range(n_segs_per_chain[ci], n_rows):
            axes[r, ci].set_axis_off()

    if title:
        fig.suptitle(title, fontsize=13)

    return fig, axes
