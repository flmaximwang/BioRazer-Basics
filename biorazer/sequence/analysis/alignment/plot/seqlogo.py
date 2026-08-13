import argparse
import bisect
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from biotite.sequence import SequenceProfile
from biotite.sequence.align import Alignment
from biotite.sequence.graphics import plot_sequence_logo
from biotite.sequence.graphics.colorschemes import get_color_scheme
from biotite.visualize import plot_scaled_text
from matplotlib.transforms import ScaledTranslation

from biorazer.sequence.io import Fasta_Profile


def _pad_profile_left(profile, n_empty):
    """在 profile 左侧补 n_empty 个全空列 (符号计数与 gap 均为 0)。

    绝对网格模式下把窗口内第一个字母右移到它的 tick 位置: 字母由
    plot_sequence_logo 从切片第 0 列画起, 只有把前导空列补进切片,
    字母才会落在绝对网格上。全空列不画任何字母 (biotite 内部对
    全空列的高度为 NaN, 调用处用 np.errstate 屏蔽警告)。
    """
    if n_empty == 0:
        return profile
    n_sym = profile.symbols.shape[1]
    symbols = np.concatenate(
        [
            np.zeros((n_empty, n_sym), dtype=profile.symbols.dtype),
            profile.symbols,
        ],
        axis=0,
    )
    gaps = np.concatenate(
        [np.zeros(n_empty, dtype=profile.gaps.dtype), profile.gaps], axis=0
    )
    return SequenceProfile(symbols, gaps, profile.alphabet)


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


def _column_entropies(profile):
    """每列 Shannon 熵 (bits), 基于全字母表残基频率。

    H_i = -sum_j f_ij * log2(f_ij), f_ij = 第 i 列符号 j 计数 / 第 i 列总计数,
    与 info 模式堆高同源 (堆高 = log2(|alphabet|) - H_i)。全 0 空列记为 NaN
    (不参与任何标注)。
    """
    counts = profile.symbols.astype(float)
    totals = counts.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        freq = counts / totals[:, None]
        # 0*log2(0) 按惯例取 0; 空列 freq 为 NaN, 下面统一盖掉
        safe = np.where(freq > 0, freq, 1.0)
        ent = -(safe * np.log2(safe)).sum(axis=1)
    ent[totals == 0] = np.nan
    return ent


def _resolve_renumber_res(renumber_res, chain_ids):
    """归一化 renumber_res -> {chain_id: {position: res_id}} 锚点表。

    只接受 dict: {chain_id: int | {position: res_id}}。纯整数 / 列表形式
    已废除 (CLI 上单链用 P_ID, 多链用 C_P_ID, 均解析为 dict)。
    编号规则: 第 p 个残基 (0-based) 的编号 = 最近一个 <= p 的锚点 (pos, id)
    的 id + (p - pos); 第一个锚点 (P) 之前的残基不重编号, 保持默认从
    1 起连续编号。同一链多个锚点
    可使编号跳变 (gap), 如 {0: {0: 1, 20: 31}} -> 编号 1..20, 31, 32, ...
    """
    if renumber_res is None:
        return {}
    if not isinstance(renumber_res, dict):
        raise ValueError(
            "renumber_res must be a dict of {chain_id: int | {position: "
            f"res_id}}, got {renumber_res!r}"
        )
    anchors_by_chain = {}
    for cid, spec in renumber_res.items():
        if not isinstance(cid, int):
            raise ValueError(
                f"renumber_res dict keys must be int chain ids, got {cid!r}"
            )
        if isinstance(spec, int):
            anchors_by_chain[cid] = {0: spec}
        elif isinstance(spec, dict):
            anchors = {}
            for pos, rid in spec.items():
                if not isinstance(pos, int) or pos < 0:
                    raise ValueError(
                        f"renumber_res anchor positions must be "
                        f"non-negative ints, got {pos!r}"
                    )
                if not isinstance(rid, int):
                    raise ValueError(
                        f"renumber_res anchor ids must be ints, got {rid!r}"
                    )
                anchors[pos] = rid
            if not anchors:
                raise ValueError(
                    f"renumber_res for chain {cid} has no anchors"
                )
            anchors_by_chain[cid] = anchors
        else:
            raise ValueError(
                f"renumber_res for chain {cid} must be an int or a "
                f"dict of {{position: res_id}}, got {spec!r}"
            )
    unknown = set(anchors_by_chain) - set(chain_ids)
    if unknown:
        raise ValueError(
            f"renumber_res references unknown chain ids {sorted(unknown)}; "
            f"known: {list(chain_ids)}"
        )
    return anchors_by_chain


def _resolve_first_tick_id(first_tick_id, chain_ids):
    """归一化 first_tick_id -> {chain_id: int} (绝对网格起始编号)。

    None -> {}; int 只允许单链 (多链校验在 plot_seqlogo 主体);
    dict {chain_id: int} 逐链指定 (多链必须用此形式)。
    """
    if first_tick_id is None:
        return {}
    if isinstance(first_tick_id, int):
        return {cid: first_tick_id for cid in chain_ids}
    if isinstance(first_tick_id, dict):
        ftids = {}
        for cid, n in first_tick_id.items():
            if not isinstance(cid, int):
                raise ValueError(
                    f"first_tick_id dict keys must be int chain ids, "
                    f"got {cid!r}"
                )
            if not isinstance(n, int):
                raise ValueError(
                    f"first_tick_id for chain {cid} must be an int, "
                    f"got {n!r}"
                )
            ftids[cid] = n
        unknown = set(ftids) - set(chain_ids)
        if unknown:
            raise ValueError(
                f"first_tick_id references unknown chain ids "
                f"{sorted(unknown)}; known: {list(chain_ids)}"
            )
        return ftids
    raise ValueError(
        "first_tick_id must be an int or a dict of {chain_id: int}, "
        f"got {first_tick_id!r}"
    )


def _resolve_mark_res_ids(mark_res_ids, chain_ids):
    """归一化 mark_res_ids -> {chain_id: set[int]} (逐链标记)。

    None -> {}; int / list[int] 只允许单链 (作用于唯一链);
    dict {chain_id: int | list[int]} 逐链指定 (多链必须用此形式)。
    """
    if mark_res_ids is None:
        return {}
    if isinstance(mark_res_ids, dict):
        marks = {}
        for cid, m in mark_res_ids.items():
            if not isinstance(cid, int):
                raise ValueError(
                    f"mark_res_ids dict keys must be int chain ids, "
                    f"got {cid!r}"
                )
            if isinstance(m, int):
                vals = [m]
            elif isinstance(m, (list, tuple)):
                vals = list(m)
            else:
                raise ValueError(
                    f"mark_res_ids for chain {cid} must be an int or a "
                    f"list of int, got {m!r}"
                )
            if any(not isinstance(v, int) for v in vals):
                raise ValueError(
                    f"mark_res_ids for chain {cid} entries must be ints, "
                    f"got {vals}"
                )
            marks[cid] = set(vals)
        unknown = set(marks) - set(chain_ids)
        if unknown:
            raise ValueError(
                f"mark_res_ids references unknown chain ids "
                f"{sorted(unknown)}; known: {list(chain_ids)}"
            )
        return marks
    if isinstance(mark_res_ids, int):
        mark_res_ids = [mark_res_ids]
    if isinstance(mark_res_ids, (list, tuple)):
        vals = list(mark_res_ids)
        if any(not isinstance(v, int) for v in vals):
            raise ValueError(
                f"mark_res_ids entries must be ints, got {vals}"
            )
        # 无链号形式只允许单链 (多链校验在 plot_seqlogo 主体)
        return {chain_ids[0]: set(vals)}
    raise ValueError(
        "mark_res_ids must be an int, a list of int, or a dict of "
        f"{{chain_id: int | list[int]}}, got {mark_res_ids!r}"
    )


def _resolve_res_id_range(res_id_range, chain_ids):
    """归一化 res_id_range -> {chain_id: (start, end)} (含首尾)。"""
    if res_id_range is None:
        return {}
    if isinstance(res_id_range, (tuple, list)) and len(res_id_range) == 2:
        start, end = res_id_range
        if not (isinstance(start, int) and isinstance(end, int)):
            raise ValueError(
                f"res_id_range bounds must be ints, got {res_id_range}"
            )
        if start > end:
            raise ValueError(
                f"res_id_range start must be <= end, got {res_id_range}"
            )
        return {cid: (start, end) for cid in chain_ids}
    if isinstance(res_id_range, dict):
        ranges = {}
        for cid, r in res_id_range.items():
            if not isinstance(cid, int):
                raise ValueError(
                    f"res_id_range dict keys must be int chain ids, got {cid!r}"
                )
            if not (isinstance(r, (tuple, list)) and len(r) == 2):
                raise ValueError(
                    f"res_id_range for chain {cid} must be a (start, end) "
                    f"pair, got {r!r}"
                )
            start, end = r
            if not (isinstance(start, int) and isinstance(end, int)):
                raise ValueError(f"res_id_range bounds must be ints, got {r}")
            if start > end:
                raise ValueError(f"res_id_range start must be <= end, got {r}")
            ranges[cid] = (start, end)
        unknown = set(ranges) - set(chain_ids)
        if unknown:
            raise ValueError(
                f"res_id_range references unknown chain ids "
                f"{sorted(unknown)}; known: {list(chain_ids)}"
            )
        return ranges
    raise ValueError(
        "res_id_range must be a (start, end) pair or a dict of "
        f"(start, end), got {res_id_range!r}"
    )


def _res_ids_for_chain(anchors, length):
    """锚点表 -> 每位置编号列表 (0-based 位置, 长度 length)。

    编号 = 最近一个 <= p 的锚点 id + (p - pos); 第一个锚点 (P) 之前的
    残基不重编号, 保持默认从 1 起连续编号 (p + 1)。锚点间及锚点衔接处
    的编号必须递增 (顺延), 回退视为非法输入。
    """
    if not anchors:
        return list(range(1, length + 1))
    positions = sorted(anchors)
    first = positions[0]
    # 第一个锚点之前的残基保持默认编号 (1..first), 故锚点编号必须
    # >= first + 1, 否则锚点处编号回退 (如 {5: 4} -> 1,2,3,4,5,4,5,...)
    if anchors[first] < first + 1:
        raise ValueError(
            f"renumber_res 锚点编号必须沿位置递增: 锚点前默认编号到 "
            f"{first}, 锚点 {first}->{anchors[first]} 回退"
        )
    for a, b in zip(positions, positions[1:]):
        if anchors[b] < anchors[a] + (b - a):
            raise ValueError(
                f"renumber_res 锚点编号必须沿位置递增: "
                f"{a}->{anchors[a]} 到 {b}->{anchors[b]} 回退"
            )
    res_ids = []
    for p in range(length):
        if p < first:
            res_ids.append(p + 1)
            continue
        i = bisect.bisect_right(positions, p) - 1
        ap = positions[i]
        res_ids.append(anchors[ap] + (p - ap))
    return res_ids


def plot_seqlogo(
    profiles: (
        SequenceProfile
        | Alignment
        | dict[int, SequenceProfile]
        | list[SequenceProfile]
    ),
    mode: str = "info",
    per_line: int = 50,
    renumber_res: dict[int, int | dict[int, int]] | None = None,
    res_id_range: tuple[int, int] | dict[int, tuple[int, int]] | None = None,
    first_tick_id: int | dict[int, int] | None = None,
    mark_res_ids: (
        int | list[int] | dict[int, int | list[int]] | None
    ) = None,
    mark_conservative_res: bool = True,
    entropy_cutoff: float = 0.0,
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
    renumber_res : dict, optional
        Real residue numbering per chain, expressed as an anchor table
        ``{chain_id: int | {position: res_id}}`` (only dict is accepted;
        plain int / list forms were removed). A value ``int`` is
        shorthand for ``{0: int}`` (the number of position 0). For an
        anchor table ``{position: res_id}``, position ``p`` gets the
        ``res_id`` of the nearest anchor at or before ``p`` plus the
        offset ``p - anchor_position``; residues before the first
        anchor keep the default numbering starting at 1 (e.g.
        ``{5: 10}`` numbers positions 0-4 as 1-5, position 5 as 10,
        onwards 11, 12, ...). Multiple anchors per chain
        produce gapped numbering, e.g. ``{0: {0: 1, 20: 31}}``
        numbers positions 0-19 as 1-20, position 20 as 31, onwards
        32, 33, ... (default: every chain starts at 1).
    res_id_range : tuple of int or dict, optional
        Inclusive ``(start, end)`` range of residue numbers to draw
        per chain; positions whose ``renumber_res`` number falls
        outside the range are not drawn. The ``(start, end)`` pair
        form is only allowed for a single chain; with multiple chains
        a dict mapping chain ids to pairs is required (default: draw
        all).
    first_tick_id : int or dict, optional
        Anchor of the absolute x-axis tick grid: a residue numbered
        ``R`` is drawn at tick ``R - first_tick_id + 1``, so the
        leftmost tick shows ``first_tick_id`` and lower-numbered ticks
        stay empty (e.g. ``first_tick_id=100`` with first residue 107
        leaves 7 empty positions). Rows then cover a fixed window of
        ``per_line`` ticks each, aligned to ``first_tick_id``, so
        chains and segments line up by residue number even when the
        numbering differs; a row may hold fewer than ``per_line``
        letters. A plain int is only allowed for a single chain; with
        multiple chains a dict ``{chain_id: int}`` is required (each
        chain may anchor its grid at a different number). Residues
        numbered below a chain's ``first_tick_id`` raise ValueError.
        Default ``None`` keeps the packed layout (rows of up to
        ``per_line`` consecutive residues starting at tick 1).
    mark_res_ids : int, list of int or dict, optional
        Residue numbers (after ``renumber_res`` numbering) to mark
        with a red triangle above the stack. The plain int/list form
        is only allowed for a single chain; with multiple chains a
        dict ``{chain_id: int | list[int]}`` marks each chain with
        only its own ids (default: no marks).
    mark_conservative_res : bool, optional
        Mark positions whose column Shannon entropy is at most
        ``entropy_cutoff`` with a red underline below their tick label
        (residue number) (default: True).
    entropy_cutoff : float, optional
        Entropy threshold in bits; a position is underlined when its
        column entropy is ``<= entropy_cutoff`` (default 0.0 = only
        fully conserved columns). Entropy is computed from the column
        symbol frequencies over the full alphabet — the same quantity
        that determines stack height in ``'info'`` mode — so a
        negative value is invalid.
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
        ``entropy_cutoff`` is negative,
        ``renumber_res``/``res_id_range``/``first_tick_id``/
        ``mark_res_ids`` contain invalid values, a chain-id-less form
        (``res_id_range`` pair, ``first_tick_id`` int,
        ``mark_res_ids`` int/list) is used with multiple chains, the
        anchor numbering goes backwards, a chain is missing from any
        per-chain dict (unknown chain id), ``res_id_range`` has no
        intersection with a chain's numbering, or a chain's smallest
        residue number is below its ``first_tick_id``.

    Examples
    --------
    >>> profiles = FASTA_PROFILE(input_io="design.fa").read(chains=[0, 1])
    >>> fig, axes = plot_seqlogo(
    ...     profiles, renumber_res={0: {0: 1}, 1: {0: 95}},
    ...     mark_res_ids={0: [9, 31], 1: [95, 96]}
    ... )
    >>> fig, axes = plot_seqlogo(
    ...     profiles, renumber_res={0: {0: 8, 20: 31}},
    ...     res_id_range={0: (28, 106)}
    ... )
    >>> fig.savefig("design_logo.png", dpi=300)
    """
    if mode not in ("info", "freq"):
        raise ValueError(
            f"Invalid mode '{mode}'. Valid options are: 'info', 'freq'"
        )
    if per_line <= 0:
        raise ValueError(f"per_line must be positive, got {per_line}")
    if entropy_cutoff < 0:
        raise ValueError(
            f"entropy_cutoff must be non-negative, got {entropy_cutoff}"
        )

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
    chain_ids = [cid for cid, _ in chain_items]
    chain_lens = [len(profile) for _, profile in chain_items]
    if any(length == 0 for length in chain_lens):
        raise ValueError("profiles contain an empty chain")

    # 无链号形式 (单个值作用于全部链) 只允许单链: 多链必须显式带链号
    if n_chains > 1:
        if isinstance(res_id_range, (tuple, list)) and len(res_id_range) == 2:
            raise ValueError(
                "res_id_range 的 (start, end) 形式只支持单链; "
                "多链请用 {chain: (start, end)}"
            )
        if isinstance(first_tick_id, int):
            raise ValueError(
                "first_tick_id 单整数只支持单链; 多链请用 {chain: int}"
            )
        if isinstance(mark_res_ids, (int, list, tuple)):
            raise ValueError(
                "mark_res_ids 无链号形式只支持单链; 多链请用 "
                "{chain: int | [int, ...]}"
            )

    # 各链残基编号 (锚点表, 见 _resolve_renumber_res)
    anchors_by_chain = _resolve_renumber_res(renumber_res, chain_ids)
    # 各链 res_id 过滤范围 (含首尾), 不在范围内不画
    ranges_by_chain = _resolve_res_id_range(res_id_range, chain_ids)
    # 各链绝对网格起始编号 (见 _resolve_first_tick_id)
    ftids_by_chain = _resolve_first_tick_id(first_tick_id, chain_ids)

    # 逐链: 每个位置的编号 + 保留位置 (编号落在范围内的)
    res_ids_per_chain = []  # [(res_ids, kept_positions)]
    for cid, profile in chain_items:
        res_ids = _res_ids_for_chain(
            anchors_by_chain.get(cid, {}), len(profile)
        )
        r = ranges_by_chain.get(cid)
        if r is None:
            kept = list(range(len(profile)))
        else:
            start, end = r
            kept = [p for p in range(len(profile)) if start <= res_ids[p] <= end]
            if not kept:
                raise ValueError(
                    f"res_id_range {r} 与 chain {cid} 的编号 "
                    f"{res_ids[0]}..{res_ids[-1]} 无交集"
                )
        res_ids_per_chain.append((res_ids, kept))

    # 要标记的残基编号 (renumber_res 编号后, 逐链, 见 _resolve_mark_res_ids)
    mark_sets_by_chain = _resolve_mark_res_ids(mark_res_ids, chain_ids)

    # 分段: 列 = 链, 行 = 段。first_tick_id 未给出时按保留残基个数紧凑分段
    # (每行最多 per_line 个, 行号 = 段序); 给出时 x 轴为绝对 tick 网格:
    # 每行固定覆盖 per_line 个 tick、窗口起点对齐 first_tick_id, 行号 =
    # 窗口号相对全局最小窗口 (跨链/跨段按编号严格对齐, 空 tick 留空)。
    segments = []  # (chain_index, row_index, kept_start, kept_end, base_id)
    rows_per_chain = [set() for _ in range(n_chains)]
    if not ftids_by_chain:
        for ci, (res_ids, kept) in enumerate(res_ids_per_chain):
            for s in range(0, len(kept), per_line):
                e = min(s + per_line, len(kept))
                row = len(rows_per_chain[ci])
                segments.append((ci, row, s, e, res_ids[kept[s]]))
                rows_per_chain[ci].add(row)
    else:
        # 每链的 (窗口号 w, kept 起止) 列表; kept 内编号严格递增,
        # 同一窗口内的保留位置在 kept 中连续, 可直接切片
        chain_windows = []
        w_min = None
        for ci, (res_ids, kept) in enumerate(res_ids_per_chain):
            ftid = ftids_by_chain[chain_ids[ci]]
            # 校验: 编号低于该链 first_tick_id 的残基画不出 (落在轴左侧)
            if res_ids[kept[0]] < ftid:
                raise ValueError(
                    f"first_tick_id {ftid} 大于 chain {chain_ids[ci]} "
                    f"的最小编号 {res_ids[kept[0]]} (tick 从 first_tick_id "
                    f"起, 更小的编号无法画出)"
                )
            ws = []
            cur_w, s = None, 0
            for i, p in enumerate(kept):
                w = (res_ids[p] - ftid) // per_line
                if cur_w is None:
                    cur_w = w
                elif w != cur_w:
                    ws.append((cur_w, s, i))
                    cur_w, s = w, i
            ws.append((cur_w, s, len(kept)))
            chain_windows.append(ws)
            for w, _, _ in ws:
                if w_min is None or w < w_min:
                    w_min = w
        for ci, ws in enumerate(chain_windows):
            ftid = ftids_by_chain[chain_ids[ci]]
            for w, s, e in ws:
                row = w - w_min
                base = ftid + w * per_line  # 本窗口第一个 tick 的编号
                segments.append((ci, row, s, e, base))
                rows_per_chain[ci].add(row)
    n_rows = max(max(rows) for rows in rows_per_chain) + 1
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
        res_ids, kept = res_ids_per_chain[ci]
        entropies = _column_entropies(profile)
        rows_ci = rows_per_chain[ci]
        col_segs = sorted(
            (r, s, e, base) for (cc, r, s, e, base) in segments if cc == ci
        )
        for r, s, e, base in col_segs:
            ax = axes[r, ci]
            p0, p1 = kept[s], kept[e - 1]  # 段首/段尾在原始 profile 中的位置
            n = e - s  # 本段残基数
            # 字母 (i) 左下角在 i+0.5、宽 1 -> 中心 = i+1。局部下标 = res_id - base:
            # 紧凑模式 base = 段首 res_id (连续排列); 网格模式 base = 窗口第一个
            # tick 的编号 (空 tick 留空, 跨链按编号对齐)
            local = [res_ids[kept[s + i]] - base for i in range(n)]
            sub = profile[p0 : p1 + 1]  # SequenceProfile 切片,与整链字母宽度一致
            # 网格模式下前导空列必须补进切片再画, 否则字母仍从列 0 画起,
            # 与 tick 标签错位 (标签在网格位置, 字母在左侧)
            sub = _pad_profile_left(sub, local[0])
            if mode == "freq":
                # 全空列 (网格模式前导空位) 频率为 NaN, 不画字母, 屏蔽警告
                with np.errstate(divide="ignore", invalid="ignore"):
                    _plot_sequence_logo_freq(ax, sub)
            else:
                with np.errstate(divide="ignore", invalid="ignore"):
                    plot_sequence_logo(ax, sub)  # 默认配色 flower
            # 不足 per_line 的窗口/末段也按满宽画,字母宽度与其它段一致
            ax.set_xlim(0.5, per_line + 0.5)
            ax.set_title(
                f"Chain {chain_id} ({res_ids[p0]}-{res_ids[p1]})",
                fontsize=9,
            )
            ax.set_xticks([i + 1.0 for i in local])
            ax.set_xticklabels(
                [str(res_ids[kept[s + i]]) for i in range(n)], fontsize=7
            )
            if mode == "freq":
                ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
            else:
                max_entropy = np.log2(len(profile.alphabet))
                ax.set_yticks(range(0, int(np.ceil(max_entropy)) + 1))
            ax.tick_params(axis="y", labelsize=7)
            ax.spines[["top", "right"]].set_visible(False)
            # 标记指定残基 (renumber_res 编号后, 逐链): 字母上方红色倒三角
            mark_set = mark_sets_by_chain.get(chain_id, set())
            if mark_set:
                if mode == "freq":
                    top_y = 1.0 + 0.08
                    tri_y = top_y - 0.05
                else:
                    top_y = np.log2(len(profile.alphabet)) + 0.4
                    tri_y = top_y - 0.25
                ax.set_ylim(0, top_y)
                for i in range(n):
                    if res_ids[kept[s + i]] in mark_set:
                        ax.plot(
                            local[i] + 1.0,
                            tri_y,
                            marker="v",
                            color="#d62728",
                            markersize=7,
                            zorder=5,
                            clip_on=False,
                        )
            # 保守位点 (熵 <= entropy_cutoff, 默认 0 = 完全保守):
            # ticklabel (残基编号) 下方的红色下划线, 横跨整列字母宽度。
            # 用 x 轴混合变换 (x 为数据坐标, y 为轴盒子相对位置) 叠加固定
            # 英寸偏移: 线画在坐标轴盒子下方 ~0.24 in, 正好压在编号文字下方
            # (编号底边距轴盒子约 14 pt), 不触碰字母区也不扩展 ylim。
            if mark_conservative_res:
                cons_ids = [
                    i
                    for i in range(n)
                    if entropies[kept[s + i]] <= entropy_cutoff
                ]
                if cons_ids:
                    under_trans = ax.get_xaxis_transform() + ScaledTranslation(
                        0, -0.24, fig.dpi_scale_trans
                    )
                    for i in cons_ids:
                        ax.plot(
                            [local[i] + 0.5, local[i] + 1.5],
                            [0.0, 0.0],
                            color="#d62728",
                            lw=1.2,
                            zorder=5,
                            clip_on=False,
                            transform=under_trans,
                        )
            if not ftids_by_chain:
                top_row, bot_row = 0, n_rows - 1
            else:
                top_row, bot_row = min(rows_ci), max(rows_ci)
            if r == top_row:
                ax.set_ylabel(
                    "Frequency" if mode == "freq" else "Information (bits)"
                )
            if r == bot_row:
                ax.set_xlabel("Position")

    # 关闭无数据的格子 (网格模式下行号是绝对窗口号, 链间可能有空行)
    for ci in range(n_chains):
        for r in range(n_rows):
            if r not in rows_per_chain[ci]:
                axes[r, ci].set_axis_off()

    if title:
        fig.suptitle(title, fontsize=13)

    return fig, axes


# --- plot-seqlogo 子命令 ---

def _parse_int_list(spec, what):
    """'1,2,3' -> list[int]; None -> None; 非法输入直接报错退出"""
    if spec is None:
        return None
    try:
        return [int(t) for t in spec.split(",")]
    except ValueError:
        sys.exit(f"error: {what} 需为逗号分隔的整数 (如 1,2,3), got {spec!r}")


def _parse_renumber_res(spec, what, n_chains):
    """解析 --renumber-res (纯整数 N 形式已废除):
    - P_ID 两段: 单链锚点, 如 "0_8" (位置 0 编号 8)、"20_31";
      只支持单链, 多链报错
    - C_P_ID 三段: 逐链锚点, 如 "0_0_1,0_20_31,1_0_1"
      (同一链可给多个锚点以产生 gap 编号)
    两种形式不能混用, 均返回 {chain: {position: res_id}}。
    """
    if spec is None:
        return None
    tokens = [t.strip() for t in spec.split(",")]
    n_parts = {len(t.split("_")) for t in tokens}
    if len(n_parts) != 1:
        sys.exit(
            f"error: {what} 的 P_ID (两段) 与 C_P_ID (三段) 形式不能混用,"
            f" got {spec!r}"
        )
    n_parts = n_parts.pop()
    if n_parts == 2:
        if n_chains > 1:
            sys.exit(
                f"error: {what} 的 P_ID 形式 (如 0_8) 只支持单链; "
                f"多链请用 C_P_ID (如 0_0_8,1_0_34)"
            )
        anchors = {}
        for t in tokens:
            try:
                p, rid = (int(x) for x in t.split("_"))
            except ValueError:
                sys.exit(f"error: {what} 非法 P_ID {t!r}")
            if p < 0:
                sys.exit(f"error: {what} 锚点位置须 >= 0, got {t!r}")
            anchors[p] = rid
        return {0: anchors}
    if n_parts == 3:
        anchors = {}
        for t in tokens:
            try:
                c, p, rid = (int(x) for x in t.split("_"))
            except ValueError:
                sys.exit(f"error: {what} 非法 C_P_ID {t!r}")
            if p < 0:
                sys.exit(f"error: {what} 锚点位置须 >= 0, got {t!r}")
            anchors.setdefault(c, {})[p] = rid
        return anchors
    sys.exit(
        f"error: {what} 需为 P_ID (如 0_8) 或 C_P_ID (如 0_0_8); "
        f"纯整数形式已废除, got {spec!r}"
    )


def _parse_res_id_range(spec, what, n_chains):
    """解析 --res-id-range:
    - start_end 两段: 单链范围 (如 "28_106", 含首尾), 多链报错
    - chain_start_end 三段: 逐链范围 (如 "0_28_106,1_34_45")
    两种形式不能混用。
    """
    if spec is None:
        return None
    tokens = [t.strip() for t in spec.split(",")]
    global_range = None
    per_chain = {}
    for t in tokens:
        parts = t.split("_")
        try:
            vals = [int(x) for x in parts]
        except ValueError:
            sys.exit(
                f"error: {what} 每段需为整数 (start_end 或 chain_start_end),"
                f" got {t!r}"
            )
        if len(vals) == 2:
            if global_range is not None:
                sys.exit(f"error: {what} 全局范围只能出现一次")
            global_range = (vals[0], vals[1])
        elif len(vals) == 3:
            per_chain[vals[0]] = (vals[1], vals[2])
        else:
            sys.exit(
                f"error: {what} 每段需 2 或 3 个整数"
                f" (start_end 或 chain_start_end), got {t!r}"
            )
    if global_range is not None and per_chain:
        sys.exit(f"error: {what} 不能同时指定全局范围与逐链范围")
    if global_range is not None:
        if n_chains > 1:
            sys.exit(
                f"error: {what} 的 start_end 形式 (如 28_106) 只支持单链; "
                f"多链请用 chain_start_end (如 0_28_106,1_34_45)"
            )
        return global_range
    return per_chain


def _parse_first_tick_id(spec, what, n_chains):
    """解析 --first-tick-id:
    - N 单整数: 单链的网格起始编号 (如 "101"), 多链报错
    - C_N 两段: 逐链起始编号 (如 "0_101,1_105")
    两种形式不能混用。拒绝 '_' 数字分隔符 (防止手误)。
    """
    if spec is None:
        return None
    tokens = [t.strip() for t in spec.split(",")]
    plain = []
    for t in tokens:
        if "_" in t:
            plain = None
            break
        try:
            plain.append(int(t))
        except ValueError:
            plain = None
            break
    if plain is not None:
        if len(plain) != 1:
            sys.exit(
                f"error: {what} 纯整数形式只接受单个值 (单链), got {spec!r}"
            )
        if n_chains > 1:
            sys.exit(
                f"error: {what} 单整数只支持单链; 多链请用 C_N "
                f"(如 0_101,1_105)"
            )
        return plain[0]
    ftids = {}
    for t in tokens:
        parts = t.split("_")
        if len(parts) != 2:
            sys.exit(
                f"error: {what} 需为单整数 (N) 或 chain_N 两段 "
                f"(如 0_101,1_105), got {t!r}"
            )
        try:
            c, n = (int(x) for x in parts)
        except ValueError:
            sys.exit(f"error: {what} 非法 C_N {t!r}")
        ftids[c] = n
    return ftids


def _parse_mark_res_ids(spec, what, n_chains):
    """解析 --mark-res-ids:
    - N[,N...] 纯整数: 单链的标记编号列表 (如 "109,131,137"), 多链报错
    - C_N 两段: 逐链标记 (如 "0_109,0_131,1_164")
    两种形式不能混用。
    """
    if spec is None:
        return None
    tokens = [t.strip() for t in spec.split(",")]
    plain = []
    for t in tokens:
        if "_" in t:
            plain = None
            break
        try:
            plain.append(int(t))
        except ValueError:
            plain = None
            break
    if plain is not None:
        if n_chains > 1:
            sys.exit(
                f"error: {what} 纯整数形式只支持单链; 多链请用 C_N "
                f"(如 0_109,1_164)"
            )
        return plain
    marks = {}
    for t in tokens:
        parts = t.split("_")
        if len(parts) != 2:
            sys.exit(
                f"error: {what} 需为逗号分隔的整数 (单链) 或 chain_N "
                f"两段 (如 0_109,1_164), got {t!r}"
            )
        try:
            c, n = (int(x) for x in parts)
        except ValueError:
            sys.exit(f"error: {what} 非法 C_N {t!r}")
        marks.setdefault(c, []).append(n)
    return marks


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
    p.add_argument("--renumber-res", metavar="P_ID|C_P_ID[,...]", default=None,
                   help="残基重编号 (纯整数 N 形式已废除): P_ID 两段"
                        " (如 5_10 = 位置 5 编号 10) 只支持单链;"
                        " C_P_ID 三段 (如 0_0_1,0_20_31,1_0_1) 逐链指定,"
                        " 同一链多个锚点可产生 gap 编号; 第一个锚点之前"
                        " 的残基不重编号, 保持 1 起连续编号"
                        " (默认不指定 = 从 1 起)")
    p.add_argument("--res-id-range", metavar="S_E|C_S_E[,...]", default=None,
                   help="只画编号在 [起,止] (含首尾) 内的残基:"
                        " 两段 start_end (如 28_106) 只支持单链,"
                        " 三段 chain_start_end (如 0_28_106,1_34_45)"
                        " 逐链, 两种形式不能混用")
    p.add_argument("--first-tick-id", metavar="N|C_N[,...]", default=None,
                   help="x 轴绝对网格的起始编号: 编号 R 的残基画在 tick"
                        " (R-N+1) 处, 更小的 tick 留空 (如 N=100 且首残基"
                        " 107 时开头空 7 格); 每行固定 per_line 个 tick、"
                        " 窗口起点对齐 N, 跨链按编号对齐。单整数只支持"
                        " 单链, 多链用 C_N 逐链 (如 0_101,1_105)"
                        " (默认不指定 = 紧凑排列, 从 1 起)")
    p.add_argument("--mark-res-ids", metavar="N[,N...]|C_N[,...]", default=None,
                   help="标记这些残基 (renumber_res 编号后): 纯整数列表"
                        " 只支持单链; 多链用 C_N 逐链 (如 0_109,0_131,1_164)")
    p.add_argument("--mark-conservative-res",
                   action=argparse.BooleanOptionalAction, default=True,
                   help="在残基编号 (ticklabel) 下方用红色下划线标注保守位点"
                        " (列熵 <= --entropy-cutoff): 默认开启,"
                        " 传 --no-mark-conservative-res 关闭")
    p.add_argument("--entropy-cutoff", type=float, default=0.0,
                   metavar="FLOAT",
                   help="熵阈值 (bits): 列熵 <= 该值的位点被下划线标注"
                        " (默认 0 = 仅完全保守位点)")
    p.set_defaults(func=_run_seqlogo)
    return p


def _run_seqlogo(args) -> None:
    """执行 plot-seqlogo 子命令"""
    ignore = _parse_int_list(args.ignore_seqs, "--ignore-seqs")
    chains = _parse_int_list(args.chains, "--chains")
    try:
        profiles = Fasta_Profile(input_io=str(args.input)).read(
            ignore_seqs=ignore, chains=chains
        )
    except (OSError, ValueError) as e:
        sys.exit(f"error: 读取 {args.input} 失败: {e}")
    # 单链 FASTA (无 ':') read() 直接返回 SequenceProfile, 多链返回 dict
    n_chains = 1 if isinstance(profiles, SequenceProfile) else len(profiles)
    if n_chains == 0:
        sys.exit("error: 没有可用的序列/链 (检查 --ignore-seqs / --chains)")

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
            renumber_res=_parse_renumber_res(
                args.renumber_res, "--renumber-res", n_chains
            ),
            res_id_range=_parse_res_id_range(
                args.res_id_range, "--res-id-range", n_chains
            ),
            first_tick_id=_parse_first_tick_id(
                args.first_tick_id, "--first-tick-id", n_chains
            ),
            mark_res_ids=_parse_mark_res_ids(
                args.mark_res_ids, "--mark-res-ids", n_chains
            ),
            mark_conservative_res=args.mark_conservative_res,
            entropy_cutoff=args.entropy_cutoff,
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
