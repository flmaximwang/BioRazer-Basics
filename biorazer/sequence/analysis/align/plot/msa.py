import argparse
import sys

import matplotlib.pyplot as plt

from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment, align_multiple, SubstitutionMatrix
from biotite.sequence.graphics import (
    plot_alignment_type_based,
    plot_alignment_similarity_based,
    plot_alignment,
)

from biorazer.sequence.io import Fasta_StrDict


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
    if sequences is None:
        formatted_sequences = list(msa.sequences)
    else:
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
    if labels is None:
        labels = [f"{i+1}" for i in range(len(formatted_sequences))]
    plot_params["labels"] = labels
    plot_params["symbols_per_line"] = symbols_per_line
    plot_params["symbol_spacing"] = symbol_spacing
    plot_params["show_numbers"] = True
    plot_params["number_functions"] = number_functions

    # Get MSA
    if msa is None:
        msa, order, tree, distance_matrix = align_multiple(
            formatted_sequences,
            SubstitutionMatrix.std_protein_matrix(),
            **(msa_params or {}),
        )

    plot_func = method_map[plot_method]
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    plot_func(ax, msa, **plot_params)

    return fig, ax


# --- plot-msa 子命令 ---

def _parse_figsize(spec):
    """'WxH' -> (W, H) inches; None -> None"""
    if spec is None:
        return None
    try:
        w, h = spec.lower().split("x")
        return (float(w), float(h))
    except ValueError:
        sys.exit(f"error: --figsize 需为 WxH 格式 (如 10x5), got {spec!r}")


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
