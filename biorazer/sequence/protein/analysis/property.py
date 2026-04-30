from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt
from biotite.sequence import ProteinSequence


STANDARD_AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def plot_amino_acid_distribution(
    sequence: str | ProteinSequence,
    normalize: bool = True,
    ax: plt.Axes | None = None,
    figsize: tuple[float, float] = (10, 4),
    color: str = "tab:blue",
) -> plt.Axes:
    """
    Plot the distribution of a protein sequence across the 20 standard amino acids.

    Parameters
    ----------
    sequence : str or ProteinSequence
            Protein sequence to analyze.
    normalize : bool, default=True
            If True, show residue fractions. Otherwise show raw counts.
    ax : matplotlib.axes.Axes, optional
            Existing axes to draw on. If None, a new figure is created.
    figsize : tuple of float, default=(10, 4)
            Figure size used only when ``ax`` is None.
    color : str, default="tab:blue"
            Bar color.

    Returns
    -------
    matplotlib.axes.Axes
            Axes containing the amino-acid distribution plot.
    """
    sequence_str = str(sequence).strip().upper()
    if not sequence_str:
        raise ValueError("sequence must not be empty")

    invalid_residues = sorted(set(sequence_str) - set(STANDARD_AMINO_ACIDS))
    if invalid_residues:
        raise ValueError(
            "sequence contains non-standard amino acids: " + ", ".join(invalid_residues)
        )

    residue_counts = Counter(sequence_str)
    values = [residue_counts.get(amino_acid, 0) for amino_acid in STANDARD_AMINO_ACIDS]
    if normalize:
        sequence_length = len(sequence_str)
        values = [count / sequence_length for count in values]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.bar(STANDARD_AMINO_ACIDS, values, color=color)
    ax.set_xlabel("Amino acid")
    ax.set_ylabel("Fraction" if normalize else "Count")
    ax.set_title("Amino-acid distribution")
    ymax = max(values) if values else 0.0
    ax.set_ylim(0.0, ymax * 1.1 if ymax > 0 else 1.0)

    return ax


__all__ = ["plot_amino_acid_distribution", "STANDARD_AMINO_ACIDS"]
