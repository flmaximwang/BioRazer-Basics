# BioRazer

A Python package for analyzing various biological information, built from practical lab experience. Install with pip and you're ready for common bioinformatics analysis work.

## Features

- **Protein sequences** — Translation, reverse translation, codon/protein dictionaries
- **Multiple Sequence Alignment (MSA)** — Generate MSA via ColabFold MMseqs2 API, visualize coverage, analyze amino acid frequencies
- **Structure analysis** — Static analysis (contacts, hydrogen bonds, surface selection), dynamic trajectory analysis (MD trajectory view, XVG/XPM plots)
- **Database access** — Query AFDB, RCSB PDB, UniProt, Ensembl
- **Protein design** — Sequence design, library generation, single test entries

## Installation

```bash
pip install biorazer
```

### Dependencies

- Python >= 3.11
- biotite, numpy, scipy, matplotlib, hydride, umap-learn, rcsb-api
- tabulate (for formatted output)

### Development

```bash
# Install with dev and test dependencies
poetry install --with dev,test

# Or using pip with test dependencies
pip install biorazer
pip install pytest pytest-cov
```

## Usage

### ColabFold MSA via MMseqs2 API

Generate protein MSA by calling the ColabFold public API — zero additional dependencies, pure stdlib:

```python
from biorazer.sequence.protein.analysis.align.query import run_search

# Single-chain MSA (unpaired, default)
files, _ = run_search(
    ["MTSENLYFQGAMG..."],
    out_dir="msa_out/",
)

# Multi-chain paired MSA (for AF3 multimers)
files, _ = run_search(
    ["CHAIN1_SEQUENCE", "CHAIN2_SEQUENCE"],
    out_dir="msa_out/",
    pair_mode="paired",
    pair_strategy="greedy",   # or "complete"
)
```

Output: A3M files (`uniref.a3m`, `bfd.mgnify30.*.a3m`, `pair.a3m`) ready for downstream folding pipelines. Supports template search (`--templates`) and custom MMseqs2 server URLs.

### MSA Visualization

```python
from biorazer.sequence.protein.analysis.align import plot_msa

fig, ax = plot_msa(
    sequences=["MTSENLYFQG", "MTSENLXFQG"],
    labels=["Wild-type", "Mutant"],
)
fig.savefig("msa_plot.png")
```

## Testing

```bash
pytest tests/ -v
```

Tests cover: FASTA parsing, sequence validation, A3M merging, and module constants. All tests are pure (no network required).

## Project Structure

```
biorazer/
├── access/         # External database APIs (AFDB, RCSB, UniProt, Ensembl)
├── design/         # Protein design tools
├── sequence/       # Sequence analysis
│   ├── nucleotide/
│   ├── protein/
│   │   ├── analysis/align/   # MSA generation & analysis
│   │   │   └── query/        # ColabFold MMseqs2 API
│   │   └── scripts/          # MSA visualizer, analyzer
│   └── translation/
├── structure/      # Structure analysis & I/O
└── util/           # Utility modules (dictionaries)
```

## License

This project is intended for academic and research use.
