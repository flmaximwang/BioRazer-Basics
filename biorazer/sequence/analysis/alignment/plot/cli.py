"""biorazer plot 子命令门面: plot-seqlogo / plot-msa / plot-msa-coverage。

三个子命令的 parser/runner 分别位于同目录的 seqlogo / msa / msa_coverage 模块,
本模块负责组合注册,并向上层模块 (biorazer.cli) 暴露 register_subcommand 与绘图 API。
"""
from .msa import _add_msa_parser, plot_msa
from .msa_coverage import _add_msa_coverage_parser, plot_msa_coverage
from .seqlogo import _add_seqlogo_parser, plot_seqlogo

__all__ = [
    "register_subcommand",
    "plot_msa",
    "plot_msa_coverage",
    "plot_seqlogo",
]


def register_subcommand(sub) -> None:
    """在 argparse subparsers 上注册 plot 相关子命令"""
    _add_seqlogo_parser(sub)
    _add_msa_parser(sub)
    _add_msa_coverage_parser(sub)
