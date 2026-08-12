"""biorazer CLI — 顶层多子命令入口"""
import argparse
import sys

from biorazer.access import cli as access_cli
from biorazer.sequence.analysis.align.query import colabfold_api
from biorazer.sequence.analysis.align.plot import cli as plot_cli


def main():
    parser = argparse.ArgumentParser(
        prog="biorazer",
        description="BioRazer 生物信息分析平台",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", help="子命令")

    # 注册子命令
    access_cli.register_subcommand(sub)
    colabfold_api.register_subcommand(sub)
    plot_cli.register_subcommand(sub)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
