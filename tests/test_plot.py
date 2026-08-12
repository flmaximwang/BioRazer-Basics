"""Tests for biorazer.sequence.analysis.align.plot 模块拆分。

plot 包已拆分为 msa / msa_coverage / seqlogo 三个子模块:
每个模块自包含绘图 API 与对应 CLI (parser/runner);
__init__.py 仅做绘图 API 再导出,cli.py 为 CLI 门面 (register_subcommand + API 再导出)。
本测试锁定公共 API 的导入路径、身份与 CLI 端到端出图。
"""

import argparse

import matplotlib

matplotlib.use("Agg")  # 无显示环境, 必须先于 pyplot 导入

import pytest

from biorazer.sequence.analysis.align.plot import (
    plot_msa,
    plot_msa_coverage,
    plot_seqlogo,
)
from biorazer.sequence.analysis.align.plot.msa import plot_msa as _sub_plot_msa
from biorazer.sequence.analysis.align.plot.msa_coverage import (
    plot_msa_coverage as _sub_plot_msa_coverage,
)
from biorazer.sequence.analysis.align.plot.seqlogo import (
    plot_seqlogo as _sub_plot_seqlogo,
)
from biorazer.sequence.analysis.align.plot.seqlogo import _plot_sequence_logo_freq


def _build_plot_parser():
    """构建注册了三个 plot 子命令的 argparse parser,返回 (parser, subparsers, plot_cli)。"""
    from biorazer.sequence.analysis.align.plot import cli as plot_cli

    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    plot_cli.register_subcommand(sub)
    return parser, sub, plot_cli


class TestModuleSplit:
    """拆分后各函数在子模块与包路径下身份一致。"""

    def test_all_exports(self):
        """__all__ 列出且仅列出三个公共函数。"""
        import biorazer.sequence.analysis.align.plot as plot_pkg

        assert plot_pkg.__all__ == ["plot_msa_coverage", "plot_msa", "plot_seqlogo"]
        assert set(plot_pkg.__all__) <= set(dir(plot_pkg))

    def test_plot_msa_identity(self):
        """plot_msa 定义于 msa 子模块,经包路径再导出。"""
        assert plot_msa is _sub_plot_msa

    def test_plot_msa_coverage_identity(self):
        """plot_msa_coverage 定义于 msa_coverage 子模块,经包路径再导出。"""
        assert plot_msa_coverage is _sub_plot_msa_coverage

    def test_plot_seqlogo_identity(self):
        """plot_seqlogo 及其私有辅助定义于 seqlogo 子模块。"""
        assert plot_seqlogo is _sub_plot_seqlogo
        assert callable(_plot_sequence_logo_freq)

    def test_star_import_via_align(self):
        """align 包 `from .plot import *` 仍能拿到三个函数。"""
        from biorazer.sequence.analysis.align import plot_msa as _align_plot_msa

        assert _align_plot_msa is plot_msa


class TestCliSplit:
    """cli 门面: 仍向上层暴露 register_subcommand 与绘图 API,实现已下沉到各模块。"""

    def test_cli_reexports_api(self):
        """上层 (biorazer.cli) 通过 plot_cli 拿到的 API 保持可用。"""
        _, _, plot_cli = _build_plot_parser()
        assert callable(plot_cli.register_subcommand)
        assert plot_cli.plot_msa is plot_msa
        assert plot_cli.plot_msa_coverage is plot_msa_coverage
        assert plot_cli.plot_seqlogo is plot_seqlogo

    def test_register_subcommand_registers_three(self):
        """register_subcommand 注册三个子命令。"""
        _, sub, _ = _build_plot_parser()
        assert {"plot-seqlogo", "plot-msa", "plot-msa-coverage"} <= set(sub.choices)

    def test_subcommand_impls_live_in_plot_modules(self):
        """三个子命令的 parser/runner 定义在对应绘图模块。"""
        from biorazer.sequence.analysis.align.plot import msa, msa_coverage, seqlogo

        assert callable(seqlogo._add_seqlogo_parser) and callable(seqlogo._run_seqlogo)
        assert callable(msa._add_msa_parser) and callable(msa._run_msa)
        assert callable(msa_coverage._add_msa_coverage_parser)
        assert callable(msa_coverage._run_msa_coverage)

    def test_parser_wired_to_module_runner(self):
        """每个子命令的 func 指向对应模块的 runner。"""
        from biorazer.sequence.analysis.align.plot import msa, msa_coverage, seqlogo

        _, sub, _ = _build_plot_parser()
        assert sub.choices["plot-seqlogo"].get_default("func") is seqlogo._run_seqlogo
        assert sub.choices["plot-msa"].get_default("func") is msa._run_msa
        assert (
            sub.choices["plot-msa-coverage"].get_default("func")
            is msa_coverage._run_msa_coverage
        )


class TestCliEndToEnd:
    """三个子命令经 argparse 全流程可出图 (Agg 后端, tmp 目录)。"""

    def _run(self, argv, tmp_path, content, name="in.fa"):
        infa = tmp_path / name
        infa.write_text(content)
        outpng = tmp_path / "out.png"
        parser, _, _ = _build_plot_parser()
        args = parser.parse_args([*argv, "-i", str(infa), "-o", str(outpng)])
        args.func(args)
        assert outpng.exists() and outpng.stat().st_size > 0

    def test_plot_seqlogo(self, tmp_path):
        self._run(
            ["plot-seqlogo", "--mode", "freq"],
            tmp_path,
            ">q\nACDEFGHIKL\n>s1\nACDEFGHIKL\n>s2\nACDEEGHIKL\n",
        )

    def test_plot_seqlogo_numbering_flags(self, tmp_path):
        """新参数 (P_ID/C_P_ID 锚点 + 范围过滤 + 绝对网格) 经 CLI 全流程可出图。"""
        self._run(
            ["plot-seqlogo", "--mode", "freq",
             "--renumber-res", "0_0_1,0_20_31",
             "--res-id-range", "0_1_32",
             "--first-tick-id", "1"],
            tmp_path,
            ">q\nACDEFGHIKL\n>s1\nACDEFGHIKL\n>s2\nACDEEGHIKL\n",
        )

    def test_plot_seqlogo_multi_chain_flags(self, tmp_path):
        """多链必须用 C_ 前缀形式 (C_P_ID / C_S_E / C_N), 全流程可出图。"""
        self._run(
            ["plot-seqlogo", "--mode", "freq",
             "--renumber-res", "0_0_1,1_0_34",
             "--res-id-range", "0_1_3,1_34_36",
             "--first-tick-id", "0_1,1_34",
             "--mark-res-ids", "0_2,1_35"],
            tmp_path,
            ">s\nAAA:BBB\n",
        )

    def test_plot_msa(self, tmp_path):
        self._run(
            ["plot-msa"],
            tmp_path,
            ">a\nACDEFGHIK\n>b\nACDEFGHIK\n>c\nACDEGHIK\n",
        )

    def test_plot_msa_coverage(self, tmp_path):
        self._run(
            ["plot-msa-coverage", "--part-lengths", "10"],
            tmp_path,
            ">q\nACDEFGHIKL\n>s1\nACD-FGHIKL\n>s2\nACDEFGHI-L\n",
        )


class TestSeqlogoNumbering:
    """renumber_res / res_id_range / first_tick_id / mark_res_ids 的编号与过滤
    (渲染后检查坐标轴标签与字母位置)。无链号形式只允许单链, 多链必须 dict。"""

    @staticmethod
    def _profiles(n_chains=1, length=20, char="A"):
        """构造 n_chains 条等长全 A 的 profile dict (走 FASTA 路径, 链以 ':' 分隔)。"""
        import io

        from biorazer.sequence.io import Fasta_Profile

        fasta = f">s\n{':'.join([char * length] * n_chains)}\n"
        return Fasta_Profile(input_io=io.StringIO(fasta)).read()

    @staticmethod
    def _labels(axes, row=0, col=0):
        return [t.get_text() for t in axes[row, col].get_xticklabels()]

    def test_renumber_res_single_chain(self):
        """单链 dict {0: 8}: 从 8 起顺延。"""
        fig, axes = plot_seqlogo(self._profiles(length=5), renumber_res={0: 8})
        assert self._labels(axes) == ["8", "9", "10", "11", "12"]
        assert axes[0, 0].get_title() == "Chain 0 (8-12)"

    def test_renumber_res_dict_shorthand(self):
        """dict 值取 int = {0: 编号} 简写, 逐链。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=3),
            renumber_res={0: 8, 1: 34},
        )
        assert self._labels(axes, col=0) == ["8", "9", "10"]
        assert self._labels(axes, col=1) == ["34", "35", "36"]

    def test_renumber_res_anchor_gap(self):
        """多锚点产生 gap 编号: 位置 0-4 编号 1-5, 位置 5 起编号 31。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=8), renumber_res={0: {0: 1, 5: 31}}
        )
        assert self._labels(axes) == ["1", "2", "3", "4", "5", "31", "32", "33"]

    def test_renumber_res_anchor_before_first(self):
        """第一个锚点之前按它逆推: {20: 31} -> 位置 0 编号 11。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=8), renumber_res={0: {20: 31}}
        )
        assert self._labels(axes) == [str(v) for v in range(11, 19)]

    def test_renumber_res_plain_int_rejected(self):
        """纯整数 / 列表形式已废除 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(self._profiles(length=5), renumber_res=8)
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(n_chains=2, length=3), renumber_res=[8, 34]
            )

    def test_renumber_res_backward_anchor(self):
        """锚点编号回退 (不满足顺延) 报错。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=8), renumber_res={0: {0: 100, 5: 1}}
            )

    def test_renumber_res_unknown_chain(self):
        """renumber_res 引用不存在的链 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=8), renumber_res={3: {0: 1}}
            )

    def test_res_id_range_single_chain(self):
        """单链 tuple 范围 (含首尾), 标题显示实际编号。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=10),
            renumber_res={0: 8},
            res_id_range=(10, 12),
        )
        assert self._labels(axes) == ["10", "11", "12"]
        assert axes[0, 0].get_title() == "Chain 0 (10-12)"

    def test_res_id_range_per_chain(self):
        """dict 逐链范围, 未指定的链不过滤。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=10),
            renumber_res={0: 1, 1: 100},
            res_id_range={0: (3, 5)},
        )
        assert self._labels(axes, col=0) == ["3", "4", "5"]
        assert self._labels(axes, col=1) == [
            str(v) for v in range(100, 110)
        ]

    def test_res_id_range_multi_chain_tuple_rejected(self):
        """多链 + (start, end) 元组 -> ValueError (必须 dict)。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(n_chains=2, length=10),
                renumber_res={0: 1, 1: 100},
                res_id_range=(3, 5),
            )

    def test_res_id_range_no_intersection(self):
        """范围与某链编号无交集 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=10),
                renumber_res={0: 8},
                res_id_range=(20, 30),
            )

    def test_res_id_range_bad_bounds(self):
        """start > end -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(self._profiles(length=8), res_id_range=(5, 2))

    def test_first_tick_id_offset(self):
        """first_tick_id=100 且首残基 107 -> 开头空 7 格, 字母中心从 8.0 起。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=5), renumber_res={0: 107}, first_tick_id=100
        )
        assert list(axes[0, 0].get_xticks()) == [8.0, 9.0, 10.0, 11.0, 12.0]
        assert self._labels(axes) == ["107", "108", "109", "110", "111"]

    def test_first_tick_id_per_chain(self):
        """多链 dict {chain: int} 逐链网格起点。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=5),
            renumber_res={0: 107, 1: 34},
            first_tick_id={0: 100, 1: 30},
        )
        assert list(axes[0, 0].get_xticks()) == [8.0, 9.0, 10.0, 11.0, 12.0]
        assert list(axes[0, 1].get_xticks()) == [5.0, 6.0, 7.0, 8.0, 9.0]

    def test_first_tick_id_multi_chain_int_rejected(self):
        """多链 + 单整数 -> ValueError (必须 dict)。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(n_chains=2, length=5),
                renumber_res={0: 107, 1: 34},
                first_tick_id=100,
            )

    def test_first_tick_id_cross_chain_alignment(self):
        """不同编号起点的两条链: 同编号残基落在同一 x 位置。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=50),
            renumber_res={0: 8, 1: 34},
            first_tick_id={0: 8, 1: 8},
        )
        assert list(axes[0, 0].get_xticks())[0] == 1.0  # 链 0 首残基 8 在 tick 1
        assert list(axes[0, 1].get_xticks())[0] == 27.0  # 链 1 首残基 34 在 tick 27
        # 链 1 开头 26 格空出 (tick 1-26 对应编号 8-33)
        x_by_id = {
            col: dict(zip(self._labels(axes, col=col), axes[0, col].get_xticks()))
            for col in (0, 1)
        }
        assert x_by_id[0]["34"] == x_by_id[1]["34"] == 27.0

    def test_first_tick_id_windows(self):
        """长链按 per_line 个 tick 分窗, 每行窗口起点对齐 first_tick_id。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=150),
            renumber_res={0: 8},
            first_tick_id=8,
            per_line=50,
        )
        assert axes.shape == (3, 1)
        assert self._labels(axes, row=0) == [str(v) for v in range(8, 58)]
        assert self._labels(axes, row=1) == [str(v) for v in range(58, 108)]
        assert self._labels(axes, row=2) == [str(v) for v in range(108, 158)]

    def test_first_tick_id_window_rows_absolute(self):
        """行号 = 窗口号 (相对全局最小窗口): 链间错开的窗口不压缩行号。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=50),
            renumber_res={0: 8, 1: 158},
            first_tick_id={0: 8, 1: 8},
            per_line=50,
        )
        assert axes.shape == (4, 2)
        assert self._labels(axes, row=0, col=0)[0] == "8"
        assert self._labels(axes, row=3, col=1)[0] == "158"
        # 无数据的格子不画标题 (有数据的格子一定有 "Chain ..." 标题)
        for r, c in [(1, 0), (2, 0), (1, 1), (2, 1)]:
            assert axes[r, c].get_title() == ""

    def test_first_tick_id_with_range(self):
        """res_id_range 过滤后仍按绝对网格定位 (编号 10 在 tick 3)。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=10),
            renumber_res={0: 8},
            res_id_range=(10, 12),
            first_tick_id=8,
        )
        assert self._labels(axes) == ["10", "11", "12"]
        assert list(axes[0, 0].get_xticks()) == [3.0, 4.0, 5.0]

    def test_first_tick_id_below_min(self):
        """first_tick_id 大于链的最小编号 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=5), renumber_res={0: 8}, first_tick_id=100
            )

    def test_first_tick_id_bad_type(self):
        """first_tick_id 非 int/dict -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(self._profiles(length=5), first_tick_id="100")

    def test_first_tick_id_unknown_chain(self):
        """first_tick_id dict 引用不存在的链 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=5), first_tick_id={3: 100}
            )

    def test_first_tick_id_letter_positions(self):
        """字母本身也要落在绝对网格上 (用户场景: 锚点编号 + 范围 + 网格)。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=160),
            renumber_res={0: {0: 1, 85: 107}},
            res_id_range=(107, 179),
            first_tick_id=101,
        )
        assert axes.shape == (2, 1)
        # 窗口 0 (tick 1-50 = 编号 101-150): 编号 107 的字母局部下标 6,
        # 左下角 x = 6.5 -> 中心 7.0, 与标签 107 所在 tick 一致
        ax0 = axes[0, 0]
        # TextPath 有 ~1e-5 的固有 x 偏移, 用 1e-3 容差
        x0 = ax0.patches[0].get_path().get_extents().x0
        assert abs(x0 - 6.5) < 1e-3, f"row0 first letter x0 = {x0}"
        assert list(ax0.get_xticks())[0] == 7.0
        assert self._labels(axes, row=0)[0] == "107"
        # 窗口 1 (tick 51-100 = 编号 151-200): 编号 151 在局部下标 0, 无需 padding
        ax1 = axes[1, 0]
        x1 = ax1.patches[0].get_path().get_extents().x0
        assert abs(x1 - 0.5) < 1e-3, f"row1 first letter x0 = {x1}"
        assert list(ax1.get_xticks())[0] == 1.0
        assert self._labels(axes, row=1)[0] == "151"

    def test_mark_res_ids_single_chain(self):
        """单链纯整数列表: 标记这些编号的残基。"""
        fig, axes = plot_seqlogo(
            self._profiles(length=20),
            renumber_res={0: 8},
            mark_res_ids=[9, 12],
        )
        ax = axes[0, 0]
        mark_x = sorted(l.get_xdata()[0] for l in ax.lines)
        assert mark_x == [2.0, 5.0]  # 9-8+1=2, 12-8+1=5

    def test_mark_res_ids_per_chain(self):
        """多链 dict {chain: [...]} 逐链标记, 互不串链。"""
        fig, axes = plot_seqlogo(
            self._profiles(n_chains=2, length=10),
            renumber_res={0: 8, 1: 34},
            mark_res_ids={0: [9], 1: [35]},
        )
        assert sorted(l.get_xdata()[0] for l in axes[0, 0].lines) == [2.0]
        assert sorted(l.get_xdata()[0] for l in axes[0, 1].lines) == [2.0]
        # 链 0 的标记 9 不应出现在链 1
        assert len(axes[0, 0].lines) == 1 and len(axes[0, 1].lines) == 1

    def test_mark_res_ids_multi_chain_plain_rejected(self):
        """多链 + 纯整数列表 -> ValueError (必须 dict)。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(n_chains=2, length=10),
                renumber_res={0: 8, 1: 34},
                mark_res_ids=[9, 35],
            )

    def test_mark_res_ids_unknown_chain(self):
        """mark_res_ids dict 引用不存在的链 -> ValueError。"""
        with pytest.raises(ValueError):
            plot_seqlogo(
                self._profiles(length=5), mark_res_ids={3: [1]}
            )

    def test_pad_profile_left(self):
        """_pad_profile_left: 前补 n 列全空列, 原符号与 gap 计数保留。"""
        import numpy as np

        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _pad_profile_left,
        )

        prof = self._profiles(length=3)[0]
        padded = _pad_profile_left(prof, 2)
        assert padded.symbols.shape == (5, prof.symbols.shape[1])
        assert not padded.symbols[:2].any()
        assert np.array_equal(padded.symbols[2:], prof.symbols)
        assert padded.gaps.shape == (5,)
        assert not padded.gaps[:2].any()
        assert np.array_equal(padded.gaps[2:], prof.gaps)
        assert _pad_profile_left(prof, 0) is prof


class TestSeqlogoCliNumbering:
    """plot-seqlogo 的 --renumber-res / --res-id-range / --first-tick-id /
    --mark-res-ids 解析与接线 (无链号形式只允许单链, 多链必须 C_ 前缀)。"""

    @staticmethod
    def _parse(argv):
        parser, _, _ = _build_plot_parser()
        return parser.parse_args(argv)

    def test_renumber_res_p_id_single_chain(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_renumber_res,
        )

        assert _parse_renumber_res(
            "0_8,20_31", "--renumber-res", 1
        ) == {0: {0: 8, 20: 31}}

    def test_renumber_res_p_id_multi_chain_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_renumber_res,
        )

        with pytest.raises(SystemExit):
            _parse_renumber_res("0_8", "--renumber-res", 2)

    def test_renumber_res_c_p_id(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_renumber_res,
        )

        assert _parse_renumber_res(
            "0_0_1,0_20_31,1_0_1", "--renumber-res", 2
        ) == {0: {0: 1, 20: 31}, 1: {0: 1}}

    def test_renumber_res_plain_int_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_renumber_res,
        )

        # 纯整数 N 形式已废除 (单链也要写 P_ID / C_P_ID)
        with pytest.raises(SystemExit):
            _parse_renumber_res("8", "--renumber-res", 1)
        with pytest.raises(SystemExit):
            _parse_renumber_res("8,34", "--renumber-res", 2)

    def test_renumber_res_mixed_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_renumber_res,
        )

        with pytest.raises(SystemExit):
            _parse_renumber_res("0_8,0_0_1", "--renumber-res", 1)

    def test_res_id_range_single_chain(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_res_id_range,
        )

        assert _parse_res_id_range("28_106", "--res-id-range", 1) == (28, 106)

    def test_res_id_range_multi_chain_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_res_id_range,
        )

        with pytest.raises(SystemExit):
            _parse_res_id_range("28_106", "--res-id-range", 2)

    def test_res_id_range_per_chain(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_res_id_range,
        )

        assert _parse_res_id_range(
            "0_28_106,1_34_45", "--res-id-range", 2
        ) == {0: (28, 106), 1: (34, 45)}

    def test_res_id_range_mixed_global_per_chain(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_res_id_range,
        )

        with pytest.raises(SystemExit):
            _parse_res_id_range("28_106,1_34_45", "--res-id-range", 2)

    def test_first_tick_id_single(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_first_tick_id,
        )

        assert _parse_first_tick_id("100", "--first-tick-id", 1) == 100
        assert _parse_first_tick_id(None, "--first-tick-id", 1) is None

    def test_first_tick_id_cn(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_first_tick_id,
        )

        assert _parse_first_tick_id(
            "0_101,1_105", "--first-tick-id", 2
        ) == {0: 101, 1: 105}

    def test_first_tick_id_multi_chain_plain_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_first_tick_id,
        )

        with pytest.raises(SystemExit):
            _parse_first_tick_id("100", "--first-tick-id", 2)

    def test_first_tick_id_invalid(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_first_tick_id,
        )

        with pytest.raises(SystemExit):
            _parse_first_tick_id("abc", "--first-tick-id", 1)
        with pytest.raises(SystemExit):
            _parse_first_tick_id("101,102", "--first-tick-id", 1)
        with pytest.raises(SystemExit):
            _parse_first_tick_id("1_0_1", "--first-tick-id", 1)
        # "1_00" 语法上是合法 C_N (链 1 编号 0): 数字分隔符误写 (1_00=100)
        # 在单链时由 unknown chain 校验兜底报错, 不会静默生效
        assert _parse_first_tick_id("1_00", "--first-tick-id", 1) == {1: 0}

    def test_mark_res_ids_single_chain(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_mark_res_ids,
        )

        assert _parse_mark_res_ids(
            "109,131,137", "--mark-res-ids", 1
        ) == [109, 131, 137]

    def test_mark_res_ids_cn(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_mark_res_ids,
        )

        assert _parse_mark_res_ids(
            "0_109,0_131,1_164", "--mark-res-ids", 2
        ) == {0: [109, 131], 1: [164]}

    def test_mark_res_ids_multi_chain_plain_rejected(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_mark_res_ids,
        )

        with pytest.raises(SystemExit):
            _parse_mark_res_ids("109,131", "--mark-res-ids", 2)

    def test_mark_res_ids_invalid(self):
        from biorazer.sequence.analysis.align.plot.seqlogo import (
            _parse_mark_res_ids,
        )

        with pytest.raises(SystemExit):
            _parse_mark_res_ids("0_109_1", "--mark-res-ids", 2)

    def test_flag_names(self):
        """旧 flag --first-res-id / --res-id-shift 已移除, 新 flag 走字符串原样到 runner。"""
        args = self._parse(
            [
                "plot-seqlogo",
                "-i",
                "x.fa",
                "--renumber-res",
                "0_0_1,0_20_31",
                "--res-id-range",
                "0_1_32",
                "--first-tick-id",
                "100",
                "--mark-res-ids",
                "0_109,0_131",
            ]
        )
        assert args.renumber_res == "0_0_1,0_20_31"
        assert args.res_id_range == "0_1_32"
        assert args.first_tick_id == "100"
        assert args.mark_res_ids == "0_109,0_131"
        assert not hasattr(args, "first_res_id")
        assert not hasattr(args, "res_id_shift")


class TestEdgeCases:
    """拆分后保留原有的输入校验行为(不触发绘图)。"""

    def test_plot_msa_requires_input(self):
        """sequences 与 msa 必须二选一。"""
        with pytest.raises(ValueError):
            plot_msa(sequences=None, msa=None)
        with pytest.raises(ValueError):
            plot_msa(sequences=["ACDEFGHIK"], msa="anything")

    def test_plot_msa_invalid_method(self):
        """非法 plot_method 直接报错。"""
        with pytest.raises(ValueError):
            plot_msa(sequences=["ACDEFGHIK"], plot_method="bad")

    def test_plot_seqlogo_invalid_mode(self):
        """mode 仅接受 info / freq。"""
        with pytest.raises(ValueError):
            plot_seqlogo([], mode="bad")

    def test_plot_seqlogo_empty_profiles(self):
        """空 profiles 报错。"""
        with pytest.raises(ValueError):
            plot_seqlogo([])
