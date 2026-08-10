"""Tests for biorazer.sequence.io.protein A3M converters.

A3m_Alignment / Alignment_A3m were merged from
biorazer.sequence.analysis.align.io (now removed).
"""

import importlib
import io

import pytest

from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment

from biorazer.sequence.io import A3m_Alignment, Alignment_A3m
from biorazer.sequence.io.protein import (
    A3m_Alignment as ModuleA3m,
    Alignment_A3m as ModuleAlignmentA3m,
)


def _make_alignment():
    """两序列 + 一个 gap 的 biotite Alignment。"""
    seqs = [ProteinSequence("ACDEFGHIK"), ProteinSequence("ACDEGHIK")]
    traces = Alignment.trace_from_strings(["ACDEFGHIK", "ACDE-GHIK"])
    return Alignment(seqs, traces)


class TestA3mImportLocation:
    def test_exported_from_sequence_io_package(self):
        assert A3m_Alignment is ModuleA3m
        assert Alignment_A3m is ModuleAlignmentA3m

    def test_old_align_io_module_removed(self):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("biorazer.sequence.analysis.align.io")


class TestA3mRead:
    def test_read_basic(self):
        text = ">seq1\nACDEFGHIK\n>seq2\nACDE-GHIK\n"
        aln = A3m_Alignment(input_io=io.StringIO(text)).read()
        assert len(aln.sequences) == 2
        assert aln.trace.shape == (9, 2)
        # '-' 已从序列中清除
        assert str(aln.sequences[0]) == "ACDEFGHIK"
        assert str(aln.sequences[1]) == "ACDEGHIK"

    def test_read_strips_lowercase_insertions(self):
        """小写字母为插入位,剔除后各行长度须一致;gap 结构保留在 trace。"""
        text = ">x\nACDEFGHIK\n>y\nACDEg-GHIK\n"
        aln = A3m_Alignment(input_io=io.StringIO(text)).read()
        assert str(aln.sequences[1]) == "ACDEGHIK"  # 小写 'g' 被剔除
        assert aln.trace[4, 1] == -1  # gap 列保留
        assert aln.trace.shape == (9, 2)

    def test_read_labels(self):
        text = ">seq1\nACDEFGHIK\n>seq2 with description\nACDE-GHIK\n"
        labels = A3m_Alignment(input_io=io.StringIO(text)).read_labels()
        assert labels == ["seq1", "seq2 with description"]

    def test_read_empty_raises(self):
        """空输入/单序列: biorazer 底层依赖 biotite,要求至少两条序列。"""
        with pytest.raises(ValueError, match="at least two sequences"):
            A3m_Alignment(input_io=io.StringIO("")).read()
        with pytest.raises(ValueError, match="at least two sequences"):
            A3m_Alignment(input_io=io.StringIO(">only\nACDEFGHIK\n")).read()


class TestAlignmentA3mWrite:
    def test_write_with_labels(self):
        out = io.StringIO()
        Alignment_A3m(output_io=out).write(_make_alignment(), labels=["a", "b"])
        assert out.getvalue() == ">a\nACDEFGHIK\n>b\nACDE-GHIK\n"

    def test_write_without_labels_uses_index(self):
        out = io.StringIO()
        Alignment_A3m(output_io=out).write(_make_alignment())
        assert out.getvalue() == ">0\nACDEFGHIK\n>1\nACDE-GHIK\n"

    def test_roundtrip(self):
        """Alignment_A3m 写出后可由 A3m_Alignment 读回等价对齐。"""
        out = io.StringIO()
        Alignment_A3m(output_io=out).write(_make_alignment(), labels=["a", "b"])
        out.seek(0)
        aln = A3m_Alignment(input_io=out).read()
        assert len(aln.sequences) == 2
        assert aln.trace.shape == (9, 2)
        assert str(aln.sequences[0]) == "ACDEFGHIK"
        assert str(aln.sequences[1]) == "ACDEGHIK"
