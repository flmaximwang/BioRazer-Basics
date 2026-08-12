"""Tests for biorazer.access.server.colabfold_msa (重构自 sequence/.../query/colabfold_api.py)。"""

import argparse
import importlib
import os
import subprocess
import sys
import tempfile
import pytest

from biorazer.access.server.colabfold_msa.http import DEFAULT_HOST, DEFAULT_UA
from biorazer.access.server.colabfold_msa.io import merge_a3m, parse_fasta, validate
from biorazer.access.server.colabfold_msa.pipeline import SearchResult, SeqResult


class TestParseFasta:
    def test_single_sequence(self):
        text = ">seq1\nMTSENLYFQG\n"
        assert parse_fasta(text) == [("seq1", "MTSENLYFQG")]

    def test_multiple_sequences(self):
        text = ">seq1\nMTSENLYFQG\n>seq2\nACDEFGHIK\n>seq3\nWPKL\n"
        assert parse_fasta(text) == [
            ("seq1", "MTSENLYFQG"),
            ("seq2", "ACDEFGHIK"),
            ("seq3", "WPKL"),
        ]

    def test_sequence_with_description(self):
        text = ">seq1 some description\nMTSENLYFQG\n"
        assert parse_fasta(text) == [("seq1", "MTSENLYFQG")]

    def test_multi_line_fasta(self):
        text = ">seq1\nMTSEN\nLYFQG\n"
        assert parse_fasta(text) == [("seq1", "MTSENLYFQG")]

    def test_empty(self):
        assert parse_fasta("") == []
        assert parse_fasta("   \n\n") == []

    def test_no_header(self):
        assert parse_fasta("MTSENLYFQG") == [("default", "MTSENLYFQG")]

    def test_uppercase_conversion(self):
        text = ">seq1\nmtsenlyfqg\n"
        assert parse_fasta(text) == [("seq1", "MTSENLYFQG")]

    def test_empty_sequence_after_header(self):
        text = ">seq1\n>seq2\nMTS\n"
        assert parse_fasta(text) == [("seq2", "MTS")]

    def test_trailing_newline(self):
        text = ">seq1\nMTSENLYFQG\n\n"
        assert parse_fasta(text) == [("seq1", "MTSENLYFQG")]

    def test_header_with_spaces(self):
        """取 header 第一个空白前的部分"""
        text = ">my_protein something else\nMTSEN\n"
        assert parse_fasta(text) == [("my_protein", "MTSEN")]

    def test_empty_header(self):
        text = ">\nMTSEN\n"
        assert parse_fasta(text) == [("default", "MTSEN")]


class TestValidate:
    def test_valid_standard(self):
        validate([("seq1", "MTSENLYFQG"), ("seq2", "ACDEFGHIK")])

    def test_valid_single_aa(self):
        validate([("a", "A")])
        validate([("c", "C")])
        validate([("y", "Y")])

    def test_valid_with_all_standard_aas(self):
        all_aas = "ACDEFGHIKLMNPQRSTVWY"
        validate([("all", all_aas)])

    def test_lowercase_is_valid(self):
        """内部做了 .upper()，小写字母合法。"""
        validate([("low", "mtsenlyfqg")])

    def test_invalid_character_exits(self):
        with pytest.raises(SystemExit):
            validate([("bad", "MTSENLXFQG")])  # X 不是标准氨基酸

    def test_empty_string_is_valid(self):
        validate([("empty", "")])

    def test_mixed_valid_and_invalid(self):
        with pytest.raises(SystemExit):
            validate([("good", "MTSENLYFQG"), ("bad", "ACD*FGHIK")])  # * 非法

    def test_multiple_invalid_sequences(self):
        with pytest.raises(SystemExit):
            validate([("a", "MTS"), ("b", "ABC"), ("c", "BAD1SEQ")])

    def test_numeric_character_exits(self):
        with pytest.raises(SystemExit):
            validate([("num", "MTS3NLYFQG")])


class TestMergeA3M:
    def test_merge_single_file(self, sample_a3m_files):
        result = merge_a3m([sample_a3m_files[0]])
        assert ">seq1" in result
        assert "MTSENLYFQG" in result
        assert result.endswith("\n")

    def test_merge_multiple_files(self, sample_a3m_files):
        result = merge_a3m(sample_a3m_files)
        assert result.count(">seq") == 2
        assert "MTSENLYFQG" in result
        assert "WPKL" in result

    def test_merge_order_preserved(self, sample_a3m_files):
        result = merge_a3m(sample_a3m_files)
        pos1 = result.index("MTSENLYFQG")
        pos2 = result.index("WPKL")
        assert pos1 < pos2

    def test_merge_empty_file(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".a3m", delete=False) as f:
            f.write("")
            fname = f.name
        result = merge_a3m([fname])
        # 空文件 .rstrip() → "" → "" + "\n" → "\n"
        assert result == "\n"
        os.unlink(fname)


class TestConstants:
    def test_default_host(self):
        assert DEFAULT_HOST == "https://api.colabfold.com"

    def test_default_ua(self):
        assert "colabfold_msa/2.0" in DEFAULT_UA
        assert "sokrypton/ColabFold" in DEFAULT_UA


# run_search 需要 mock HTTP 请求，纯函数测试见上。


class TestImportLocation:
    """重构身份验证: 新路径的模块可导入, 旧路径 (alignment/query) 已删除。"""

    def test_old_align_query_module_removed(self):
        """旧路径 biorazer.sequence.analysis.alignment.query.colabfold_api 必须不复存在。"""
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("biorazer.sequence.analysis.alignment.query.colabfold_api")
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("biorazer.sequence.analysis.alignment.query")

    def test_new_modules_importable(self):
        """拆分后的各 module 均可独立导入。"""
        importlib.import_module("biorazer.access.server.colabfold_msa.api")
        importlib.import_module("biorazer.access.server.colabfold_msa.cli")
        importlib.import_module("biorazer.access.server.colabfold_msa.http")
        importlib.import_module("biorazer.access.server.colabfold_msa.io")
        importlib.import_module("biorazer.access.server.colabfold_msa.paired")
        importlib.import_module("biorazer.access.server.colabfold_msa.pipeline")
        importlib.import_module("biorazer.access.server.colabfold_msa.template")
        importlib.import_module("biorazer.access.server.colabfold_msa.unpaired")


class TestCliWiring:
    """CLI 接线: colabfold-msa 经新位置注册到 biorazer 顶层入口。"""

    def test_subcommand_registered(self):
        from biorazer.access.server.colabfold_msa import cli

        parser = argparse.ArgumentParser(prog="biorazer")
        sub = parser.add_subparsers(dest="command")
        cli.register_subcommand(sub)
        assert "colabfold-msa" in sub.choices
        assert sub.choices["colabfold-msa"].get_default("func") is cli._run_colabfold

    def test_top_level_cli_help_lists_colabfold_msa(self):
        """经 biorazer.cli 主入口 (最上层消费路径) 可见该子命令。"""
        result = subprocess.run(
            [sys.executable, "-m", "biorazer.cli", "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "colabfold-msa" in result.stdout
