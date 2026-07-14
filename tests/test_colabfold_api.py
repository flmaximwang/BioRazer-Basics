"""Tests for biorazer.sequence.protein.analysis.align.query.colabfold_api."""

import os
import tempfile
import pytest

from biorazer.sequence.protein.analysis.align.query.colabfold_api import (
    DEFAULT_HOST,
    DEFAULT_UA,
    parse_fasta,
    validate,
    merge_a3m,
)


class TestParseFasta:
    def test_single_sequence(self):
        text = ">seq1\nMTSENLYFQG\n"
        assert parse_fasta(text) == ["MTSENLYFQG"]

    def test_multiple_sequences(self):
        text = ">seq1\nMTSENLYFQG\n>seq2\nACDEFGHIK\n>seq3\nWPKL\n"
        assert parse_fasta(text) == ["MTSENLYFQG", "ACDEFGHIK", "WPKL"]

    def test_sequence_with_description(self):
        text = ">seq1 some description\nMTSENLYFQG\n"
        assert parse_fasta(text) == ["MTSENLYFQG"]

    def test_multi_line_fasta(self):
        text = ">seq1\nMTSEN\nLYFQG\n"
        assert parse_fasta(text) == ["MTSENLYFQG"]

    def test_empty(self):
        assert parse_fasta("") == []
        assert parse_fasta("   \n\n") == []

    def test_no_header(self):
        assert parse_fasta("MTSENLYFQG") == ["MTSENLYFQG"]

    def test_uppercase_conversion(self):
        text = ">seq1\nmtsenlyfqg\n"
        assert parse_fasta(text) == ["MTSENLYFQG"]

    def test_empty_sequence_after_header(self):
        text = ">seq1\n>seq2\nMTS\n"
        assert parse_fasta(text) == ["MTS"]

    def test_trailing_newline(self):
        text = ">seq1\nMTSENLYFQG\n\n"
        assert parse_fasta(text) == ["MTSENLYFQG"]


class TestValidate:
    def test_valid_standard(self):
        validate(["MTSENLYFQG", "ACDEFGHIK"])  # should not raise/exit

    def test_valid_single_aa(self):
        validate(["A"])
        validate(["C"])
        validate(["Y"])

    def test_valid_with_all_standard_aas(self):
        all_aas = "ACDEFGHIKLMNPQRSTVWY"
        validate([all_aas])

    def test_lowercase_is_valid(self):
        """内部做了 .upper()，小写字母合法。"""
        validate(["mtsenlyfqg"])

    def test_invalid_character_exits(self):
        with pytest.raises(SystemExit):
            validate(["MTSENLXFQG"])  # X 不是标准氨基酸

    def test_empty_string_is_valid(self):
        validate([""])

    def test_mixed_valid_and_invalid(self):
        with pytest.raises(SystemExit):
            validate(["MTSENLYFQG", "ACD*FGHIK"])  # * 非法

    def test_multiple_invalid_sequences(self):
        with pytest.raises(SystemExit):
            validate(["MTS", "ABC", "BAD1SEQ"])

    def test_numeric_character_exits(self):
        with pytest.raises(SystemExit):
            validate(["MTS3NLYFQG"])


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
