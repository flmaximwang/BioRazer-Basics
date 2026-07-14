"""pytest fixtures for the biorazer test suite."""

import tempfile
import os
import pytest


@pytest.fixture
def tmp_msa_dir():
    """提供临时目录用于 MSA 输出测试，自动清理。"""
    with tempfile.TemporaryDirectory(prefix="biorazer_test_msa_") as d:
        yield d


@pytest.fixture
def sample_a3m_files():
    """创建两个临时的 A3M 文件用于 merge_a3m 测试。"""
    files = []
    contents = [
        (">seq1\nMTSENLYFQG\n", "a3m_1.a3m"),
        (">seq2\nWPKL\n", "a3m_2.a3m"),
    ]
    for content, name in contents:
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=f"_{name}", delete=False
        )
        f.write(content)
        f.close()
        files.append(f.name)
    yield files
    for p in files:
        os.unlink(p)
