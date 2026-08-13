"""Tests for colabfold_msa output layout (unpaired/ vs paired/ isolation).

不联网: mock _submit_and_download 生成假 tar.gz, 验证:
- 各模式结果写入 unpaired/ 与 paired/ 子目录;
- a3m 均为每链拆分文件 (unpaired_N.a3m / paired_N.a3m), 不产合并形态;
- 单链时无 _N 后缀 (unpaired.a3m / paired.a3m);
- 未请求的模式目录仍建, 仅含 query 序列;
- logo/coverage 对拆分后的每链文件绘制;
- 模板 (pdb70_N.m8) 按链拆分并归 unpaired/。
"""

import io
import os
import tarfile

import pytest

from biorazer.access.server.colabfold_msa import paired as paired_mod
from biorazer.access.server.colabfold_msa import unpaired as unpaired_mod
from biorazer.access.server.colabfold_msa.pipeline import run_search


def _fake_tar_gz(path, files):
    """写一个包含 files {name: content} 的 tar.gz"""
    with tarfile.open(path, "w:gz") as tar:
        for name, content in files.items():
            data = content.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def _make_submitter(tar_files):
    def _submit_and_download(seqs, out_dir, endpoint, mode, host, ua,
                             tar_name="out.tar.gz"):
        path = os.path.join(out_dir, tar_name)
        _fake_tar_gz(path, tar_files)
        return path, list(range(101, 101 + len(seqs)))
    return _submit_and_download


# 假服务器返回: 2 链, 每链 1 行 query + 1 行 hit
PAIRED_TAR = {
    "pair.a3m": ">101\nAAA\n>|uniref|hit1\nAAA\n"
                ">102\nBBB\n>|uniref|hit2\nBBB\n",
}
UNPAIRED_TAR = {
    "uniref.a3m": ">101\nAAA\n>|uniref|hit1\nAAA\n"
                  ">102\nBBB\n>|uniref|hit2\nBBB\n",
}
SINGLE_UNPAIRED_TAR = {
    "uniref.a3m": ">101\nAAA\n>|uniref|hit1\nAAA\n",
}


@pytest.fixture
def mock_pair(monkeypatch):
    monkeypatch.setattr(paired_mod, "_submit_and_download",
                        _make_submitter(PAIRED_TAR))
    yield


@pytest.fixture
def mock_unpaired(monkeypatch):
    monkeypatch.setattr(unpaired_mod, "_submit_and_download",
                        _make_submitter(UNPAIRED_TAR))
    yield


class TestLayoutPairedPlusUnpaired:
    def test_separate_mode_dirs_no_merged(self, tmp_path, mock_pair, mock_unpaired):
        run_search([("complex", "AAA:BBB")], str(tmp_path),
                   pair_mode="paired+unpaired",
                   use_env=False, use_filter=False,
                   host="mock", ua="mock")
        base = tmp_path / "complex"
        # 每链拆分文件
        assert (base / "unpaired" / "unpaired_0.a3m").is_file()
        assert (base / "unpaired" / "unpaired_1.a3m").is_file()
        assert (base / "paired" / "paired_0.a3m").is_file()
        assert (base / "paired" / "paired_1.a3m").is_file()
        # 服务器返回的原始文件原样保留
        assert (base / "unpaired" / "uniref.a3m").is_file()
        assert (base / "unpaired" / "uniref.a3m").read_text() \
            == ">101\nAAA\n>|uniref|hit1\nAAA\n>102\nBBB\n>|uniref|hit2\nBBB\n"
        # 不产 padded 的整复合物合并形态
        assert not (base / "unpaired" / "unpaired.a3m").exists()
        assert not (base / "merged.a3m").exists()
        # paired 保留整复合物 merged 版本 (真实配对结果)
        assert (base / "paired" / "paired.a3m").is_file()

    def test_per_chain_files_content(self, tmp_path, mock_pair, mock_unpaired):
        run_search([("complex", "AAA:BBB")], str(tmp_path),
                   pair_mode="paired+unpaired",
                   use_env=False, use_filter=False,
                   host="mock", ua="mock")
        base = tmp_path / "complex"
        # 每链段: 各自 query + hit (不 padding)
        assert (base / "unpaired" / "unpaired_0.a3m").read_text() \
            == ">101\nAAA\n>|uniref|hit1\nAAA\n"
        assert (base / "unpaired" / "unpaired_1.a3m").read_text() \
            == ">102\nBBB\n>|uniref|hit2\nBBB\n"
        assert (base / "paired" / "paired_0.a3m").read_text() \
            == ">101\nAAA\n>|uniref|hit1\nAAA\n"
        assert (base / "paired" / "paired_1.a3m").read_text() \
            == ">102\nBBB\n>|uniref|hit2\nBBB\n"

    def test_plots_per_split_file(self, tmp_path, mock_pair, mock_unpaired):
        result = run_search([("complex", "AAA:BBB")], str(tmp_path),
                            pair_mode="paired+unpaired",
                            use_env=False, use_filter=False,
                            host="mock", ua="mock")
        base = tmp_path / "complex"
        # logo/coverage 对拆分后的每链文件绘制; 整复合物 paired.a3m 也出图
        assert (base / "unpaired" / "logo_0.png").is_file()
        assert (base / "unpaired" / "coverage_0.png").is_file()
        assert (base / "unpaired" / "logo_1.png").is_file()
        assert (base / "paired" / "logo_0.png").is_file()
        assert (base / "paired" / "logo_1.png").is_file()
        assert (base / "paired" / "logo.png").is_file()
        plots = result.per_seq["complex"].plots
        assert any(p.endswith("unpaired/logo_0.png") for p in plots)
        assert any(p.endswith("paired/logo_1.png") for p in plots)
        assert any(p.endswith("paired/logo.png") for p in plots)

    def test_files_list_uses_mode_dirs(self, tmp_path, mock_pair, mock_unpaired):
        result = run_search([("complex", "AAA:BBB")], str(tmp_path),
                            pair_mode="paired+unpaired",
                            use_env=False, use_filter=False,
                            host="mock", ua="mock")
        files = result.per_seq["complex"].files
        assert any(f.endswith("unpaired/unpaired_0.a3m") for f in files)
        assert any(f.endswith("unpaired/unpaired_1.a3m") for f in files)
        assert any(f.endswith("paired/paired.a3m") for f in files)
        assert any(f.endswith("paired/paired_0.a3m") for f in files)
        assert any(f.endswith("paired/paired_1.a3m") for f in files)


class TestLayoutSingleMode:
    def test_unpaired_only(self, tmp_path, mock_unpaired):
        run_search([("complex", "AAA:BBB")], str(tmp_path), pair_mode="unpaired",
                   use_env=False, use_filter=False, host="mock", ua="mock")
        base = tmp_path / "complex"
        assert (base / "unpaired" / "unpaired_0.a3m").is_file()
        # paired 未请求: 目录仍建, 仅 query 序列
        assert (base / "paired" / "paired.a3m").is_file()
        assert (base / "paired" / "paired_0.a3m").read_text() == ">101\nAAA\n"
        assert (base / "paired" / "paired_1.a3m").read_text() == ">102\nBBB\n"

    def test_paired_only(self, tmp_path, mock_pair):
        run_search([("complex", "AAA:BBB")], str(tmp_path), pair_mode="paired",
                   use_env=False, use_filter=False, host="mock", ua="mock")
        base = tmp_path / "complex"
        assert (base / "paired" / "paired.a3m").is_file()
        assert (base / "paired" / "paired_0.a3m").is_file()
        # unpaired 未请求: 目录仍建, 仅 query 序列, 不出图
        assert (base / "unpaired" / "unpaired_0.a3m").read_text() == ">101\nAAA\n"
        assert (base / "unpaired" / "unpaired_1.a3m").read_text() == ">102\nBBB\n"
        assert not (base / "unpaired" / "logo_0.png").exists()

    def test_single_chain_no_suffix(self, tmp_path, monkeypatch):
        """单链: a3m 无 _N 后缀, 图名无序号"""
        monkeypatch.setattr(unpaired_mod, "_submit_and_download",
                            _make_submitter(SINGLE_UNPAIRED_TAR))
        result = run_search([("prot", "AAA")], str(tmp_path), pair_mode="unpaired",
                            use_env=False, use_filter=False, host="mock", ua="mock")
        base = tmp_path / "prot"
        assert (base / "unpaired" / "unpaired.a3m").is_file()
        assert not (base / "unpaired" / "unpaired_0.a3m").exists()
        assert (base / "unpaired" / "logo.png").is_file()
        assert not (base / "unpaired" / "logo_0.png").exists()
        assert any(f.endswith("unpaired/unpaired.a3m") for f in
                   result.per_seq["prot"].files)


class TestTemplatesPerChain:
    def test_pdb70_split_by_chain(self, tmp_path, mock_unpaired):
        """原始 pdb70.m8 保留, 同时按链拆分出 pdb70_N.m8, 归 unpaired/"""
        monkeypatch = pytest.MonkeyPatch()
        tar_files = dict(UNPAIRED_TAR)
        tar_files["pdb70.m8"] = (
            "101\t1abc_A\t1.0\t0.5\t10\t5\t1\t2\t3\t4\t5\t6\n"
            "102\t2def_A\t1.0\t0.5\t10\t5\t1\t2\t3\t4\t5\t6\n"
        )
        monkeypatch.setattr(unpaired_mod, "_submit_and_download",
                            _make_submitter(tar_files))
        try:
            run_search([("complex", "AAA:BBB")], str(tmp_path),
                       pair_mode="unpaired",
                       use_env=False, use_filter=False,
                       host="mock", ua="mock")
        finally:
            monkeypatch.undo()
        base = tmp_path / "complex"
        # 服务器返回的原始 pdb70.m8 原样保留
        assert (base / "unpaired" / "pdb70.m8").is_file()
        assert "101\t" in (base / "unpaired" / "pdb70.m8").read_text()
        assert "102\t" in (base / "unpaired" / "pdb70.m8").read_text()
        # 按链拆分文件
        assert (base / "unpaired" / "pdb70_0.m8").is_file()
        assert (base / "unpaired" / "pdb70_1.m8").is_file()
        assert not (base / "pdb70.m8").exists()
        # 各链文件只含对应链的 hit
        assert "101\t" in (base / "unpaired" / "pdb70_0.m8").read_text()
        assert "102\t" in (base / "unpaired" / "pdb70_1.m8").read_text()

    def test_pdb70_single_chain_no_suffix(self, tmp_path, monkeypatch):
        tar_files = dict(SINGLE_UNPAIRED_TAR)
        tar_files["pdb70.m8"] = "101\t1abc_A\t1.0\t0.5\t10\t5\t1\t2\t3\t4\t5\t6\n"
        monkeypatch.setattr(unpaired_mod, "_submit_and_download",
                            _make_submitter(tar_files))
        run_search([("prot", "AAA")], str(tmp_path), pair_mode="unpaired",
                   use_env=False, use_filter=False, host="mock", ua="mock")
        base = tmp_path / "prot"
        # 单链: 原始 pdb70.m8 保留, 不另拆 pdb70_0.m8
        assert (base / "unpaired" / "pdb70.m8").is_file()
        assert not (base / "unpaired" / "pdb70_0.m8").exists()
