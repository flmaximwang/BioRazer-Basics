"""Tests for biorazer.access 的 fetch CLI 子命令。

fetch 从 RCSB / UniProt / AlphaFold DB 下载条目文件:
- biorazer/access/cli.py 为 CLI 门面 (register_subcommand + runner);
- 各源 fetch 实现位于 biorazer/access/<SOURCE>/files.py。
所有测试断网运行: 通过 monkeypatch 替换源 fetch 函数与网络层。
"""

import argparse
from pathlib import Path

import pytest

from biorazer.access import cli
from biorazer.access.database.AFDB import files as afdb_files
from biorazer.access.database.AFDB.info import AFDBEntry


def _build_parser():
    """构建注册了 fetch 子命令的 argparse parser, 返回 (parser, subparsers)。"""
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command")
    cli.register_subcommand(sub)
    return parser, sub


class TestCliRegistration:
    """fetch 子命令的注册与参数校验。"""

    def test_register_fetch(self):
        """register_subcommand 注册 fetch, func 指向 runner。"""
        _, sub = _build_parser()
        assert "fetch" in sub.choices
        assert sub.choices["fetch"].get_default("func") is cli._run_fetch

    def test_source_choices(self):
        """--source 仅接受三个数据源, 大小写不敏感。"""
        parser, _ = _build_parser()
        args = parser.parse_args(["fetch", "3e2c", "--source", "rcsb", "--fmt", "pdb"])
        assert args.source == "RCSB"

    def test_source_invalid(self):
        """未知数据源直接报 argparse 错误。"""
        parser, _ = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["fetch", "3e2c", "--source", "ENSEMBL"])

    def test_source_required(self):
        """--source 必填。"""
        parser, _ = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["fetch", "3e2c"])

    def test_biorazer_cli_wires_access_cli(self):
        """顶层 biorazer.cli 引入 access_cli 注册 fetch。"""
        import biorazer.cli

        assert biorazer.cli.access_cli is cli


class TestResolveTarget:
    """-o 参数三种语义: 省略 / 目录 / 文件。"""

    def test_default_current_dir(self):
        """省略 -o -> 当前目录 ./{id}.{fmt}。"""
        assert cli._resolve_target("3e2c", "pdb", None) == Path("3e2c.pdb")

    def test_existing_dir(self, tmp_path):
        """已存在目录 -> 目录内 {id}.{fmt}。"""
        assert (
            cli._resolve_target("3e2c", "pdb", str(tmp_path))
            == tmp_path / "3e2c.pdb"
        )

    def test_trailing_slash(self):
        """以 / 结尾 (含未创建目录) -> 目录内 {id}.{fmt}。"""
        assert cli._resolve_target("3e2c", "pdb", "sub/") == Path("sub/3e2c.pdb")

    def test_file_path(self):
        """其他情况视为完整文件路径。"""
        assert cli._resolve_target("3e2c", "pdb", "out.pdb") == Path("out.pdb")


class TestRunFetch:
    """runner 全流程: 断网 (monkeypatch 源 fetch), 验证落盘与命名规则。"""

    def _run(self, argv, tmp_path, monkeypatch, source="RCSB"):
        calls = []

        def fake_fetch(identifier, fmt="pdb", download_dir=".", overwrite=False, logger=None):
            calls.append((identifier, fmt, download_dir, overwrite))
            Path(download_dir, f"{identifier}.{fmt}").write_bytes(b"DATA")

        monkeypatch.setitem(cli.SOURCES, source, (fake_fetch, cli.SOURCES[source][1]))
        monkeypatch.chdir(tmp_path)
        parser, _ = _build_parser()
        args = parser.parse_args(argv)
        args.func(args)
        return calls

    def test_default_target_in_cwd(self, tmp_path, monkeypatch):
        """缺省 -o: 文件落在当前目录, 名为 {id}.{fmt}。"""
        calls = self._run(["fetch", "3e2c", "--source", "RCSB"], tmp_path, monkeypatch)
        assert (tmp_path / "3e2c.pdb").read_bytes() == b"DATA"
        # 源 fetch 收到的是临时目录 + 强制覆盖, 临时文件随后移动到目标
        assert len(calls) == 1
        assert calls[0][0] == "3e2c" and calls[0][1] == "pdb" and calls[0][3] is True

    def test_output_file_name(self, tmp_path, monkeypatch):
        """-o 为文件名: 下载到该路径, 默认名 {id}.{fmt} 不残留。"""
        self._run(["fetch", "3e2c", "--source", "RCSB", "-o", "out.pdb"], tmp_path, monkeypatch)
        assert (tmp_path / "out.pdb").read_bytes() == b"DATA"
        assert not (tmp_path / "3e2c.pdb").exists()

    def test_output_dir(self, tmp_path, monkeypatch):
        """-o 为目录: 下载到目录内, 名为 {id}.{fmt}。"""
        sub = tmp_path / "sub"
        sub.mkdir()
        self._run(
            ["fetch", "3e2c", "--source", "RCSB", "-o", str(sub)],
            tmp_path, monkeypatch,
        )
        assert (sub / "3e2c.pdb").read_bytes() == b"DATA"

    def test_output_dir_trailing_slash(self, tmp_path, monkeypatch):
        """-o 以 / 结尾 (目录未创建): 视为目录并创建。"""
        self._run(
            ["fetch", "3e2c", "--source", "RCSB", "-o", "sub/"],
            tmp_path, monkeypatch,
        )
        assert (tmp_path / "sub" / "3e2c.pdb").read_bytes() == b"DATA"

    def test_skip_if_exists(self, tmp_path, monkeypatch):
        """目标已存在且未 --overwrite: 跳过, 不调用源 fetch。"""
        (tmp_path / "3e2c.pdb").write_bytes(b"OLD")
        calls = self._run(["fetch", "3e2c", "--source", "RCSB"], tmp_path, monkeypatch)
        assert calls == []
        assert (tmp_path / "3e2c.pdb").read_bytes() == b"OLD"

    def test_overwrite(self, tmp_path, monkeypatch):
        """--overwrite: 重新下载覆盖。"""
        (tmp_path / "3e2c.pdb").write_bytes(b"OLD")
        calls = self._run(
            ["fetch", "3e2c", "--source", "RCSB", "--overwrite"],
            tmp_path, monkeypatch,
        )
        assert len(calls) == 1
        assert (tmp_path / "3e2c.pdb").read_bytes() == b"DATA"

    def test_uniprot_default_fmt(self, tmp_path, monkeypatch):
        """UniProt 缺省格式为 fasta。"""
        calls = self._run(
            ["fetch", "P0DP23", "--source", "UNIPROT"], tmp_path, monkeypatch,
            source="UNIPROT",
        )
        assert calls[0][1] == "fasta"
        assert (tmp_path / "P0DP23.fasta").read_bytes() == b"DATA"

    def test_afdb_default_fmt(self, tmp_path, monkeypatch):
        """AFDB 缺省格式为 cif。"""
        calls = self._run(
            ["fetch", "P0DP23", "--source", "AFDB"], tmp_path, monkeypatch,
            source="AFDB",
        )
        assert calls[0][1] == "cif"
        assert (tmp_path / "P0DP23.cif").read_bytes() == b"DATA"

    def test_fmt_explicit(self, tmp_path, monkeypatch):
        """显式 --fmt 覆盖缺省格式。"""
        calls = self._run(
            ["fetch", "3e2c", "--source", "RCSB", "--fmt", "cif"],
            tmp_path, monkeypatch,
        )
        assert calls[0][1] == "cif"
        assert (tmp_path / "3e2c.cif").read_bytes() == b"DATA"

    def test_valueerror_exits(self, tmp_path, monkeypatch, capsys):
        """源 fetch 抛 ValueError (如不支持格式): 打印错误并以状态码 1 退出。"""

        def bad_fetch(identifier, fmt="pdb", download_dir=".", overwrite=False, logger=None):
            raise ValueError(f"Unsupported format: {fmt}")

        monkeypatch.setitem(cli.SOURCES, "RCSB", (bad_fetch, "pdb"))
        monkeypatch.chdir(tmp_path)
        parser, _ = _build_parser()
        args = parser.parse_args(["fetch", "3e2c", "--source", "RCSB"])
        with pytest.raises(SystemExit) as exc:
            args.func(args)
        assert exc.value.code == 1
        assert "错误: Unsupported format: pdb" in capsys.readouterr().err


class TestAfdbFetch:
    """AFDB fetch: 断网 (monkeypatch uniprot_to_entries 与 AFDBEntry.download)。"""

    def _entry(self, uid):
        return AFDBEntry(data={"entryId": f"AF-{uid}-F1-model_v4", "cifUrl": "http://x.cif"})

    def test_fetch_writes_entry_file(self, tmp_path, monkeypatch):
        """按 entryId 命名落盘并返回路径。"""
        monkeypatch.setattr(
            afdb_files, "uniprot_to_entries", lambda uid: [self._entry("P0DP23")]
        )

        def fake_download(self_, file_type, folder_dir):
            Path(folder_dir, f"{self_.id}.{file_type}").write_bytes(b"CIF")

        monkeypatch.setattr(AFDBEntry, "download", fake_download)
        out = afdb_files.fetch("P0DP23", fmt="cif", download_dir=str(tmp_path))
        assert out == tmp_path / "AF-P0DP23-F1-model_v4.cif"
        assert out.read_bytes() == b"CIF"

    def test_skip_if_exists(self, tmp_path, monkeypatch):
        """已存在且未 overwrite: 不调用 download, 直接返回。"""
        monkeypatch.setattr(
            afdb_files, "uniprot_to_entries", lambda uid: [self._entry("P0DP23")]
        )
        target = tmp_path / "AF-P0DP23-F1-model_v4.cif"
        target.write_bytes(b"OLD")

        def fake_download(self_, file_type, folder_dir):
            raise AssertionError("不应调用 download")

        monkeypatch.setattr(AFDBEntry, "download", fake_download)
        out = afdb_files.fetch("P0DP23", download_dir=str(tmp_path))
        assert out == target and target.read_bytes() == b"OLD"

    def test_unsupported_fmt(self, tmp_path):
        """不支持的格式抛 ValueError。"""
        with pytest.raises(ValueError, match="Unsupported format"):
            afdb_files.fetch("P0DP23", fmt="xyz", download_dir=str(tmp_path))

    def test_no_entry(self, tmp_path, monkeypatch):
        """API 无条目 (或返回非条目数据) 抛 ValueError。"""
        monkeypatch.setattr(afdb_files, "uniprot_to_entries", lambda uid: [])
        with pytest.raises(ValueError, match="No AFDB entry"):
            afdb_files.fetch("P0DP23", download_dir=str(tmp_path))

    def test_apierror_dict_not_entry(self, tmp_path, monkeypatch):
        """API 报错返回 dict (如 {"detail": ...}) 时同样报无条目, 而非 TypeError。"""
        monkeypatch.setattr(afdb_files, "uniprot_to_entries", lambda uid: {"detail": "Not found"})
        with pytest.raises(ValueError, match="No AFDB entry"):
            afdb_files.fetch("P0DP23", download_dir=str(tmp_path))
