"""colabfold_msa.template - 模板搜索结果的整理: pdb70.m8 分发 + 模板结构。

模板结构来自两种来源:
- 网络: 从 files.rcsb.org 下载 CIF, 按链拆分 (_download_and_split_templates);
- 本地: 从 RCSB divided PDB 镜像读取 PDB → biotite 转 CIF, 按链拆分
  (_split_templates_from_local, 由 local_rcsb_database 触发)。
两种模式输出格式一致 (CIF), 进入各序列 unpaired/ 子目录的 templates/。
pdb70.m8 只由 msa (unpaired) ticket 附带, 故模板整体归 unpaired/。
"""

import gzip
import io
import os
import shutil
import sys
import time
from typing import Dict, List, Optional, Tuple
from urllib.request import Request, urlopen

from biotite.structure import AtomArray

from biorazer.structure.io.protein import AtomArray_Cif, Pdb_AtomArray


# ── 模板下载与链拆分 ────────────────────────────────
def _extract_chain_from_cif(cif_path: str, chain: str, out_path: str) -> bool:
    """从 mmCIF 文件中提取指定链/Entity 的原子记录，写入新 CIF。

    chain 参数可以是:
    - auth_asym_id (如 "C", "a") — 匹配 auth_asym_id 列
    - entity_id (数字, 如 "1") — 提取该 entity 的所有 chain

    自动处理大小写不匹配。
    """
    try:
        with open(cif_path) as f:
            lines = f.readlines()

        # === 第一遍：扫描结构 ===
        # 1) 解析 entity → strand 映射
        entity_strands: dict = {}  # entity_id → [auth_asym_id, ...]
        current_eid = None
        for line in lines:
            if line.startswith("_entity_poly.entity_id"):
                val = line.split(None, 1)
                if len(val) >= 2:
                    current_eid = val[1].strip().rstrip(";")
            elif line.startswith("_entity_poly.pdbx_strand_id") and current_eid is not None:
                val = line.split(None, 1)
                if len(val) >= 2:
                    strand_str = val[1].strip().rstrip(";")
                    entity_strands[current_eid] = [s.strip() for s in strand_str.split(",")]

        # 2) 解析 atom_site 列索引
        label_col = auth_col = -1
        in_atom_loop = False
        cols: List[str] = []
        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; cols = []; continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    col_name = line.split()[0]
                    cols.append(col_name)
                    if col_name == "_atom_site.label_asym_id": label_col = len(cols) - 1
                    if col_name == "_atom_site.auth_asym_id":  auth_col = len(cols) - 1
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False

        # 3) 收集所有可用的 auth_asym_id
        auth_set: set = set()
        in_atom_loop = False; idx = 0
        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; idx = 0; continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    idx += 1
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False
                elif auth_col >= 0 and line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    if len(parts) > auth_col:
                        auth_set.add(parts[auth_col])

        # === 匹配目标 ===
        # chain 参数可能是 auth_asym_id (字母) 或 entity_id (数字)
        target_chains: list = []

        # 先尝试作为 entity_id
        if chain in entity_strands:
            target_chains = entity_strands[chain]
        else:
            # 尝试作为 auth_asym_id (精确匹配)
            if chain in auth_set:
                target_chains = [chain]
            else:
                # 尝试忽略大小写匹配 auth_asym_id
                for ac in auth_set:
                    if ac.lower() == chain.lower():
                        target_chains.append(ac)
                # 最后尝试 label_asym_id 忽略大小写
                if not target_chains and label_col >= 0:
                    label_set: set = set()
                    in_atom_loop = False; idx = 0
                    for line in lines:
                        if line.startswith("loop_"):
                            in_atom_loop = True; idx = 0; continue
                        if in_atom_loop:
                            if line.startswith("_atom_site."): idx += 1
                            elif line.strip()=="" or line.startswith("#"): in_atom_loop=False
                            elif line.startswith(("ATOM","HETATM")):
                                parts=line.split()
                                if len(parts)>label_col: label_set.add(parts[label_col])
                    for lc in label_set:
                        if lc.lower() == chain.lower() and lc not in target_chains:
                            target_chains.append(lc)

        if not target_chains:
            print(f"  [!] CIF {cif_path} 中未找到 chain/entity '{chain}' "
                  f"(auth: {sorted(auth_set)})", file=sys.stderr)
            return False

        # === 第二遍：筛选原子行 ===
        in_atom_loop = False; idx = 0
        auth_col2 = label_col2 = -1
        cols2: List[str] = []
        out: List[str] = []

        def _should_keep(parts: list) -> bool:
            for tc in target_chains:
                if auth_col2 >= 0 and len(parts) > auth_col2 and parts[auth_col2] == tc:
                    return True
                if label_col2 >= 0 and len(parts) > label_col2 and parts[label_col2] == tc:
                    return True
            return False

        for line in lines:
            if line.startswith("loop_"):
                in_atom_loop = True; idx = 0; cols2 = []; auth_col2 = label_col2 = -1
                out.append(line); continue
            if in_atom_loop:
                if line.startswith("_atom_site."):
                    col_name = line.split()[0]; cols2.append(col_name)
                    if col_name == "_atom_site.label_asym_id": label_col2 = len(cols2) - 1
                    if col_name == "_atom_site.auth_asym_id":  auth_col2 = len(cols2) - 1
                    out.append(line)
                elif line.strip() == "" or line.startswith("#"):
                    in_atom_loop = False; out.append(line)
                else:
                    parts = line.split()
                    if _should_keep(parts):
                        out.append(line)
            else:
                out.append(line)

        atom_count = sum(1 for l in out if l.startswith(("ATOM", "HETATM")))
        if atom_count == 0:
            print(f"  [!] CIF {cif_path}: 匹配 chain(s) {target_chains} 后无原子记录",
                  file=sys.stderr)
            return False

        with open(out_path, "w") as f:
            f.writelines(out)
        return True
    except Exception as e:
        print(f"  [!] CIF 链过滤失败 {cif_path} [{chain}]: {e}", file=sys.stderr)
        return False


# ── 本地 RCSB 镜像读取 ──────────────────────────────
def _find_local_pdb(db_root: str, pdb_id: str) -> Optional[str]:
    """在本地 RCSB rsync 镜像中定位 PDB 文件。

    标准 RCSB divided 布局: {db_root}/{id[1:3]}/pdb{id}.ent.gz
    (由 rsync.rcsb.org:33444 ftp_data/structures/divided/pdb/ 同步而来)。
    """
    if len(pdb_id) < 3:
        return None
    candidates = [
        os.path.join(db_root, pdb_id[1:3], f"pdb{pdb_id}.ent.gz"),
        os.path.join(db_root, pdb_id[1:3], f"pdb{pdb_id}.ent"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def _split_templates_from_local(
    chain_templates: Dict[int, List[Tuple[str, Optional[str]]]],
    Ms: List[int],
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    db_root: str,
) -> None:
    """从本地 RCSB 镜像读取 PDB → 转 CIF → 按链拆分到各序列的 templates/ 目录。

    与下载模式输出格式一致 (CIF)，只把网络下载替换为本地读取。
    """
    # 收集所有不重复的 PDB ID
    all_pdbs: set = set()
    for entries in chain_templates.values():
        for pdb_id, _ in entries:
            all_pdbs.add(pdb_id.lower())
    if not all_pdbs:
        return

    pdb_list = sorted(all_pdbs)
    total = len(pdb_list)
    print(f"  [→] 本地模板: {total} 个 PDB (镜像 {db_root})", file=sys.stderr)

    cache = os.path.join(out_dir, ".template_cache")
    os.makedirs(cache, exist_ok=True)

    # 逐个从本地镜像读取 PDB → AtomArray
    arrays: Dict[str, AtomArray] = {}
    missing: List[str] = []
    for pdb_id in pdb_list:
        src = _find_local_pdb(db_root, pdb_id)
        if src is None:
            missing.append(pdb_id)
            continue
        try:
            if src.endswith(".gz"):
                with gzip.open(src, "rt") as f:
                    text = f.read()
            else:
                with open(src) as f:
                    text = f.read()
            buf = io.StringIO(text)
            arrays[pdb_id] = Pdb_AtomArray(buf, io.StringIO()).read()
        except Exception as e:
            print(f"  [!] {pdb_id} 读取失败: {e}", file=sys.stderr)
            missing.append(pdb_id)
    if missing:
        print(f"  [!] 本地镜像缺少: {', '.join(missing)}", file=sys.stderr)

    # 为每条序列从 AtomArray 中提取对应链，转 CIF
    chain_downloaded = 0
    for i, (name, _) in enumerate(named_seqs):
        M = Ms[i] if i < len(Ms) else -1
        if M not in chain_templates:
            continue
        seq_dir = os.path.join(out_dir, name, "unpaired")
        tpl_dir = os.path.join(seq_dir, "templates")
        os.makedirs(tpl_dir, exist_ok=True)
        chain_count = 0
        for pdb_id, chain in chain_templates[M]:
            pdb_id = pdb_id.lower()
            arr = arrays.get(pdb_id)
            if arr is None:
                continue
            try:
                if chain:
                    sub = arr[arr.chain_id == chain]
                    if sub.array_length() == 0:
                        print(f"  [!] {pdb_id} 无链 {chain}，跳过", file=sys.stderr)
                        continue
                    dst = os.path.join(tpl_dir, f"{pdb_id}_{chain}.cif")
                else:
                    sub = arr
                    dst = os.path.join(tpl_dir, f"{pdb_id}.cif")
                buf = io.StringIO()
                AtomArray_Cif(io.StringIO(), buf).write(sub)
                with open(dst, "w") as f:
                    f.write(buf.getvalue())
                chain_count += 1
            except Exception as e:
                print(f"  [!] {pdb_id} CIF 写出失败: {e}", file=sys.stderr)
        chain_downloaded += chain_count
        if chain_count > 0:
            print(f"    [✓] {name}: {chain_count} 个模板 -> {tpl_dir}", file=sys.stderr)

    # 清理缓存（本次未用，保留占位目录以防其它流程引用）
    if chain_downloaded > 0:
        print(f"  [✓] 模板: {chain_downloaded} 个 CIF (按链拆分, 本地镜像)", file=sys.stderr)


def _download_and_split_templates(
    chain_templates: Dict[int, List[Tuple[str, Optional[str]]]],
    Ms: List[int],
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    host: str, ua: str,
    local_rcsb_database: Optional[str] = None,
) -> None:
    """从 RCSB PDB 下载或本地镜像复制 CIF → 按链拆分到各序列的 templates/ 目录。

    当 local_rcsb_database 给出时，跳过网络下载，直接从本地 RCSB 镜像读取
    (divided PDB 布局: {db}/{id[1:3]}/pdb{id}.ent.gz)，输出格式仍为 CIF。
    """
    # 收集所有不重复的 PDB ID
    all_pdbs: set = set()
    for entries in chain_templates.values():
        for pdb_id, _ in entries:
            all_pdbs.add(pdb_id.lower())
    if not all_pdbs:
        return

    if local_rcsb_database:
        _split_templates_from_local(
            chain_templates, Ms, named_seqs, out_dir, local_rcsb_database
        )
        return

    pdb_list = sorted(all_pdbs)
    total = len(pdb_list)
    print(f"  [→] 下载模板: {total} 个 PDB (RCSB)", file=sys.stderr)

    cache = os.path.join(out_dir, ".template_cache")
    os.makedirs(cache, exist_ok=True)

    # 逐个从 RCSB 下载 CIF
    downloaded = 0
    failed = 0
    rcsb_url = "https://files.rcsb.org/download"
    for pdb_id in pdb_list:
        dst = os.path.join(cache, f"{pdb_id}.cif")
        if os.path.isfile(dst):
            downloaded += 1
            continue
        url = f"{rcsb_url}/{pdb_id}.cif"
        ok = False
        for attempt in range(3):
            try:
                req = Request(url, headers={"User-Agent": ua})
                with urlopen(req, timeout=120) as resp:
                    with open(dst, "wb") as f:
                        while True:
                            chunk = resp.read(8192)
                            if not chunk:
                                break
                            f.write(chunk)
                # 验证是有效 CIF（含 data_ 行）
                with open(dst) as f:
                    first = f.read(200)
                if "data_" in first:
                    ok = True
                    break
                else:
                    os.remove(dst)
                    raise ValueError("不是有效 CIF 文件")
            except Exception as e:
                if attempt < 2:
                    print(f"    [!] {pdb_id} 下载失败 (尝试 {attempt + 1}/3): {e}",
                          file=sys.stderr)
                    time.sleep(3)
        if ok:
            downloaded += 1
        else:
            failed += 1
        # 进度
        done = downloaded + failed
        if done % max(1, total // 10) == 0 or done == total:
            print(f"    [{done}/{total}] 已下载 {downloaded}，失败 {failed}",
                  file=sys.stderr)

    # 为每条序列从 CIF 中提取对应链
    chain_downloaded = 0
    for i, (name, _) in enumerate(named_seqs):
        M = Ms[i] if i < len(Ms) else -1
        if M not in chain_templates:
            continue
        seq_dir = os.path.join(out_dir, name, "unpaired")
        tpl_dir = os.path.join(seq_dir, "templates")
        os.makedirs(tpl_dir, exist_ok=True)
        chain_count = 0
        for pdb_id, chain in chain_templates[M]:
            cif_name = f"{pdb_id}.cif"
            src = os.path.join(cache, cif_name)
            if not os.path.isfile(src):
                continue
            if chain:
                dst = os.path.join(tpl_dir, f"{pdb_id}_{chain}.cif")
                if _extract_chain_from_cif(src, chain, dst):
                    chain_count += 1
            else:
                dst = os.path.join(tpl_dir, cif_name)
                # copy2 会保留元数据，但 CIF 内容更重要
                with open(src) as f_in:
                    content = f_in.read()
                with open(dst, "w") as f_out:
                    f_out.write(content)
                chain_count += 1
        chain_downloaded += chain_count
        if chain_count > 0:
            print(f"    [✓] {name}: {chain_count} 个模板 -> {tpl_dir}", file=sys.stderr)

    # 清理缓存
    shutil.rmtree(cache, ignore_errors=True)
    if chain_downloaded > 0:
        print(f"  [✓] 模板: {chain_downloaded} 个 CIF (按链拆分)", file=sys.stderr)


# ── 模板 (pdb70.m8 → 各序列目录 + 下载/本地读取) ──────
def _handle_templates(tmp_dir: str, out_dir: str, Ms: List[int],
                      named_seqs: List[Tuple[str, str]],
                      host: str, ua: str,
                      local_rcsb_database: Optional[str]) -> None:
    """解析 tar 内 pdb70.m8 → 各序列目录 pdb70.m8 + 模板结构。

    注意: 仅 unpaired 的 msa ticket 附带 pdb70.m8; pair ticket 不做
    模板搜索 (与 colabfold 一致: use_pairing 时 use_templates=False)。
    """
    m8_path = os.path.join(tmp_dir, "pdb70.m8")
    if not os.path.isfile(m8_path):
        return
    chain_templates: Dict[int, List[Tuple[str, Optional[str]]]] = {}
    chain_m8_lines: Dict[int, List[str]] = {}
    with open(m8_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                M = int(parts[0])
                target = parts[1]
                chain: Optional[str] = None
                pdb_id = target
                if "_" in target:
                    pdb_id, chain = target.split("_", 1)
                chain_templates.setdefault(M, []).append((pdb_id, chain))
                chain_m8_lines.setdefault(M, []).append(line)
    print(f"  [→] 模板搜索: {sum(len(v) for v in chain_templates.values())} 个 hit",
          file=sys.stderr)
    # 服务器返回的原始 pdb70.m8 原样保留, 同时按链拆分出 pdb70_N.m8
    # (多链 → pdb70_0.m8, pdb70_1.m8, ...; 单链 → 仅原始 pdb70.m8)
    for name, _ in named_seqs:
        seq_dir = os.path.join(out_dir, name, "unpaired")
        os.makedirs(seq_dir, exist_ok=True)
        if os.path.isfile(m8_path):
            shutil.copy2(m8_path, os.path.join(seq_dir, "pdb70.m8"))
        if len(Ms) > 1:
            for i, M in enumerate(Ms):
                lines_out = chain_m8_lines.get(M, [])
                if not lines_out:
                    continue
                m8_out = os.path.join(seq_dir, f"pdb70_{i}.m8")
                with open(m8_out, "w") as f:
                    for raw_line in lines_out:
                        f.write(raw_line + "\n")
    _download_and_split_templates(chain_templates, Ms, named_seqs, out_dir, host, ua,
                                  local_rcsb_database=local_rcsb_database)