# ColabFold MSA 使用指南

调用 ColabFold MMseqs2 公共 API (`https://api.colabfold.com`) 生成蛋白质多序列比对 (MSA)，保存为 A3M 格式。

**零第三方依赖** — 纯 Python stdlib 实现，无需 CUDA / JAX / AlphaFold。

---

## 快速开始

### 单链 MSA

```python
from biorazer.access.server.colabfold_msa.pipeline import run_search

result = run_search(
    [("my_protein", "MTSENLYFQGAMGSMTSENLYFQGAMG")],
    "msa_out/",
)
# result.per_seq["my_protein"].files   -> ['msa_out/my_protein/unpaired/uniref.a3m', ...]
# result.per_seq["my_protein"].plots   -> ['msa_out/my_protein/unpaired/logo.png']
```

### 多链配对 MSA（用于 AF3 / OpenDDE 多聚体）

ColabFold 多链约定: **一条记录内以 `:` 分隔各链** (如 `>complex\nAAA:BBB`),
提交时每条链拆成独立记录 (`>101`, `>102`, ...) 到配对端点, 服务器端配对,
返回后合并为整复合物的一份 a3m。

```python
result = run_search(
    [("complex", "EVQLVESGGGLVQPGGSLRLSCAASGFTFS:DIQMTQSPSSLSASVGDRVTITCRASQGIR")],
    "ab_msa/",
    pair_mode="paired",
)
# 输出: ab_msa/complex/paired/pair.a3m  (整复合物配对 MSA, 每行 = 各链拼接)
#       ab_msa/complex/paired/pair_0.a3m, pair_1.a3m ... (每链原始段)
```

注意: 多条记录 = 多个独立任务 (各自提交, 互不配对)。同一复合物的多链
必须写在**一条记录**里用 `:` 分隔, 不要拆成多条记录。

### 从 FASTA 文件读取

```python
from biorazer.access.server.colabfold_msa.io import parse_fasta
from biorazer.access.server.colabfold_msa.pipeline import run_search

with open("input.fasta") as f:
    named_seqs = parse_fasta(f.read())

result = run_search(named_seqs, "msa_out/")
```

---

## Python API 参考

### `run_search()`

核心函数，提交任务到 MMseqs2 服务器，轮询等待，下载并解压结果。

```python
run_search(
    named_seqs: List[Tuple[str, str]],
    out_dir: str,
    pair_mode: str = "unpaired",
    use_env: bool = True,
    use_filter: bool = True,
    host: str = "https://api.colabfold.com",
    ua: str = "colabfold_msa/2.0 ...",
    pair_strategy: str = "greedy",
    local_rcsb_database: str | None = None,
) -> SearchResult
```

**参数说明**

| 参数 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `named_seqs` | `List[Tuple[str, str]]` | 必填 | (名称, 序列) 列表 |
| `out_dir` | `str` | 必填 | 输出目录（自动创建） |
| `pair_mode` | `str` | `"unpaired"` | `"unpaired"` / `"paired"` / `"paired+unpaired"` |
| `use_env` | `bool` | `True` | 是否使用环境数据库 (BFD/MGnify) |
| `use_filter` | `bool` | `True` | 是否过滤 MSA 结果 |
| `host` | `str` | `https://api.colabfold.com` | MMseqs2 服务器地址 |
| `ua` | `str` | ... | User-Agent 标识 |
| `pair_strategy` | `str` | `"greedy"` | 配对策略: `"greedy"` (快) / `"complete"` (全) |
| `local_rcsb_database` | `str \| None` | `None` | 本地 RCSB 镜像根目录（divided PDB 布局，如 `/mnt/data/public/RCSB`）。给出时模板结构直接从本地读取，不从 RCSB 下载 |

**返回值**

```python
SearchResult(
    per_seq = {
        "my_protein": SeqResult(
            files=["msa_out/my_protein/unpaired/uniref.a3m",
                   "msa_out/my_protein/unpaired/unpaired.a3m"],
            plots=["msa_out/my_protein/unpaired/logo.png"],
            report=""
        ),
    },
    merged="",                          # 预留: 顶层合并 (当前不产出)
    templates=None,                     # PDB 模板映射
)
```

### `parse_fasta()`

解析 FASTA 文本，返回 `[(名称, 序列), ...]`。

```python
parse_fasta(text: str) -> List[Tuple[str, str]]
```

- 跳过空行和注释
- 自动大写
- 支持多条序列
- 裸字符串输入 → `[("default", 序列)]`
- header 取第一个空白前的部分（`>my_protein something` → `"my_protein"`）

### `validate()`

验证序列只含 20 种标准氨基酸字符 (`ACDEFGHIKLMNPQRSTVWY`)。多链序列以 `:` 分隔各链（ColabFold 多链约定，如 `>complex\nAAA:BBB`），逐链验证。

```python
validate(seqs: List[Tuple[str, str]]) -> None
```

### `merge_a3m()`

合并多个 A3M 文件为一个字符串。

```python
merge_a3m(file_list: List[str]) -> str
```

---

## 输出目录结构

结果按模式分目录; 模式未请求或服务器无结果时, 该模式目录仍会创建,
但只包含 query 序列 (保证下游路径恒定)。处理顺序: **先对服务器返回的
原始文件 (pdb70/uniref/bfd) 按链段拆分, 再逐链合并** — 不产合并形态的
整复合物 a3m。

```
msa_out/
├── my_protein/
│   ├── unpaired/                          # unpaired 模式结果
│   │   ├── uniref.a3m                     # 服务器原始返回 (原样保留)
│   │   ├── bfd.mgnify30.metaeuk30.smag30.a3m  # 环境数据库原始返回 (use_env=True)
│   │   ├── unpaired_0.a3m, unpaired_1.a3m, ...  # 每链拆分 (各库段拼接, 不 padding; 单链时即 unpaired.a3m)
│   │   ├── logo_0.png / coverage_0.png, ...      # 每链各一张 (单链时 logo.png / coverage.png)
│   │   ├── pdb70.m8                     # 服务器原始返回 (原样保留)
│   │   ├── pdb70_0.m8, pdb70_1.m8, ...          # 模板 hit 按链拆分 (单链时无)
│   │   └── templates/                     # 模板 CIF, 按链命名 {pdb}_{chain}.cif
│   ├── paired/                            # paired 模式结果
│   │   ├── paired.a3m                     # 整复合物配对 MSA (各链段逐行拼接, 真实配对结果)
│   │   ├── paired_0.a3m, paired_1.a3m, ...      # 每链原始段 (单链时不拆分)
│   │   ├── logo.png / coverage.png              # 整复合物 (对比各链 MSA 差异)
│   │   ├── logo_0.png / coverage_0.png, ...      # 每链各一张 (单链时无 _N)
│   │   └── ...
│   └── msa.sh                             # 本次提交的配置脚本
└── .template_cache/                       # 模板下载缓存 (自动删除)
```

说明:
- 多链 complex (一条记录, `:` 分隔) 的各链结果都在**同一个** uniref.a3m /
  pair.a3m 里 (服务器按段返回, 段头 `>101`, `>102`, ...), 因此只有一个
  unpaired/ 和 paired/ 目录。本地按段头把原始文件拆成每链文件
  `unpaired_N.a3m` / `paired_N.a3m` (N 为链序号, 从 0 开始)。
- unpaired 的服务器原始文件 (uniref.a3m / bfd.*.a3m / pdb70.m8) 原样保留;
  拆分文件 = 原始文件按段头切出 (unpaired_N.a3m = uniref + bfd 对应链段
  拼接, 先拆分后合并; pdb70_N.m8 = pdb70.m8 对应链的 hit 行)。
- logo/coverage 对拆分后的每链文件绘制; 整复合物 paired.a3m 也出图
  (便于对比各链 MSA 差异)。
- 模板只由 msa (unpaired) ticket 附带, 故 `pdb70.m8` / `pdb70_N.m8` 与
  `templates/` 归 unpaired/; paired 模式不产出模板。

---

## 速率限制与重试

- 自动检测 `RATELIMIT` / `UNKNOWN` 状态，等待 5-10s 后重新提交
- HTTP 错误（连接超时、5xx）自动指数退避重试（最大 10 次，2s → 60s）
- 任务 3 次超时后放弃，递增 seed 重新提交
- 公共服务器 `api.colabfold.com` 有匿名速率限制，建议请求间隔至少 5s

---

## 与项目现有 MSA 工具的关系

```
biorazer/access/server/colabfold_msa/   ← 本模块: MSA 生成 (调用 ColabFold API)
├── http.py        ← HTTP 传输层 (POST/GET/下载, 纯 stdlib urllib)
├── api.py         ← MMseqs2 API 任务生命周期: 提交 → 轮询 → 下载
├── io.py          ← FASTA/A3M 解析、合并 (parse_fasta/validate/merge_a3m)
├── paired.py      ← 多链配对 (paired) 提交与结果合并
├── unpaired.py    ← 非配对 (unpaired) 提交与结果整理
├── template.py    ← 模板搜索: pdb70.m8 分发 + CIF 按链拆分 (网络/本地镜像)
├── pipeline.py    ← 流水线编排: run_search / SeqResult / SearchResult
└── cli.py         ← colabfold-msa 子命令 (register_subcommand)

biorazer/sequence/analysis/alignment/
├── plot.py        ← 可视化 MSA (coverage + 序列比对图)
├── report.py      ← MSA 分析报告
└── util.py        ← 对齐工具函数
```

其他 MSA 相关脚本：

| 文件 | 作用 |
|------|------|
| `biorazer/sequence/protein/scripts/MSA_visualizer.py` | 读 A3M → 绘制 coverage 图 |
| `biorazer/sequence/protein/scripts/msa_analyzer.py` | 统计每个位置的氨基酸频率 |

---

## 注意

1. **网络要求**：需要访问 `https://api.colabfold.com`（境外服务），建议在校园网 / 代理下使用
2. **无 CUDA 要求**：所有计算在远端完成，本地只需 Python stdlib + matplotlib + biotite
3. **A3M 格式**：结果为大写字母（匹配）+ 小写字母（插入）+ `-`（缺失），读取时注意小写 → `-` 转换
4. **序列名称**：从 FASTA header 第一个空白前提取，裸序列命名为 `"default"`
5. **tar.gz**：解压后自动删除
6. **Sequence logo**：使用 biotite 的 `SequenceProfile` + `plot_sequence_logo` 自动生成
7. **配对策略**：`greedy` 快但可能不是最优配对，`complete` 慢但完整
8. **本地 RCSB 镜像**：`--local-rcsb-database DIR`（或 `run_search(local_rcsb_database=...)`）可跳过模板网络下载——从 divided PDB 镜像（`DIR/{id[1:3]}/pdb{id}.ent.gz`）直接读取，经 biotite 转 CIF 后按链拆分，输出与下载模式完全一致
