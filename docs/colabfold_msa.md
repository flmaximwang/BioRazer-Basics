# ColabFold MSA 使用指南

调用 ColabFold MMseqs2 公共 API (`https://api.colabfold.com`) 生成蛋白质多序列比对 (MSA)，保存为 A3M 格式。

**零第三方依赖** — 纯 Python stdlib 实现，无需 CUDA / JAX / AlphaFold。

---

## 快速开始

### 单链 MSA

```python
from biorazer.sequence.protein.analysis.align.query import run_search

files, tmpl_map = run_search(
    ["MTSENLYFQGAMGSMTSENLYFQGAMG"],
    "msa_out/",
    paired=False,
)
print(files)  # ['msa_out/uniref.a3m', ...]
```

### 多链配对 MSA（用于 AF3 / OpenDDE 多聚体）

```python
files, tmpl_map = run_search(
    ["EVQLVESGGGLVQPGGSLRLSCAASGFTFS",   # 重链
     "DIQMTQSPSSLSASVGDRVTITCRASQGIR"],   # 轻链
    "msa_out/",
    paired=True,
)
# 输出包含 unpaired.a3m + paired.a3m
```

### 从 FASTA 文件读取

```python
from biorazer.sequence.protein.analysis.align.query import parse_fasta, run_search

with open("input.fasta") as f:
    seqs = parse_fasta(f.read())

files, tmpl_map = run_search(seqs, "msa_out/", paired=len(seqs) > 1)
```

---

## Python API 参考

### `run_search()`

核心函数，提交任务到 MMseqs2 服务器，轮询等待，下载并解压结果。

```python
run_search(
    seqs: List[str],
    out_dir: str,
    paired: bool,
    use_env: bool = True,
    use_filter: bool = True,
    host: str = "https://api.colabfold.com",
    ua: str = "colabfold_msa/2.0 ...",
    pair_strategy: str = "greedy",
) -> Tuple[List[str], Optional[dict]]
```

**参数说明**

| 参数 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `seqs` | `List[str]` | 必填 | 蛋白序列列表（多链时每链一条） |
| `out_dir` | `str` | 必填 | 输出目录（自动创建） |
| `paired` | `bool` | 必填 | 多链时是否做配对搜索 |
| `use_env` | `bool` | `True` | 是否使用环境数据库 (BFD/MGnify) |
| `use_filter` | `bool` | `True` | 是否过滤 MSA 结果 |
| `host` | `str` | `https://api.colabfold.com` | MMseqs2 服务器地址 |
| `ua` | `str` | ... | User-Agent 标识 |
| `pair_strategy` | `str` | `"greedy"` | 配对策略: `"greedy"` (快) / `"complete"` (全) |

**返回值**

- `files`: 生成的 A3M 文件路径列表
  - 单链: `[uniref.a3m]` (+ 可选 `bfd.mgnify30.*.a3m`)
  - 配对: `[pair.a3m]`
- `tmpl_map`: PDB 模板映射 `{序列索引: [pdb_id, ...]}`，无模板则返回 `None`

**mode 对应关系**

| paired | use_filter | use_env | 实际 mode |
|--------|------------|---------|-----------|
| False  | True       | True    | `env` |
| False  | True       | False   | `all` |
| False  | False      | True    | `env-nofilter` |
| False  | False      | False   | `nofilter` |
| True   | —          | True    | `pairgreedy-env` / `paircomplete-env` |
| True   | —          | False   | `pairgreedy` / `paircomplete` |

### `parse_fasta()`

解析 FASTA 文本，返回序列列表。

```python
parse_fasta(text: str) -> List[str]
```

- 跳过空行和注释
- 自动大写
- 支持多条序列

### `validate()`

验证序列只含 20 种标准氨基酸字符 (`ACDEFGHIKLMNPQRSTVWY`)。

```python
validate(seqs: List[str]) -> None
```

- 非法字符时打印错误到 stderr 并 `sys.exit(1)`
- 通常在提交前调用

### `merge_a3m()`

合并多个 A3M 文件为一个字符串。

```python
merge_a3m(file_list: List[str]) -> str
```

---

## 完整流程示例：单链 MSA + 可视化

```python
from biorazer.sequence.protein.analysis.align.query import run_search, validate
from biorazer.sequence.protein.scripts.MSA_visualizer import read_a3m, plot_msa

# 1. 定义序列
seq = "MTSENLYFQGAMGSMTSENLYFQGAMG"
validate([seq])

# 2. 生成 MSA
files, _ = run_search([seq], "msa_out/", paired=False)

# 3. 读取生成的 A3M
msa = read_a3m(files[0])

# 4. 可视化 coverage
plt = plot_msa(msa, seq_len_list=[len(seq)])
plt.savefig("coverage.png")
```

---

## 多链流程示例：配对抗体 MSA

```python
from biorazer.sequence.protein.analysis.align.query import run_search, merge_a3m, parse_fasta

heavy = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
light = "DIQMTQSPSSLSASVGDRVTITCRASQGIR"

files, tmpl_map = run_search(
    [heavy, light],
    "ab_msa/",
    paired=True,
    use_env=True,
)

# 多链时默认生成: ab_msa/unpaired/  +  ab_msa/paired/
# 自动合并为 unpaired.a3m + paired.a3m

with open("ab_unpaired.a3m", "w") as f:
    f.write(merge_a3m(files))
```

---

## 模板下载

```python
files, tmpl_map = run_search(
    [seq],
    "msa_with_templates/",
    paired=False,
    use_env=True,
)

if tmpl_map:
    from biorazer.sequence.protein.analysis.align.query.colabfold_api import _write_templates
    tpl_dir = _write_templates(tmpl_map, "msa_with_templates/",
                                host="https://api.colabfold.com",
                                ua="colabfold_msa/2.0")
    print(f"Templates: {tpl_dir}")
```

> `_write_templates()` 是私有函数，当前版本需直接导入。后续版本将开放为公共 API。

---

## 自建 MMseqs2 服务器

如果部署了自己的 MMseqs2 API 服务器，通过 `host` 参数指定：

```python
files, _ = run_search([seq], "msa_out/", paired=False,
                       host="http://localhost:8080")
```

---

## 输出目录结构

### 单链

```
msa_out/
├── uniref.a3m                    # UniRef30 搜索结果
├── bfd.mgnify30.metaeuk30.smag30.a3m  # 环境数据库 (use_env=True)
├── out.tar.gz                    # 原始下载缓存
└── pdb70.m8                      # 模板搜索结果（如有）
```

### 多链配对

```
msa_out/
├── unpaired/                     # unpaired MSA 子目录
│   ├── uniref.a3m
│   └── out.tar.gz
├── paired/                       # paired MSA 子目录
│   └── out.tar.gz
├── unpaired.a3m                  # 合并后的 unpaired
├── paired.a3m                    # 合并后的 paired
└── templates/                    # 模板文件（--templates 时）
```

---

## 速率限制与重试

- 自动检测 `RATELIMIT` / `UNKNOWN` 状态，等待 5-10s 后重新提交
- HTTP 错误（连接超时、5xx）自动指数退避重试（最大 10 次，2s → 60s）
- 任务 3 次超时后放弃，递增 seed 重新提交
- 公共服务器 `api.colabfold.com` 有匿名速率限制，建议请求间隔至少 5s
- 缓存检测：如果 `out.tar.gz` 已存在则跳过提交

---

## 与项目现有 MSA 工具的关系

```
biorazer/sequence/protein/analysis/align/
├── query/               ← 本模块: MSA 生成（调用 ColabFold API）
├── io.py                ← 读取 A3M 文件
├── plot.py              ← 可视化 MSA（coverage + 序列比对图）
├── report.py            ← MSA 分析报告
└── util.py              ← 对齐工具函数
```

其他 MSA 相关脚本：

| 文件 | 作用 |
|------|------|
| `biorazer/sequence/protein/scripts/MSA_visualizer.py` | 读 A3M → 绘制 coverage 图 |
| `biorazer/sequence/protein/scripts/msa_analyzer.py` | 统计每个位置的氨基酸频率 |

---

## 注意

1. **网络要求**：需要访问 `https://api.colabfold.com`（境外服务），建议在校园网 / 代理下使用
2. **无 CUDA 要求**：所有计算在远端完成，本地只需 Python stdlib
3. **A3M 格式**：结果为大写字母（匹配）+ 小写字母（插入）+ `-`（缺失），读取时注意小写 → `-` 转换
4. **序列去重**：内部已做 `dict.fromkeys()` 去重，重复序列只提交一次
5. **配对策略**：`greedy` 快但可能不是最优配对，`complete` 慢但完整
