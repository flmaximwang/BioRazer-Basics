# ColabFold MSA 使用指南

调用 ColabFold MMseqs2 公共 API (`https://api.colabfold.com`) 生成蛋白质多序列比对 (MSA)，保存为 A3M 格式。

**零第三方依赖** — 纯 Python stdlib 实现，无需 CUDA / JAX / AlphaFold。

---

## 快速开始

### 单链 MSA

```python
from biorazer.sequence.protein.analysis.align.query import run_search

result = run_search(
    [("my_protein", "MTSENLYFQGAMGSMTSENLYFQGAMG")],
    "msa_out/",
)
# result.per_seq["my_protein"].files   -> ['msa_out/my_protein/uniref.a3m', ...]
# result.per_seq["my_protein"].plots   -> ['msa_out/my_protein/logo.png']
```

### 多链配对 MSA（用于 AF3 / OpenDDE 多聚体）

```python
result = run_search(
    [("heavy", "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"),   # 重链
     ("light", "DIQMTQSPSSLSASVGDRVTITCRASQGIR")],   # 轻链
    "ab_msa/",
    pair_mode="paired",
)
# 输出: ab_msa/heavy/pair.a3m, ab_msa/light/pair.a3m
```

### 从 FASTA 文件读取

```python
from biorazer.sequence.protein.analysis.align.query import parse_fasta, run_search

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

**返回值**

```python
SearchResult(
    per_seq = {
        "my_protein": SeqResult(
            files=["msa_out/my_protein/uniref.a3m", "msa_out/my_protein/merged.a3m"],
            plots=["msa_out/my_protein/logo.png"],
            report=""
        ),
    },
    merged="msa_out/merged.a3m",       # 所有链×数据库合并
    templates={101: ["1abc", ...]},     # PDB 模板映射
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

验证序列只含 20 种标准氨基酸字符 (`ACDEFGHIKLMNPQRSTVWY`)。

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

### 单链

```
msa_out/
├── my_protein/
│   ├── uniref.a3m                    # UniRef30 搜索结果
│   ├── bfd.mgnify30.metaeuk30.smag30.a3m  # 环境数据库 (use_env=True)
│   ├── merged.a3m                    # 该链所有数据库合并
│   ├── logo.png                      # sequence logo (50 AA/行, 10 AA/刻度)
│   ├── coverage.png                  # coverage 热图
│   ├── pdb70.m8                      # 该链的模板搜索结果
│   └── templates/                    # 该链的模板 PDB 文件
├── merged.a3m                        # 所有链×所有数据库合并
└── pdb70.m8                          # 模板搜索结果（如有）
```

### 多链（unpaired）

```
msa_out/
├── chain_A/
│   ├── uniref.a3m
│   ├── bfd.mgnify30.*.a3m
│   ├── merged.a3m
│   ├── logo.png
│   ├── coverage.png
│   ├── pdb70.m8
│   └── templates/
├── chain_B/
│   ├── uniref.a3m
│   ├── bfd.mgnify30.*.a3m
│   ├── merged.a3m
│   ├── logo.png
│   ├── coverage.png
│   ├── pdb70.m8
│   └── templates/
├── merged.a3m
├── templates/                        # 模板文件
└── pdb70.m8
```

---

## 速率限制与重试

- 自动检测 `RATELIMIT` / `UNKNOWN` 状态，等待 5-10s 后重新提交
- HTTP 错误（连接超时、5xx）自动指数退避重试（最大 10 次，2s → 60s）
- 任务 3 次超时后放弃，递增 seed 重新提交
- 公共服务器 `api.colabfold.com` 有匿名速率限制，建议请求间隔至少 5s

---

## 与项目现有 MSA 工具的关系

```
biorazer/sequence/protein/analysis/align/
├── query/               ← 本模块: MSA 生成（调用 ColabFold API）
├── io.py                ← 读取 A3M 文件（A3M2ALIGN → biotite Alignment）
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
2. **无 CUDA 要求**：所有计算在远端完成，本地只需 Python stdlib + matplotlib + biotite
3. **A3M 格式**：结果为大写字母（匹配）+ 小写字母（插入）+ `-`（缺失），读取时注意小写 → `-` 转换
4. **序列名称**：从 FASTA header 第一个空白前提取，裸序列命名为 `"default"`
5. **tar.gz**：解压后自动删除
6. **Sequence logo**：使用 biotite 的 `SequenceProfile` + `plot_sequence_logo` 自动生成
7. **配对策略**：`greedy` 快但可能不是最优配对，`complete` 慢但完整
