"""colabfold_msa.http - ColabFold MMseqs2 API 的 HTTP 传输层 (纯 stdlib urllib)。

只包含裸 HTTP 请求 (POST / GET / 下载字节流) 与服务器默认值,
不涉及任何 API 语义或重试逻辑 (见 api.py)。
"""

import json
from urllib.request import Request, urlopen

# ── 默认值 ──────────────────────────────────────────
DEFAULT_HOST = "https://api.colabfold.com"
DEFAULT_UA = "colabfold_msa/2.0 (standalone; https://github.com/sokrypton/ColabFold)"


def _post(url: str, data: dict, ua: str) -> dict:
    """POST application/x-www-form-urlencoded，返回 JSON 响应字典。"""
    body = "&".join(f"{k}={v}" for k, v in data.items()).encode()
    req = Request(url, data=body, method="POST",
                  headers={"User-Agent": ua,
                           "Content-Type": "application/x-www-form-urlencoded"})
    with urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _get(url: str, ua: str) -> dict:
    """GET JSON 响应字典。"""
    req = Request(url, headers={"User-Agent": ua})
    with urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _download(url: str, out_path: str, ua: str) -> None:
    """GET 下载字节流到文件。"""
    req = Request(url, headers={"User-Agent": ua})
    with urlopen(req, timeout=300) as resp:
        with open(out_path, "wb") as f:
            while True:
                chunk = resp.read(8192)
                if not chunk:
                    break
                f.write(chunk)