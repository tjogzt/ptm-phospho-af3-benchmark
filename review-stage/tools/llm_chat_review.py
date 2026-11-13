#!/usr/bin/env python3
"""llm-chat HTTP reviewer transport for the ARIS `auto-review-loop` skill.

中文说明：
    ARIS 的 auto-review-loop 默认通过 `mcp__llm-chat__review` 这个 MCP 服务把评审
    请求发到 OpenAI 兼容的 HTTP 端点。本机没有安装该 MCP 服务，且无法在会话中途注册，
    因此本脚本在本地实现了同一份 HTTP 契约：把评审 prompt 与 `files:` 中列出的
    主证据文件原文，一起 POST 到 OpenAI 兼容的 /chat/completions 端点；随后记录
    服务方回报的真实模型 id，供 ARIS 停止门 `review_gate.py --round-backend llm-chat`
    校验「executor 与 reviewer 跨家族」这一独立性条件。

English:
    The ARIS auto-review-loop normally reaches its reviewer through the
    `mcp__llm-chat__review` MCP server. That server is not installed on this host
    and cannot be registered mid-session, so this script implements the same HTTP
    contract locally: it POSTs the review prompt together with the verbatim
    contents of every `files:` primary artifact to an OpenAI-compatible
    /chat/completions endpoint, then records the provider-reported model id so the
    ARIS stop gate (`review_gate.py --round-backend llm-chat`) can verify
    cross-family independence between executor and reviewer.

Usage:
    LLM_BASE_URL / LLM_API_KEY / LLM_MODEL are read from the environment, or
    passed explicitly with --base-url / --api-key / --model.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time

# ---------------------------------------------------------------------------
# 模型家族判定 / Model-family derivation.
# 与 tools/review_gate.py 的 derive_model_family 保持一致，用于本地自检。
# Kept in sync with tools/review_gate.py so the local self-check matches the gate.
# ---------------------------------------------------------------------------
_FAMILY_PATTERNS = (
    ("openai", r"(^|[^a-z0-9])(gpt|chatgpt|codex|oracle|o1|o3|o4)([^a-z0-9]|$)"),
    ("anthropic", r"(^|[^a-z0-9])(claude|sonnet|opus|haiku|anthropic)([^a-z0-9]|$)"),
    ("google", r"(^|[^a-z0-9])(gemini|google)([^a-z0-9]|$)"),
    ("deepseek", r"(^|[^a-z0-9])(deepseek)[0-9.]*([^a-z0-9]|$)"),
    ("moonshot", r"(^|[^a-z0-9])(kimi|moonshot)[0-9.]*([^a-z0-9]|$)"),
    ("qwen", r"(^|[^a-z0-9])(qwen|tongyi)[0-9.]*([^a-z0-9]|$)"),
    ("zhipu", r"(^|[^a-z0-9])(glm|zhipu)[0-9.]*([^a-z0-9]|$)"),
    ("minimax", r"(^|[^a-z0-9])(minimax|abab)[0-9.]*([^a-z0-9]|$)"),
    ("xiaomi", r"(^|[^a-z0-9])(mimo|xiaomi)[0-9.]*([^a-z0-9]|$)"),
    ("bytedance", r"(^|[^a-z0-9])(doubao|bytedance|volcengine)[0-9.]*([^a-z0-9]|$)"),
    ("xai", r"(^|[^a-z0-9])(grok)[0-9.]*([^a-z0-9]|$)"),
    ("meta", r"(^|[^a-z0-9])(llama)[0-9.]*([^a-z0-9]|$)"),
    ("mistral", r"(^|[^a-z0-9])(mistral|mixtral)[0-9.]*([^a-z0-9]|$)"),
)


def derive_model_family(model: str) -> str:
    """返回模型家族；无法唯一判定时返回 'unknown'。
    Return the model family, or 'unknown' when it cannot be decided uniquely."""
    name = (model or "").strip().lower()
    hits = {fam for fam, pat in _FAMILY_PATTERNS if re.search(pat, name)}
    return next(iter(hits)) if len(hits) == 1 else "unknown"


def _read(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def build_user_message(prompt: str, files: list[str], label_map: dict[str, str]) -> str:
    """拼装 user message：先放 prompt 原文，再逐文件附上原文内容。
    Assemble the user message: the verbatim prompt first, then each artifact."""
    parts = [prompt.rstrip(), "", "=== PRIMARY ARTIFACTS (verbatim; read these yourself) ==="]
    for path in files:
        label = label_map.get(path, path)
        parts.append(f"\n----- BEGIN FILE: {label} -----")
        parts.append(_read(path))
        parts.append(f"----- END FILE: {label} -----")
    return "\n".join(parts)


def call_endpoint(base_url: str, api_key: str, model: str, system: str,
                  user: str, max_tokens: int, timeout: int) -> dict:
    """调用 OpenAI 兼容端点 / Call the OpenAI-compatible endpoint.

    中文说明：本机 Python 的 TLS 握手会被网络中间设备重置（SSLEOFError），而 curl
    可以正常建链，因此这里通过 subprocess 调用 curl 作为传输层。
    English: this host's Python TLS handshake is reset by a middlebox
    (SSLEOFError) while curl connects fine, so curl is used as the transport.
    """
    url = base_url.rstrip("/") + "/chat/completions"
    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        "max_tokens": max_tokens,
        "temperature": 0.2,
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False)
        payload_path = fh.name
    try:
        started = time.time()
        proc = subprocess.run(
            [
                "curl", "-sS", "--fail-with-body", "-m", str(timeout),
                url,
                "-H", "Content-Type: application/json",
                "-H", f"Authorization: Bearer {api_key}",
                "--data-binary", f"@{payload_path}",
            ],
            capture_output=True, text=True,
        )
        elapsed = round(time.time() - started, 2)
    finally:
        os.unlink(payload_path)

    if proc.returncode != 0:
        raise SystemExit(
            f"REVIEW_UNAVAILABLE: curl exit {proc.returncode} for {url}: "
            f"{(proc.stderr or proc.stdout)[:2000]}"
        )
    try:
        body = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"REVIEW_UNAVAILABLE: non-JSON response from {url}: {proc.stdout[:2000]}"
        ) from exc
    if isinstance(body, dict) and body.get("error"):
        raise SystemExit(f"REVIEW_UNAVAILABLE: provider error: {json.dumps(body['error'])[:2000]}")
    return {"raw": body, "status": 200, "elapsed_seconds": elapsed}


def main() -> int:
    ap = argparse.ArgumentParser(description="llm-chat HTTP reviewer transport")
    ap.add_argument("--prompt-file", required=True)
    ap.add_argument("--files", nargs="*", default=[])
    ap.add_argument("--response-out", required=True)
    ap.add_argument("--meta-out", required=True)
    ap.add_argument("--base-url", default=os.environ.get("LLM_BASE_URL", ""))
    ap.add_argument("--api-key", default=os.environ.get("LLM_API_KEY", ""))
    ap.add_argument("--model", default=os.environ.get("LLM_MODEL", ""))
    ap.add_argument("--system", default=(
        "You are a senior machine-learning and computational-biology reviewer "
        "for a top venue (NeurIPS / ICML / Bioinformatics). Be adversarial and "
        "verify claims against the supplied artifacts rather than trusting the "
        "author's framing. Answer in English."
    ))
    ap.add_argument("--executor-model", default="")
    ap.add_argument("--max-tokens", type=int, default=12000)
    ap.add_argument("--timeout", type=int, default=900)
    args = ap.parse_args()

    missing = [n for n, v in (("LLM_BASE_URL", args.base_url), ("LLM_API_KEY", args.api_key),
                              ("LLM_MODEL", args.model)) if not v]
    if missing:
        raise SystemExit("REVIEW_UNAVAILABLE: missing config: " + ", ".join(missing))

    prompt = _read(args.prompt_file)
    label_map = {p: os.path.relpath(p) for p in args.files}
    user = build_user_message(prompt, args.files, label_map)

    result = call_endpoint(args.base_url, args.api_key, args.model, args.system,
                           user, args.max_tokens, args.timeout)
    body = result["raw"]

    try:
        content = body["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError):
        raise SystemExit(
            "REVIEW_UNAVAILABLE: malformed completion payload: "
            + json.dumps(body)[:2000]
        )
    reported_model = body.get("model") or args.model

    with open(args.response_out, "w", encoding="utf-8") as fh:
        fh.write(content)

    reviewer_family = derive_model_family(reported_model)
    executor_family = derive_model_family(args.executor_model) if args.executor_model else "unknown"
    independence_ok = (
        reviewer_family != "unknown"
        and executor_family != "unknown"
        and reviewer_family != executor_family
    )
    meta = {
        "transport": "llm-chat (local HTTP client; MCP server unavailable on host)",
        "endpoint": args.base_url,
        "requested_model": args.model,
        "reviewer_model": reported_model,
        "reviewer_model_source": "backend-reported",
        "reviewer_family": reviewer_family,
        "executor_model": args.executor_model,
        "executor_model_source": "caller-declared",
        "executor_family": executor_family,
        "family_relation": (
            "unknown"
            if "unknown" in (reviewer_family, executor_family)
            else ("same" if reviewer_family == executor_family else "different")
        ),
        "independence_verified": independence_ok,
        "n_files": len(args.files),
        "user_message_chars": len(user),
        "http_status": result["status"],
        "elapsed_seconds": result["elapsed_seconds"],
        "usage": body.get("usage", {}),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    with open(args.meta_out, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2, ensure_ascii=False)

    print(json.dumps({k: meta[k] for k in (
        "reviewer_model", "reviewer_family", "executor_family",
        "family_relation", "independence_verified", "usage", "elapsed_seconds")},
        ensure_ascii=False))
    if not independence_ok:
        print("WARN: cross-family independence NOT verified", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
