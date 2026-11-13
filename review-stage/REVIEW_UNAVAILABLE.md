# AUTO REVIEW LOOP — REVIEW_UNAVAILABLE

- **run_id**: `run_20261111_7f3a9c21`
- **date**: 2026-11-11
- **requested skill**: `auto-review-loop`
- **subject**: PTM 磷酸化结构效应基准（`paper/main.tex` + `results/`）
- **terminal status**: **REVIEW_UNAVAILABLE**（fail-closed，未执行任何评审轮次）

## 1. 后端可用性探针（逐项取证）

| 后端 | 探针命令 / 依赖 | 结果 |
|---|---|---|
| `copilot-native` | 需 `copilot_native_evidence.py` + Copilot root session | ❌ helper 缺失；`copilot` CLI 缺失 |
| `copilot`（兼容驱动） | 需 `copilot --agent` (须支持 `--model/--effort/--allow-tool`) | ❌ `command -v copilot` → MISSING |
| `codex`（MCP） | 需 `mcp__codex__codex` 工具 | ❌ ToolSearch 返回 0 |
| `codex exec`（nightmare 传输） | `codex exec` over CLI | ⚠️ CLI 存在（0.160.0）但**不可用且同族**（见 §2） |
| `manual` | 需 `mcp__manual_review__review` 工具 | ❌ ToolSearch 返回 0 |
| **stop gate** | `review_gate.py`（helper 链解析） | ❌ 三处候选路径均未解析到，`shared-references/` 目录不存在 |

> 技能目录 `/Users/taozhu/.trae/skills/auto-review-loop/` 内**仅有 `SKILL.md`**，其引用的 `../shared-references/*.md`、`tools/review_gate.py`、`save_trace.sh`、`copilot_native_evidence.py` 均不存在。项目内无 `.aris/`，无 `~/.aris/repo`。

## 2. `codex exec` 探针的完整证据

```
$ codex exec "Reply with exactly the token OK ..." --skip-git-repo-check
ERROR: unexpected status 502 Bad Gateway: CC Switch local proxy failed while handling
Codex endpoint /responses. Provider: default; model: deepseek-v4-pro;
cause: 转发失败: 连接失败: error sending request for url (http://127.0.0.1:8000/v1/responses),
url: http://127.0.0.1:15721/v1/responses
```

```
$ lsof -nP -iTCP:8000 -sTCP:LISTEN     → 无监听进程
$ lsof -nP -iTCP:15721 -sTCP:LISTEN    → cc-switch (PID 741) 在监听
$ grep model ~/.codex/config.toml      → model_provider = "deepseek" / model = "deepseek-v4-pro"
```

**双重不合格**：
1. **链路故障**：`cc-switch` 代理的上游 `127.0.0.1:8000` 无服务 → 502，评审请求无法送达。
2. **同族违规**：即便链路恢复，其配置模型 `deepseek-v4-pro` 与执行模型 `DeepSeek-V4.1-Flash` **同属 DeepSeek 家族**，不满足技能强制的「cross-family / different-family」独立评审要求。同族 verdict 不得作为判定依据。

## 3. 根因

**环境缺失，非稿件问题。** 本机安装了 ARIS 技能的 `SKILL.md`，但没有配套的 ARIS toolchain（`review_gate.py`、`save_trace.sh`、`copilot_native_evidence.py`、`shared-references/`），也没有任何**异族**评审模型通道（`copilot` CLI 缺失；Codex 被指向同族 DeepSeek；manual MCP 未注册）。

按技能 fail-closed 条款：
> "If the helper is missing, native acceptance is unavailable; use Codex only if that external backend is positively available, otherwise emit `REVIEW_UNAVAILABLE`."
> "If the helper cannot be resolved or executed, emit `REVIEW_UNAVAILABLE`; do not improvise a transition."

故终止为 `REVIEW_UNAVAILABLE`，未创建 `REVIEW_STATE.json`、未写入 `ACQUITTAL_LOG.jsonl`、未调用任何替代模型冒充评审。

## 4. 恢复路径（任选其一即可重新启用本循环）

- **A. 装齐 ARIS toolchain**：将 `tools/review_gate.py`、`scripts/save_trace.sh`、`scripts/copilot_native_evidence.py`、`shared-references/` 安装到 `/Users/taozhu/.trae/skills/`（或设 `ARIS_REPO` 指向 ARIS 仓库根）。
- **B. 提供异族评审模型**：`copilot` CLI 缺失，可改为指定 `— reviewer: manual`，但需注册 `mcp__manual_review__review` MCP 服务，并保证其模型来自**非 DeepSeek** 家族（OpenAI / Anthropic / Google / Moonshot / Qwen）。
- **C. 修复 Codex 通道 + 换异族模型**：启动 `127.0.0.1:8000` 上游服务，并在 `~/.codex/config.toml` 将 `model` 指向异族模型（如 `gpt-*`），再以 nightmare 传输（`codex exec`）运行。

## 5. 未做的事（诚实声明）

- 未以同模型冒充独立评审。
- 未伪造 Score / Verdict。
- 未创建 `CLAIMS_FROM_RESULTS.md`。
