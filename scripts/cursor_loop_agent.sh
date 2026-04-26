#!/usr/bin/env bash
# Bridge stdin → Cursor Agent CLI for scripts/autonomous_loop.py (same idea as
# `codex exec ... -`). Installs: https://cursor.com/install
#
# Usage (from repo root):
#   python3 scripts/autonomous_loop.py --iterations 1 --agent-command "bash scripts/cursor_loop_agent.sh"
#
# Model override:
#   CURSOR_LOOP_MODEL=composer-2         bash scripts/cursor_loop_agent.sh
#   CURSOR_LOOP_MODEL=composer-2-fast    (default; matches in-product "fast" tier)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN="${CURSOR_LOOP_AGENT:-}"
if [ -z "$BIN" ]; then
  if command -v agent >/dev/null 2>&1; then
    BIN=agent
  elif command -v cursor-agent >/dev/null 2>&1; then
    BIN=cursor-agent
  else
    echo "ERROR: neither 'agent' nor 'cursor-agent' on PATH." >&2
    echo "Install the Cursor Agent CLI: https://cursor.com/install" >&2
    exit 127
  fi
fi

MODEL="${CURSOR_LOOP_MODEL:-composer-2-fast}"
PROMPT=$(cat) || true

# Non-interactive print mode: https://cursor.com/docs/cli/using
# --force: allow tools without prompt (only in a trusted, intentionally automated repo)
exec "$BIN" -p --force --model "$MODEL" --workspace "$ROOT" -- "$PROMPT"
