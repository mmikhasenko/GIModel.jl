# Autonomous Loop

This repository has a guarded autonomous loop inspired by
`karpathy/autoresearch`: an agent receives a Markdown program, makes one bounded
change, verification runs, and the attempt is either committed or preserved as a
failed attempt.

Unlike metric-only ML autoresearch, this project should ratchet on auditability:
formula provenance, normalization consistency, reproducibility, and clear
classification of residuals.

## Files

- `docs/autonomous_program.md`: the program/instructions for the autonomous
  agent.
- `scripts/autonomous_loop.py`: the outer loop runner.
- `docs/autonomous_runs/`: per-run prompts, outputs, diffs, verification logs,
  and summaries.

## Chat-driven series (IDE)

The same ratchet can be done **without** `autonomous_loop.py`: one bounded
change, run the full verification gate, commit with
`autonomous: iteration N`, repeat from a clean `main`. This matches multi-step
work in the editor; the Python loop is for unattended **terminal** sub-agents
(Codex, Cursor `agent` CLI) with the same program file.

## Quick Start

Run one guarded iteration:

```bash
python3 scripts/autonomous_loop.py --iterations 1
```

Run for several iterations:

```bash
python3 scripts/autonomous_loop.py --iterations 6 --sleep-seconds 20
```

Run for a time budget:

```bash
python3 scripts/autonomous_loop.py --iterations 100 --max-minutes 240
```

The default agent command is:

```bash
codex exec --full-auto --cd . -m gpt-5.2 -
```

It reads the generated prompt from stdin. The explicit model avoids inheriting a
local Codex CLI default that may require a newer CLI. Override it with
`--agent-command` if you want a different coding agent.

### Cursor Agent with Composer 2 (or Composer 2 Fast)

Yes. The same runner works: point `--agent-command` at a helper that reads
**stdin** (the loop always pipes the program text on stdin) and calls the
Cursor **Agent** CLI. Install the CLI from [Cursor’s install
page](https://cursor.com/install) so `agent` (or `cursor-agent`) is on your
`PATH`—it is a separate install from the `cursor` editor launcher.

This repo provides:

```bash
python3 scripts/autonomous_loop.py --iterations 1 \
  --agent-command "bash scripts/cursor_loop_agent.sh"
```

The wrapper [scripts/cursor_loop_agent.sh](../scripts/cursor_loop_agent.sh) runs
the agent in non-interactive **print** mode (`-p`) with `--force` so file/shell
tools can run without manual approval. **Only use that on a trusted check-out**
(the same risk profile as `codex exec --full-auto` in an open terminal).

- Default model: **Composer 2 Fast** (`composer-2-fast`, overridable with
  `CURSOR_LOOP_MODEL`, e.g. `composer-2` for the standard tier).
- Override the binary with `CURSOR_LOOP_AGENT=agent` if you need to be
  explicit.

The Cursor Agent CLI is still evolving; if a flag in the wrapper is wrong for
your version, run `agent --help` and adjust
[scripts/cursor_loop_agent.sh](../scripts/cursor_loop_agent.sh) locally (model
id, or whether the subcommand is `agent` vs `agent chat`).

**Note:** This is the **terminal Agent** (same family as the editor’s agent), not
the in-IDE “Composer 2 in chat” window. You get a comparable model when the CLI
`--model` matches the Composer 2 / Composer 2 Fast ids your account exposes.

## Ratchet Behavior

Each iteration starts from the last accepted git commit.

If the agent makes changes and verification passes, the loop commits them with:

```text
autonomous: iteration N
```

If verification fails, the loop writes the diff to the iteration log and stashes
the changes by default. Failed attempts can be recovered with `git stash list`.

Use `--on-fail keep` to leave a failed worktree dirty for manual inspection.

If the child agent reports a Codex usage limit, the loop marks the iteration
with `usage_limit` and stops the current batch early instead of burning the
remaining iterations.

## Verification Gate

By default, every dirty successful agent iteration must pass:

```bash
python3 scripts/validate_seed.py
julia --project=. test/runtests.jl
julia --project=. scripts/analyze_heavy_quarkonium.jl
julia --project=. scripts/run_all_spectrum_checks.jl
```

Override or shorten the gate with repeated `--verify-command` flags, for
example:

```bash
python3 scripts/autonomous_loop.py \
  --iterations 1 \
  --verify-command "python3 scripts/validate_seed.py" \
  --verify-command "julia --project=. test/runtests.jl"
```

## Safety Notes

- The runner refuses to start from a dirty tree unless `--allow-dirty-start` is
  passed.
- Failed changes are stashed, not deleted.
- The autonomous program tells the agent not to ask the user for input and not
  to silently refit parameters.
- Do not use `--dangerously-bypass-approvals-and-sandbox` unless the whole
  environment is externally sandboxed.
