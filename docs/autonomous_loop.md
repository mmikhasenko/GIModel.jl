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

## Ratchet Behavior

Each iteration starts from the last accepted git commit.

If the agent makes changes and verification passes, the loop commits them with:

```text
autonomous: iteration N
```

If verification fails, the loop writes the diff to the iteration log and stashes
the changes by default. Failed attempts can be recovered with `git stash list`.

Use `--on-fail keep` to leave a failed worktree dirty for manual inspection.

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
