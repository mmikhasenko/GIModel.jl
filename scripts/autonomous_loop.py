#!/usr/bin/env python3
"""Run a guarded autonomous agent loop for this repository.

The loop is inspired by autoresearch, but the ratchet is project-specific:
iterations are accepted only when the repository verification gate passes.
Failed attempts are preserved in logs and optionally stashed, so the next
iteration starts from the last accepted state.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROGRAM = ROOT / "docs" / "autonomous_program.md"
DEFAULT_LOG_DIR = ROOT / "docs" / "autonomous_runs"
DEFAULT_VERIFY = [
    "python3 scripts/validate_seed.py",
    "python3 scripts/verify_table_ii_toml.py",
    "python3 scripts/validate_reference_spectra.py",
    "julia --project=. test/runtests.jl",
    "julia scripts/compare_central_pointwise_vs_a7a8.jl",
    "julia scripts/analyze_heavy_quarkonium.jl",
    "julia scripts/run_all_spectrum_checks.jl",
]


def run(
    command: str,
    *,
    cwd: Path = ROOT,
    timeout: int | None = None,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        shell=True,
        text=True,
        input=input_text,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
    )


def timestamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d-%H%M%S")


def git_status_porcelain() -> str:
    return run("git status --porcelain").stdout.strip()


def git_current_head() -> str:
    return run("git rev-parse --short HEAD").stdout.strip()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def append_jsonl(path: Path, record: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def compose_prompt(program: str, iteration: int, max_minutes: int | None) -> str:
    budget = (
        f"The outer loop has a total budget of {max_minutes} minutes."
        if max_minutes is not None
        else "The outer loop has no explicit time budget."
    )
    return f"""You are running autonomous iteration {iteration} for this repository.

{budget}

Read and follow the program below. Make one bounded, reviewable improvement.
Do not ask the user for input. If you are blocked, document the blocker in a
repo file or your final answer and leave the tree clean when possible.

Autonomous program:

{program}
"""


def run_agent(args: argparse.Namespace, prompt: str, log_dir: Path) -> int:
    command = args.agent_command
    output_file = log_dir / "agent.out"
    completed = run(
        command,
        timeout=args.agent_timeout_seconds,
        input_text=prompt,
    )
    write_text(output_file, completed.stdout)
    return completed.returncode


def output_contains_usage_limit(log_dir: Path) -> bool:
    output_file = log_dir / "agent.out"
    if not output_file.exists():
        return False
    output = output_file.read_text(encoding="utf-8", errors="replace").lower()
    return "usage limit" in output or "purchase more credits" in output


def run_verify(args: argparse.Namespace, log_dir: Path) -> tuple[bool, list[dict]]:
    results: list[dict] = []
    ok = True
    for index, command in enumerate(args.verify_command, start=1):
        started = time.time()
        completed = run(command, timeout=args.verify_timeout_seconds)
        elapsed = round(time.time() - started, 3)
        log_path = log_dir / f"verify-{index}.out"
        write_text(log_path, completed.stdout)
        result = {
            "command": command,
            "returncode": completed.returncode,
            "elapsed_seconds": elapsed,
            "log": str(log_path.relative_to(ROOT)),
        }
        results.append(result)
        if completed.returncode != 0:
            ok = False
            break
    return ok, results


def save_diff(log_dir: Path) -> None:
    diff = run("git diff --binary").stdout
    staged = run("git diff --cached --binary").stdout
    write_text(log_dir / "worktree.diff", diff)
    write_text(log_dir / "staged.diff", staged)


def accept_iteration(log_dir: Path, iteration: int) -> str:
    run("git add .")
    diff_stat = run("git diff --cached --stat").stdout.strip()
    if not diff_stat:
        return "no_changes"
    write_text(log_dir / "accepted.diffstat", diff_stat + "\n")
    message = f"autonomous: iteration {iteration}"
    completed = run(f"git commit -m {shlex.quote(message)}")
    write_text(log_dir / "commit.out", completed.stdout)
    if completed.returncode != 0:
        return "commit_failed"
    return "accepted"


def reject_iteration(log_dir: Path, iteration: int, mode: str) -> str:
    save_diff(log_dir)
    if mode == "keep":
        return "failed_kept_dirty"
    if mode == "stash":
        message = f"autonomous failed iteration {iteration}"
        completed = run(f"git stash push -u -m {shlex.quote(message)}")
        write_text(log_dir / "stash.out", completed.stdout)
        return "failed_stashed" if completed.returncode == 0 else "stash_failed"
    raise ValueError(f"unknown reject mode: {mode}")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--max-minutes", type=int)
    parser.add_argument("--sleep-seconds", type=int, default=0)
    parser.add_argument("--program", type=Path, default=DEFAULT_PROGRAM)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument(
        "--agent-command",
        default="codex exec --full-auto --cd . -m gpt-5.2 -",
        help="Command that reads the prompt on stdin.",
    )
    parser.add_argument(
        "--agent-timeout-seconds",
        type=int,
        default=1800,
        help="Per-iteration agent timeout.",
    )
    parser.add_argument(
        "--verify-timeout-seconds",
        type=int,
        default=1800,
        help="Timeout for each verification command.",
    )
    parser.add_argument(
        "--verify-command",
        action="append",
        default=None,
        help="Verification command. May be repeated. Defaults to the project gate.",
    )
    parser.add_argument(
        "--on-fail",
        choices=("stash", "keep"),
        default="stash",
        help="What to do with a failed dirty iteration.",
    )
    parser.add_argument(
        "--allow-dirty-start",
        action="store_true",
        help="Allow starting with a dirty tree. Use only for deliberate manual runs.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_command is None:
        args.verify_command = DEFAULT_VERIFY
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    os.chdir(ROOT)

    dirty = git_status_porcelain()
    if dirty and not args.allow_dirty_start:
        print("Refusing to start from a dirty tree. Commit, stash, or pass --allow-dirty-start.")
        print(dirty)
        return 2
    if dirty and args.allow_dirty_start and not args.dry_run:
        print("Refusing a non-dry-run autonomous loop from a dirty tree.")
        print("This protects pre-existing work from being committed or stashed by the ratchet.")
        print(dirty)
        return 2

    program = args.program.read_text(encoding="utf-8")
    run_id = timestamp()
    run_root = args.log_dir / run_id
    run_root.mkdir(parents=True, exist_ok=True)
    write_text(run_root / "program.md", program)

    start = time.time()
    accepted = 0
    for iteration in range(1, args.iterations + 1):
        if args.max_minutes is not None:
            elapsed_minutes = (time.time() - start) / 60.0
            if elapsed_minutes >= args.max_minutes:
                break

        iter_dir = run_root / f"iteration-{iteration:03d}"
        iter_dir.mkdir(parents=True, exist_ok=True)
        head = git_current_head()
        prompt = compose_prompt(program, iteration, args.max_minutes)
        write_text(iter_dir / "prompt.md", prompt)

        if args.dry_run:
            append_jsonl(run_root / "events.jsonl", {
                "iteration": iteration,
                "status": "dry_run",
                "head": head,
            })
            continue

        agent_code = run_agent(args, prompt, iter_dir)
        save_diff(iter_dir)
        dirty_after_agent = git_status_porcelain()

        if agent_code != 0:
            hit_usage_limit = output_contains_usage_limit(iter_dir)
            status = reject_iteration(iter_dir, iteration, args.on_fail) if dirty_after_agent else "agent_failed_clean"
            if hit_usage_limit:
                status = f"{status}_usage_limit"
            verify_results: list[dict] = []
        elif not dirty_after_agent:
            status = "no_changes"
            verify_results = []
        else:
            verify_ok, verify_results = run_verify(args, iter_dir)
            status = accept_iteration(iter_dir, iteration) if verify_ok else reject_iteration(iter_dir, iteration, args.on_fail)
            if status == "accepted":
                accepted += 1

        append_jsonl(run_root / "events.jsonl", {
            "iteration": iteration,
            "head": head,
            "agent_returncode": agent_code,
            "status": status,
            "verify": verify_results,
        })

        print(f"iteration {iteration}: {status}")
        if "usage_limit" in status:
            print("stopping batch early: agent reported a usage limit")
            break
        if args.sleep_seconds and iteration != args.iterations:
            time.sleep(args.sleep_seconds)

    write_text(run_root / "summary.json", json.dumps({
        "run_id": run_id,
        "accepted_iterations": accepted,
        "final_head": git_current_head(),
    }, indent=2, sort_keys=True) + "\n")
    print(f"logs: {run_root.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
