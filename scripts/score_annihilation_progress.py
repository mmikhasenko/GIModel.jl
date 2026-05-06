#!/usr/bin/env python3
"""Generate the step-4 isoscalar annihilation progress scorecard."""

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CLEAN = DATA / "clean"
DOCS = ROOT / "docs"
REPORTS = DOCS / "residual_reports"
SCORECARD = REPORTS / "annihilation_model_scorecard.md"

FORMULA_ITEMS = (
    "Eq. (16)",
    "Eq. (17)",
    "Eq. (18a)",
    "Eq. (18b)",
    "Table III",
)

MODELS = (":none", ":calibrated_p1", ":paper_p1", ":paper_p2")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def formula_points() -> tuple[int, list[str]]:
    text = read_text(DOCS / "formula_map.md")
    points = 0
    notes = []
    for item in FORMULA_ITEMS:
        lines = [line for line in text.splitlines() if item in line]
        audited = any("image-audited" in line.lower() or "pdf/image-audited" in line.lower() for line in lines)
        if audited:
            points += 4
            notes.append(f"{item}: audited")
        else:
            notes.append(f"{item}: pending")
    return points, notes


def clean_target_points() -> tuple[int, list[str]]:
    points = 0
    notes = []
    masses = read_csv(CLEAN / "masses.csv")
    pseudoscalars = [
        row for row in masses
        if row.get("sector") == "isoscalar"
        and row.get("L") == "S"
        and row.get("S") == "0"
        and row.get("J") == "0"
        and row.get("n") in {"1", "2"}
        and row.get("mass_model_MeV")
    ]
    if len(pseudoscalars) >= 4:
        points += 10
        notes.append(f"clean isoscalar pseudoscalar masses: {len(pseudoscalars)} rows")
    else:
        notes.append(f"clean isoscalar pseudoscalar masses: {len(pseudoscalars)} rows")

    mixings = read_csv(CLEAN / "mixings.csv")
    audited_mixings = [
        row for row in mixings
        if row.get("model") in {"P1", "P2"}
        and row.get("confidence") in {"high", "medium"}
        and row.get("provenance_file")
    ]
    if audited_mixings:
        models = sorted({row.get("model") for row in audited_mixings})
        if {"P1", "P2"} <= set(models):
            points += 10
        else:
            points += 5
        notes.append(f"audited Table III mixing rows: {len(audited_mixings)}")
    else:
        notes.append("audited Table III mixing rows: 0")
    return points, notes


def implementation_points() -> tuple[int, list[str]]:
    source = read_text(ROOT / "src" / "sector_comparison.jl")
    tests = read_text(ROOT / "test" / "runtests.jl")
    points = 0
    notes = []

    checks = {
        ":none": ("scheme == :none", "plain = compare("),
        ":calibrated_p1": (":calibrated_p1", "isoscalar_pseudoscalar_annihilation = :calibrated_p1"),
        ":paper_p1": (":paper_p1", "isoscalar_pseudoscalar_annihilation = :paper_p1"),
        ":paper_p2": (":paper_p2", "isoscalar_pseudoscalar_annihilation = :paper_p2"),
    }
    for mode, (dispatch_token, test_token) in checks.items():
        implemented = dispatch_token in source
        tested = test_token in tests
        if implemented and tested:
            points += 5
            notes.append(f"{mode}: implemented/tested")
        elif implemented:
            notes.append(f"{mode}: implemented without regression test")
        else:
            notes.append(f"{mode}: pending")
    return points, notes


def calibrated_spectral_score() -> tuple[int, str, str]:
    text = read_text(REPORTS / "isoscalar_residuals.md")
    residuals = []
    for line in text.splitlines():
        if not line.startswith("| `"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 4:
            continue
        state = cells[0].strip("`")
        if state not in {"1^1S_0", "2^1S_0"}:
            continue
        try:
            residuals.append(abs(float(cells[3].replace("+", ""))))
        except ValueError:
            continue
        if len(residuals) == 4:
            break
    if len(residuals) < 4:
        return 0, "n/a", "isoscalar residual report lacks four pseudoscalar rows"
    mean_abs = sum(residuals) / len(residuals)
    if mean_abs <= 25:
        points = 20
    elif mean_abs <= 50:
        points = 10
    else:
        points = 0
    return points, f"{mean_abs:.1f}", "from docs/residual_reports/isoscalar_residuals.md"


def mixing_score(model: str) -> tuple[int, str, str]:
    mixings = read_csv(CLEAN / "mixings.csv")
    literal = model.replace(":paper_", "").upper()
    audited = [
        row for row in mixings
        if row.get("model") == literal and row.get("confidence") in {"high", "medium"}
    ]
    if not audited:
        return 0, "n/a", "no audited clean Table III amplitudes"
    return 0, "n/a", "audited amplitudes present; eigenvector comparison pending"


def write_report(rows: list[dict[str, str]], notes: list[str]) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    with SCORECARD.open("w", encoding="utf-8") as f:
        f.write("# Annihilation Model Scorecard\n\n")
        f.write("Generated by `python3 scripts/score_annihilation_progress.py`.\n\n")
        f.write(
            "| model | formula points | clean target points | implementation points | "
            "spectral points | mixing points | total points | mean abs mass residual MeV | mixing RMS error |\n"
        )
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            f.write(
                "| {model} | {formula_points} | {clean_target_points} | "
                "{implementation_points} | {spectral_points} | {mixing_points} | "
                "{total_points} | {mean_abs_mass_residual_MeV} | {mixing_rms_error} |\n".format(
                    **row
                )
            )
        f.write("\n## Notes\n\n")
        for note in notes:
            f.write(f"- {note}\n")


def main() -> int:
    formula, formula_notes = formula_points()
    clean, clean_notes = clean_target_points()
    implementation, implementation_notes = implementation_points()
    spectral, mean_abs, spectral_note = calibrated_spectral_score()

    rows = []
    notes = (
        [f"Formula provenance: {formula}/20"]
        + formula_notes
        + [f"Clean targets: {clean}/20"]
        + clean_notes
        + [f"Implementation modes: {implementation}/20"]
        + implementation_notes
        + [f"Calibrated spectral fidelity: {spectral}/20 ({spectral_note})"]
    )

    for model in MODELS:
        model_spectral = spectral if model == ":calibrated_p1" else 0
        model_mean_abs = mean_abs if model == ":calibrated_p1" else "n/a"
        mixing, rms, mixing_note = mixing_score(model)
        total = formula + clean + implementation + model_spectral + mixing
        rows.append({
            "model": model,
            "formula_points": str(formula),
            "clean_target_points": str(clean),
            "implementation_points": str(implementation),
            "spectral_points": str(model_spectral),
            "mixing_points": str(mixing),
            "total_points": str(total),
            "mean_abs_mass_residual_MeV": model_mean_abs,
            "mixing_rms_error": rms,
        })
        notes.append(f"{model} mixing fidelity: {mixing}/20 ({mixing_note})")

    write_report(rows, notes)
    print(f"wrote {SCORECARD.relative_to(ROOT)}")
    print("top score:", max(int(row["total_points"]) for row in rows), "/ 100")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
