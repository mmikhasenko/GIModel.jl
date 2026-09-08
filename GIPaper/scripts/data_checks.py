#!/usr/bin/env python3
"""Repository data promotion and progress scorecards.

Package data validation lives in the Julia tests. This script is limited to
write-oriented CSV promotion and report scorecard maintenance.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CLEAN = DATA / "clean"
DOCS = ROOT / "docs"
REPORTS = DOCS / "residual_reports"
MASS_FIELDS = [
    "clean_id",
    "source",
    "page",
    "figure",
    "panel",
    "panel_label",
    "sector",
    "quark_content",
    "jp_or_jpc",
    "state_label",
    "n",
    "L",
    "S",
    "J",
    "P",
    "C",
    "mass_model_MeV",
    "mass_exp_MeV",
    "assignment",
    "source_status",
    "confidence",
    "provenance_file",
    "provenance_row",
    "extraction_method",
    "notes",
]

MIXING_FIELDS = [
    "clean_id",
    "source",
    "page",
    "table",
    "state_group",
    "state_name",
    "model",
    "basis",
    "amplitude",
    "mass_shift_MeV",
    "confidence",
    "provenance_file",
    "provenance_row",
    "notes",
]

TABLE_III_PSEUDOSCALAR_ROWS = [
    ("1 S0", "eta(548)", "P1", 370, 2, [("1 ns", +0.67), ("1 ss", -0.73), ("1 cc", +0.001), ("1 bb", +2e-4), ("2 ns", +0.11), ("2 ss", +0.042), ("2 cc", -5e-4)]),
    ("1 S0", "eta(548)", "P2", 340, 3, [("1 ns", +0.68), ("1 ss", -0.73), ("1 cc", -0.005), ("1 bb", +3e-4), ("2 ns", +0.09), ("2 ss", +0.051), ("2 cc", +0.002)]),
    ("1 S0", "eta_prime(958)", "P1", 810, 4, [("1 ns", +0.58), ("1 ss", +0.62), ("1 cc", +0.004), ("1 bb", +5e-4), ("2 ns", +0.47), ("2 ss", +0.13), ("2 cc", -0.002)]),
    ("1 S0", "eta_prime(958)", "P2", 780, 5, [("1 ns", +0.48), ("1 ss", +0.78), ("1 cc", -0.008), ("1 bb", +9e-4), ("2 ns", +0.34), ("2 ss", +0.14), ("2 cc", +0.004)]),
    ("1 S0", "eta_c(2980)", "P1", "", 6, [("1 ns", -0.008), ("1 ss", -0.005), ("1 cc", +1.000), ("1 bb", +2e-4), ("2 ns", +0.004), ("2 ss", +0.003), ("2 cc", -0.002)]),
    ("1 S0", "eta_c(2980)", "P2", "", 7, [("1 ns", -0.002), ("1 ss", -0.001), ("1 cc", +1.000), ("1 bb", +2e-4), ("2 ns", +0.001), ("2 ss", +0.001), ("2 cc", -0.002)]),
    ("2 S0", "eta_r(?)", "P1", "", 8, [("1 ns", -0.26), ("1 ss", -0.17), ("1 cc", -0.003), ("1 bb", -3e-4), ("2 ns", +0.79), ("2 ss", -0.44), ("2 cc", +0.001)]),
    ("2 S0", "eta_r(?)", "P2", "", 9, [("1 ns", +0.07), ("1 ss", +0.08), ("1 cc", -0.005), ("1 bb", -2e-4), ("2 ns", +0.99), ("2 ss", +0.07), ("2 cc", +0.002)]),
    ("2 S0", "eta_r_prime(?)", "P1", "", 10, [("1 ns", -0.17), ("1 ss", -0.10), ("1 cc", -0.003), ("1 bb", -3e-4), ("2 ns", +0.26), ("2 ss", +0.86), ("2 cc", +0.001)]),
    ("2 S0", "eta_r_prime(?)", "P2", "", 11, [("1 ns", +0.09), ("1 ss", +0.08), ("1 cc", -0.006), ("1 bb", -1e-4), ("2 ns", -0.16), ("2 ss", +0.97), ("2 cc", +0.003)]),
]

TABLE_III_NONPSEUDOSCALAR_ROWS = [
    ("1 3S1", "omega(783)", "", 10, 50, [("1 ns", +0.999), ("1 ss", -0.020), ("1 cc", -0.001), ("1 bb", -1e-4)]),
    ("1 3S1", "phi(1019)", "", 250, 51, [("1 ns", +0.020), ("1 ss", +0.999), ("1 cc", -0.0006), ("1 bb", -7e-5)]),
    ("1 3S1", "J/psi(3097)", "", -40, 52, [("1 ns", +0.0009), ("1 ss", +0.0006), ("1 cc", +1.000), ("1 bb", -3e-5)]),
    ("1 3P2", "f2(1270)", "", 215, 53, [("1 ns", +0.997), ("1 ss", +0.060), ("1 cc", +0.007), ("1 bb", +9e-4)]),
    ("1 3P2", "f2_prime(1515)", "", "", 54, [("1 ns", -0.070), ("1 ss", +0.997), ("1 cc", +0.004), ("1 bb", +4e-4)]),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def spin_from_multiplicity(value: str) -> str:
    try:
        mult = int(value)
    except (TypeError, ValueError):
        return ""
    return str((mult - 1) // 2) if (mult - 1) % 2 == 0 else f"{(mult - 1)}/2"


def parity_charge(jpc: str) -> tuple[str, str]:
    match = re.match(r"^\d+([+-])([+-])?$", (jpc or "").strip())
    return ((match.group(1) or ""), (match.group(2) or "")) if match else ("", "")


def promoted_table_iii_mixings() -> list[dict[str, str]]:
    rows = []
    provenance = "data/raw/digitized_tables/table_iii_isoscalar_mixings/table_iii_visible_rows.provisional.csv"
    clean_index = 1
    for state_group, state_name, model, shift, raw_row, amplitudes in TABLE_III_PSEUDOSCALAR_ROWS:
        for basis, amplitude in amplitudes:
            rows.append({
                "clean_id": f"GI1985-X{clean_index:04d}",
                "source": "Godfrey-Isgur-1985.pdf",
                "page": "11",
                "table": "Table III",
                "state_group": state_group,
                "state_name": state_name,
                "model": model,
                "basis": basis,
                "amplitude": f"{amplitude:.6g}",
                "mass_shift_MeV": str(shift),
                "confidence": "medium",
                "provenance_file": provenance,
                "provenance_row": str(raw_row),
                "notes": "Image-audited visible pseudoscalar row from page-011 markdown.",
            })
            clean_index += 1
    for state_group, state_name, model, shift, raw_row, amplitudes in TABLE_III_NONPSEUDOSCALAR_ROWS:
        for basis, amplitude in amplitudes:
            rows.append({
                "clean_id": f"GI1985-X{clean_index:04d}",
                "source": "Godfrey-Isgur-1985.pdf",
                "page": "11",
                "table": "Table III",
                "state_group": state_group,
                "state_name": state_name,
                "model": model,
                "basis": basis,
                "amplitude": f"{amplitude:.6g}",
                "mass_shift_MeV": str(shift),
                "confidence": "medium",
                "provenance_file": "paper/vision_ocr/pages/page-011.md",
                "provenance_row": str(raw_row),
                "notes": "Image-audited visible non-pseudoscalar row from page-011 markdown; general Eq. (16) reproduction is pending.",
            })
            clean_index += 1
    return rows


def promote_clean(_args: argparse.Namespace) -> int:
    CLEAN.mkdir(parents=True, exist_ok=True)
    rows = []
    clean_index = 1
    for path in sorted(DATA.glob("reference_spectrum_*.csv")):
        with path.open(newline="", encoding="utf-8") as f:
            for raw_index, row in enumerate(csv.DictReader(f), start=2):
                jpc = row.get("jp_or_jpc") or row.get("jpc") or ""
                p, c = parity_charge(jpc)
                rows.append({
                    "clean_id": f"GI1985-M{clean_index:04d}",
                    "source": row.get("source", ""),
                    "page": row.get("page", ""),
                    "figure": row.get("figure", ""),
                    "panel": row.get("panel", ""),
                    "panel_label": row.get("panel_label", ""),
                    "sector": row.get("sector", ""),
                    "quark_content": row.get("quark_content", ""),
                    "jp_or_jpc": jpc,
                    "state_label": row.get("composition_raw", ""),
                    "n": row.get("n", ""),
                    "L": row.get("L", ""),
                    "S": spin_from_multiplicity(row.get("multiplicity", "")),
                    "J": row.get("J", ""),
                    "P": p,
                    "C": c,
                    "mass_model_MeV": f"{1000 * float(row['mass_GeV']):.1f}",
                    "mass_exp_MeV": "",
                    "assignment": row.get("composition_raw", ""),
                    "source_status": row.get("label_type", "model_label"),
                    "confidence": row.get("confidence", ""),
                    "provenance_file": str(path.relative_to(ROOT)),
                    "provenance_row": str(raw_index),
                    "extraction_method": row.get("extraction_method", ""),
                    "notes": row.get("notes", ""),
                })
                clean_index += 1

    with (CLEAN / "masses.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MASS_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    mixings = promoted_table_iii_mixings()
    with (CLEAN / "mixings.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MIXING_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(mixings)
    print(f"wrote data/clean/masses.csv ({len(rows)} rows)")
    print(f"wrote data/clean/mixings.csv ({len(mixings)} rows)")
    return 0


def annihilation_score(_args: argparse.Namespace) -> int:
    formula_items = ("Eq. (16)", "Eq. (17)", "Eq. (18a)", "Eq. (18b)", "Table III")
    formula_text = read_text(DOCS / "formula_map.md")
    formula_points = sum(
        4
        for item in formula_items
        if any(
            ("image-audited" in line.lower() or "implemented" in line.lower())
            for line in formula_text.splitlines()
            if item in line
        )
    )

    masses = read_csv(CLEAN / "masses.csv")
    clean_points = 10 if len([
        row for row in masses
        if row.get("sector") == "isoscalar" and row.get("L") == "S"
        and row.get("S") == "0" and row.get("J") == "0" and row.get("n") in {"1", "2"}
    ]) >= 4 else 0
    mixings = read_csv(CLEAN / "mixings.csv")
    if {"P1", "P2"} <= {row.get("model") for row in mixings if row.get("confidence") in {"high", "medium"}}:
        clean_points += 10

    source = read_text(ROOT / "src" / "comparison.jl")
    tests = "\n".join(
        read_text(path)
        for path in sorted((ROOT / "test").glob("*.jl"))
    )
    modes = {
        ":none": ("scheme in (:none", "plain = compare_reference("),
        ":calibrated_p1": (":calibrated_p1", "isoscalar_pseudoscalar_annihilation = :calibrated_p1"),
        ":paper_p1": ("PaperP1Annihilation", "isoscalar_pseudoscalar_annihilation = :paper_p1"),
        ":paper_p2": ("PaperP2Annihilation", "isoscalar_pseudoscalar_annihilation = :paper_p2"),
    }
    implementation_points = sum(5 for dispatch, test in modes.values() if dispatch in source and test in tests)

    residuals = []
    for line in read_text(REPORTS / "isoscalar_residuals.md").splitlines():
        if not line.startswith("| `"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) >= 4 and cells[0].strip("`") in {"1^1S_0", "2^1S_0"}:
            try:
                residuals.append(abs(float(cells[3].replace("+", ""))))
            except ValueError:
                pass
        if len(residuals) == 4:
            break
    mean_abs = sum(residuals) / len(residuals) if len(residuals) == 4 else None
    spectral_points = 20 if mean_abs is not None and mean_abs <= 25 else 10 if mean_abs is not None and mean_abs <= 50 else 0

    # Table III eigenvector fidelity from the mixing audit: mean amplitude RMS
    # per literal model ("Mean vector RMS for `P1`: `0.180`; ...").
    mixing_rms: dict[str, float] = {}
    literal_mass_residual: dict[str, float] = {}
    audit_path = REPORTS / "table_iii_mixing_audit.md"
    if audit_path.exists():
        for line in read_text(audit_path).splitlines():
            m = re.search(r"Mean vector RMS for `(P[12])`: `([0-9.]+)`", line)
            if m:
                mixing_rms[m.group(1)] = float(m.group(2))
            m = re.search(r"Mean abs mass residual for `(P[12])`: `([0-9.]+) MeV`", line)
            if m:
                literal_mass_residual[m.group(1)] = float(m.group(2))

    def mixing_points_for(rms: float | None) -> int:
        if rms is None:
            return 0
        if rms <= 0.20:
            return 20
        if rms <= 0.35:
            return 10
        return 0

    model_rms = {":paper_p1": mixing_rms.get("P1"), ":paper_p2": mixing_rms.get("P2")}
    model_mass = {":paper_p1": literal_mass_residual.get("P1"), ":paper_p2": literal_mass_residual.get("P2")}
    if mean_abs is not None:
        model_mass[":calibrated_p1"] = mean_abs

    def spectral_points_for(residual: float | None) -> int:
        if residual is None:
            return 0
        if residual <= 25:
            return 20
        if residual <= 50:
            return 10
        return 0

    rows = []
    for model in modes:
        model_spectral = spectral_points_for(model_mass.get(model))
        model_mixing = mixing_points_for(model_rms.get(model))
        total = formula_points + clean_points + implementation_points + model_spectral + model_mixing
        rows.append((model, formula_points, clean_points, implementation_points, model_spectral, model_mixing, total))

    REPORTS.mkdir(parents=True, exist_ok=True)
    with (REPORTS / "annihilation_model_scorecard.md").open("w", encoding="utf-8") as f:
        f.write("# Annihilation Model Scorecard\n\n")
        f.write("Generated by `python3 scripts/data_checks.py score-annihilation`.\n\n")
        f.write("| model | formula points | clean target points | implementation points | spectral points | mixing points | total points | mean abs mass residual MeV | mixing RMS error |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for model, fp, cp, ip, sp, mp, total in rows:
            residual = model_mass.get(model)
            mean = f"{residual:.1f}" if residual is not None else "n/a"
            rms = model_rms.get(model)
            rms_text = f"{rms:.3f}" if rms is not None else "n/a"
            f.write(f"| {model} | {fp} | {cp} | {ip} | {sp} | {mp} | {total} | {mean} | {rms_text} |\n")
        f.write("\n## Missing Points\n\n")
        f.write("- Mixing points come from the mean Table III amplitude RMS in `table_iii_mixing_audit.md` (20 if <= 0.20, 10 if <= 0.35), evaluated on HO-basis wavefunctions with the Φ(0)>0 annihilation phase convention.\n")
        f.write("- Spectral points use the same thresholds for every model (20 if mean abs residual <= 25 MeV, 10 if <= 50): the calibrated control is scored against the digitized Fig. 5 labels, the literal modes against the paper-model mass targets in `table_iii_mixing_audit.md`. Literal-mode mass residuals are reported there directly rather than attributed to an obsolete solver discrepancy.\n")
    print("wrote docs/residual_reports/annihilation_model_scorecard.md")
    print("top score:", max(row[-1] for row in rows), "/ 100")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("promote-clean").set_defaults(func=promote_clean)
    sub.add_parser("score-annihilation").set_defaults(func=annihilation_score)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
