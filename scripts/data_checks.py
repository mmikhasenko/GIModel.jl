#!/usr/bin/env python3
"""Repository data promotion, validation, and progress scorecards.

This is the one Python data/check entry point. Physics calculations and report
generation should stay in Julia unless they are pure CSV/provenance bookkeeping.
"""

import argparse
import csv
import re
import sys
import tomllib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CLEAN = DATA / "clean"
DOCS = ROOT / "docs"
REPORTS = DOCS / "residual_reports"

SEED_REQUIRED = {
    "sector",
    "state_label",
    "assignment",
    "mass_MeV",
    "source_status",
    "source_short",
    "source_url",
}

REFERENCE_REQUIRED = {
    "sector",
    "quark_content",
    "composition_raw",
    "n",
    "multiplicity",
    "L",
    "J",
    "mass_GeV",
    "confidence",
}

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

TABLE_II_MAP = {
    "m_ud_avg": (("masses", "m_ud_avg_MeV"), 1.0),
    "m_s": (("masses", "m_s_MeV"), 1.0),
    "m_c": (("masses", "m_c_MeV"), 1.0),
    "m_b": (("masses", "m_b_MeV"), 1.0),
    "b": (("potential", "b_GeV2"), 1.0),
    "Lambda": (("potential", "Lambda_MeV"), 1.0),
    "c": (("potential", "c_MeV"), 1.0),
    "sigma0": (("relativistic_smearing", "sigma0_GeV"), 1.0),
    "s": (("relativistic_smearing", "s"), 1.0),
    "epsilon_c": (("relativistic_factors", "epsilon_c"), 1.0),
    "epsilon_t": (("relativistic_factors", "epsilon_t"), 1.0),
    "epsilon_so_vector": (("relativistic_factors", "epsilon_so_vector"), 1.0),
    "epsilon_so_scalar": (("relativistic_factors", "epsilon_so_scalar"), 1.0),
}

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


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def fail(lines: list[str]) -> int:
    for line in lines:
        print(line, file=sys.stderr)
    return 1


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
                "notes": "Image-audited visible pseudoscalar row from page-011.png; non-pseudoscalar Table III rows remain raw.",
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
        writer = csv.DictWriter(f, fieldnames=MASS_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    mixings = promoted_table_iii_mixings()
    with (CLEAN / "mixings.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MIXING_FIELDS)
        writer.writeheader()
        writer.writerows(mixings)
    print(f"wrote data/clean/masses.csv ({len(rows)} rows)")
    print(f"wrote data/clean/mixings.csv ({len(mixings)} rows)")
    return 0


def validate_seed() -> int:
    path = DATA / "seed" / "godfrey_isgur_seed_masses.csv"
    if not path.exists():
        return fail([f"missing seed table: {path}"])
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = SEED_REQUIRED.difference(reader.fieldnames or [])
        rows = list(reader)
    if missing:
        return fail(["seed missing required columns: " + ", ".join(sorted(missing))])
    errors = []
    keys = Counter()
    for index, row in enumerate(rows, start=2):
        for column in ("sector", "assignment", "source_status"):
            if not row.get(column, "").strip():
                errors.append(f"seed row {index}: missing {column}")
        try:
            float(row.get("mass_MeV", ""))
        except ValueError:
            errors.append(f"seed row {index}: mass_MeV is not numeric")
        keys[(row.get("sector", ""), row.get("state_label", ""), row.get("assignment", ""), row.get("mass_MeV", ""))] += 1
    duplicates = {key: count for key, count in keys.items() if count > 1}
    for key, count in sorted(duplicates.items()):
        errors.append(f"seed duplicate x{count}: {key}")
    if errors:
        return fail(errors)
    print(f"seed: OK ({len(rows)} rows)")
    return 0


def validate_table_ii() -> int:
    rows = {r["parameter_key"]: float(r["value"]) for r in read_csv(DATA / "table_ii_parameters.csv")}
    with (DATA / "parameters.provisional.toml").open("rb") as f:
        toml = tomllib.load(f)
    bad = []
    for key, (path, scale) in TABLE_II_MAP.items():
        if key not in rows:
            bad.append(f"Table II missing CSV key: {key}")
            continue
        value = toml
        for part in path:
            value = value[part]
        csv_value = rows[key] * scale
        if abs(csv_value - float(value)) > 1e-9 * max(1.0, abs(float(value))):
            bad.append(f"{key}: CSV {rows[key]!r} maps to {csv_value}, TOML {path} is {value!r}")
    if bad:
        return fail(bad)
    print("Table II/TOML: OK")
    return 0


def validate_references() -> int:
    files = sorted(DATA.glob("reference_spectrum_*.csv"))
    if not files:
        return fail(["no reference_spectrum_*.csv under data/"])
    bad = []
    for path in files:
        with path.open(newline="", encoding="utf-8") as f:
            header = next(csv.reader(f), None)
        if not header:
            bad.append(f"{path.name}: empty file")
            continue
        missing = sorted(REFERENCE_REQUIRED - {column.strip() for column in header})
        if missing:
            bad.append(f"{path.name}: missing columns {missing}")
    if bad:
        return fail(bad)
    print(f"reference spectra: OK ({len(files)} files)")
    return 0


def validate_clean() -> int:
    bad = []
    mass_rows = read_csv(CLEAN / "masses.csv")
    if not mass_rows:
        bad.append("missing or empty data/clean/masses.csv")
    seen = set()
    for index, row in enumerate(mass_rows, start=2):
        cid = row.get("clean_id", "")
        if not cid or cid in seen:
            bad.append(f"masses row {index}: missing or duplicate clean_id {cid!r}")
        seen.add(cid)
        try:
            if float(row.get("mass_model_MeV", "")) <= 0:
                bad.append(f"masses row {index}: non-positive mass_model_MeV")
        except ValueError:
            bad.append(f"masses row {index}: invalid mass_model_MeV")
        if not (ROOT / row.get("provenance_file", "")).exists():
            bad.append(f"masses row {index}: missing provenance file")
        if row.get("confidence") not in {"high", "medium", "low"}:
            bad.append(f"masses row {index}: unexpected confidence {row.get('confidence')!r}")

    mixing_rows = read_csv(CLEAN / "mixings.csv")
    seen = set()
    for index, row in enumerate(mixing_rows, start=2):
        cid = row.get("clean_id", "")
        if not cid or cid in seen:
            bad.append(f"mixings row {index}: missing or duplicate clean_id {cid!r}")
        seen.add(cid)
        try:
            float(row.get("amplitude", ""))
        except ValueError:
            bad.append(f"mixings row {index}: invalid amplitude")
        if not (ROOT / row.get("provenance_file", "")).exists():
            bad.append(f"mixings row {index}: missing provenance file")
        if row.get("confidence") not in {"high", "medium", "low"}:
            bad.append(f"mixings row {index}: unexpected confidence {row.get('confidence')!r}")
    if bad:
        return fail(bad)
    print(f"clean data: OK ({len(mass_rows)} mass rows, {len(mixing_rows)} mixing rows)")
    return 0


def validate_all(_args: argparse.Namespace) -> int:
    for check in (validate_seed, validate_table_ii, validate_references, validate_clean):
        rc = check()
        if rc != 0:
            return rc
    return 0


def annihilation_score(_args: argparse.Namespace) -> int:
    formula_items = ("Eq. (16)", "Eq. (17)", "Eq. (18a)", "Eq. (18b)", "Table III")
    formula_text = read_text(DOCS / "formula_map.md")
    formula_points = sum(
        4
        for item in formula_items
        if any("image-audited" in line.lower() for line in formula_text.splitlines() if item in line)
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

    source = read_text(ROOT / "src" / "sector_comparison.jl")
    tests = read_text(ROOT / "test" / "runtests.jl")
    modes = {
        ":none": ("scheme == :none", "plain = compare("),
        ":calibrated_p1": (":calibrated_p1", "isoscalar_pseudoscalar_annihilation = :calibrated_p1"),
        ":paper_p1": (":paper_p1", "isoscalar_pseudoscalar_annihilation = :paper_p1"),
        ":paper_p2": (":paper_p2", "isoscalar_pseudoscalar_annihilation = :paper_p2"),
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

    rows = []
    for model in modes:
        model_spectral = spectral_points if model == ":calibrated_p1" else 0
        total = formula_points + clean_points + implementation_points + model_spectral
        rows.append((model, formula_points, clean_points, implementation_points, model_spectral, 0, total))

    REPORTS.mkdir(parents=True, exist_ok=True)
    with (REPORTS / "annihilation_model_scorecard.md").open("w", encoding="utf-8") as f:
        f.write("# Annihilation Model Scorecard\n\n")
        f.write("Generated by `python3 scripts/data_checks.py score-annihilation`.\n\n")
        f.write("| model | formula points | clean target points | implementation points | spectral points | mixing points | total points | mean abs mass residual MeV | mixing RMS error |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for model, fp, cp, ip, sp, mp, total in rows:
            mean = f"{mean_abs:.1f}" if model == ":calibrated_p1" and mean_abs is not None else "n/a"
            f.write(f"| {model} | {fp} | {cp} | {ip} | {sp} | {mp} | {total} | {mean} | n/a |\n")
        f.write("\n## Missing Points\n\n")
        f.write("- Literal `:paper_p1` and `:paper_p2` implementation modes are not implemented yet.\n")
        f.write("- Mixing RMS scoring is blocked on comparing model eigenvectors to `data/clean/mixings.csv`.\n")
    print("wrote docs/residual_reports/annihilation_model_scorecard.md")
    print("top score:", max(row[-1] for row in rows), "/ 100")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("promote-clean").set_defaults(func=promote_clean)
    sub.add_parser("validate").set_defaults(func=validate_all)
    sub.add_parser("score-annihilation").set_defaults(func=annihilation_score)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
