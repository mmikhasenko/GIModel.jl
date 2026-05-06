#!/usr/bin/env python3
"""Promote audited working reference spectra into data/clean/masses.csv.

This script is intentionally conservative: it preserves row provenance back to
the top-level reference_spectrum CSVs and carries confidence flags forward
instead of pretending all rows are equally audited.
"""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CLEAN = DATA / "clean"

OUT_FIELDS = [
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


def spin_from_multiplicity(value: str) -> str:
    try:
        mult = int(value)
    except (TypeError, ValueError):
        return ""
    if (mult - 1) % 2 == 0:
        return str((mult - 1) // 2)
    return f"{(mult - 1)}/2"


def parity_charge(jpc: str) -> tuple[str, str]:
    text = (jpc or "").strip()
    match = re.match(r"^\d+([+-])([+-])?$", text)
    if not match:
        return "", ""
    return match.group(1) or "", match.group(2) or ""


def promoted_table_iii_mixings() -> list[dict[str, str]]:
    rows = []
    clean_index = 1
    provenance = "data/raw/digitized_tables/table_iii_isoscalar_mixings/table_iii_visible_rows.provisional.csv"
    for state_group, state_name, model, shift, raw_row, amplitudes in TABLE_III_PSEUDOSCALAR_ROWS:
        for basis, amplitude in amplitudes:
            rows.append(
                {
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
                }
            )
            clean_index += 1
    return rows


def main() -> int:
    CLEAN.mkdir(parents=True, exist_ok=True)
    rows = []
    clean_index = 1
    for path in sorted(DATA.glob("reference_spectrum_*.csv")):
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for raw_index, row in enumerate(reader, start=2):
                jpc = row.get("jp_or_jpc") or row.get("jpc") or ""
                p, c = parity_charge(jpc)
                rows.append(
                    {
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
                    }
                )
                clean_index += 1

    out = CLEAN / "masses.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    mixing_rows = promoted_table_iii_mixings()
    mixings = CLEAN / "mixings.csv"
    with mixings.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MIXING_FIELDS)
        writer.writeheader()
        writer.writerows(mixing_rows)

    print(f"wrote {out} ({len(rows)} rows)")
    print(f"wrote {mixings} ({len(mixing_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
