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

    mixings = CLEAN / "mixings.csv"
    if not mixings.exists():
        with mixings.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
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
            )

    print(f"wrote {out} ({len(rows)} rows)")
    print(f"ensured {mixings} schema")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
