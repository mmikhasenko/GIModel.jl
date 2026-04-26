#!/usr/bin/env python3
"""Require every `data/reference_spectrum_*.csv` to have columns the Julia loader uses."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

REQUIRED = {
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


def main() -> int:
    files = sorted(DATA.glob("reference_spectrum_*.csv"))
    if not files:
        print("No reference_spectrum_*.csv under data/", file=sys.stderr)
        return 1
    bad = []
    for path in files:
        with path.open(newline="", encoding="utf-8") as f:
            row0 = next(csv.reader(f), None)
        if not row0:
            bad.append(f"{path.name}: empty file")
            continue
        names = {c.strip() for c in row0}
        miss = sorted(REQUIRED - names)
        if miss:
            bad.append(f"{path.name}: missing columns {miss}")
    if bad:
        for b in bad:
            print(b, file=sys.stderr)
        return 1
    print(f"reference_spectrum: OK ({len(files)} files, required columns present).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
