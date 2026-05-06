#!/usr/bin/env python3
"""Validate Phase-2 clean data promotion files."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "data" / "clean"

REQUIRED_MASSES = {
    "clean_id",
    "source",
    "page",
    "sector",
    "quark_content",
    "state_label",
    "n",
    "L",
    "S",
    "J",
    "mass_model_MeV",
    "source_status",
    "confidence",
    "provenance_file",
    "provenance_row",
}

REQUIRED_MIXINGS = {
    "clean_id",
    "source",
    "page",
    "table",
    "state_group",
    "state_name",
    "model",
    "basis",
    "amplitude",
    "confidence",
    "provenance_file",
    "provenance_row",
}


def main() -> int:
    path = CLEAN / "masses.csv"
    if not path.exists():
        print(f"missing {path}", file=sys.stderr)
        return 1

    bad = []
    ids = set()
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED_MASSES - set(reader.fieldnames or [])
        if missing:
            bad.append(f"masses.csv missing columns: {sorted(missing)}")
        count = 0
        for row_num, row in enumerate(reader, start=2):
            count += 1
            cid = row.get("clean_id", "")
            if not cid:
                bad.append(f"row {row_num}: missing clean_id")
            elif cid in ids:
                bad.append(f"row {row_num}: duplicate clean_id {cid}")
            ids.add(cid)
            try:
                mass = float(row.get("mass_model_MeV", ""))
                if mass <= 0:
                    bad.append(f"row {row_num}: non-positive mass_model_MeV")
            except ValueError:
                bad.append(f"row {row_num}: invalid mass_model_MeV")
            prov = ROOT / row.get("provenance_file", "")
            if not prov.exists():
                bad.append(f"row {row_num}: missing provenance file {prov}")
            if row.get("confidence") not in {"high", "medium", "low"}:
                bad.append(f"row {row_num}: unexpected confidence {row.get('confidence')!r}")

    mixings = CLEAN / "mixings.csv"
    if not mixings.exists():
        bad.append("missing data/clean/mixings.csv")
    else:
        mixing_ids = set()
        with mixings.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            missing = REQUIRED_MIXINGS - set(reader.fieldnames or [])
            if missing:
                bad.append(f"mixings.csv missing columns: {sorted(missing)}")
            for row_num, row in enumerate(reader, start=2):
                cid = row.get("clean_id", "")
                if not cid:
                    bad.append(f"mixings row {row_num}: missing clean_id")
                elif cid in mixing_ids:
                    bad.append(f"mixings row {row_num}: duplicate clean_id {cid}")
                mixing_ids.add(cid)
                try:
                    float(row.get("amplitude", ""))
                except ValueError:
                    bad.append(f"mixings row {row_num}: invalid amplitude")
                prov = ROOT / row.get("provenance_file", "")
                if not prov.exists():
                    bad.append(f"mixings row {row_num}: missing provenance file {prov}")
                if row.get("confidence") not in {"high", "medium", "low"}:
                    bad.append(
                        f"mixings row {row_num}: unexpected confidence {row.get('confidence')!r}"
                    )

    if bad:
        for item in bad:
            print(item, file=sys.stderr)
        return 1
    print(f"clean data: OK ({count} promoted mass rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
