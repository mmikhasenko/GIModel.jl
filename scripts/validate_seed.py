#!/usr/bin/env python3
import csv
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "data" / "seed" / "godfrey_isgur_seed_masses.csv"

REQUIRED = {
    "sector",
    "state_label",
    "assignment",
    "mass_MeV",
    "source_status",
    "source_short",
    "source_url",
}


def main() -> int:
    if not SEED.exists():
        print(f"missing seed table: {SEED}", file=sys.stderr)
        return 1

    with SEED.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing_columns = REQUIRED.difference(reader.fieldnames or [])
        if missing_columns:
            print("missing required columns: " + ", ".join(sorted(missing_columns)))
            return 1

        rows = list(reader)

    errors = []
    keys = Counter()
    for index, row in enumerate(rows, start=2):
        for column in ("sector", "assignment", "source_status"):
            if not row.get(column, "").strip():
                errors.append(f"row {index}: missing {column}")

        try:
            float(row.get("mass_MeV", ""))
        except ValueError:
            errors.append(f"row {index}: mass_MeV is not numeric")

        key = (
            row.get("sector", "").strip(),
            row.get("state_label", "").strip(),
            row.get("assignment", "").strip(),
            row.get("mass_MeV", "").strip(),
        )
        keys[key] += 1

    duplicates = {key: count for key, count in keys.items() if count > 1}

    print(f"loaded rows: {len(rows)}")
    print(f"duplicate exact rows: {len(duplicates)}")

    if duplicates:
        for key, count in sorted(duplicates.items()):
            print(f"duplicate x{count}: {key}")

    if errors:
        print("schema errors:")
        for error in errors:
            print(f"- {error}")
        return 1

    print("seed schema validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

