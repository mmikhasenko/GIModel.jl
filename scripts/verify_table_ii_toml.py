#!/usr/bin/env python3
"""Fail if `data/parameters.provisional.toml` disagrees with `data/table_ii_parameters.csv`
for all Table II entries that the solver takes from the digitization (not diagnostic k_*)."""

import sys
import tomllib
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "data" / "table_ii_parameters.csv"
TOML_PATH = ROOT / "data" / "parameters.provisional.toml"

# (parameter_key from CSV) -> (tuple path in TOML, optional scale: MeV->GeV etc.)
MAP = {
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


def get_nested(d, path: tuple) -> float:
    for p in path:
        d = d[p]
    return float(d)


def main() -> int:
    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        rows = {r["parameter_key"]: float(r["value"]) for r in csv.DictReader(f)}

    with TOML_PATH.open("rb") as f:
        toml = tomllib.load(f)

    bad = []
    for key, (tpath, _scale) in MAP.items():
        if key not in rows:
            bad.append(f"missing key in CSV: {key}")
            continue
        v_csv = rows[key] * _scale
        v_t = get_nested(toml, tpath)
        if abs(v_csv - v_t) > 1e-9 * max(1.0, abs(v_t)):
            bad.append(f"{key}: CSV {rows[key]!r} (mapped {v_csv}) vs TOML {tpath} = {v_t!r}")

    if bad:
        print("Table II / TOML mismatch:", file=sys.stderr)
        for b in bad:
            print(" ", b, file=sys.stderr)
        return 1
    print("Table II digitization and parameters.provisional.toml agree on mapped entries.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
