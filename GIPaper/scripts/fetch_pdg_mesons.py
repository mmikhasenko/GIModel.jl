#!/usr/bin/env python3
"""Download the official PDG SQLite API snapshot and export conventional mesons.

The output contains one row per isospin/charge member, including the implied
partners of I=1/2 doublets and I=1 triplets.  ``observed`` is deliberately
conservative: it is true only when the PDG records contain charge-specific
measurement metadata or an unambiguous exclusive decay with that total charge.
When a family has several summary mass properties, an explicitly named
``... MASS`` property is preferred, then the first stable property ID.

Usage:
    python3 GIPaper/scripts/fetch_pdg_mesons.py
    python3 GIPaper/scripts/fetch_pdg_mesons.py --database /path/to/pdg-2026.0.sqlite
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
from pathlib import Path
import re
import sqlite3
import tempfile
import urllib.request


EDITION = "2026"
DATABASE_URL = f"https://pdg.lbl.gov/{EDITION}/api/pdg-{EDITION}.0.sqlite"
DEFAULT_OUTPUT = (
    Path(__file__).resolve().parents[1] / "data" / f"pdg_mesons_{EDITION}.csv"
)
SUPPLEMENT = Path(__file__).resolve().parents[1] / "data" / "meson_reference_supplement.csv"

# A few very precisely measured ground-state charge partners have separate PDG
# family identifiers.  They are one isospin multiplet for this visualization.
MERGED_MULTIPLETS = (
    ("S008", "S009"),       # pi+/- and pi0
    ("S010", "S011"),       # K+/- and K0
    ("S031", "S032"),       # D+/- and D0
    ("M061", "M062"),       # D*(2007)0 and D*(2010)+/-
    ("S041", "S042"),       # B+/- and B0
)


def particle_id(db: sqlite3.Connection, pdgid: str) -> int:
    row = db.execute(
        "select min(id) from pdgparticle where pdgid = ?", (pdgid,)
    ).fetchone()
    if row is None or row[0] is None:
        raise ValueError(f"PDG anchor {pdgid!r} is absent")
    return int(row[0])


def class_boundaries(db: sqlite3.Connection) -> list[tuple[int, str]]:
    """Resolve stable PDG family anchors instead of hard-coding row numbers."""
    anchors = [
        ("S008", "light"),
        ("S010", "kaons"),
        ("S031", "D"),
        ("S034", "Ds"),
        ("S041", "B"),
        ("S086", "Bs"),
        ("S091", "Bc"),
        ("M026", "charmonium"),
        ("M171", "bottomonium"),
    ]
    return [(particle_id(db, pdgid), name) for pdgid, name in anchors]


def meson_class(particle_row: sqlite3.Row, boundaries: list[tuple[int, str]]) -> str:
    group = boundaries[0][1]
    for start, candidate in boundaries:
        if particle_row["id"] < start:
            break
        group = candidate
    if group == "light":
        return "light isovector" if particle_row["quantum_i"] == "1" else "light isoscalar"
    return group


def descendants(db: sqlite3.Connection, family_pdgid: str) -> list[str]:
    return [
        row[0]
        for row in db.execute(
            """
            with recursive tree(pdgid) as (
                values (?)
                union all
                select child.pdgid from pdgid child join tree parent
                  on child.parent_pdgid = parent.pdgid
            ) select pdgid from tree
            """,
            (family_pdgid,),
        )
    ]


def normalize_charge_evidence(charge_text: str, isospin: str) -> set[int]:
    charges: set[int] = set()
    text = charge_text.strip()
    if "0" in text:
        charges.add(0)
    if "+" in text:
        charges.add(1)
    if "-" in text:
        charges.add(-1)
    # For an I=1/2 q-qbar family we draw the particle doublet (Q=0,+1), not
    # its charge-conjugate doublet.  Evidence for an antiparticle is evidence
    # for the corresponding particle by charge conjugation.
    if isospin == "1/2" and -1 in charges:
        charges.remove(-1)
        charges.add(1)
    return charges


def item_charge(db: sqlite3.Connection, item_id: int, name: str) -> float | None:
    values = [
        row[0]
        for row in db.execute(
            "select distinct charge from pdgparticle where pdgitem_id = ?",
            (item_id,),
        )
    ]
    if len(values) == 1:
        return float(values[0])
    # Common neutral generic aliases are not particle rows themselves.
    neutral_aliases = {
        "gamma", "eta", "eta'", "omega", "phi", "J/psi", "psi",
        "Upsilon", "pi0", "K(S)0", "K(L)0",
    }
    if name in neutral_aliases or name.endswith("0"):
        return 0.0
    return None


def observed_charge_evidence(
    db: sqlite3.Connection, family_pdgid: str, isospin: str
) -> dict[int, set[str]]:
    evidence: dict[int, set[str]] = {}
    nodes = descendants(db, family_pdgid)
    placeholders = ",".join("?" for _ in nodes)

    for row in db.execute(
        f"select distinct charge from pdgmeasurement where pdgid in ({placeholders}) "
        "and charge is not null and trim(charge) != ''",
        nodes,
    ):
        for charge in normalize_charge_evidence(row[0], isospin):
            evidence.setdefault(charge, set()).add("PDG charge-tagged measurement")

    for row in db.execute(
        f"select description from pdgid where pdgid in ({placeholders})", nodes
    ):
        description = row[0].upper()
        if "NEUTRAL ONLY" in description:
            evidence.setdefault(0, set()).add("PDG neutral-only property")
        if "CHARGED ONLY" in description:
            targets = (1,) if isospin == "1/2" else (-1, 1)
            for charge in targets:
                evidence.setdefault(charge, set()).add("PDG charged-only property")
        if "+-" in description and re.search(r"\+\-.*\b(MASS|WIDTH)\b", description):
            targets = (1,) if isospin == "1/2" else (-1, 1)
            for charge in targets:
                evidence.setdefault(charge, set()).add("charge-specific PDG property")
        if re.search(r"0\s+(MASS|WIDTH)\b", description):
            evidence.setdefault(0, set()).add("charge-specific PDG property")
        if re.search(r"\+\s+(MASS|WIDTH)\b", description):
            evidence.setdefault(1, set()).add("charge-specific PDG property")
        if re.search(r"-\s+(MASS|WIDTH)\b", description):
            target = 1 if isospin == "1/2" else -1
            evidence.setdefault(target, set()).add("charge-specific PDG property")

    decay_ids = [
        row[0]
        for row in db.execute(
            f"select pdgid from pdgid where pdgid in ({placeholders}) "
            "and data_type like 'BFX%'",
            nodes,
        )
    ]
    for decay_id in decay_ids:
        daughters = list(
            db.execute(
                "select pdgitem_id, name, multiplier from pdgdecay "
                "where pdgid = ? and is_outgoing = 1 order by sort",
                (decay_id,),
            )
        )
        charges = [item_charge(db, row[0], row[1]) for row in daughters]
        if daughters and all(charge is not None for charge in charges):
            total = int(round(sum(charge * row[2] for charge, row in zip(charges, daughters))))
            normalized = normalize_charge_evidence(
                "0" if total == 0 else "+" if total == 1 else "-" if total == -1 else "",
                isospin,
            )
            for charge in normalized:
                evidence.setdefault(charge, set()).add(f"exclusive decay {decay_id}")
    return evidence


def expected_charges(isospin: str, particles: list[sqlite3.Row]) -> list[int]:
    if isospin == "1":
        return [-1, 0, 1]
    if isospin == "1/2":
        return [0, 1]
    preferred = [row for row in particles if row["cc_type"] != "A"]
    representative = (preferred or particles)[0]
    return [int(round(representative["charge"]))]


def member_name(particles: list[sqlite3.Row], charge: int, fallback: str) -> str:
    candidates = [
        row for row in particles
        if int(round(row["charge"])) == charge
        and (row["cc_type"] != "A" or charge < 0)
    ]
    if candidates:
        return str(candidates[0]["name"])
    suffix = "0" if charge == 0 else "+" if charge > 0 else "-"
    stem = fallback
    if stem.endswith(("+", "-", "0")):
        stem = stem[:-1]
    return stem + suffix


def summary_masses(db: sqlite3.Connection, family_pdgid: str) -> list[sqlite3.Row]:
    return list(
        db.execute(
            """
            with recursive descendants(pdgid) as (
                values (?)
                union all
                select child.pdgid
                from pdgid child join descendants parent
                  on child.parent_pdgid = parent.pdgid
            )
            select d.pdgid, d.value, d.error_positive, d.error_negative,
                   d.display_value_text, d.value_type, d.comment,
                   node.description
            from descendants tree
            join pdgid node on node.pdgid = tree.pdgid
            join pdgdata d on d.pdgid = node.pdgid
            where node.data_type = 'M' and d.edition = ?
              and d.in_summary_table = 1 and d.unit_text = 'MeV'
              and d.value is not null and d.limit_type is null
            order by (upper(node.description) like '%MASS%') desc, d.pdgid
            """,
            (family_pdgid, EDITION),
        )
    )


def value_or_blank(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{float(value):.9g}"


def export_mesons(database: Path, output: Path) -> int:
    db = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    info = dict(db.execute("select name, value from pdginfo"))
    if info.get("edition") != EDITION or info.get("status") != "production":
        raise ValueError(f"expected the production PDG {EDITION} database")

    bounds = class_boundaries(db)
    first_id = bounds[0][0]
    # The last conventional bottomonium family is Upsilon(11020).  The next
    # PDG block consists of tetraquark candidates and is deliberately outside
    # the ten q-qbar classes requested by this plot.
    last_id = max(
        row[0]
        for row in db.execute("select id from pdgparticle where pdgid = 'M093'")
    )
    particles = list(
        db.execute(
            """
            select p.* from pdgparticle p where p.id between ? and ?
            order by p.id
            """,
            (first_id, last_id),
        )
    )

    families: dict[str, list[sqlite3.Row]] = {}
    for particle in particles:
        families.setdefault(particle["pdgid"], []).append(particle)

    merged_by_family = {
        pdgid: group for group in MERGED_MULTIPLETS for pdgid in group
    }
    units: list[tuple[str, ...]] = []
    seen_families: set[str] = set()
    for family_pdgid in families:
        if family_pdgid in seen_families:
            continue
        unit = tuple(
            pdgid for pdgid in merged_by_family.get(family_pdgid, (family_pdgid,))
            if pdgid in families
        )
        units.append(unit)
        seen_families.update(unit)

    supplements = list(csv.DictReader(SUPPLEMENT.open())) if SUPPLEMENT.is_file() else []
    overrides = {
        (row["family_pdgid"], int(row["charge_e"])): row
        for row in supplements if row["operation"] == "override"
    }

    fieldnames = [
        "class", "particle", "family_pdgid", "multiplet_pdgids", "member_index",
        "mcid", "charge_e",
        "I", "J", "P", "C", "mass_MeV", "error_plus_MeV",
        "error_minus_MeV", "error_stat_MeV", "error_syst_MeV", "error_ext_MeV",
        "selected_mass_pdgid", "n_mass_values",
        "mass_pdgids", "pdg_display",
        "value_types", "observed", "charge_evidence", "edition", "retrieved",
        "data_status", "source_url",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        output_rows = 0
        for unit in units:
            family_pdgid = unit[0]
            members = [particle for pdgid in unit for particle in families[pdgid]]
            particle = members[0]
            isospin = particle["quantum_i"] or ""
            evidence: dict[int, set[str]] = {}
            for pdgid in unit:
                for charge, reasons in observed_charge_evidence(db, pdgid, isospin).items():
                    evidence.setdefault(charge, set()).update(reasons)
            charges = expected_charges(isospin, members)
            for member_index, charge in enumerate(charges, start=1):
                matching = [
                    row for row in members
                    if int(round(row["charge"])) == charge
                    and (row["cc_type"] != "A" or isospin == "1")
                ]
                actual = matching[0] if matching else particle
                source_pdgid = actual["pdgid"] if matching else family_pdgid
                masses = summary_masses(db, source_pdgid)
                primary = masses[0] if masses else None
                observed = charge in evidence or (len(charges) == 1 and primary is not None)
                output_row = {
                    "class": meson_class(particle, bounds),
                    "particle": member_name(members, charge, particle["name"]),
                    "family_pdgid": source_pdgid,
                    "multiplet_pdgids": ";".join(unit),
                    "member_index": member_index,
                    "mcid": "" if actual["mcid"] is None else actual["mcid"],
                    "charge_e": charge,
                    "I": isospin,
                    "J": particle["quantum_j"] or "",
                    "P": particle["quantum_p"] or "",
                    "C": particle["quantum_c"] or "",
                    "mass_MeV": value_or_blank(None if primary is None else primary["value"]),
                    "error_plus_MeV": value_or_blank(
                        None if primary is None else primary["error_positive"]
                    ),
                    "error_minus_MeV": value_or_blank(
                        None if primary is None else primary["error_negative"]
                    ),
                    "error_stat_MeV": "",
                    "error_syst_MeV": "",
                    "error_ext_MeV": "",
                    "selected_mass_pdgid": "" if primary is None else primary["pdgid"],
                    "n_mass_values": len(masses),
                    "mass_pdgids": ";".join(row["pdgid"] for row in masses),
                    "pdg_display": " | ".join(
                        row["display_value_text"] for row in masses
                    ),
                    "value_types": ";".join(sorted({row["value_type"] for row in masses})),
                    "observed": "true" if observed else "false",
                    "charge_evidence": "; ".join(sorted(evidence.get(charge, set()))),
                    "edition": EDITION,
                    "retrieved": dt.date.today().isoformat(),
                    "data_status": "PDG 2026",
                    "source_url": f"https://pdglive.lbl.gov/Particle.action?node={particle['pdgid']}",
                }
                override = overrides.get((source_pdgid, charge))
                if override is not None:
                    for field in fieldnames:
                        if override.get(field, "") != "":
                            output_row[field] = override[field]
                writer.writerow(output_row)
                output_rows += 1
        for addition in (row for row in supplements if row["operation"] == "add"):
            writer.writerow({field: addition.get(field, "") for field in fieldnames})
            output_rows += 1
    db.close()
    return output_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, help="reuse a downloaded PDG SQLite file")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.database is not None:
        database = args.database
        if not database.is_file():
            parser.error(f"database does not exist: {database}")
        count = export_mesons(database, args.output)
    else:
        with tempfile.TemporaryDirectory(prefix="pdg-mesons-") as tempdir:
            database = Path(tempdir) / f"pdg-{EDITION}.0.sqlite"
            print(f"downloading {DATABASE_URL}")
            urllib.request.urlretrieve(DATABASE_URL, database)
            count = export_mesons(database, args.output)
    print(f"wrote {count} PDG meson charge states to {args.output}")


if __name__ == "__main__":
    main()
