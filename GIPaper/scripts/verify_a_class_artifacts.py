#!/usr/bin/env python3
"""Verify the committed trace → coefficient → plot ledger, without Julia."""
import csv
import gzip
import math
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "docs" / "a_class_universality"


def rows(name):
    # Prefer committed archives, so a stale uncompressed file cannot mask them.
    path = OUT / f"{name}.csv.gz"
    with gzip.open(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, stem):
    return complex(float(row[stem + "_re"]), float(row[stem + "_im"]))


def close(a, b):
    assert abs(a - b) <= 2e-10 * max(abs(a), abs(b), 1e-8), (a, b)


def main():
    traces = rows("traces")
    comparisons = rows("comparisons")
    memberships = rows("plot_membership")
    certificates = {r["certificate_id"]: r for r in rows("certificates")}
    wave_ids = {r["wave_id"] for r in rows("waves")}
    assert len(traces) == 51840
    raw = {r["trace_id"]: r for r in traces}
    assert len(raw) == len(traces)
    assert len(comparisons) == 9720
    assert len(certificates) == 3240
    assert len(memberships) == 7290
    assert {r["wave_backend"] for r in traces} == {"HO", "SHO", "FD1200", "FD1800"}
    assert {r["transition"] for r in traces} == {f"A{i}" for i in range(1, 10)}
    assert {r["family"] for r in traces} == {"A", "Aprime", "Adoubleprime"}
    for row in traces:
        assert row["parent_mixing"] == row["daughter_mixing"] == row["emitted_mixing"] == "1.0"
        for role in ("parent", "daughter"):
            assert row[role + "_wave_id"].startswith("SHO|") or row[role + "_wave_id"] in wave_ids
        assert row["emitter"] in {"u", "d", "s"}
        assert float(row["q_GeV"]) > 0

    def reconstruct(ids, L):
        samples = [raw[x] for x in ids.split(";")]
        if len(samples) == 1:
            return number(samples[0], "coefficient")
        assert len(samples) == 5
        qa, qb = (float(s["q_GeV"]) for s in samples[-2:])
        a, b = (number(s, "coefficient") / float(s["q_GeV"]) ** L for s in samples[-2:])
        return (qa * qa * b - qb * qb * a) / (qa * qa - qb * qb)

    points = {r["point_id"]: r for r in comparisons}
    assert len(points) == len(comparisons)
    for r in comparisons:
        assert r["certificate_id"] in certificates
        assert r["certified"] == certificates[r["certificate_id"]]["certified"]
        for backend, prefix in (("native", ""), ("sho", "SHO_")):
            for piece in ("g", "h"):
                actual = reconstruct(r[f"{backend}_{piece}_trace_ids"], int(r["relative_L"]))
                close(actual, number(r, prefix + "C" + piece))
        if r["scalar_stable"] == "true":
            cg, ch = number(r, "Cg"), number(r, "Ch")
            sg, sh = number(r, "SHO_Cg"), number(r, "SHO_Ch")
            close(cg / sg, number(r, "Rg"))
            close(ch / sh, number(r, "Rh"))
            close((ch / cg) / (sh / sg), number(r, "U_hg"))
            assert abs(number(r, "U_hg").imag) < 1e-8

    expected = {r["point_id"] for r in comparisons if r["backend"] == "HO"}
    for figure in ("histograms", "transition_panels", "threshold_histogram"):
        entries = [m for m in memberships if m["figure"] == figure]
        selected = expected if figure != "threshold_histogram" else {k for k in expected if float(points[k]["q_GeV"]) == 0}
        assert len(entries) == len(selected)
        assert {m["point_id"] for m in entries} == selected
        for m in entries:
            r = points[m["point_id"]]
            for key in ("native_g_trace_ids", "native_h_trace_ids", "sho_g_trace_ids", "sho_h_trace_ids"):
                assert m[key] == r[key]
    with (OUT / "threshold_bins.csv").open(newline="") as stream:
        edges = {int(r["bin"]): (float(r["left"]), float(r["right"])) for r in csv.DictReader(stream)}
    publication = [m for m in memberships if m["figure"] == "threshold_histogram"]
    assert len(publication) == 810
    for m in publication:
        if m["scalar_stable"] == "true":
            left, right = edges[int(m["bin"])]
            assert left <= float(m["x"]) <= right
    print(f"Verified {len(traces)} raw pieces, {len(comparisons)} comparisons, "
          f"{len(certificates)} certificates and {len(memberships)} plot contributions.")


if __name__ == "__main__":
    main()
