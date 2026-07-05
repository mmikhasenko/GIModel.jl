#!/usr/bin/env python3
import csv
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "data" / "raw" / "digitized_figures" / "fig03_isovector_mesons"
CSV_PATH = FIG_DIR / "figure_03_labels.csv"
SVG_PATH = FIG_DIR / "figure_03_replot.svg"


JPC_ORDER = [
    "0-+",
    "1--",
    "1+-",
    "0++",
    "1++",
    "2++",
    "2-+",
    "2--",
    "3--",
    "3+-",
    "3++",
    "4++",
    "4-+",
    "4--",
    "5--",
]

WIDTH = 1900
HEIGHT = 1050
LEFT = 150
RIGHT = 80
TOP = 70
BOTTOM = 115
Y_MIN = 0.0
Y_MAX = 2.55
BAR_HALF_WIDTH = 38


def y_pos(mass: float) -> float:
    plot_height = HEIGHT - TOP - BOTTOM
    return TOP + (Y_MAX - mass) / (Y_MAX - Y_MIN) * plot_height


def x_pos(jpc: str) -> float:
    first_column = LEFT + 72
    plot_width = WIDTH - RIGHT - first_column
    index = JPC_ORDER.index(jpc)
    return first_column + index / (len(JPC_ORDER) - 1) * plot_width


def spectroscopic_label(row: dict[str, str]) -> str:
    return f"{row['n']}^{row['multiplicity']}{row['L']}_{row['J']}({float(row['mass_GeV']):.2f})"


def label_svg(row: dict[str, str]) -> str:
    n = escape(row["n"])
    mult = escape(row["multiplicity"])
    L = escape(row["L"])
    J = escape(row["J"])
    mass = f"{float(row['mass_GeV']):.2f}"
    return (
        f"<tspan>{n}</tspan>"
        f"<tspan baseline-shift=\"super\" font-size=\"70%\">{mult}</tspan>"
        f"<tspan>{L}</tspan>"
        f"<tspan baseline-shift=\"sub\" font-size=\"70%\">{J}</tspan>"
        f"<tspan>({mass})</tspan>"
    )


def label_offsets(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], tuple[float, float]]:
    offsets: dict[tuple[str, str, str], tuple[float, float]] = {}
    by_jpc: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_jpc.setdefault(row["jpc"], []).append(row)

    for jpc, group in by_jpc.items():
        group.sort(key=lambda row: float(row["mass_GeV"]))
        clusters: list[list[dict[str, str]]] = []
        for row in group:
            if not clusters or float(row["mass_GeV"]) - float(clusters[-1][-1]["mass_GeV"]) > 0.06:
                clusters.append([row])
            else:
                clusters[-1].append(row)

        for cluster in clusters:
            for i, row in enumerate(cluster):
                key = (row["jpc"], row["composition_raw"], row["mass_GeV"])
                centered = i - (len(cluster) - 1) / 2
                offsets[key] = (centered * 24, -8 - i * 10)

    manual = {
        ("0-+", "1^1S_0", "0.15"): (22, -10),
        ("0-+", "2^1S_0", "1.30"): (24, -10),
        ("0-+", "3^1S_0", "1.88"): (34, -10),
        ("1--", "1^3S_1", "0.77"): (0, -10),
        ("1--", "2^3S_1", "1.45"): (0, -10),
        ("1--", "1^3D_1", "1.66"): (0, -10),
        ("1--", "3^3S_1", "2.00"): (0, -10),
        ("1--", "2^3D_1", "2.15"): (0, -10),
        ("2++", "1^3P_2", "1.31"): (18, -10),
        ("2++", "2^3P_2", "1.82"): (38, -10),
        ("2++", "1^3F_2", "2.05"): (0, -10),
        ("2-+", "1^1D_2", "1.68"): (0, -10),
        ("2-+", "2^1D_2", "2.13"): (0, -10),
        ("2--", "1^3D_2", "1.70"): (0, -10),
        ("2--", "2^3D_2", "2.15"): (0, -10),
        ("3--", "1^3D_3", "1.68"): (0, -10),
        ("3--", "2^3D_3", "2.13"): (0, -10),
        ("3--", "1^3G_3", "2.37"): (0, -10),
        ("4-+", "1^1G_4", "2.33"): (-24, -10),
        ("5--", "1^3G_5", "2.30"): (18, -10),
    }
    offsets.update(manual)
    return offsets


def make_svg(rows: list[dict[str, str]]) -> str:
    offsets = label_offsets(rows)
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family: Times, "Times New Roman", serif; fill: #000;} .axis{stroke:#000;stroke-width:3;} .tick{stroke:#000;stroke-width:2;} .bar{stroke:#000;stroke-width:4;stroke-linecap:square;} .label{font-size:24px;font-weight:700;} .xtick{font-size:30px;} .ytick{font-size:32px;} .unit{font-size:34px;}</style>',
    ]

    # Axes
    x0 = LEFT
    y0 = y_pos(Y_MIN)
    parts.append(f'<line class="axis" x1="{x0}" y1="{TOP}" x2="{x0}" y2="{y0}" />')
    parts.append(f'<line class="axis" x1="{x0}" y1="{y0}" x2="{WIDTH - RIGHT}" y2="{y0}" />')
    parts.append(f'<path d="M {x0} {TOP - 22} L {x0 - 9} {TOP + 4} L {x0 + 9} {TOP + 4} Z" fill="#000"/>')
    parts.append(f'<text class="unit" x="{x0 - 78}" y="{TOP + 28}">GeV</text>')

    for tick in [0.40, 0.80, 1.20, 1.60, 2.00, 2.40]:
        y = y_pos(tick)
        parts.append(f'<line class="tick" x1="{x0}" y1="{y:.1f}" x2="{x0 + 18}" y2="{y:.1f}" />')
        parts.append(f'<text class="ytick" x="{x0 - 70}" y="{y + 11:.1f}">{tick:.2f}</text>')

    for jpc in JPC_ORDER:
        x = x_pos(jpc)
        parts.append(f'<line class="tick" x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0 - 22}" />')
        parts.append(f'<text class="xtick" x="{x:.1f}" y="{y0 + 46}" text-anchor="middle">{escape(jpc)}</text>')

    # Bars and labels
    for row in sorted(rows, key=lambda r: (JPC_ORDER.index(r["jpc"]), float(r["mass_GeV"]))):
        x = x_pos(row["jpc"])
        y = y_pos(float(row["mass_GeV"]))
        parts.append(f'<line class="bar" x1="{x - BAR_HALF_WIDTH:.1f}" y1="{y:.1f}" x2="{x + BAR_HALF_WIDTH:.1f}" y2="{y:.1f}" />')
        key = (row["jpc"], row["composition_raw"], row["mass_GeV"])
        dx, dy = offsets.get(key, (0, -10))
        parts.append(
            f'<text class="label" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="middle">'
            f"{label_svg(row)}</text>"
        )

    parts.append(
        f'<text x="{WIDTH / 2:.1f}" y="{HEIGHT - 25}" text-anchor="middle" font-size="24">'
        "Fig. 3 replot from digitized Godfrey-Isgur isovector model labels only"
        "</text>"
    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    with CSV_PATH.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    SVG_PATH.write_text(make_svg(rows))
    print(f"wrote {SVG_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
