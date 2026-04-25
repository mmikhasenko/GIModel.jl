#!/usr/bin/env python3
import csv
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
FIG_ROOT = ROOT / "data" / "raw" / "digitized_figures"

HIDDEN_ORDER = [
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

CONFIGS = {
    "04": {
        "folder": "fig04_strange_mesons",
        "title": "Fig. 4 strange-meson replot from digitized model labels only",
        "order": ["0-", "1-", "0+", "1+", "2+", "2-", "3-", "3+", "4+", "4-", "5-"],
        "y_min": 0.30,
        "y_max": 2.90,
        "ticks": [0.40, 0.80, 1.20, 1.60, 2.00, 2.40, 2.80],
    },
    "05": {
        "folder": "fig05_isoscalar_mesons",
        "title": "Fig. 5 isoscalar-meson replot from digitized model labels only",
        "order": HIDDEN_ORDER,
        "y_min": 0.30,
        "y_max": 2.90,
        "ticks": [0.40, 0.80, 1.20, 1.60, 2.00, 2.40, 2.80],
    },
    "06": {
        "folder": "fig06_charmonia",
        "title": "Fig. 6 charmonium replot from digitized model labels only",
        "order": HIDDEN_ORDER[:12],
        "y_min": 2.75,
        "y_max": 4.90,
        "ticks": [2.80, 3.20, 3.60, 4.00, 4.40, 4.80],
    },
    "07": {
        "folder": "fig07_charmed_mesons",
        "title": "Fig. 7 charmed-meson replots from digitized model labels only",
        "order": ["0-", "1-", "0+", "1+", "2+", "3-", "4+"],
        "y_min": 1.70,
        "y_max": 3.45,
        "ticks": [1.80, 2.20, 2.60, 3.00, 3.40],
    },
    "08": {
        "folder": "fig08_bottomonia",
        "title": "Fig. 8 bottomonium replot from digitized model labels only",
        "order": HIDDEN_ORDER[:12],
        "y_min": 9.20,
        "y_max": 11.25,
        "ticks": [9.20, 9.60, 10.00, 10.40, 10.80, 11.20],
    },
    "09": {
        "folder": "fig09_b_flavored_mesons",
        "title": "Fig. 9 b-flavored meson replots from digitized model labels only",
        "order": ["0-", "1-", "2+", "3-", "4+"],
        "y_min": 5.10,
        "y_max": 7.55,
        "ticks": [5.20, 5.80, 6.20, 6.60, 7.00, 7.40],
    },
}

PANEL_ORDER = {
    "07": ["a", "b"],
    "09": ["a", "b", "c"],
}

WIDTH = 2600
PANEL_WIDTH = 620
HEIGHT = 1050
LEFT = 125
RIGHT = 60
TOP = 70
BOTTOM = 120
BAR_HALF_WIDTH = 34


def spectroscopic_label(row: dict[str, str]) -> str:
    return f"{row['n']}^{row['multiplicity']}{row['L']}_{row['J']}({float(row['mass_GeV']):.2f})"


def label_svg(row: dict[str, str]) -> str:
    n = escape(row["n"])
    mult = escape(row["multiplicity"])
    term_l = escape(row["L"])
    j = escape(row["J"])
    mass = f"{float(row['mass_GeV']):.2f}"
    return (
        f"<tspan>{n}</tspan>"
        f"<tspan baseline-shift=\"super\" font-size=\"70%\">{mult}</tspan>"
        f"<tspan>{term_l}</tspan>"
        f"<tspan baseline-shift=\"sub\" font-size=\"70%\">{j}</tspan>"
        f"<tspan>({mass})</tspan>"
    )


def y_pos(mass: float, cfg: dict[str, object]) -> float:
    plot_height = HEIGHT - TOP - BOTTOM
    return TOP + (float(cfg["y_max"]) - mass) / (float(cfg["y_max"]) - float(cfg["y_min"])) * plot_height


def x_pos(jpc: str, order: list[str], x0: float, x1: float) -> float:
    if len(order) == 1:
        return (x0 + x1) / 2
    return x0 + order.index(jpc) / (len(order) - 1) * (x1 - x0)


def cluster_offsets(rows: list[dict[str, str]]) -> dict[tuple[str, str, str, str], tuple[float, float, str]]:
    offsets = {}
    by_col: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        by_col.setdefault((row["panel"], row["jp_or_jpc"]), []).append(row)

    for (panel, jpc), group in by_col.items():
        group.sort(key=lambda row: float(row["mass_GeV"]))
        clusters: list[list[dict[str, str]]] = []
        for row in group:
            if not clusters or float(row["mass_GeV"]) - float(clusters[-1][-1]["mass_GeV"]) > 0.065:
                clusters.append([row])
            else:
                clusters[-1].append(row)

        for cluster in clusters:
            for i, row in enumerate(cluster):
                centered = i - (len(cluster) - 1) / 2
                anchor = "middle"
                dx = centered * 58
                dy = -8 - i * 22
                if len(cluster) >= 3 and i == len(cluster) - 1:
                    dx += 22
                offsets[(panel, jpc, row["composition_raw"], row["mass_GeV"])] = (dx, dy, anchor)
    return offsets


def make_svg(fig_id: str, cfg: dict[str, object], rows: list[dict[str, str]]) -> str:
    panels = PANEL_ORDER.get(fig_id, ["main"])
    panel_width = PANEL_WIDTH if len(panels) > 1 else WIDTH - LEFT - RIGHT
    width = LEFT + RIGHT + len(panels) * panel_width + max(0, len(panels) - 1) * 60
    order = cfg["order"]
    y0 = y_pos(float(cfg["y_min"]), cfg)
    offsets = cluster_offsets(rows)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{HEIGHT}" viewBox="0 0 {width} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family: Times, "Times New Roman", serif; fill:#000;} .axis{stroke:#000;stroke-width:3;} .tick{stroke:#000;stroke-width:2;} .bar{stroke:#000;stroke-width:4;stroke-linecap:square;} .label{font-size:22px;font-weight:700;} .xtick{font-size:30px;} .ytick{font-size:30px;} .unit{font-size:34px;} .panel{font-size:30px;font-weight:700;} .caption{font-size:24px;}</style>',
    ]

    for panel_index, panel in enumerate(panels):
        panel_x = LEFT + panel_index * (panel_width + 60)
        x_axis0 = panel_x
        x_axis1 = panel_x + panel_width
        first_column = x_axis0 + 72
        last_column = x_axis1 - 38
        panel_rows = [row for row in rows if row["panel"] == panel]

        parts.append(f'<line class="axis" x1="{x_axis0}" y1="{TOP}" x2="{x_axis0}" y2="{y0}" />')
        parts.append(f'<line class="axis" x1="{x_axis0}" y1="{y0}" x2="{x_axis1}" y2="{y0}" />')
        parts.append(f'<path d="M {x_axis0} {TOP - 22} L {x_axis0 - 9} {TOP + 4} L {x_axis0 + 9} {TOP + 4} Z" fill="#000"/>')
        parts.append(f'<text class="unit" x="{x_axis0 - 78}" y="{TOP + 28}">GeV</text>')

        if len(panels) > 1:
            panel_label = panel_rows[0]["panel_label"] if panel_rows else panel
            parts.append(f'<text class="panel" x="{(x_axis0 + x_axis1) / 2:.1f}" y="{TOP - 26}" text-anchor="middle">{escape(panel_label)}</text>')

        for tick in cfg["ticks"]:
            y = y_pos(float(tick), cfg)
            parts.append(f'<line class="tick" x1="{x_axis0}" y1="{y:.1f}" x2="{x_axis0 + 18}" y2="{y:.1f}" />')
            parts.append(f'<text class="ytick" x="{x_axis0 - 72}" y="{y + 10:.1f}">{float(tick):.2f}</text>')

        for jpc in order:
            x = x_pos(jpc, order, first_column, last_column)
            parts.append(f'<line class="tick" x1="{x:.1f}" y1="{y0}" x2="{x:.1f}" y2="{y0 - 22}" />')
            parts.append(f'<text class="xtick" x="{x:.1f}" y="{y0 + 46}" text-anchor="middle">{escape(jpc)}</text>')

        for row in sorted(panel_rows, key=lambda r: (order.index(r["jp_or_jpc"]), float(r["mass_GeV"]))):
            x = x_pos(row["jp_or_jpc"], order, first_column, last_column)
            y = y_pos(float(row["mass_GeV"]), cfg)
            parts.append(f'<line class="bar" x1="{x - BAR_HALF_WIDTH:.1f}" y1="{y:.1f}" x2="{x + BAR_HALF_WIDTH:.1f}" y2="{y:.1f}" />')
            key = (row["panel"], row["jp_or_jpc"], row["composition_raw"], row["mass_GeV"])
            dx, dy, anchor = offsets.get(key, (0, -10, "middle"))
            parts.append(f'<text class="label" x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}">{label_svg(row)}</text>')

    parts.append(f'<text class="caption" x="{width / 2:.1f}" y="{HEIGHT - 26}" text-anchor="middle">{escape(str(cfg["title"]))}</text>')
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def render_png(svg_path: Path) -> None:
    png_path = svg_path.with_suffix(".png")
    try:
        subprocess.run(["rsvg-convert", "-o", str(png_path), str(svg_path)], check=True)
    except FileNotFoundError:
        print(f"warning: rsvg-convert not found; skipped {png_path}")


def main() -> int:
    for fig_id, cfg in CONFIGS.items():
        fig_dir = FIG_ROOT / str(cfg["folder"])
        csv_path = fig_dir / f"figure_{fig_id}_labels.csv"
        svg_path = fig_dir / f"figure_{fig_id}_replot.svg"
        with csv_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        svg_path.write_text(make_svg(fig_id, cfg, rows))
        render_png(svg_path)
        print(f"wrote {svg_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
