#!/usr/bin/env python3
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "paper" / "text" / "godfrey_isgur_1985.md"
OUTPUT = ROOT / "paper" / "text" / "godfrey_isgur_1985_polished.md"


HEADER_NOTE = """# Godfrey-Isgur 1985 Polished Working Markdown

Derived from `paper/text/godfrey_isgur_1985.md`.

This file applies conservative OCR cleanup for smoother searching and reading.
It is not an authority for equations, signs, symbols, state labels, or numerical
tables. Check those against `paper/Godfrey-Isgur-1985.pdf`,
`paper/text/pdftotext_layout.txt`, and rendered page images.

"""


TEXT_REPLACEMENTS = {
    "\x02": "",
    "Uniuersity": "University",
    "hive revolutionized": "have revolutionized",
    "urilikely": "unlikely",
    "HAMII. TONIAN": "HAMILTONIAN",
    "HAMII.TONIAN": "HAMILTONIAN",
    "Schrodinger": "Schrodinger",
    "FICx.": "FIG.",
    "FICs.": "FIG.",
    "FICi.": "FIG.",
    "FIGr.": "FIG.",
    "FIGr.": "FIG.",
    "Fig.": "Fig.",
    "Grev": "GeV",
    "GreV": "GeV",
    "Ge V": "GeV",
    "Me V": "MeV",
    "Mev": "MeV",
    "Ph. D.": "Ph.D.",
    "Ph. D": "Ph.D.",
    "m-italic": "m",
    "ensuremath": "",
    "%'hile": "While",
    "%'e": "We",
    "%'idth": "Width",
    "%'idths": "Widths",
    "%'ave": "Wave",
    "%'here": "Where",
    "%'ith": "With",
    "We .substantiate": "We substantiate",
    "nonperturbatiUe": "nonperturbative",
    "@CD": "QCD",
    "FinaHy": "Finally",
    "potent1als": "potentials",
    "approxi-": "approxi-",
    "Coulomb- plus-linear": "Coulomb-plus-linear",
    "one-gluon- exchange": "one-gluon-exchange",
    "linear- confinement": "linear-confinement",
    "momentum- dependent": "momentum-dependent",
    "spin- independent": "spin-independent",
    "spin- orbit": "spin-orbit",
    "color- magnetic": "color-magnetic",
    "Thomas- precession": "Thomas-precession",
    "non- relativistic": "nonrelativistic",
    "nonre- lativistic": "nonrelativistic",
    "relativ- istic": "relativistic",
    "phenomenologi- cal": "phenomenological",
    "experi- ment": "experiment",
    "experi- mental": "experimental",
    "configura- tions": "configurations",
    "modifi- cations": "modifications",
    "approxi- mately": "approximately",
    "simultane- ously": "simultaneously",
    "con- finement": "confinement",
}


REGEX_REPLACEMENTS = [
    (re.compile(r"mesons\s+from the m to the Y"), "mesons from the pi to the upsilon"),
    (re.compile(r"mesons\s+from the pion to the upsilon\s+\"an"), "mesons from the pion to the upsilon can"),
    (re.compile(r"\bTABLE\s+([IVXLCDM]+)\."), r"TABLE \1."),
    (re.compile(r"\bFig\.\s+([0-9]+)"), r"Fig. \1"),
    (re.compile(r"\bFIG\.\s*([0-9]+)"), r"FIG. \1"),
    (re.compile(r"[ \t]+$"), ""),
]


def polish(text: str) -> str:
    for old, new in TEXT_REPLACEMENTS.items():
        text = text.replace(old, new)

    lines = []
    for line in text.splitlines():
        for pattern, replacement in REGEX_REPLACEMENTS:
            line = pattern.sub(replacement, line)
        lines.append(line.rstrip())

    return "\n".join(lines) + "\n"


def main() -> int:
    if not INPUT.exists():
        raise SystemExit(f"missing input markdown: {INPUT}")

    source = INPUT.read_text(errors="replace")
    first_page = source.find("## PDF Page 1")
    if first_page == -1:
        raise SystemExit("input markdown does not contain page markers")

    body = source[first_page:]
    OUTPUT.write_text(HEADER_NOTE + polish(body))
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
