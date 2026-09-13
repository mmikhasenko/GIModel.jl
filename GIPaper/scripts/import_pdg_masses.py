#!/usr/bin/env python3
"""Refresh the reviewed PDG property selections from an official SQLite snapshot.

Usage: python3 GIPaper/scripts/import_pdg_masses.py /path/to/pdg-2026.0.sqlite
Download: https://pdg.lbl.gov/2026/api/pdg-2026.0.sqlite
No dependency installation or network access is required by this importer.
Assignments and historical inputs are reviewed separately, never overwritten.
"""
import csv
import datetime
import hashlib
import json
from pathlib import Path
import sqlite3
import sys


def main(database):
    directory = Path(__file__).resolve().parents[1] / 'data' / 'mass_inputs'
    db = sqlite3.connect(f'file:{database.resolve()}?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    info = dict(db.execute('select name,value from pdginfo'))
    if info['edition'] != '2026' or info['status'] != 'production':
        raise ValueError('Expected production PDG 2026 database; review selections for a new edition')
    rows = []
    for selection in csv.DictReader((directory / 'selection.csv').open()):
        pid = selection['pdg_id']
        values = list(db.execute('''select d.*, i.description from pdgdata d
            join pdgid i on i.pdgid=d.pdgid where d.pdgid=? and d.edition='2026'
            and d.in_summary_table=1 and i.data_type='M' ''', (pid,)))
        if len(values) != 1:
            raise ValueError(f'{pid}: ambiguous or absent summary mass')
        r = values[0]
        if r['unit_text'] != 'MeV' or r['value'] is None or r['limit_type'] is not None:
            raise ValueError(f'{pid}: not a finite central mass in MeV')
        rows.append(dict(key=selection['key'], particle=r['description'],
            mass_GeV=r['value']/1000, pdg_id=pid, pdg_display=r['display_value_text'],
            edition='2026', status={'E':'PDG estimate', 'FC':'PDG fit', 'AC':'PDG average', 'V':'PDG quoted value'}.get(r['value_type'], 'PDG input'), note=r['comment'] or '',
            source_url='https://pdglive.lbl.gov/DataBlock.action?node='+pid))
    with (directory / 'pdg-2026.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    (directory / 'source.json').write_text(json.dumps(dict(edition='2026',
        database_url='https://pdg.lbl.gov/2026/api/pdg-2026.0.sqlite',
        sha256=hashlib.sha256(database.read_bytes()).hexdigest(),
        retrieved=str(datetime.date.today()),license=info['license'],citation=info['citation']),indent=2)+'\n')
    print(f'Imported {len(rows)} reviewed PDG mass properties')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
