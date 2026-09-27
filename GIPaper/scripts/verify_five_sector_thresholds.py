#!/usr/bin/env python3
"""Independent reconstruction from the earlier Python flavor census."""
from pathlib import Path
from collections import defaultdict,Counter
import gzip,json,tomllib
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/pseudoscalar_census'
old=json.loads(gzip.decompress((OUT/'census.json.gz').read_bytes()))
new=tomllib.loads(gzip.decompress((OUT/'five_sector_thresholds.toml.gz').read_bytes()).decode())
mt=tomllib.loads((OUT/'five_sector_masses_refined.toml').read_text())
masses={(r['sector'],(r['n'],'SPDFG'.index(r['L']),(r['multiplicity']-1)//2,r['J'])):r['unmixed_GeV'] for r in mt['states']}
import csv
pdg={r['key']:float(r['mass_GeV']) for r in csv.DictReader((ROOT/'data/mass_inputs/pdg-2026.csv').open())}
alias={'pi0':'pi0','pi+':'pi+','pi-':'pi+','K+':'K+','K-':'K+','K0':'K0','K0bar':'K0','eta':'eta','eta-prime':'eta_prime'}
mp={'pi+':'pi','pi0':'pi','pi-':'pi','eta':'eta','eta-prime':'eta-prime','K+':'K','K0':'K','K-':'Kbar','K0bar':'Kbar',
    'cubar':'D','cdbar':'D','csbar':'Ds','ucbar':'D_anti','dcbar':'D_anti','scbar':'Ds_anti',
    'bubar':'B','bdbar':'B','bsbar':'Bs','ubbar':'B_anti','dbbar':'B_anti','sbbar':'Bs_anti'}
cc={'pi+':'pi-','pi-':'pi+','pi0':'pi0','eta':'eta','eta-prime':'eta-prime','K+':'K-','K-':'K+','K0':'K0bar','K0bar':'K0'}
for h in 'cb':
 for q in 'uds':cc[h+q+'bar']=q+h+'bar';cc[q+h+'bar']=h+q+'bar'
def spec(t):return f'{t[0]}^{2*t[2]+1}{"SPDFG"[t[1]]}_{t[3]}'
def key(ti,tf,p,d,e,iso=False,conjugate=False):
 names=[cc[x] if conjugate else x for x in (p,d,e)]
 if iso:names=[mp[x] for x in names]
 if tf==(1,0,0,0):names[1:]=sorted(names[1:])
 return '|'.join([spec(ti),spec(tf)]+names)
def mass(t,name):
 if t==(1,0,0,0) and name in alias:return pdg[alias[name]]
 if name.startswith('pi') or name=='eta':sec='nn'
 elif name=='eta-prime':sec='ss'
 elif name.startswith('K'):sec='K'
 else:sec=mp[name].replace('_anti','')
 return masses[sec,t]
byid={r['id']:r for r in old['catalog']}
expected={}
for r in old['cases']['ideal/nonet']['records']:
 t=byid[r['transition']];ti,tf=tuple(t['initial']),tuple(t['final'])
 p,d,e=r['parent'],r['daughter'],r['emitted'];k=key(ti,tf,p,d,e)
 delta=mass(ti,p)-mass(tf,d)-pdg[alias[e]]
 iso=key(ti,tf,p,d,e,True);isoC=min(iso,key(ti,tf,p,d,e,True,True))
 if k not in expected:expected[k]=dict(sector=r['sector'],delta=delta,waves=set(),iso=iso,isoC=isoC)
 assert abs(expected[k]['delta']-delta)<1e-12
 expected[k]['waves'].add(t['ell'])
actual={r['key']:r for r in new['channels']}
assert actual.keys()==expected.keys()
for k,e in expected.items():
 r=actual[k]
 assert abs(r['excess_GeV']-e['delta'])<1e-12
 assert r['open']==(e['delta']>0)
 assert set(r['waves'])==e['waves']
 assert r['isospin_key']==e['iso'] and r['isospin_C_key']==e['isoC']
# Expected D doublet and Ds singlet really are different groups.
assert mp['cubar']==mp['cdbar'] and mp['csbar']!=mp['cdbar']
for s in new['summary']:
 allrows=[r for r in expected.values() if r['sector']==s['sector']]
 rows=[r for r in allrows if r['delta']>0]
 assert len(allrows)==s['before_channels'] and len(rows)==s['open_channels']
 assert sum(len(r['waves']) for r in rows)==s['open_partial_waves']
 assert len({r['iso'] for r in rows})==s['isospin_channels']
 assert len({r['isoC'] for r in rows})==s['isospin_C_channels']
 assert len({(r['iso'],l) for r in rows for l in r['waves']})==s['isospin_partial_waves']
 assert len({(r['isoC'],l) for r in rows for l in r['waves']})==s['isospin_C_partial_waves']
print(f'Independent Python reconstruction matches all {len(expected)} final-channel keys, mass cuts, partial waves and isospin/C groups.')
print('Open totals:',sum(r['open_channels'] for r in new['summary']),sum(r['isospin_channels'] for r in new['summary']),sum(r['isospin_C_channels'] for r in new['summary']))
