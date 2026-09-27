#!/usr/bin/env python3
"""Strict on-shell mass cut on the previously enumerated charm basis channels."""
from pathlib import Path
import csv,json,tomllib,hashlib
from collections import Counter
from audit_charmed_census import GROUPS,ROUTES,meson,EMIT,WAVE,S0,T1,B1,D3
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/pseudoscalar_census'

def main():
    coarse_source=OUT/'charmed_threshold_masses.toml'
    source=OUT/'charmed_threshold_masses_refined.toml'
    coarse=tomllib.loads(coarse_source.read_text())
    masses=tomllib.loads(source.read_text())
    pdgpath=ROOT/'data/mass_inputs/pdg-2026.csv'
    pdg={r['key']:r for r in csv.DictReader(pdgpath.open())}
    keymap={'pi0':'pi0','pi+':'pi+','pi-':'pi+','K+':'K+','K-':'K+','K0':'K0','K0bar':'K0','eta':'eta','eta-prime':'eta_prime'}
    emass={e:float(pdg[k]['mass_GeV']) for e,k in keymap.items()}
    table={(r['sector'],(r['n'],'SPDFG'.index(r['L']),(r['multiplicity']-1)//2,r['J'])):r for r in masses['states']}
    coarse_table={(r['sector'],(r['n'],'SPDFG'.index(r['L']),(r['multiplicity']-1)//2,r['J'])):r for r in coarse['states']}
    def coarse_mass(level,flavor):return coarse_table['Ds' if 's' in flavor else 'D',level]['unmixed_GeV']
    def mass(level,flavor,mode): return table['Ds' if 's' in flavor else 'D',level][mode+'_GeV']
    records=[]
    for tid,i,f,waves in GROUPS:
        for n,(p,d,e) in enumerate(ROUTES,1):
            mi,mf=mass(i,p,'unmixed'),mass(f,d,'unmixed')
            delta=mi-mf-emass[e]
            mixed_delta=mass(i,p,'mixed')-mass(f,d,'mixed')-emass[e]
            records.append(dict(id=f'{tid}/R{n:02}',template=tid,parent=meson(i,p),daughter=meson(f,d),emitted=EMIT[e],
                parent_flavor=p,daughter_flavor=d,emission_key=e,relative_L=sorted(l for _,l in waves),
                parent_GeV=mi,daughter_GeV=mf,emitted_GeV=emass[e],excess_GeV=delta,open=delta>0,
                mixed_label_excess_GeV=mixed_delta,mixed_label_open=mixed_delta>0,
                coarse_excess_GeV=coarse_mass(i,p)-coarse_mass(f,d)-emass[e]))
    assert len(records)==len({r['id'] for r in records})==448
    for r in records:
        p=r['parent_flavor'];d=r['daughter_flavor']
        ce={'pi+':'pi-','pi-':'pi+','K+':'K-','K-':'K+','K0':'K0bar','K0bar':'K0'}.get(r['emission_key'],r['emission_key'])
        conjugate=next(x for x in records if x['template']==r['template'] and x['parent_flavor']==p[::-1] and x['daughter_flavor']==d[::-1] and x['emission_key']==ce)
        assert r['excess_GeV']==conjugate['excess_GeV']
    def counts(rows):return dict(final_channels=len(rows),partial_waves=sum(len(r['relative_L']) for r in rows))
    opened=[r for r in records if r['open']];closed=[r for r in records if not r['open']]
    summary=dict(before=counts(records),open=counts(opened),closed=counts(closed))
    assert summary['open']['partial_waves']+summary['closed']['partial_waves']==504
    payload=dict(summary=summary,criterion='M_parent > M_daughter + m_emitted (strict, unrounded)',
        charm_mass_prescription='unmixed fixed-sector computed masses; matches the basis census',
        emitted_mass_prescription='stored physical pseudoscalar masses, PDG snapshot; no new data lookup',
        inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,coarse_source,pdgpath)},records=records)
    (OUT/'charmed_thresholds.json').write_text(json.dumps(payload,indent=2)+'\n')
    lines=['# Charmed census after the phase-space cut','',
       'Reproduce: run `julia --project=GIPaper/scripts GIPaper/scripts/charmed_threshold_masses.jl` once without flags and once with `--refined`, then `python3 GIPaper/scripts/filter_charmed_thresholds.py`.','',
       f'**{len(opened)} of 448 basis final channels are open; {summary["open"]["partial_waves"]} of 504 partial-wave entries remain.**',
       f'The cut removes {len(closed)} final channels and {summary["closed"]["partial_waves"]} partial-wave entries. Charge conjugates are included; positive charm alone gives half these totals.','',
       'Criterion: **M(parent) - M(daughter) - m(P) > 0**. Equality is closed. The decision uses unrounded values, not the displayed MeV numbers.',
       'Parent/daughter masses are freshly calculated GI fixed-sector eigenvalues including contact and fine structure, before inter-state mixing. This matches the unmixed-basis 504/448 inventory. They are not the central-only P/D masses used in the earlier ratio calculation.',
       'For emitted π, K, η, η′ this run uses the stored physical masses. Thus this is explicitly a calculated-charm/physical-pseudoscalar threshold test, not an all-model-mass test. No threshold fitting, finite-width smearing, or amplitude calculation is performed.','',
       f'Numerics: `{masses["solver"]}`. Level inventory: {masses["level_inventory"]}. Parameter SHA-256: `{masses["parameter_sha256"]}`.','',
       '## Charm mass inputs (MeV)','',
       '| State label | D before mixing | Ds before mixing | D after mixing | Ds after mixing |','|---|---:|---:|---:|---:|']
    for level in sorted({i for _,i,f,rs in GROUPS}|{f for _,i,f,rs in GROUPS}):
        dn,ds=table['D',level],table['Ds',level]
        lines.append(f'| {dn["label"]} | {1000*dn["unmixed_GeV"]:.6f} | {1000*ds["unmixed_GeV"]:.6f} | {1000*dn["mixed_GeV"]:.6f} | {1000*ds["mixed_GeV"]:.6f} |')
    lines+=['','## Emitted mass inputs (MeV)','','| Meson | Mass |','|---|---:|']
    for k in ('pi0','pi+','K+','K0','eta','eta_prime'):lines.append(f'| {k} | {1000*float(pdg[k]["mass_GeV"]):.6f} |')
    lines+=['','## Per-template selection','','| ID | Transition | Waves | Open final / 28 | Open partial waves |','|---|---|---|---:|---:|']
    for tid,i,f,rs in GROUPS:
        good=[r for r in opened if r['template']==tid]
        lines.append(f'| {tid} | {meson(i,"cu")} → {meson(f,"cu")} + P | {", ".join(WAVE[l] for _,l in sorted(rs,key=lambda x:x[1]))} | {len(good)} | {sum(len(r["relative_L"]) for r in good)} |')
    lines+=['','## Near-threshold decisions (within 25 MeV)','','| Channel | Excess MeV | Status |','|---|---:|---|']
    for r in sorted(records,key=lambda r:abs(r['excess_GeV'])):
        if abs(r['excess_GeV'])<.025:
            lines.append(f'| {r["id"]}: {r["parent"]} → {r["daughter"]} + {r["emitted"]} | {1000*r["excess_GeV"]:+.6f} | {"open" if r["open"] else "closed"} |')
    mixed=[r for r in records if r['mixed_label_open']]
    coarse_open=[r for r in records if r['coarse_excess_GeV']>0]
    changed=[r for r in records if (r['coarse_excess_GeV']>0)!=r['open']]
    lines+=['','## Grid-refinement check','',
       f'The 450-point grid gives {len(coarse_open)} open final channels / {sum(len(r["relative_L"]) for r in coarse_open)} partial-wave entries. The main result uses 900 points at the same rmax=24 GeV⁻¹. {len(changed)} channel classifications change. This is a finite-grid check, not a model-mass uncertainty estimate.']
    for r in changed:
        lines.append(f'- {r["id"]}: excess changes from {1000*r["coarse_excess_GeV"]:+.6f} to {1000*r["excess_GeV"]:+.6f} MeV.')
    lines+=['','## Mixing distinction','',
       f'As a mass-assignment sensitivity check, using final mixed eigenvalues under the same 16 catalog labels leaves **{len(mixed)} final channels / {sum(len(r["relative_L"]) for r in mixed)} partial waves**. This is not a complete physical-channel rotation; in particular the isolated axial daughter needs both physical D₁/D₁′ projections, as explained in the scratch pad.',
       'The main answer above deliberately filters the previously counted 448/504 basis entries. It does not silently replace them by a larger, differently mixed channel inventory.','',
       '## Selected channels: complete open list','',
       'Every retained final channel appears once below. S/D waves share one line. Closed channels remain in `charmed_thresholds.json`, with their masses and negative threshold excess.','',
       '| ID | Open channel | Waves | Excess MeV |','|---|---|---|---:|']
    for r in opened:
        lines.append(f'| {r["id"]} | {r["parent"]} → {r["daughter"]} + {r["emitted"]} | {", ".join(WAVE[l] for l in r["relative_L"])} | {1000*r["excess_GeV"]:+.6f} |')
    (OUT/'charmed_thresholds.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
