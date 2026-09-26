#!/usr/bin/env python3
"""Independent, explicit reconstruction of the charm subset of the census."""
from pathlib import Path
from collections import Counter
import gzip
import json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/pseudoscalar_census/charmed_scratchpad.md'
# Explicit flavors, independent of the census's route-generation algorithm.
ROUTES=[
 ('cu','cu','pi0'),('cu','cu','eta'),('cu','cu','eta-prime'),('cu','cd','pi-'),('cu','cs','K-'),
 ('cd','cu','pi+'),('cd','cd','pi0'),('cd','cd','eta'),('cd','cd','eta-prime'),('cd','cs','K0bar'),
 ('cs','cu','K+'),('cs','cd','K0'),('cs','cs','eta'),('cs','cs','eta-prime'),
 ('uc','uc','pi0'),('uc','uc','eta'),('uc','uc','eta-prime'),('uc','dc','pi+'),('uc','sc','K+'),
 ('dc','uc','pi-'),('dc','dc','pi0'),('dc','dc','eta'),('dc','dc','eta-prime'),('dc','sc','K0'),
 ('sc','uc','K-'),('sc','dc','K0bar'),('sc','sc','eta'),('sc','sc','eta-prime'),
]
# Group the table by parent/daughter, explicitly combining its two S/D pairs.
# label, initial spectroscopic tuple, final tuple, [(old row id, relative L)]
S0=(1,0,0,0); S1=(1,0,1,1); P0=(1,1,1,0); T1=(1,1,1,1)
B1=(1,1,0,1); P2=(1,1,1,2); Dv=(1,2,1,1); D3=(1,2,1,3)
GROUPS=[
 ('T01',S1,S0,[('A1',1)]),
 ('T02',P2,S0,[('A2',2)]),('T03',P2,S1,[('A3',2)]),
 ('T04',T1,S1,[('A4',2),('S1',0)]),('T05',B1,S1,[('A5',2),('S2',0)]),
 ('T06',D3,S0,[('A6',3)]),('T07',D3,S1,[('A7',3)]),
 ('T08',B1,P0,[('Aprime1',1)]),('T09',D3,B1,[('Adoubleprime1',2)]),
 ('T10',T1,P0,[('A01',1)]),('T11',P0,S0,[('S3',0)]),
 ('T12',Dv,S0,[('D1',1)]),('T13',Dv,S1,[('D2',1)]),
 ('T14',(2,0,0,0),S1,[('P1',1)]),('T15',(2,0,1,1),S0,[('P2',1)]),
 ('T16',(2,0,1,1),S1,[('P3',1)]),
]
EMIT={'pi0':'π⁰','pi+':'π⁺','pi-':'π⁻','eta':'η','eta-prime':'η′','K+':'K⁺','K-':'K⁻','K0':'K⁰','K0bar':'K̄⁰'}
FLAVOR={'cu':('D','⁰'),'cd':('D','⁺'),'cs':('D_s','⁺'),'uc':('D̄','⁰'),'dc':('D','⁻'),'sc':('D_s','⁻')}
SUFFIX={S0:'',S1:'*',P0:'₀*',T1:'₁[t]',B1:'₁′[s]',P2:'₂*',Dv:'(1³D₁)',D3:'(1³D₃)',(2,0,0,0):'(2¹S₀)',(2,0,1,1):'(2³S₁)'}
WAVE={0:'S',1:'P',2:'D',3:'F'}
def meson(level,flavor):
    base,charge=FLAVOR[flavor]
    return base+SUFFIX[level]+charge

def main():
    assert len(ROUTES)==len(set(ROUTES))==28
    q3={'u':2,'d':-1,'s':-1,'c':2}
    pq={'pi0':0,'pi+':3,'pi-':-3,'eta':0,'eta-prime':0,'K+':3,'K-':-3,'K0':0,'K0bar':0}
    for parent,daughter,p in ROUTES:
        assert q3[parent[0]]-q3[parent[1]]==q3[daughter[0]]-q3[daughter[1]]+pq[p]
        assert ('c'==parent[0])==('c'==daughter[0])
    # Build full keys without importing the original counting program.
    keys={(row,p+'bar',d+'bar',e) for _,i,f,rs in GROUPS for row,l in rs for p,d,e in ROUTES}
    data=json.loads(gzip.decompress((ROOT/'docs/pseudoscalar_census/census.json.gz').read_bytes()))
    old=[x for x in data['cases']['mixed/nonet']['records'] if x['sector']=='charmed']
    oldkeys={(x['transition'],x['parent'],x['daughter'],x['emitted']) for x in old}
    assert len(old)==len(oldkeys)==len(keys)==504
    assert keys==oldkeys
    catalog={x['id']:x for x in data['catalog']}
    for _,i,f,rs in GROUPS:
        for row,l in rs:
            r=catalog[row]
            assert (tuple(r['initial']),tuple(r['final']),r['ell'])==(i,f,l)
    distinct={(i,f,p,d,e) for _,i,f,rs in GROUPS for p,d,e in ROUTES}
    assert len(distinct)==448
    assert Counter(len(rs) for _,i,f,rs in GROUPS)=={1:14,2:2}
    # Generic axial rotation closes the lone singlet daughter onto both slots.
    physical_pairs={(i,f) for _,i,f,rs in GROUPS} | {(D3,T1)}
    physical_waves={(i,f,l) for _,i,f,rs in GROUPS for _,l in rs} | {(D3,T1,2)}
    assert len(physical_pairs)*len(ROUTES)==476
    assert len(physical_waves)*len(ROUTES)==532
    lines=['# Charmed-meson counting scratch pad','',
      'Generated and checked by `python3 GIPaper/scripts/audit_charmed_census.py`.','',
      'The follow-up [phase-space selection](charmed_thresholds.md) applies calculated charm masses and physical emitted-meson masses; its JSON ledger retains open and closed decisions for every entry.','',
      '## Verdict','',
      '**504 = 18 partial-wave templates × 28 charge/flavor routes.**',
      '**448 = 16 spectroscopic parent/daughter templates × 28 charge/flavor routes.**',
      'The difference is **56 = 2 × 28**: the two axial-basis -> D* templates each have both an S wave and a D wave. Those are two matrix elements of one final channel, not two decays.',
      'An independently written list of 28 routes and 16 templates reproduces every one of the original 504 keys exactly; there are no duplicated keys. Charge conservation and preservation of the heavy spectator are checked separately.','',
      '**Important correction to the interpretation:** these are counts in the unmixed spectroscopic basis, not a completed inventory of physical D1/D1′ channels. The prior description of 448 as physical final channels was too strong. Mixing expands the isolated A-double-prime daughter into both physical axial states (see below). No numerical widths were computed in this census.','',
      '## Names and finite scope','',
      '| Spectroscopic state | Scratch-pad name |',
      '|---|---|',
      '| 1¹S₀ | D |','| 1³S₁ | D* |','| 1³P₀ | D₀* |','| 1³P₂ | D₂* |',
      '| 1³P₁ | D₁[t], a triplet **basis** state |','| 1¹P₁ | D₁′[s], a singlet **basis** state |',
      '| 1³D₁, 1³D₃ | D(1³D₁), D(1³D₃) |',
      '| 2¹S₀, 2³S₁ | D(2¹S₀), D(2³S₁) |','',
      'The primes above are bookkeeping labels, not assignments of measured resonances to pure spins. Physical D₁ and D₁′ are mixtures of [t] and [s]. Replace D by D_s for c sbar, and use the indicated antiparticle and charge for cbar q. No resonance masses are assigned.',
      'The “charmed” sector includes D and D_s, particles and antiparticles. It includes the three 2S templates printed in Table IV. It does not include all possible D-wave parents: 1¹D₂ and 1³D₂ are absent from the printed template list. Thresholds are deliberately not imposed.','',
      '## The 28 flavor routes for any one template X -> Y + P','',
      'For positive charm there are 5 + 5 + 4 = 14; their distinct charge conjugates give another 14. X and Y stand for whichever spectroscopic parent and daughter are selected below.','',
      '| Route | Parent flavor slot | Daughter flavor slot | Emitted meson |',
      '|---|---|---|---|']
    for n,(p,d,e) in enumerate(ROUTES,1):
        lines.append(f'| R{n:02} | {meson(S0,p)} | {meson(S0,d)} | {EMIT[e]} |')
    lines += ['', 'For example, the first five routes mean X⁰ -> Y⁰π⁰, Y⁰η, Y⁰η′, Y⁺π⁻, Y_s⁺K⁻. For X_s⁺ there are Y⁰K⁺, Y⁺K⁰, Y_s⁺η, Y_s⁺η′. The channel X_s⁺ -> Y_s⁺π⁰ vanishes in exact isospin and is not counted.',
      'Only the light antiquark emits for c qbar; only the light quark emits for q cbar. There is no additional factor of two for emitter topology. Neutral π⁰, η and η′ are each counted once as a meson, not once per flavor component. A heavy daughter and a light emitted pseudoscalar are distinct, so their interchange contributes no extra charm channel.','',
      '## All 16 templates (18 partial waves)','',
      '| ID | Representative neutral channel | Relative waves | Table-IV entries | Final channels | Partial-wave channels |',
      '|---|---|---|---|---:|---:|']
    for tid,i,f,rs in GROUPS:
        lines.append(f'| {tid} | {meson(i,"cu")} -> {meson(f,"cu")} + P | {", ".join(WAVE[l] for _,l in sorted(rs,key=lambda x:x[1]))} | {", ".join(r for r,_ in rs)} | 28 | {28*len(rs)} |')
    lines += ['| Total | | | | **448** | **504** |','',
      'Rows T04 and T05 are the only double-wave rows. T08 and T10 have different axial basis parents; T06 and T12 have J=3 and J=1 parents, respectively. None are duplicates.',
      'The original table also prints two 1³D₃ -> 1S entries with relative L=1; both fail the J triangle and neither is counted here. The retained 1³D₃ -> D and D* entries have L=3.','',
      '## Useful reductions of the large numbers','',
      '| Scope | Routes per template | Partial-wave entries | Distinct basis final channels |',
      '|---|---:|---:|---:|',
      '| Full nonet, charm and anticharm | 28 | 504 | 448 |',
      '| Full nonet, positive charm only | 14 | 252 | 224 |',
      '| π, K, η emission, both conjugates | 22 | 396 | 352 |',
      '| π, η emission, both conjugates | 14 | 252 | 224 |',
      '| Full nonet, both conjugates, omit the three 2S templates | 28 | 420 | 364 |','',
      'Charge-conjugate channels are separately named channels, not independent reduced amplitudes. Exact isospin further relates charge channels. Neither 504 nor 448 counts independent model parameters or independent experimental constraints.','',
      '## What changes when D₁ and D₁′ mean physical mesons?','',
      'Write D₁ = cos(θ) D₁[t] + sin(θ) D₁′[s], and D₁′ = -sin(θ) D₁[t] + cos(θ) D₁′[s], separately in each light-flavor sector. The phase convention is illustrative; only the two-dimensional span matters here.',
      'T04/T05 already cover both axial basis parents -> D* with S and D waves; rotation leaves two physical parents with these wave possibilities. T08/T10 likewise already cover both axial parents -> D₀*.',
      '**T09 is different:** A-double-prime includes only D(1³D₃) -> D₁′[s] + P with relative L=2. A generic nonzero mixing angle projects this daughter onto both physical D₁ and D₁′. Thus the represented subspace supplies two physical daughter channels per route, not one.',
      'Consequently a generic mixing completion of this restricted template set gives **19 × 28 = 532 partial-wave channels** and **17 × 28 = 476 distinct final channels**. This is a projection/counting statement, not a numerical prediction: a complete amplitude can also require omitted basis matrix elements. Special mixing angles or dynamical cancellations can remove channels. It is still not an exhaustive through-D survey.',
      'For a paper, label 504/448 explicitly as an **unmixed-basis census**. Do not attach those totals unqualified to named physical D₁/D₁′ resonances.','',
      '## Fully expanded ledger: all 448 basis final channels','',
      'Each line below is one distinct final channel. The bracket lists every retained relative partial wave; two-wave lines account for the extra 56 matrix elements.']
    total=0
    for tid,i,f,rs in GROUPS:
        waves=', '.join(WAVE[l] for _,l in sorted(rs,key=lambda x:x[1]))
        lines += ['',f'### {tid}: {meson(i,"cu")} -> {meson(f,"cu")} + P','']
        for n,(p,d,e) in enumerate(ROUTES,1):
            total+=1
            lines.append(f'- {tid}/R{n:02}: {meson(i,p)} -> {meson(f,d)} + {EMIT[e]} [{waves}]')
    assert total==448
    OUT.write_text('\n'.join(lines)+'\n')
    print(f'Validated 504 unique partial-wave keys and 448 unique basis final channels; wrote {OUT}')

if __name__=='__main__': main()
