#!/usr/bin/env python3
"""Exact integer flavor census of the printed Table-IV partial waves.
No masses, widths, fits or numerical wavefunctions are used. See generated report.
"""
from collections import Counter, defaultdict
from pathlib import Path
import json
import gzip

OUT = Path(__file__).resolve().parents[1] / 'docs/pseudoscalar_census'
# n, internal L, constituent spin S, total J.
S0=(1,0,0,0); S1=(1,0,1,1); P0=(1,1,1,0)
P1=(1,1,1,1); B1=(1,1,0,1); P2=(1,1,1,2)
D1=(1,2,1,1); D3=(1,2,1,3)
CATALOG=[]
def add(family, initial, final, ell):
    CATALOG.append(dict(id=f'{family}{1+sum(r["family"]==family for r in CATALOG)}',
                        family=family,initial=initial,final=final,ell=ell))
for i,f,l in [(S1,S0,1),(P2,S0,2),(P2,S1,2),(P1,S1,2),(B1,S1,2),(D3,S0,3),(D3,S1,3)]: add('A',i,f,l)
add('Aprime',B1,P0,1); add('Adoubleprime',D3,B1,2); add('A0',P1,P0,1)
for i,f in [(P1,S1),(B1,S1),(P0,S0)]: add('S',i,f,0)
# Literal print, not a silently inferred correction: its last two D rows fail J.
for i in (D1,D3):
    for f in (S0,S1): add('D',i,f,1)
for i,f in [((2,0,0,0),S1),((2,0,1,1),S0),((2,0,1,1),S1)]: add('P',i,f,1)

def allowed(r):
    ni,li,si,ji=r['initial']; nf,lf,sf,jf=r['final']; ell=r['ell']
    return abs(ji-jf)<=ell<=ji+jf and (li-lf-ell-1)%2==0 and (si,sf)!=(0,0)

def pole(name,sector,components): return dict(name=name,sector=sector,components=components)
def flavors(mixed=True):
    states=[pole('pi+','isovector',{('u','d'):1}),pole('pi0','isovector',{('u','u'):1,('d','d'):-1}),pole('pi-','isovector',{('d','u'):1})]
    # Generic orthogonal isoscalar mixtures, normalization irrelevant for zeros.
    # Symbols eta/eta-prime denote the two isoscalar flavor slots; for other
    # spectroscopic levels they are the corresponding two isoscalar poles.
    states += [pole('eta','isoscalar',{('u','u'):1,('d','d'):1,('s','s'):1} if mixed else {('u','u'):1,('d','d'):1}),
               pole('eta-prime','isoscalar',{('u','u'):1,('d','d'):1,('s','s'):-2} if mixed else {('s','s'):1})]
    for name,pair in [('K+',('u','s')),('K0',('d','s')),('K-',('s','u')),('K0bar',('s','d'))]: states.append(pole(name,'strange',{pair:1}))
    for heavy,sector in [('c','charmed'),('b','bottom-flavored')]:
        for light in 'uds':
            for pair in [(heavy,light),(light,heavy)]: states.append(pole(pair[0]+pair[1]+'bar',sector,{pair:1}))
    return states

charge3=dict(u=2,d=-1,s=-1,c=2,b=-1)
def charge(state):
    values={charge3[a]-charge3[b] for a,b in state['components']}
    assert len(values)==1
    return values.pop()
def mass_class(f): return 'n' if f in 'ud' else f

def flavor_blocks(initial,final,emitted,sign):
    """Sum isospin/C-related paths, keeping different radial integrals independent.

    Charge conjugation relates the anti-emitter coefficient to the quark one
    by (-1)^(Li+Si+Lf+Sf). Blocks retain emitter/spectator/produced mass classes;
    hence SU(3)-breaking overlaps cannot spuriously cancel against each other.
    """
    blocks=defaultdict(int)
    for (a,b),ci in initial['components'].items():
        for z in 'uds':
            if a in 'uds':
                cf=final['components'].get((z,b),0); cp=emitted['components'].get((a,z),0)
                blocks[tuple(map(mass_class,(a,b,z)))]+=ci*cf*cp
            if b in 'uds':
                cf=final['components'].get((a,z),0); cp=emitted['components'].get((z,b),0)
                blocks[tuple(map(mass_class,(b,a,z)))]+=sign*ci*cf*cp
    return {','.join(k):v for k,v in blocks.items() if v}

def census(mixed=True, emission_set='nonet'):
    states=flavors(mixed); fields=flavors(True)[:9]
    if emission_set=='no_eta_prime': fields=[p for p in fields if p['name']!='eta-prime']
    if emission_set=='pi_eta': fields=[p for p in fields if p['sector']=='isovector' or p['name']=='eta']
    records=[]
    for r in CATALOG:
        if not allowed(r): continue
        sign=(-1)**(sum(r['initial'][1:3])+sum(r['final'][1:3]))
        for i in states:
            for f in (flavors(True) if r['final']==S0 else states):
                for p in fields:
                    blocks=flavor_blocks(i,f,p,sign)
                    if not blocks: continue
                    assert charge(i)==charge(f)+charge(p)
                    # Same ground pseudoscalar multiplet: interchange final
                    # mesons is a second Eq19 role assignment, not another decay.
                    final_key=tuple(sorted((f['name'],p['name']))) if r['final']==S0 else (f['name'],p['name'])
                    records.append(dict(transition=r['id'],sector=i['sector'],parent=i['name'],daughter=f['name'],emitted=p['name'],
                                        final_key=final_key,blocks=blocks))
    totals={}
    for sector in ['isovector','isoscalar','strange','charmed','bottom-flavored']:
        rows=[x for x in records if x['sector']==sector]
        pw={(x['transition'],x['parent'],tuple(x['final_key'])) for x in rows}
        by_id={r['id']:r for r in CATALOG}
        decays={(by_id[x['transition']]['initial'],by_id[x['transition']]['final'],x['parent'],tuple(x['final_key'])) for x in rows}
        totals[sector]=dict(ordered_kernel_channels=len(rows),partial_wave_channels=len(pw),distinct_final_channels=len(decays))
    return dict(totals=totals,records=records)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert len(CATALOG)==20 and sum(map(allowed,CATALOG))==18
    assert len({(r['initial'],r['final']) for r in CATALOG if allowed(r)})==16
    # Symmetry controls: rho0->pi0 pi0 forbidden; rho+->pi+ pi0 allowed;
    # isoscalar->isovector eta forbidden, and heavy spectator never changes.
    st={x['name']:x for x in flavors()}
    assert not flavor_blocks(st['pi0'],st['pi0'],st['pi0'],-1)
    assert flavor_blocks(st['pi+'],st['pi+'],st['pi0'],-1)
    assert not flavor_blocks(st['eta'],st['pi0'],st['eta'],1)
    # Every retained channel has its charge-conjugate counterpart.
    def conjugate_name(name):
        state=st[name]
        flipped={(b,a):v for (a,b),v in state['components'].items()}
        return next(k for k,v in st.items() if v['components']==flipped)
    rows=census()['records']
    keys={(x['transition'],x['parent'],x['daughter'],x['emitted']) for x in rows}
    assert all((t,conjugate_name(i),conjugate_name(f),conjugate_name(p)) in keys for t,i,f,p in keys)
    cases={f'{basis}/{emissions}':census(basis=='mixed',emissions)
           for basis in ('mixed','ideal') for emissions in ('nonet','no_eta_prime','pi_eta')}
    payload=dict(catalog=[dict(r,angular_allowed=allowed(r)) for r in CATALOG],cases=cases)
    (OUT/'census.json.gz').write_bytes(gzip.compress((json.dumps(payload,indent=2)+'\n').encode(),mtime=0))
    lines=['# Pseudoscalar-emission channel census','',
           'Generated by `python3 GIPaper/scripts/count_pseudoscalar_channels.py`.','',
           'For the charm subset, see the [named-channel scratch pad](charmed_scratchpad.md) and the [mass-threshold selection](charmed_thresholds.md).','',
           'The [five-sector threshold/isospin scratch pad](five_sector_scratchpad.md) adds model masses, threshold cuts, and explicit isospin and charge-conjugation grouping. Its excited isoscalars use ideal nn/ss basis states; do not mix its totals with the generic-flavor-mixing estimates below.','',
           '## Printed Table IV','',
           '7 A + 1 A-prime + 1 A-double-prime + 1 A0 + 3 S + 4 D + 3 P = **20 partial-wave entries**.',
           'Ac and Sc are charm-specific formulas for existing transitions, not additional transitions.',
           'Two printed D entries (3D3 -> 1S0 or 3S1 with relative L=1) violate the angular-momentum triangle. They are excluded, leaving **18 allowed entries and 16 distinct spectroscopic parent/daughter pairs**. No conjectural erratum is applied.',
           'The table includes 2S -> 1S. It is a representative list, not an exhaustive list of transitions through D waves.','',
           '## Same counting convention as the histogram','',
           '| Parent sector | Flavor/emitter routes | Nine plotted rows | 20 printed rows | 18 allowed rows |',
           '|---|---:|---:|---:|---:|']
    basis_parents={'isovector':['ud','du'],'isoscalar':['uu','dd','ss'],
                   'strange':['us','ds','su','sd'],
                   'charmed':[p for q in 'uds' for p in ('c'+q,q+'c')],
                   'bottom-flavored':[p for q in 'uds' for p in ('b'+q,q+'b')]}
    basis_counts={sector:sum(3 for pair in pairs for emitter in pair if emitter in 'uds')
                  for sector,pairs in basis_parents.items()}
    assert sum(basis_counts.values())==90
    for sector,n in basis_counts.items():
        lines.append(f'| {sector} | {n} | {9*n} | {20*n} | {18*n} |')
    lines+=['| Total | 90 | 810 | 1800 | 1620 |','',
        'These are individual ordered flavor/emitter basis contributions, not physical decays. The old histogram calls uu/dd/ss components isoscalar and does not separately include a coherent neutral isovector. This convention must not be reinterpreted as physical isospin multiplets.','',
        '## Charge-resolved spectroscopic-basis channel census','',
        'The spin labels are unmixed spectroscopic basis states, including the two axial states in open-flavor sectors. Isoscalar flavor mixing below does not include axial spin mixing. See [the explicit charm audit](charmed_scratchpad.md): 504/448 are basis counts; generic axial mixing expands the represented charm set to 532/476. These are not exhaustive physical-resonance counts.',
        'Each spectroscopic level has a full isovector triplet, two isoscalar flavor states, four strange states, six open-charm states and six open-bottom states. Heavy quarks are spectators. Both charge conjugates are counted. Isospin is exact (mu=md). Flavor blocks with different emitter/spectator/produced masses remain independent; no SU(3)-symmetric radial cancellation is assumed.',
        'The mixed case uses two generic isoscalar mixtures for parent and surviving meson; ideal uses pure nn and ss except that ground-pseudoscalar daughters always use the same mixed eta/eta-prime states as the emitted fields. Mixing coefficients are algebraic representatives, not fitted angles. Special mixing angles can reduce counts.',
        'For each flavor block the algebra is Fq + epsilon Fanti, where epsilon=(-1)^(Li+Si+Lf+Sf). Spin interchange and the spatial/operator parity give this relative sign. Nine retained partial waves have epsilon=+1 and nine have epsilon=-1. In the mixed/nonet case the per-row counts are respectively (20,22,32,28,28) and (14,8,32,28,28) in sector order, so the ordered total is 9*(130+110)=2160.',
        'Quark and antiquark emission paths are summed with the relative charge-conjugation sign before counting nonzero flavor blocks. Charge conservation is asserted for every retained entry. No phase-space cut is made. Accidental dynamical zeros are not predicted.',
        'Ordered kernels distinguish surviving and emitted pseudoscalars. Partial-wave channels identify their interchange when both are ground pseudoscalars. Distinct final channels additionally merge the S/D waves for the same decay.','']
    for name,case in cases.items():
        lines += [f'### {name}','', '| Sector | Ordered kernels | Partial-wave channels | Distinct basis final channels |','|---|---:|---:|---:|']
        for sec,v in case['totals'].items(): lines.append(f'| {sec} | {v["ordered_kernel_channels"]} | {v["partial_wave_channels"]} | {v["distinct_final_channels"]} |')
        sums=[sum(v[k] for v in case['totals'].values()) for k in ('ordered_kernel_channels','partial_wave_channels','distinct_final_channels')]
        lines += ['| Total | '+' | '.join(map(str,sums))+' |','']
    lines += ['## Interpretation and scope','',
              'These are algebraically permitted candidates, not computed widths or a count of PDG measurements. The original nine-row calculation remains the only native-wave numerical survey here. A PDG comparison additionally needs masses and thresholds, physical mixing, coupling calibration, final-state symmetrization and uncertainties.',
              'Going beyond the table requires an explicit finite level inventory (radial range as well as maximum orbital L). There is no finite universal count from "up to D" alone: radial excitations are unbounded. Additional D2 transitions and other P/D daughter channels are absent from this representative table. The operator API has no Table-IV row restriction; numerical validation of an extension is a separate step.',
              'Suggested paper wording: Even the angular-momentum-consistent representative transitions of Table IV generate about two thousand charge-resolved partial-wave candidates in the spectroscopic basis across the light, strange, charm and bottom sectors when the full pseudoscalar nonet is included. This is an algebraic inventory before threshold cuts, not a claim that all corresponding widths have been evaluated or measured.',
              'Source: original GI (1985), Table IV, printed p. 201 (PDF p. 13), checked against the page image rather than the OCR.','']
    (OUT/'README.md').write_text('\n'.join(lines))
    print('\n'.join(lines))
if __name__=='__main__':main()
