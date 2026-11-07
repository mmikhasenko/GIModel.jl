"""Audit committed numerical products and capture focused-test outcomes."""
import csv,json,re,math
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def rows(n):return list(csv.DictReader((root/'results'/n).open()))
a=rows('rates_ho48.csv');assert len(a)==22
for r in a:
 assert math.isfinite(float(r['width_GeV'])) and float(r['width_GeV'])>0
 assert 0<=float(r['FL'])<=1
 assert 0<float(r['bf'])<1
s=json.loads((root/'results/summary.json').read_text())
assert max(s['ho32_relative_width_difference'].values())<.002
assert max(s['fd600_relative_width_difference'].values())<.002
assert max(s['quadrature_relative_difference'].values())<1e-4
r=rows('resonances.csv');assert len(r)==20
assert all(abs(float(x['ratio'])-1)<1e-4 for x in r if float(x['width_scale'])<.001)
lf=rows('lightfront.csv')
for name in {r['channel'] for r in lf}:
 a=[r for r in lf if r['channel']==name and r['Q2']=='0.0' and r['spin']=='true']
 assert max(float(r['fp']) for r in a)/min(float(r['fp']) for r in a)-1<1e-4
extended=rows('extended_convergence.csv')
assert len(extended)==6 and max(float(r['relative_difference']) for r in extended)<.001
logs={}
for name in ['tests','lightfront','rescattering']:
 text=(root/'results'/f'{name}.log').read_text()
 assert 'ERROR:' not in text and 'Test Failed' not in text
 lines=[l.strip() for l in text.splitlines() if '|' in l and not l.startswith('Test Summary')]
 logs[name]=lines
(root/'results/validation.json').write_text(json.dumps(dict(date='2026-11-07',status='passed',focused_test_summaries=logs,checks=['22 positive ground-state rates; physical FL','HO32/HO48 and FD600/HO48 differences below 0.2 percent','selected quadrature differences below 0.01 percent','four narrow-width limits within 0.01 percent','LF zero-transfer convergence below 0.01 percent','refined finite-width and spacelike LF integration differences below 0.1 percent'],caveat='These checks validate implementation/numerics, not model completeness or accuracy.'),indent=2)+'\n')
print('Output audits passed; focused test summaries saved.')
