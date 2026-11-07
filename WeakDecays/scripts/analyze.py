"""Generate research figures and numerical summaries from immutable CSV outputs."""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
root=Path(__file__).resolve().parents[1]; out=root/'report/figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight','figure.dpi':150})
def rows(name):return list(csv.DictReader((root/'results'/name).open()))
def num(r,k):return float(r[k].split(' + ')[0])
def save(fig,name):fig.tight_layout();fig.savefig(out/(name+'.pdf'));plt.close(fig)
rates=rows('rates_ho48.csv'); channels=json.loads((root/'data/channels.json').read_text()); lut={r['name']:r for r in channels}
colors={'galilean':'#176c9e','boosted_spectator':'#c16a24'}
labels={'galilean':'Galilean impulse','boosted_spectator':'Scalar spectator boost'}
# Numerical diagnostics: do not turn model spread into a confidence interval.
summary={}
for method in ['ho32','fd600']:
 other={(r['channel'],r['scheme']):r for r in rows('rates_'+method+'.csv')}
 diffs={r['channel']+'/'+r['scheme']:abs(num(other[r['channel'],r['scheme']],'width_GeV')/num(r,'width_GeV')-1) for r in rates}
 summary[method+'_relative_width_difference']=diffs
conv=rows('convergence.csv'); diffs={}
for name in {r['channel'] for r in conv}:
 for scheme in colors:
  a=[r for r in conv if r['channel']==name and r['scheme']==scheme]
  diffs[name+'/'+scheme]=abs(num(a[1],'width_GeV')/num(a[2],'width_GeV')-1)
summary['quadrature_relative_difference']=diffs
# Comparative branching fractions.
fig,ax=plt.subplots(figsize=(8,4.2));names=[r['name'] for r in channels if r['bf_ref']>0]
for scheme,color in colors.items():
 vals=[num(next(r for r in rates if r['channel']==n and r['scheme']==scheme),'ratio') for n in names]
 ax.plot(np.arange(len(names)),vals,'o-',label=labels[scheme],color=color)
ax.axhline(1,color='black',lw=1);ax.set_yscale('log');ax.set_xticks(range(len(names)),names,rotation=20);ax.set_ylabel('Calculated BF / reference BF');ax.legend();ax.set_title('Unfitted currents: agreement is channel dependent');save(fig,'benchmark_ratios')
# Decay constants.
lep=rows('leptonic_fd600.csv');fig,ax=plt.subplots(figsize=(7.5,3.5));ix=np.arange(6)
ax.bar(ix-.16,[num(r,'f_MeV') for r in lep[:6]],.32,label='GI current (FD600)',color='#176c9e');ax.bar(ix+.16,[num(r,'reference_MeV') for r in lep[:6]],.32,label='PDG 2025 preferred theory',color='#66796e');ax.set_xticks(ix,[r['meson'] for r in lep[:6]]);ax.set_ylabel(r'$f_P$ [MeV]');ax.legend();save(fig,'decay_constants')
curves=rows('curves_ho48.csv');fig,axes=plt.subplots(2,2,figsize=(8,6))
for ax,name in zip(axes.flat,['D0_K','D0_pi','B0_D','B0_pi']):
 for scheme,color in colors.items():
  a=[r for r in curves if r['channel']==name and r['scheme']==scheme]
  ax.plot([num(r,'q2') for r in a],[num(r,'fp') for r in a],label=labels[scheme],color=color)
 ax.set_title(name);ax.set_xlabel(r'$q^2$ [GeV$^2$]');ax.set_ylabel(r'$f_+(q^2)$')
axes[0,0].legend(fontsize=8);save(fig,'form_factors')
# Light front comparison and controlled interpolation/extrapolation experiment.
lf=rows('lightfront.csv');fig,ax=plt.subplots(figsize=(8,4));names=[r['name'] for r in channels if not r['vector']]
for scheme,color in colors.items():
 vals=[num(next(r for r in rates if r['channel']==n and r['scheme']==scheme),'fp0') for n in names]
 ax.plot(range(len(names)),vals,'o-',label=labels[scheme],color=color)
vals=[num(next(r for r in lf if r['channel']==n and r['Q2']=='0.0' and r['n']=='160' and r['spin']=='true'),'fp') for n in names]
ax.plot(range(len(names)),vals,'s--',label='Light-front plus current',color='#5b7d32');ax.set_xticks(range(len(names)),names);ax.set_ylabel(r'$f_+(0)$');ax.legend(fontsize=8);save(fig,'lightfront')
GF=1.1663787e-5;HBAR=6.582119569e-25
# Pole values are explicit phenomenological inputs, not solved transition poles.
poles={'D0_K':2.1122,'D0_pi':2.01027,'B0_D':6.33,'B0_pi':5.3247,'Bs_Ds':6.33}
zrows=[]
for name in names:
 r=lut[name];M,m=r['M'],r['m'];tp=(M+m)**2;tm=(M-m)**2;t0=tp*(1-np.sqrt(1-tm/tp));pole=poles[name]
 def z(q):return (np.sqrt(tp-q)-np.sqrt(tp-t0))/(np.sqrt(tp-q)+np.sqrt(tp-t0))
 a=[s for s in lf if s['channel']==name and s['n']=='100' and s['spin']=='true']
 for order in (1,2):
  for maxQ in (2.,4.):
   train=[s for s in a if num(s,'Q2')<=maxQ];q=-np.array([num(s,'Q2') for s in train]);f=np.array([num(s,'fp') for s in train])
   co=np.polynomial.polynomial.polyfit(z(q)-z(0),f*(1-q/pole**2),order)
   def ff(q):return np.polynomial.polynomial.polyval(z(q)-z(0),co)/(1-q/pole**2)
   def dwidth(q):
    p=np.sqrt(max(0,((M+m)**2-q)*((M-m)**2-q)))/(2*M)
    return GF**2*r['ckm']**2*p**3*ff(q)**2/(24*np.pi**3)
   width=quad(dwidth,0,tm,epsabs=1e-24)[0]
   residual=max(abs(ff(qi)/fi-1) for qi,fi in zip(q,f))
   zrows.append(dict(channel=name,order=order,max_Q2=maxQ,pole_GeV=pole,bf_massless=width*r['tau']/HBAR,fp_endpoint=ff(tm),max_training_relative_residual=residual,coefficients=';'.join(map(str,co))))
with (root/'results/z_extrapolation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=zrows[0].keys(),lineterminator="\n");w.writeheader();w.writerows(zrows)
summary['z_extrapolation_bf_ranges']={n:[min(r['bf_massless'] for r in zrows if r['channel']==n),max(r['bf_massless'] for r in zrows if r['channel']==n)] for n in names}
# Heavy limit.
h=rows('heavy_limit.csv');fig,axes=plt.subplots(1,2,figsize=(8,3.6))
for scheme,color in colors.items():
 a=[r for r in h if r['scheme']==scheme and r['w']=='1.2'];axes[0].plot([num(r,'heavy_scale') for r in a],[abs(num(r,'hminus')) for r in a],'o-',color=color,label=labels[scheme]);axes[1].plot([num(r,'heavy_scale') for r in a],[num(r,'slope') for r in a],'o-',color=color)
axes[0].set_xscale('log');axes[0].set_yscale('log');axes[0].set_ylabel(r'$|h_-(1.2)|$');axes[0].set_xlabel('Heavy masses / physical values');axes[0].legend(fontsize=7)
axes[1].set_xscale('log');axes[1].axhline(.75,ls='--',color='black',label='HQ sum-rule lower bound');axes[1].set_ylabel(r'Near-zero-recoil slope $\rho^2$');axes[1].set_xlabel('Heavy masses / physical values');axes[1].legend(fontsize=7);save(fig,'heavy_limit')
# Width folding.
rr=rows('resonances.csv');fig,ax=plt.subplots(figsize=(7.5,3.6))
for name,marker in [('D0_Kstar','o'),('B0_rho','s')]:
 for scheme,color in colors.items():
  a=[r for r in rr if r['channel']==name and r['scheme']==scheme and r['running']=='true']
  ax.plot([num(r,'R') for r in a],[num(r,'ratio') for r in a],marker+'-',label=name+', '+labels[scheme],color=color,alpha=1 if marker=='o' else .6)
ax.axhline(1,color='black',lw=1);ax.set_xlabel(r'Barrier radius [GeV$^{-1}$]');ax.set_ylabel('Finite width / stable rate');ax.legend(fontsize=7);save(fig,'finite_width')
# Coherence and elastic unitarization. Arbitrary common vertical normalization.
a=rows('rescattering.csv');fig,axes=plt.subplots(1,2,figsize=(8,3.6))
for ax,sgn in zip(axes,[-1,1]):
 b=[r for r in a if int(float(r['sign']))==sgn];x=np.array([num(r,'mass') for r in b]); norm=max(num(r,'incoherent') for r in b)
 for key,label,col in [('incoherent','Incoherent poles','#999999'),('coherent_bw','Coherent Breit-Wigners','#176c9e'),('kmatrix','Elastic K-matrix','#c16a24')]:ax.plot(x,[num(r,key)/norm for r in b],label=label,color=col)
 ax.set_title('Relative source sign '+str(sgn));ax.set_xlabel(r'$m_{K\pi}$ [GeV]');ax.set_ylabel('Illustrative intensity');ax.set_yscale('log');ax.set_ylim(1e-6,2.)
axes[0].legend(fontsize=7);save(fig,'rescattering')
(root/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:max(v.values()) if isinstance(v,dict) and v and isinstance(next(iter(v.values())),float) else v for k,v in summary.items()},indent=2))
