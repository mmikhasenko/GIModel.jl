"""Read-only PDG.NN bridge. Freeze inputs locally; never modify or train PDG.NN."""
import sys,json,hashlib,platform
from pathlib import Path
source=Path(sys.argv[1]).expanduser() if len(sys.argv)>1 else Path.home()/'Documents/PDG.CAT/PDG.NN'
sys.path.insert(0,str(source))
from pdgnn.features import particle
from pdgnn.flavor import v
from pdgnn.v5 import prior
from pdgnn.topology import topology
import importlib.metadata
root=Path(__file__).resolve().parents[1]
rows=json.loads((source/'data/decays.json').read_text())
# Pure flavor, S-wave channels, with charge-conjugate rates identified at tree level.
channels=[
 ('D0_K',421,321,'cu_1_0','su_1_0','c','s','u',False,[-321,-11,12]),
 ('D0_pi',421,211,'cu_1_0','ud_1_0','c','u','u',False,[-211,-11,12]),
 ('D0_Kstar',421,323,'cu_1_0','su_1_1','c','s','u',True,[-323,-11,12]),
 ('Dplus_Kstar',411,313,'cu_1_0','su_1_1','c','s','u',True,[-313,-11,12]),
 ('Ds_phi',431,333,'cs_1_0','ss_1_1','c','s','s',True,[333,-11,12]),
 ('B0_D',511,411,'bu_1_0','cu_1_0','b','c','u',False,None),
 ('B0_Dstar',511,413,'bu_1_0','cu_1_1','b','c','u',True,None),
 ('B0_pi',511,211,'bu_1_0','ud_1_0','b','u','u',False,None),
 ('B0_rho',511,213,'bu_1_0','ud_1_1','b','u','u',True,None),
 ('Bs_Ds',531,431,'bs_1_0','cs_1_0','b','c','s',False,[-431,-13,14]),
 ('Bs_Dsstar',531,433,'bs_1_0','cs_1_1','b','c','s',True,None),
]
ckms={('c','s'):abs(v('c','s')),('c','u'):abs(v('c','d')),('b','c'):abs(v('c','b')),('b','u'):abs(v('u','b'))}
# c->u here denotes equal-mass u/d spatial wave; physical charged current is c->d.
out=[]
for name,parent,daughter,wi,wf,mi,mf,ms,vec,ds in channels:
 p,d=particle(parent),particle(daughter)
 ref=next((r for r in rows if r['parent']==parent and ds and sorted(r['daughters'])==sorted(ds)),None)
 row=dict(name=name,parent=parent,daughter=daughter,wi=wi,wf=wf,mi=mi,mf=mf,ms=ms,vector=vec,M=p.mass/1000,m=d.mass/1000,tau=p.lifetime*1e-9,ckm=ckms[mi,mf],ml=0.1056583755 if name=='Bs_Ds' else 0.00051099895)
 if ref: row.update(bf_ref=ref['bf'],error=ref['error_positive'],ref_id=ref['id'],ref_url=ref['url'],description=ref['description'],prior_bf=10**prior(ref))
 else: row.update(bf_ref=-1.,error=-1.,ref_id='not_in_flat_extraction',ref_url='',description=name,prior_bf=-1.)
 out.append(row)
# JSON is archival. TOML is read by Julia without an extra dependency.
(root/'data/channels.json').write_text(json.dumps(out,indent=2)+'\n')
with (root/'data/channels.toml').open('w') as f:
 for r in out:
  f.write('[[channel]]\n')
  for k,x in r.items(): f.write(k+' = '+(json.dumps(x) if isinstance(x,(str,bool)) else str(x))+'\n')
  f.write('\n')
leps=[]
for name,pid,key,ckm,ref in [('pi',211,'ud_1_0',abs(v('u','d')),130.2),('K',321,'su_1_0',abs(v('u','s')),155.7),('D',411,'cu_1_0',abs(v('c','d')),212.6),('Ds',431,'cs_1_0',abs(v('c','s')),249.9),('B',521,'bu_1_0',abs(v('u','b')),189.5),('Bs',531,'bs_1_0',0.,230.3),('Bc',541,'bc_1_0',abs(v('c','b')),-1.)]:
 p=particle(pid);leps.append(dict(name=name,pid=pid,key=key,M=p.mass/1000,tau=p.lifetime*1e-9,ckm=ckm,fref_MeV=ref))
with (root/'data/leptonic.toml').open('w') as f:
 for r in leps:
  f.write('[[meson]]\n')
  for k,x in r.items():f.write(k+' = '+(json.dumps(x) if isinstance(x,str) else str(x))+'\n')
non=[]
for name,pid,ds in [('B0_D_pi',511,[-411,211]),('D0_K_pi',421,[-321,211])]:
 r=next(r for r in rows if r['parent']==pid and sorted(r['daughters'])==sorted(ds)); non.append(r)
(root/'data/nonleptonic_references.json').write_text(json.dumps(non,indent=2)+'\n')
manifest=dict(source=str(source),particle_version=importlib.metadata.version('particle'),date='2026-11-07',note='PDG.NN 2026 branching fractions; installed particle masses/lifetimes are separately pinned here. No assumption that editions coincide.',files={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in ['data/decays.json','data/manifest.json','pdgnn/v5.py','pdgnn/flavor.py']},decay_source_manifest=json.loads((source/'data/manifest.json').read_text()))
(root/'data/provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Frozen',len(out),'semileptonic channels and',len(leps),'decay constants')
