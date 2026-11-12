"""Independent full-column symmetry and finite polar-preparation checks."""
import json,hashlib,random
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from model.coherent_defect_wave.chiral_joint.oriented_links.model import reverse_label
from model.openwave_neutral_stationary.source import full_column
ROOT=Path(__file__).resolve().parents[3]
def change(key,axis):
 mask,links,scalars=key;bit=4>>axis;result=list(links)
 region=sum(((1<<30)-1)<<(30*site)for site in range(8)if site&bit)
 for site in range(8):result[3*(site^bit)+axis]=reverse_label(links[3*site+axis])
 return (mask,tuple(result),scalars),(-1)**((mask&region).bit_count())
def translate(key,axis):
 mask,links,scalars=key;bit=4>>axis;occupied=[i for i in range(240)if mask>>i&1]
 perm=[30*((i//30)^bit)+i%30 for i in occupied]
 sign=(-1)**sum(a>b for j,a in enumerate(perm)for b in perm[j+1:])
 outlinks=[0]*24;outscalar=[0]*8
 for site in range(8):
  outscalar[site^bit]=scalars[site]
  for direction in range(3):outlinks[3*(site^bit)+direction]=links[3*site+direction]
 return (sum(1<<i for i in perm),tuple(outlinks),tuple(outscalar)),sign

def run():
 import sympy as s
 rng=random.Random(1254);checks=0
 for _ in range(40):
  key=(rng.getrandbits(240),tuple(rng.randrange(119)for _ in range(24)),tuple(rng.randrange(5)for _ in range(8)))
  for a in range(3):
   u,sign=change(key,a);v,sign2=change(u,a);assert v==key and sign*sign2==1
   for b in range(3):
    v,sign2=translate(u,b);z,sg=translate(key,b);z,sg2=change(z,a)
    assert v==z and sign*sign2==sg*sg2*((-1)**key[0].bit_count()if a==b else 1);checks+=1
 keys=[((1<<5)|(1<<41)|(1<<209),(0,)*24,(0,)*8),
       ((1<<11)|(1<<54)|(1<<188),tuple((7*j+3)%119 for j in range(24)),(0,1,0,0,2,0,0,3))]
 column_checks=0
 for key in keys:
  column=full_column(key)
  for axis in range(3):
   left={}
   for target,value in column.items():
    out,sign=change(target,axis);left[out]=left.get(out,0)+sign*value
   target,sign=change(key,axis);right={k:sign*v for k,v in full_column(target).items()}
   assert all(s.simplify(left.get(k,0)-right.get(k,0))==0 for k in left.keys()|right.keys());column_checks+=1
 # Independent finite geometry witnesses of polar density/stationarity bounds.
 rng=np.random.default_rng(162);trials=0;worst=0.
 for eta in(.001,.01,.08):
  rank=3;dim=8;theta=np.arcsin(eta);Pcols=np.eye(dim)[:,:rank]
  for j in range(rank):Pcols[j,j]=np.cos(theta);Pcols[j+rank,j]=np.sin(theta)
  E=np.diag([1.]*rank+[0.]*(dim-rank));P=Pcols@Pcols.T;epsilon=1e-6
  M=(1-epsilon/2)*E+(epsilon/3)*(np.eye(dim)-E)
  A=Pcols.T@M@M@Pcols;q=(1-eta*eta)*(1-epsilon)**2;alpha=1-q;tau=alpha*alpha/(1-alpha)
  L=np.sqrt(q)*M@Pcols@(3*np.eye(rank)-A)/2@Pcols.T
  out=L@P@L.T/rank;success=np.trace(out);rho=out/success
  bound=4*tau/(1-tau)**2+2*epsilon/q
  error=np.linalg.norm(rho-E/rank,ord='nuc');assert error<=bound
  assert q*(1-tau)**2<=success<=q
  H=np.diag(rng.uniform(-1,1,dim));U=expm(-2.7j*H)
  assert np.linalg.norm(U@rho@U.conj().T-rho,ord='nuc')<=2*bound
  worst=max(worst,error/bound);trials+=1
 refs=['model/openwave_neutral_dispersion/'+f for f in ['symmetry.py','flatness.py','PROOF.md','results.json']]+['model/openwave_neutral_preparation/'+f for f in ['filtering.py','local.py','records.py','PROOF.md','results.json']]
 return {'verdict':'PASS','independent_translation_parity_relations':checks,'independent_complete_native_H_columns':column_checks,
  'independent_polar_density_and_stationarity_trials':trials,'maximum_error_over_bound':worst,
  'proof_review':'Full-link inversion/plane parity gives exact odd-sector momentum flatness. Complete local algebraic seeds and a contraction polar filter replace the full-dimension preparation penalty. Only reduced-source near-stationarity follows; primitive gates and moving-particle identity are unproved.',
  'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in refs}}
if __name__=='__main__':
 v=run();Path(__file__).with_suffix('.json').write_text(json.dumps(v,indent=2,sort_keys=True)+'\n');print(v['verdict'])
