"""Independent native current contractions for the charged antiparticle."""
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from itertools import product
from model.openwave_dirac_neutral_current.neutral import ix,leptonic_words,bog
from model.openwave_current_averaging.averaging import act

def amp(d,i,j,se,sn,sbar,p,o):
 parent=1<<ix(0,p,9,1);target=(1<<ix(0,se,9,0))|(1<<ix(o,sn,10,i))|(1<<ix(-o,sbar,15,j));changed={ix(0,p,9,1),ix(0,se,9,0),ix(o,sn,10,i),ix(-o,sbar,15,j)};out=0j
 for(a,b,c,e),coefficient in leptonic_words():
  choices=[]
  for k,dagger in((a,True),(b,True),(e,False),(c,False)):
   row=[(mode,create!=dagger,np.conj(z)if dagger else z)for mode,create,z in d['mapping'][k]if mode in changed]
   if not row:break
   choices.append(row)
  if len(choices)!=4:continue
  for term in product(*choices):
   k,sign=act(parent,tuple((mode,create)for mode,create,z in term))
   if k==target and sign:out+=coefficient*sign*np.prod([z for mode,create,z in term])
 return out

@lru_cache(None)
def exact_antiparticle_spin_contraction():
 rows=[]
 for o,p,se,sn,sbar in product((-1,1),range(2),range(2),range(2),range(2)):
  parent=1<<ix(0,p,9,1);target=(1<<ix(0,se,9,0))|(1<<ix(o,sn,10,0))|(1<<ix(-o,sbar,15,1))
  mapping={ix(0,1-p,11,1):(ix(0,p,9,1),True,1-2*p),ix(0,1-se,11,0):(ix(0,se,9,0),True,1-2*se),
   ix(o,sn,10,0):(ix(o,sn,10,0),False,1),ix(o,1-sbar,10,1):(ix(-o,sbar,15,1),True,-1)}
  coefficient=F(0)
  for(a,b,c,e),z in leptonic_words():
   if any(k not in mapping for k in(a,b,c,e)):continue
   ops=[];factor=1
   for k,dagger in((a,True),(b,True),(e,False),(c,False)):
    mode,create,sign=mapping[k];ops.append((mode,create!=dagger));factor*=sign
   output,sign=act(parent,tuple(ops))
   if output==target:coefficient+=F(z)*sign*factor/2
  expected=F(o)if sn==p and se!=sbar else F(0);assert coefficient==expected
  rows.append({'orientation':o,'parent_spin':p,'charged_spin':se,'neutral_spin':sn,'antineutral_spin':sbar,'exact_native_CAR_coefficient':str(coefficient)})
 return rows

def complete_antiparticle_vertex(d):
 labels=list(product((-1,1),range(3),range(3),range(2),range(2),range(2)))
 W=np.array([[amp(d,i,j,se,sn,sbar,p,o)for p in range(2)]for o,i,j,se,sn,sbar in labels]);expected=np.zeros_like(W)
 for k,(o,i,j,se,sn,sbar)in enumerate(labels):
  for p in range(2):
   if sn!=p or se==sbar:continue
   u=bog(o*d['q'],1-2*sn,d['neutral'][i])[1];v=bog(o*d['q'],1-2*(1-sbar),d['neutral'][j])[2]
   expected[k,p]=o*u*v*d['R'][0,i].conj()*d['R'][1,j]
 assert np.linalg.norm(W-expected)<3e-13 and np.linalg.norm(W.conj().T@W-np.eye(2))<3e-13
 return labels,W
