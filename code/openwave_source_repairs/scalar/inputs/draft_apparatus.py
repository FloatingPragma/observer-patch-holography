"""Complete native preparation/spin words and current-pointer instruments."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import importlib.util,sys,math
import numpy as np
import source as S
from model.observer.source_scalar import compiler as C
from model.native_effect_synthesis import intervals as iv,synthesis as syn
from model.native_interface import pointer
from model.native_clock import configuration
ROOT=S.ROOT
sp=importlib.util.spec_from_file_location('symmetric_spin_interval_reader',ROOT/'openwave/shared/native_current_clock/interval_audit.py');ivread=importlib.util.module_from_spec(sp);sys.modules[sp.name]=ivread;sp.loader.exec_module(ivread)
PAULI=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex))
def pi():return iv.binary_exact(iv.pi_grid(),iv.AQ)
def prep(module=0):
 p=pi();a=iv.binary_exact(iv.split_angle(F(2),F(1)),iv.AQ)
 word=[C.pulse('Y',-a,[3]),C.pulse('Y',-p/8,[2])]
 label=['I']*7;label[2]='Y';label[3]='Z';word+=C.pauli_rotation(''.join(label),-p/8)
 word +=[C.pulse('Y',-p/4,[0]),C.pulse('Y',-p/4,[1]),C.pulse('Y',p/4,[4])]
 if module:word +=[C.pulse('I',-p/2,[6]),C.pulse('X',p/2,[6])]
 return word

def qt(a):
 w,x,y,z=map(lambda v:float(v[0])+float(v[1])*(1+math.sqrt(5))/2,a)
 return w*np.eye(2)-1j*(x*PAULI[0]+y*PAULI[1]+z*PAULI[2])
def spinword(U):
 a,b=U[0,0],U[0,1]
 beta=math.atan2(abs(b),abs(a));aa=math.atan2(a.imag,a.real) if abs(a)>1e-15 else 0.;bb=math.atan2(b.imag,b.real) if abs(b)>1e-15 else 0.
 alpha=(aa+bb)/2;gamma=(aa-bb)/2
 return [C.pulse('Z',gamma,[0]),C.pulse('Y',beta,[0]),C.pulse('Z',alpha,[0])]
def qinterval(a):
 phi=iv.divide(iv.add(iv.exact(1),ivread.root(iv.exact(5))),2)
 return iv.add(iv.exact(a[0]),iv.mul(iv.exact(a[1]),phi))
def qtarget(a):
 w,x,y,z=map(qinterval,a)
 return [[(w,iv.neg(z)),(iv.neg(y),iv.neg(x))],[(y,iv.neg(x)),(w,z)]]

def execute(rows,arrays):
 allwords=[];preps=[];prep_errors=[]
 for m in (0,1):
  word=prep(m);allwords.append(word);U=C.execute(word,7);f=S.ray(m)/math.sqrt(24)
  e=float(np.linalg.norm(U[:,0]-f));assert e<1e-9;assert np.linalg.norm(U.conj().T@U-np.eye(128),'fro')<1e-9
  arrays[f'prep_{m}']=U;preps.append(U);prep_errors.append(e)
 spin_matrices=[];geometric=[];word_errors=[];spin_words=[]
 for a,g in rows:
  wanted=qt(a);word=spinword(wanted);actual=C.execute(word,1);spin_words.append(word);allwords.append(word)
  err=float(np.linalg.norm(actual-wanted,2));assert err<1e-9;word_errors.append(err);spin_matrices.append(actual)
  interval,radius=ivread.execute_interval(word,1);bound=ivread.error_upper(interval,qtarget(a))+radius
  assert bound<F(1,10**9);geometric.append(bound)
 arrays['spin_cover_words']=np.array(spin_matrices)
 # Complete native-current pointer superoperators, both C and E auxiliaries.
 local=pointer.local_matrix();superops={};pointer_error=[];effects={}
 for law in ('C','E'):
  counts=pointer.native_effects(law);sup=np.zeros((2,2,2,2,2),complex)
  for b in (0,1):
   for z in range(4):
    K=local[np.ix_([z,z+4],[pointer.READY_CODE,pointer.READY_CODE+4])]
    sup[b]+=counts[b][z]/2*np.einsum('ac,bd->abcd',K,K.conj())
  wanted=np.zeros_like(sup);wanted[0,0,0,0,0]=wanted[1,1,1,1,1]=1
  err=float(np.max(abs(sup-wanted)));assert err<1e-9;pointer_error.append(err);superops[law]=sup
  effects[law]=[]
  for b in (0,1):
   effect=np.zeros((2,2),complex)
   for i in range(2):
    for j in range(2):effect[j,i]=sum(sup[b,a,a,i,j] for a in range(2))
   effects[law].append(effect)
  arrays['pointer_superoperator_'+law]=sup
 p=pi();basis_words=[[],[C.pulse('Y',p/4,[0])],[C.pulse('X',-p/4,[0])]]
 basis=[C.execute(w,1) for w in basis_words];allwords+=basis_words
 states=[np.array([1,0]),np.array([0,1]),np.array([1,1])/math.sqrt(2),np.array([1,1j])/math.sqrt(2)]
 spin_outcomes={};joint={}
 for law in ('C','E'):
  probs=[];joined=[]
  for B,axis in zip(basis,(2,0,1)):
   assert np.linalg.norm(B.conj().T@PAULI[2]@B-PAULI[axis],2)<1e-9
   row=[]
   for v in states:
    bv=B@v;vals=[float(np.vdot(bv,E@bv).real) for E in effects[law]]
    assert abs(sum(vals)-1)<1e-9;row.append(vals)
   probs.append(row)
  arrays['spin_probability_'+law]=np.array(probs);spin_outcomes[law]=probs
  Uf=preps[0];f=S.ray()/math.sqrt(24);plus=f.copy();plus[16:28]*=-1
  for state in ('coherent','dephased','opposite_shell_sign'):
   if state=='coherent':rho=np.outer(f,f)
   elif state=='dephased':rho=np.diag(f*f)
   else:rho=np.outer(plus,plus)
   pos=np.diag(Uf.conj().T@rho@Uf).real
   vals=np.outer(pos,[float(np.vdot(states[0],E@states[0]).real) for E in effects[law]])
   assert vals.min()>-1e-9 and abs(vals.sum()-1)<1e-8;joined.append(vals)
  arrays['joint_energy_ray_spin_'+law]=np.array(joined);joint[law]=[float(x[0].sum()) for x in joined]
  assert max(abs(a-b) for a,b in zip(joint[law],[1.,1/24,0.]))<1e-8
 # Spin2 central phase, with the same complete source-current reader.
 prepare=[C.pulse('Y',-p/4,[1])]
 controlled=[C.pulse('Z',-p/2,[0])]+C.pauli_rotation('ZZ',p/2)
 loops=[]
 for count in (0,1,2):
  word=prepare+controlled*count+C.inverse(prepare);allwords.append(word);U=C.execute(word,2);arrays[f'central_word_{count}']=U
  input=np.array([1,0,0,0],complex);out=U@input;vals=[]
  for law in ('C','E'):
   vals.append([float(np.vdot(out,np.kron(E,np.eye(2))@out).real) for E in effects[law]])
  assert max(abs(v[0]-(1 if count%2==0 else 0)) for v in vals)<1e-8;loops.append(vals)
 # All actually emitted native primitive matrices are independently certified.
 unique=sorted({(p['kind'],p['angle']) for word in allwords for p in word})
 for k,a in unique:assert syn.certify_primitive(k,a)
 lo,hi=iv.pi_interval();assert max(abs(F(math.pi)-lo),abs(F(math.pi)-hi))<F(1,iv.AQ)
 maxprep=max(map(len,allwords[:2]));geom_prep=F(maxprep+10,iv.AQ)
 # Conservative complete assay: two preparations/inverses, one spin word,
 # 16 native bit pointers, all reset and proper-code complements.
 pulses=4*maxprep+16*6+32;resets=32;e=F(2,10**11);gain=(1+e)**(2*pulses)
 timeout=configuration.budget(full=True,mode='reset',target=pointer.ready_state())['failure_upper']
 complete=(gain-1)/2+gain*(resets*timeout+4*geom_prep+max(geometric))
 assert complete<F(2,10**6)
 return {'preparation_pulses':list(map(len,allwords[:2])),'preparation_state_errors':prep_errors,'unique_native_primitives_certified':len(unique),'spin_cover_words':120,'spin_word_max_error':max(word_errors),'spin_word_geometric_upper':str(max(geometric)),'preparation_geometric_upper':str(geom_prep),'pointer_matrix_unit_error':max(pointer_error),'joint_energy_ray_clicks':joint,'central_reference_laws':loops,'assay_pulse_bound':pulses,'assay_reset_bound':resets,'reset_failure_upper':str(timeout),'complete_halfdiamond_upper':str(complete),'source_hardware_precision_claimed':False}, {'preparations':[prep(0),prep(1)],'spin_cover':spin_words,'spin_basis':basis_words,'central_control':controlled}
