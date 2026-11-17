"""Full coherent two-path source, complete poles, currents and records."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,time
import numpy as np
import mpmath as mp

ROOT=Path(__file__).parent;mp.mp.dps=65
started=time.time();p=mp.mpf(1)/12;rt2=mp.sqrt(2);I=mp.eye(3)
load=lambda name:json.loads((ROOT/'inputs'/name).read_text())
memory=load('memory_results.json');nu=load('neutral_results.json')
clock=load('common_clock.json');reg=load('register.json');old=load('chain_result.json')
plan_sha=hashlib.sha256((ROOT/'PLAN.md').read_bytes()).hexdigest()
phi=(1+mp.sqrt(5))/2
raw=[[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],
 [0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]]
ports=[mp.matrix(v)for v in raw]
rot=memory['architecture']['rotation'];orbit=[memory['architecture']['seed_port']]
for j in range(4):orbit.append(rot[orbit[-1]])
assert orbit==[1,5,11,10,7]
P=[ports[j]*ports[j].T/(ports[j].T*ports[j])[0]for j in orbit]
weights=[mp.mpf(12)**(4-j)/22621 for j in range(5)]
Y=sum((weights[j]*P[j]for j in range(5)),mp.zeros(3))
# The exact input is over Q(sqrt5), evaluated without Python eval.
import sympy as sy
given=mp.matrix([[mp.mpf(str(sy.N(sy.sympify(v),70)))for v in row]for row in nu['source_family_Y']])
assert mp.norm(Y-given)<mp.mpf('1e-60')
face=[4,5,11]
A=mp.matrix(3);B=mp.matrix(3)
for j in range(3):
 for i in range(3):A[i,j]=ports[face[j]][i];B[i,j]=ports[face[(j+1)%3]][i]
R=B*A**-1;N=reg['selected_rounds']
z=mp.exp(N*mp.log(1-mp.mpf(2)/9+mp.mpf(2)/9*mp.exp(1j/mp.mpf(N))))
charged=(I+(z*R+z.conjugate()*R.T)/rt2)/mp.sqrt(27)
charged_masses,Fe=mp.eighe(charged.H*charged)
assert abs(charged_masses[0]-mp.mpf(str(clock['register_432_midpoints'][1])))<mp.mpf('1e-17')
unitfrac=Fraction(clock['register_432_source_unit_GeV_exact'])*10**9
unit=mp.mpf(unitfrac.numerator)/unitfrac.denominator

def loc(stage,side,branch):return 3+12*stage+6*side+3*branch
def add(M,a,b,X):
 for i in range(3):
  for j in range(3):M[a+i,b+j]+=X[i,j]
def build(phase):
 M=mp.zeros(63);Vs=[I+(phase-1)*q for q in P]
 for k in range(5):
  for branch in range(2):
   r=loc(k,0,branch);l=loc(k,1,branch)
   add(M,r,l,I);add(M,l,r,I)
   V=I if branch==0 else Vs[k]
   if k==0:add(M,0,r,p*V/rt2);add(M,r,0,p*V.T/rt2)
   else:
    for parent in range(2):
     prev=loc(k-1,1,parent)
     add(M,prev,r,p*V/2);add(M,r,prev,p*V.T/2)
 for a in range(2):
  for b in range(2):add(M,loc(4,1,a),loc(4,1,b),Y/2)
 return M,Vs
def realify(M):
 n=M.rows;H=mp.zeros(2*n)
 for i in range(n):
  for j in range(n):
   a=mp.re(M[i,j]);b=mp.im(M[i,j])
   H[i,j]=a;H[i,n+j]=-b;H[n+i,j]=-b;H[n+i,n+j]=-a
 return H
def enc(M):return {'real':[[mp.nstr(mp.re(x),60)for x in row]for row in M.tolist()],
 'imag':[[mp.nstr(mp.im(x),60)for x in row]for row in M.tolist()]}
def maxabs(M):return max(abs(x)for x in M)

M,Vs=build(-1j);H=realify(M)
print('source assembled',time.time()-started,flush=True)
checkpoint=ROOT/'EIGENSYSTEM.json'
if checkpoint.exists():
 saved=json.loads(checkpoint.read_text())
 assert saved['plan_sha256']==plan_sha and saved['mp_dps']==mp.mp.dps
 values=mp.matrix([mp.mpf(x)for x in saved['values']])
 frame=mp.matrix([[mp.mpf(x)for x in row]for row in saved['frame']])
else:
 values,frame=mp.eigsy(H)
 checkpoint.write_text(json.dumps({'plan_sha256':plan_sha,'mp_dps':mp.mp.dps,
  'values':[mp.nstr(x,68)for x in values],
  'frame':[[mp.nstr(x,68)for x in row]for row in frame.tolist()]})+'\n')
print('full126 eigensystem',time.time()-started,flush=True)
assert values[62]<0<values[63]
positive=[values[63+j]for j in range(63)]
U=mp.matrix(63)
for i in range(63):
 for j in range(63):U[i,j]=frame[i,63+j]+1j*frame[63+i,63+j]
orth=maxabs(U.H*U-mp.eye(63));takagi=maxabs(U.T*M*U-mp.diag(positive))
assert orth<mp.mpf('1e-50')and takagi<mp.mpf('1e-50')
assert min(positive[3:])>mp.mpf(13)/30
current=Fe.H*U[:3,:]
complete=maxabs(current*current.H-I);assert complete<mp.mpf('1e-50')
light=current[:,:3];gram=light.H*light
gv,gf=mp.eighe(gram);polar=light*gf*mp.diag([1/mp.sqrt(x)for x in gv])*gf.H
angles=[abs(polar[0,1])**2/(1-abs(polar[0,2])**2),abs(polar[0,2])**2,
 abs(polar[1,2])**2/(1-abs(polar[0,2])**2)]
heavy=[sum(abs(current[i,j])**2 for j in range(3,63))for i in range(3)]

# Zero-frequency inverse uses all sixty sterile coordinates.
D=M[3:,3:];C=M[:3,3:];tail=D**-1*C.T
schur=-C*tail;Fprod=mp.eye(3)
for V in Vs:Fprod=Fprod*(I+V)/2
predicted=p**10*Fprod*Y*Fprod.T
schur_error=maxabs(schur-predicted);assert schur_error<mp.mpf('1e-60')
# The physical real doubled problem has a six-dimensional residue. It
# retains particle/antiparticle mixing and the entire sterile complement.
active=[0,1,2,63,64,65];heavyidx=[i for i in range(126)if i not in active]
Bc=mp.matrix([[H[i,j]for j in heavyidx]for i in active])
Dc=mp.matrix([[H[i,j]for j in heavyidx]for i in heavyidx])
tailreal=mp.zeros(120,6)
invD=D**-1
Dr=mp.zeros(120)
for i in range(60):
 for j in range(60):
  a=mp.re(invD[i,j]);b=mp.im(invD[i,j])
  Dr[i,j]=a;Dr[i,60+j]=b;Dr[60+i,j]=b;Dr[60+i,60+j]=-a
tailreal=Dr*Bc.T
assert maxabs(Dc*tailreal-Bc.T)<mp.mpf('1e-55')
Sreal=-Bc*tailreal;Z=mp.eye(6)+tailreal.T*tailreal
zz,zf=mp.eigsy(Z);iz=zf*mp.diag([1/mp.sqrt(x)for x in zz])*zf.T
leading,_=mp.eigsy(iz*Sreal*iz)
leading=sorted(abs(leading[i])for i in range(6))[::2]

# All spectral channels enter every retained response record.
records=[]
times=[mp.mpf(0),1/positive[2]**2,10/positive[2]**2]
for momentum in [mp.mpf(1),mp.mpf(10)]:
 energies=[mp.sqrt(momentum**2+x*x)for x in positive]
 weights=[(e+momentum)/(2*e)for e in energies]
 K=mp.diag([mp.sqrt(w)for w in weights])*current.H
 for t in times:
  phase=mp.diag([mp.exp(-1j*(e-energies[0])*t)for e in energies])
  for initial in range(3):
   prep=K[:,initial];returned=K.H*phase*prep
   success=mp.re((prep.H*prep)[0]);detected=[abs(x)**2 for x in returned]
   probs=detected+[1-success,success-sum(detected)]
   assert min(probs)>-mp.mpf('1e-50')and abs(sum(probs)-1)<mp.mpf('1e-50')
   records.append({'momentum':str(momentum),'source_time':mp.nstr(t,60),'initial':initial,
    'complete_probabilities':[mp.nstr(x,60)for x in probs]})

# At phase one the exact symmetric/antisymmetric branch transform splits
# the full source into the real33 chain and thirty unit Dirac modes.
M0,_=build(1);Q=mp.zeros(63);Q[0,0]=Q[1,1]=Q[2,2]=1
for k in range(5):
 for side in range(2):
  for fam in range(3):
   plus=3+6*k+3*side+fam;minus=33+6*k+3*side+fam
   a=loc(k,side,0)+fam;b=loc(k,side,1)+fam
   Q[a,plus]=Q[b,plus]=1/rt2;Q[a,minus]=1/rt2;Q[b,minus]=-1/rt2
control=Q.T*M0*Q;expected=mp.zeros(63)
for j in range(10):
 for fam in range(3):expected[3*j+fam,3*(j+1)+fam]=expected[3*(j+1)+fam,3*j+fam]=p if j%2==0 else 1
for i in range(3):
 for j in range(3):expected[30+i,30+j]=Y[i,j]
for k in range(5):
 for fam in range(3):expected[33+6*k+fam,36+6*k+fam]=expected[36+6*k+fam,33+6*k+fam]=1
control_error=maxabs(control-expected);assert control_error<mp.mpf('1e-60')

ev=[x*unit for x in positive[:3]]
result={'schema':'oph.complete_coherent_neutral_paths.v1','plan_sha256':plan_sha,'mp_dps':mp.mp.dps,
 'declared_source':{'active_Weyl_fields':3,'sterile_Weyl_fields':60,'all_real_doubled_bands':126,
  'stages':5,'p':'1/12','phase':'-i','splitter':'1/sqrt(2)','terminal':'unchanged real Y',
  'orbit':orbit,'empirical_neutrino_inputs':[],'native_compiler_or_interacting_clock_transfer':False},
 'full_positive_masses':[mp.nstr(x,60)for x in positive],
 'full_signed_doubled_spectrum':[mp.nstr(x,60)for x in values],
 'light_eV_conditional':[mp.nstr(x,60)for x in ev],
 'splittings_eV2_conditional':[mp.nstr(ev[1]**2-ev[0]**2,60),mp.nstr(ev[2]**2-ev[0]**2,60)],
 'splitting_ratio':mp.nstr((positive[1]**2-positive[0]**2)/(positive[2]**2-positive[0]**2),60),
 'polar_orientation_sin2_12_13_23':[mp.nstr(x,60)for x in angles],
 'full_current':enc(current),'light_current_Gram':enc(light*light.H),
 'all_heavy_current_weights':[mp.nstr(x,60)for x in heavy],
 'Takagi_residual':mp.nstr(takagi,20),'unitarity_residual':mp.nstr(orth,20),
 'full_current_completeness_residual':mp.nstr(complete,20),
 'heavy_absolute_gap_lower_exact':'13/30','zero_frequency_Schur':enc(schur),
 'zero_frequency_product_formula_error':mp.nstr(schur_error,20),
 'complete_real_residue_metric':enc(Z),'residue_corrected_light_roots':[mp.nstr(x,60)for x in leading],
 'zero_phase_exact_control_residual':mp.nstr(control_error,20),
 'zero_phase_control':'real33 chain plus30 sterile Weyl modes forming15 unit Dirac pairs, all retained',
 'complete_records':records,'physical_lepton_poles_claimed':False,
 'elapsed_seconds':time.time()-started}
(ROOT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'FULL_TAKAGI.json').write_text(json.dumps({'source':enc(M),'takagi':enc(U)},indent=2)+'\n')
print(json.dumps({k:result[k]for k in ['light_eV_conditional','splittings_eV2_conditional','splitting_ratio',
 'polar_orientation_sin2_12_13_23','all_heavy_current_weights','elapsed_seconds']},indent=2),flush=True)
