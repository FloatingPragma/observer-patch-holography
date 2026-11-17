"""Fixed independent charged/neutral-frame candidates; no target search."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib
import json
import numpy as np
import sympy as s
import mpmath as mp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];RER=ROOT/'reverse-engineering-reality'
paths=[HERE/'PLAN.md',RER/'model/lepton_register_round/results.json',
 RER/'model/openwave_neutrino_oscillation/results.json',RER/'model/native_flavour_transport/results.json',
 ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection/REGISTER.json',
 ROOT/'openwave-push/rowfull/ROW_10_symmetry/ANGLES.md',RER/'model/common_source/slow_geometry.py']
pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
reg=json.loads(paths[1].read_text());nu=json.loads(paths[2].read_text());word=json.loads(paths[3].read_text());finite=json.loads(paths[4].read_text())
phi=(1+s.sqrt(5))/2
raw=s.Matrix([[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],
              [0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]])
face=reg['register']['geometry']['outward_order'];assert face==[4,5,11]
A=s.Matrix.hstack(*(raw[i,:].T for i in face));B=s.Matrix.hstack(*(raw[i,:].T for i in face[1:]+face[:1]))
clean=lambda M:M.applyfunc(s.simplify)
R=clean(B*A.inv());assert clean(R.T*R-s.eye(3))==s.zeros(3) and clean(R**3-s.eye(3))==s.zeros(3)
P=[clean(raw[i,:].T*raw[i,:]/(raw[i,:]*raw[i,:].T)[0]) for i in range(12)]
V=[s.eye(3)+(-s.I-1)*p for p in P]
primary=clean(V[10]*V[3]*V[2]);assert clean(primary-s.Matrix(word['primary']['transport']))==s.zeros(3)
loop=clean(V[3]*V[2]*V[3].H*V[2].H)
Y=s.Matrix(nu['source_family_Y']);yn=np.array(Y.evalf(30),complex)
mn,Fn=np.linalg.eigh(yn)
mp.mp.dps=60;N=finite['selected_rounds'];p=mp.mpf(2)/9
z=complex(mp.exp(N*mp.log(1-p+p*mp.exp(1j/mp.mpf(N)))))
rn=np.array(R.evalf(30),complex)
M=(np.eye(3)+(z*rn+z.conjugate()*rn.T)/np.sqrt(2))/np.sqrt(27)
Ye=M.conj().T@M;me,Fe=np.linalg.eigh(Ye)
expected=sorted(float(F(x)) for x in finite['mass_midpoints_exact'])
assert np.max(abs(me-expected))<1e-14
rows=[]
for name,vex in [('real_memory',s.eye(3)),('existing_primary',primary),('closed_commutator',loop),('inverse_loop_control',loop.H)]:
    assert clean(vex.H*vex-s.eye(3))==s.zeros(3)
    v=np.array(vex.evalf(30),complex);C=Fe.conj().T@v@Fn;mod=abs(C)**2
    Mnu=v.conjugate()@yn@v.conj().T
    assert np.max(abs((v@Fn).T@Mnu@(v@Fn)-np.diag(mn)))<1e-12
    z13=float(mod[0,2]);angles=[float(mod[0,1]/(1-z13)),z13,float(mod[1,2]/(1-z13))]
    probs=[]
    for momentum,time,initial in product((1.,10.),(0.,1.,10.),range(3)):
        E=np.hypot(momentum,mn);w=(E+momentum)/(2*E)
        K=np.diag(np.sqrt(w))@C.conj().T;U=np.diag(np.exp(-1j*(E-E[0])*time))
        state=K.conj().T@U@K[:,initial];prob=list(abs(state)**2)
        success=float(np.vdot(K[:,initial],K[:,initial]).real)
        prob +=[1-success,success-sum(prob)]
        assert min(prob)>-1e-13 and abs(sum(prob)-1)<1e-12
        probs.append(dict(momentum=momentum,time=time,initial=initial,complete_probabilities=prob))
    rows.append(dict(name=name,angles_sin2_12_13_23=angles,moduli_squared=mod.tolist(),
        J=float(np.imag(C[0,0]*C[1,1]*C[0,1].conjugate()*C[1,0].conjugate())),
        electron_muon_modulus_row_difference=float(np.max(abs(mod[0]-mod[1]))),
        complete_records=probs,charged_mass_coefficients=me.tolist(),neutral_mass_coefficients=mn.tolist(),
        minimum_neutral_over_electron=float(mn[0]/me[0])))
assert rows[0]['electron_muon_modulus_row_difference']<1e-12
comparison=np.array([.308,.02215,.470]);lower=np.array([.011,.00058,.013]);upper=np.array([.012,.00056,.017])
for row in rows:
    diff=np.array(row['angles_sin2_12_13_23'])-comparison
    row['comparison_pulls']=(diff/np.where(diff>=0,upper,lower)).tolist()
    row['all_three_inside_3sigma']=bool(np.max(abs(np.array(row['comparison_pulls'])))<=3)
out=dict(schema='oph.independent_neutral_frame_candidates.v1',input_sha256=pins,
    charged_frame='inherited face(4,5,11) outward-third-turn eigenframe of actual finite432 analyzer effect',
    neutral_frame='native five-state memory Y, increasing eigenvalue order',
    target_selected_parameters=False,all_declared_candidates_retained=True,
    neutral_scale_input=False,unscaled_neutrals_are_not_light=True,
    real_memory_obstruction='The electron and muon are the conjugate charged characters. Every real neutral eigenframe gives identical electron/muon modulus rows; a complex neutral source is necessary in this join.',
    candidates=rows,
    next_obligations=['independently generated neutral mass scale','dynamics giving observed mixing','common interacting source attachment'])
for path in paths:assert hashlib.sha256(path.read_bytes()).hexdigest()==pins[str(path.relative_to(ROOT))]
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([{k:r[k] for k in ('name','angles_sin2_12_13_23','comparison_pulls','all_three_inside_3sigma','minimum_neutral_over_electron')} for r in rows],indent=2))
