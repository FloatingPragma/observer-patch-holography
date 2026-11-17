"""Complete uncalibrated33-field neutral chain, exact small-root brackets and currents."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import mpmath as mp
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4];RER=ROOT/'reverse-engineering-reality'
paths=[HERE/'PLAN.md',RER/'model/openwave_neutrino_oscillation/results.json',
 RER/'model/native_flavour_transport/results.json',RER/'model/common_source/native_repair_mass_results.json',
 ROOT/'openwave-tuning/next/CHAIN3/species_reuse_connection/REGISTER.json',HERE.parent/'COMMON_CLOCK.json']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={str(p.relative_to(ROOT)):sha(p)for p in paths}
nu,transport,memory,reg,clock=[json.loads(p.read_text())for p in paths[1:]]
assert memory['architecture']['memory_states']==5
mp.mp.dps=90;p=F(1,12);pm=mp.mpf(1)/12
Y=s.Matrix(nu['source_family_Y']);x=s.Symbol('x');char=s.Poly(Y.charpoly(x).as_expr(),x,extension=s.sqrt(5))
assert char.all_coeffs()[:2]==[1,-1]
e2=F(str(char.all_coeffs()[2]));det=s.expand(Y.det())
aa=F(str(det.coeff(s.sqrt(5),0)));bb=F(str(det.coeff(s.sqrt(5),1)))
den=10**85;isq=int(mp.floor(mp.sqrt(5)*den));r5lo=F(isq,den);r5hi=F(isq+1,den)
assert r5lo*r5lo<5<r5hi*r5hi
dets=[aa+bb*r for r in [r5lo,r5hi]];dlo,dhi=min(dets),max(dets)
ym=mp.matrix([[mp.mpf(str(s.N(v,95)))for v in row]for row in Y.tolist()]);yv,yf=mp.eigsy(ym)
fmp=lambda q:mp.mpf(q.numerator)/q.denominator
brackets=[]
for j,z in enumerate(yv):
    cen=F(mp.nstr(z,80));eps=F(1,10**72);lo,hi=cen-eps,cen+eps
    evals=lambda q:(q**3-q*q+e2*q-dhi,q**3-q*q+e2*q-dlo)
    vl,vh=evals(lo),evals(hi)
    if j in [0,2]:assert vl[1]<0<vh[0]
    else:assert vl[0]>0>vh[1]
    brackets.append((lo,hi))

def chain_det(t,y):
    a,b=F(1),-t
    for row in range(1,11):
        link=p if (row-1)%2==0 else F(1)
        diagonal=y if row==10 else F(0)
        a,b=b,(diagonal-t)*b-link*link*a
    return b

# Exact zero-frequency Schur coefficient and its kinetic residue, before
# applying any low-mode approximation or normalizing away sterile weight.
ys=s.Symbol('y',real=True);heavy=s.zeros(10)
for q in range(9):heavy[q,q+1]=heavy[q+1,q]=1 if q%2==0 else s.Rational(1,12)
heavy[9,9]=ys;c=s.zeros(10,1);c[0]=s.Rational(1,12)
tail=heavy.inv()*c
z0=sum(p**(2*q)for q in range(6));z1=p**10*sum(p**(2*q)for q in range(5))
schur=s.factor(-(c.T*tail)[0]);residue=s.factor(1+(tail.T*tail)[0])
assert s.simplify(schur-s.Rational(p.numerator,p.denominator)**10*ys)==0
assert s.simplify(residue-s.Rational(z0.numerator,z0.denominator)-s.Rational(z1.numerator,z1.denominator)*ys**2)==0

families=[];all_modes=[];spectra=[]
for j,y in enumerate(yv):
    H=mp.zeros(11)
    for q in range(10):H[q,q+1]=H[q+1,q]=pm if q%2==0 else 1
    H[10,10]=y
    vals,vec=mp.eigsy(H);spectra.append((vals,vec))
    lam=vals[5];assert lam>0 and abs(vals[4])>.43 and vals[6]>.43
    center=F(mp.nstr(lam,80));eps=F(1,10**72);lo,hi=center-eps,center+eps
    assert max(chain_det(lo,yy)for yy in brackets[j])<0
    assert min(chain_det(hi,yy)for yy in brackets[j])>0
    norm0=1/sum(pm**(2*q)for q in range(6))
    family=dict(family=j,Y_eigenvalue=mp.nstr(y,60),Y_interval_exact=[str(z)for z in brackets[j]],
        light_mass=mp.nstr(lam,60),light_interval_exact=[str(lo),str(hi)],
        leading_light_mass=mp.nstr(pm**10*y*norm0,60),
        relative_error_from_leading=mp.nstr(lam/(pm**10*y*norm0)-1,30),
        active_light_weight=mp.nstr(vec[0,5]**2,60),
        full_signed_spectrum=[mp.nstr(z,60)for z in vals],
        full_active_spectral_weights=[mp.nstr(vec[0,k]**2,60)for k in range(11)])
    families.append(family)
    for q in range(11):all_modes.append((abs(vals[q]),j,q,vals[q],vec[0,q]))
all_modes.sort(key=lambda q:q[0])
assert [r[1]for r in all_modes[:3]]==[0,1,2]

# Existing charged frame, fixed mass ordering; no permutation fit.
phi=(1+s.sqrt(5))/2
raw=s.Matrix([[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],
 [0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]])
face=[4,5,11];A=s.Matrix.hstack(*(raw[i,:].T for i in face));B=s.Matrix.hstack(*(raw[i,:].T for i in face[1:]+face[:1]))
rn=np.array((B*A.inv()).evalf(25),complex)
N=reg['selected_rounds'];z=complex(mp.exp(N*mp.log(1-mp.mpf(2)/9+(mp.mpf(2)/9)*mp.exp(1j/N))))
charged=(np.eye(3)+(z*rn+z.conjugate()*rn.T)/np.sqrt(2))/np.sqrt(27)
me,Fe=np.linalg.eigh(charged.conj().T@charged)
assert abs(me[0]-clock['register_432_midpoints'][1])<1e-14
fn=np.array(yf.tolist(),float)
unit=F(clock['register_432_source_unit_GeV_exact'])*10**9
ev=[fmp(unit)*mp.mpf(q['light_mass'])for q in families]
records=[]
for name,V in [('existing_primary',s.Matrix(transport['primary']['transport'])),('real_control',s.eye(3))]:
    v=np.array(V.evalf(30),complex);W=v@fn;C=Fe.conj().T@W
    current=np.column_stack([C[:,j]*complex(a)*(1j if signed<0 else 1)for mass,j,q,signed,a in all_modes])
    assert np.max(abs(current@current.conj().T-np.eye(3)))<1e-12
    # Full33 Takagi reconstruction, using alternating family frames on links.
    full=np.zeros((33,33),complex)
    for q in range(10):full[3*q:3*q+3,3*q+3:3*q+6]=(float(p)if q%2==0 else 1)*np.eye(3)
    full=full+full.T;full[30:33,30:33]=v.conjugate()@np.array(ym.tolist(),float)@v.conj().T
    takagi=np.empty((33,33),complex)
    for col,(mass,j,q,signed,a)in enumerate(all_modes):
        vals,vec=spectra[j]
        for block in range(11):
            wb=W if block%2==0 else W.conjugate()
            takagi[3*block:3*block+3,col]=wb[:,j]*float(vec[block,q])*(1j if signed<0 else 1)
    err=float(np.max(abs(takagi.T@full@takagi-np.diag([float(r[0])for r in all_modes]))))
    assert err<1e-12
    assert np.max(abs(Fe.conj().T@takagi[:3,:]-current))<1e-12
    Bc=abs(C)**2;sin13=Bc[0,2]
    source_times=[mp.mpf(0),1/all_modes[2][0]**2,10/all_modes[2][0]**2]
    probs=[]
    for mom in [mp.mpf(1),mp.mpf(10)]:
        energies=[mp.sqrt(mom*mom+r[0]*r[0])for r in all_modes]
        weights=np.array([float((e+mom)/(2*e))for e in energies])
        K=np.sqrt(weights)[:,None]*current.conj().T
        for t in source_times:
            evolution=np.array([complex(mp.exp(-1j*(e-energies[0])*t))for e in energies])
            for initial in range(3):
                prepared=K[:,initial];returned=K.conj().T@(evolution*prepared)
                success=float(np.vdot(prepared,prepared).real)
                pp=list(abs(returned)**2)+[1-success,success-float(np.vdot(returned,returned).real)]
                assert min(pp)>-1e-12 and abs(sum(pp)-1)<1e-12
                probs.append(dict(momentum=str(mom),source_time=mp.nstr(t,50),initial=initial,
                                  complete_probabilities=[float(x)for x in pp]))
    records.append(dict(candidate=name,full_Takagi_reconstruction_error=err,
       full_current_completeness_error=float(np.max(abs(current@current.conj().T-np.eye(3)))),
       actual_current_real=current.real.tolist(),actual_current_imag=current.imag.tolist(),
       light_current_Gram=(current[:,:3]@current[:,:3].conj().T).real.tolist(),
       total_heavy_current_weights=np.sum(abs(current[:,3:])**2,axis=1).tolist(),
       unitary_orientation_sin2_12_13_23=[float(Bc[0,1]/(1-sin13)),float(sin13),float(Bc[1,2]/(1-sin13))],
       standard_three_state_PMNS_unitarity_exact=False,complete_records=probs))

out=dict(schema='oph.declared_neutral_five_stage_chain.v1',input_sha256=pins,
    architecture={'active_fields':3,'sterile_fields':30,'stages':5,'small_link_exact':'1/12',
       'unit_link_exact':'1','terminal_majorana':'existing V*YV†, primary native word and identity control',
       'source_change':'new coherent sterile chain; not the recorded reset instrument',
       'same_bare_kinetic_clock_assumed':True,'native_compiler_executed':False,
       'selected_after_hierarchy_failure':True,'empirical_neutrino_inputs':[]},
    families=families,heavy_absolute_gap_lower_exact='13/30',
    exact_zero_frequency_Schur=str(schur),exact_linearized_residue_metric=str(residue),
    canonical_small_root_requires_residue_metric=True,
    leading_zero_mode_normalization_exact=str(sum(p**(2*q)for q in range(6))),
    leading_suppression_exact=str(p**10/sum(p**(2*q)for q in range(6))),
    electron_anchored_conditional_light_eV=[mp.nstr(e,50)for e in ev],
    electron_anchored_conditional_splittings_eV2=[mp.nstr(ev[1]**2-ev[0]**2,50),mp.nstr(ev[2]**2-ev[0]**2,50)],
    mass_sum_conditional_eV=mp.nstr(sum(ev),50),
    modes_order=[dict(family=j,scalar_index=q,signed_mass=mp.nstr(signed,50),Takagi_mass=mp.nstr(mass,50))for mass,j,q,signed,a in all_modes],
    candidates=records,light_spectrum_mechanism_constructed=True,
    measured_neutrino_spectrum_and_mixing_solved=False,physical_clock_and_interacting_source_transfer=False)
assert pins=={str(q.relative_to(ROOT)):sha(q)for q in paths}
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['leading_suppression_exact','electron_anchored_conditional_light_eV','electron_anchored_conditional_splittings_eV2','mass_sum_conditional_eV']},indent=2))
for r in records:print(r['candidate'],r['unitary_orientation_sin2_12_13_23'],r['total_heavy_current_weights'])
