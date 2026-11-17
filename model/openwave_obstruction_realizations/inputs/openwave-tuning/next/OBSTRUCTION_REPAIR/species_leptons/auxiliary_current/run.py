"""Full auxiliary dressed currents, exact CAR moments, all spectral complements."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import hashlib
import json
import sys
import numpy as np
import sympy as sy

HERE=Path(__file__).resolve().parent;INPUT=HERE/'inputs';sys.path.insert(0,str(INPUT))
import local_auxiliary_interface as srcmod
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pin=json.loads((HERE/'INPUTS.json').read_text())
for name,h in pin['snapshots_sha256'].items():assert sha(INPUT/name)==h
assert sha(INPUT/'LOCAL_AUXILIARY_INTERFACE_MANIFEST.json')==pin['source_manifest_sha256']
s=srcmod.source();words=s['words'];ix={w:i for i,w in enumerate(words)};R=s['singlet_basis'];Ph=s['half_projector'];Pt=s['threehalf_projector']
assert len(words)==220

# Independent exact Gaussian-integer one-body lift, in the same frozen CAR order.
def lift(h):
    re=np.zeros((220,220),dtype=np.int64);im=np.zeros_like(re)
    nz=list(zip(*np.nonzero(h)))
    for j,w in enumerate(words):
        for a,b in nz:
            a,b=int(a),int(b)
            if not (w>>b)&1:continue
            z=w^(1<<b);sign=(-1)**((w&((1<<b)-1)).bit_count())
            if (z>>a)&1:continue
            sign*=(-1)**((z&((1<<a)-1)).bit_count());y=z|(1<<a)
            val=h[a,b];assert val.real==round(val.real) and val.imag==round(val.imag)
            re[ix[y],j]+=sign*int(val.real);im[ix[y],j]+=sign*int(val.imag)
    return re,im
def arr(a):return a[0]+1j*a[1]
def act(a,v):
    rr=a[0]@v[0]-a[1]@v[1];ii=a[0]@v[1]+a[1]@v[0]
    assert max(np.max(abs(rr)),np.max(abs(ii)))<10**7
    return rr,ii
def inner(a,b):
    # Entries are bounded above; all dot products are exact int64 sums.
    return int(a[0]@b[0]+a[1]@b[1]),int(a[0]@b[1]-a[1]@b[0])
pauli=[np.eye(2,dtype=complex),np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]).astype(complex)]
tau=np.array([[0,1],[0,0]],complex)
curr=[lift(np.kron(np.kron(a,tau),np.eye(3)))for a in pauli]
for got,want in zip(curr,[s['V'],*s['A']]):assert np.max(abs(arr(got)-want.toarray()))==0
spinminus=lift(np.kron(np.kron(tau.T,np.eye(2)),np.eye(3)))
spinz2=lift(np.kron(np.kron(pauli[3],np.eye(2)),np.eye(3)))
isospinz2=lift(np.kron(np.kron(np.eye(2),pauli[3]),np.eye(3)))

# Produce exact highest-weight vector from the source singlet sector, then
# verify its full CAR constraints rather than treating float selection as proof.
Sz=R.T@s['spin'][2]@R;Iz=R.T@((s['V']@s['V'].conj().T-s['V'].conj().T@s['V'])/2)@R
projector=Ph@(.5*np.eye(20)+Sz)@(.5*np.eye(20)-Iz)
vals,vec=np.linalg.eigh((projector+projector.conj().T)/2)
assert abs(vals[-1]-1)<1e-12 and max(abs(vals[:-1]))<1e-12
trial=R@vec[:,-1];assert max(abs(trial.imag))<1e-12
trial=trial.real;pivot=int(np.argmax(abs(trial)))
rat=[F(float(v/trial[pivot])).limit_denominator(128)for v in trial]
den=lcm(*(v.denominator for v in rat));qn=np.array([int(v*den)for v in rat],dtype=np.int64)
zero=np.zeros(220,dtype=np.int64);nup=(qn,zero);norm=inner(nup,nup)[0]
assert max(abs(qn))<100 and norm>0
qp=act(curr[0],nup);nd=act(spinminus,nup);pd=act(spinminus,qp);quad=[nup,nd,qp,pd]
assert np.array_equal(act(spinz2,nup)[0],qn) and np.array_equal(act(isospinz2,nup)[0],-qn)
assert not np.any(act((spinminus[0].T,-spinminus[1].T),nup)[0])
assert not np.any(act((curr[0][0].T,-curr[0][1].T),nup)[0])
for a in quad:assert inner(a,a)==(norm,0)
base=np.column_stack([arr(v)/np.sqrt(norm)for v in quad])
assert np.max(abs(base.conj().T@base-np.eye(4)))<1e-13

# Integer-scaled colour-spin generators: 4G except a7, where 4sqrt3 G.
T=[]
for a,b in [(0,1),(0,2),(1,2)]:
    x=np.zeros((3,3),complex);x[a,b]=x[b,a]=1;T.append(x)
    y=np.zeros((3,3),complex);y[a,b]=-1j;y[b,a]=1j;T.append(y)
T.extend([np.diag([1,-1,0]).astype(complex),np.diag([1,1,-2]).astype(complex)])
G=[];gden=[]
for a,t in enumerate(T):
    colour=lift(np.kron(np.eye(4),t))
    for v in quad:assert inner(act(colour,v),act(colour,v))==(0,0)
    for i in range(3):
        z=lift(np.kron(np.kron(pauli[i+1],np.eye(2)),t));G.append(z);gden.append(48 if a==7 else 16)
        assert np.max(abs(arr(z)/(4*np.sqrt(3)if a==7 else 4)-s['G'][3*a+i].toarray()))<1e-14

bare=[];excited=[];bare_norm=[];excited_norm=[]
for mu,J in enumerate(curr):
    bare.append([[str(F(inner(v,act(J,w))[0],norm))+'+'+str(F(inner(v,act(J,w))[1],norm))+'i' for w in quad[:2]]for v in quad[2:]])
    bnr=F(inner(act(J,nup),act(J,nup))[0],norm);bare_norm.append(bnr)
    exnr=F();exax=F();exay=F()
    for g,d in zip(G,gden):
        gn=act(g,nup);gp=act(g,qp);jgn=act(J,gn)
        exnr+=F(inner(jgn,jgn)[0],4*norm*d)
        rr,ii=inner(gp,jgn);exax+=F(rr,4*norm*d);exay+=F(ii,4*norm*d)
    excited_norm.append(exnr);excited.append((exax,exay))
assert F(inner(qp,act(curr[0],nup))[0],norm)==1
assert F(inner(qp,act(curr[3],nup))[0],norm)==F(5,3)
assert excited[0]==(F(1),F()) and excited[3]==(F(10,9),F())
assert sum(F(inner(act(g,nup),act(g,nup))[0],norm*d)for g,d in zip(G,gden))==4

# Exact Lie algebra argument behind the current compression.
sym=lambda q:sy.Matrix([[int(z.real)+sy.I*int(z.imag)for z in row]for row in q])
S=[sym(q)/2 for q in pauli[1:]]
exactT=[sym(t)/(2*sy.sqrt(3)if a==7 else 2)for a,t in enumerate(T)]
assert sy.simplify(sum((t*t for t in exactT),sy.zeros(3))-sy.Rational(4,3)*sy.eye(3))==sy.zeros(3)
for A in [sym(q)for q in pauli[1:]]:
    dc=sum((h*(h*A-A*h)-(h*A-A*h)*h for h in S),sy.zeros(2))
    assert (dc-2*A).applyfunc(sy.simplify)==sy.zeros(2)

cases=[]
for kappa in [F(1,2),F(1),F(2)]:
    Delta=3*kappa;disc=sy.sqrt(sy.Rational(Delta.numerator,Delta.denominator)**2+sy.Rational(1,4))
    sn2=sy.simplify((1-sy.Rational(Delta.numerator,Delta.denominator)/disc)/2)
    ga=sy.simplify(sy.Rational(5,3)*(1-sn2/3));sn=float(sn2)
    Dl,EH=srcmod.lower_isometry(s,float(kappa));Du=np.zeros_like(Dl)
    for lam,Pj in [(4,Ph),(2,Pt)]:
        di=np.sqrt(float(Delta)**2+lam/16);cc=np.sqrt((1+float(Delta)/di)/2);ss=np.sqrt((1-float(Delta)/di)/2)
        Du[0]+=ss*R@Pj
        for j,g in enumerate(s['G']):Du[j+1]+=cc/np.sqrt(lam)*(g@R@Pj)
    low=Dl.reshape(5500,20);high=Du.reshape(5500,20)
    assert np.max(abs(low.conj().T@high))<1e-13
    assert np.max(abs(high.conj().T@high-np.eye(20)))<1e-13
    coef=R.T@base;D=(low@coef).reshape(25,220,4)
    compressed=[];rows=[];maxgauss=0.;maxspin=0.;maxH=0.
    for j in range(4):
        dv=D[:,:,j];eng=(float(Delta)-np.sqrt(float(Delta)**2+.25))/2
        maxH=max(maxH,float(np.linalg.norm(srcmod.hamiltonian_action(s,float(kappa),dv)-eng*dv)))
        for a in range(8):
            maxgauss=max(maxgauss,float(np.linalg.norm(srcmod.total_generator(s,s['colour'][a],np.kron(s['adjoint_colour'][a],np.eye(3)),dv))))
        z=srcmod.total_generator(s,s['spin'][2],np.kron(np.eye(8),s['vector_spin'][2]),dv)
        maxspin=max(maxspin,float(np.linalg.norm(z-(.5 if j%2==0 else -.5)*dv)))
    for mu,J in enumerate([s['V'],*s['A']]):
        applied=np.stack([J@D[a,:,:2]for a in range(25)])
        block=D[:,:,2:].reshape(5500,2).conj().T@applied.reshape(5500,2)
        want=(1 if mu==0 else float(ga))*pauli[mu]
        assert np.max(abs(block-want))<1e-12;compressed.append(block)
        for initial in range(2):
            w=applied[:,:,initial].reshape(5500);lc=low.conj().T@w;uc=high.conj().T@w
            total=float(np.vdot(w,w).real);ln=float(np.vdot(lc,lc).real);un=float(np.vdot(uc,uc).real)
            dark=w-low@lc-high@uc
            assert abs(np.vdot(dark,dark).real-(total-ln-un))<1e-12
            assert np.linalg.norm(srcmod.hamiltonian_action(s,float(kappa),dark.reshape(25,220))-float(Delta)*dark.reshape(25,220))<1e-11
            probs=[float(np.vdot(Ph@lc,Ph@lc).real)/9,float(np.vdot(Pt@lc,Pt@lc).real)/9,
                   float(np.vdot(Ph@uc,Ph@uc).real)/9,float(np.vdot(Pt@uc,Pt@uc).real)/9,
                   float(np.vdot(dark,dark).real)/9,1-total/9]
            assert min(probs)>-1e-12 and abs(sum(probs)-1)<1e-12
            expectednorm=(1-sn)*float(bare_norm[mu])+sn*float(excited_norm[mu])
            assert abs(total-expectednorm)<1e-12
            rows.append(dict(current=mu,initial_spin=initial,complete_probabilities=probs,
                   total_current_norm_squared=total,low_half_coherent_amplitude_squared=float((1 if mu==0 else ga**2)),
                   all_complement_current_weight=total-(1 if mu==0 else float(ga)**2)))
    assert max(maxgauss,maxspin,maxH)<1e-11
    cases.append(dict(kappa=str(kappa),Delta=str(Delta),g='1/8',auxiliary_occupation_exact=str(sn2),
       auxiliary_occupation=float(sn2),vector_exact='1',axial_exact=str(ga),axial=float(ga),
       H_residual=maxH,total_colour_Gauss_residual=maxgauss,total_spin_z_residual=maxspin,
       current_blocks_real=np.array(compressed).real.tolist(),current_blocks_imag=np.array(compressed).imag.tolist(),
       complete_current_records=rows))

out=dict(schema='oph.full_auxiliary_dressed_current.v1',source_snapshots=pin,
    quark_CAR_sector_dimension=220,auxiliary_dimension=25,full_carrier_dimension=5500,
    physical_Gauss_sector_dimension=80,coupled_lower_dimension=20,coupled_upper_dimension=20,dark_dimension=40,
    exact_quartet_integer_norm=norm,exact_nup_integer_vector=qn.tolist(),
    exact_bare_current_blocks=bare,exact_auxiliary_sector_nup_to_pup=[list(map(str,x))for x in excited],
    exact_bare_current_norm_squared=list(map(str,bare_norm)),exact_auxiliary_current_norm_squared=list(map(str,excited_norm)),
    exact_double_commutator='sum_ai[G_ai,[G_ai,A]]=8A/3',
    axial_formula='(5/3)(1-s²/3), s²=(1-Delta/sqrt(Delta²+16g²))/2',
    positive_gap_all_coupling_axial_lower_bound_exact='25/18',
    instrument='K=J/3; resolve lowerS1/2,lowerS3/2,upperS1/2,upperS3/2,dark; retain I-K†K failure',
    cases=cases,empirical_axial_fit=False,physical_axial_target_solved=False,
    joined_to_two_cell_interacting_vacuum=False,Dirac_current_attachment=False,
    source_scope='new local quark+adjoint/vector auxiliary, every colour/spin complement retained')
for name,h in pin['snapshots_sha256'].items():assert sha(INPUT/name)==h
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k]for k in ['exact_bare_current_norm_squared','exact_auxiliary_current_norm_squared','positive_gap_all_coupling_axial_lower_bound_exact']},indent=2))
for r in cases:print({k:r[k]for k in ['kappa','auxiliary_occupation','axial','H_residual','total_colour_Gauss_residual','total_spin_z_residual']})
