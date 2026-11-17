"""Complete recorded five-memory response, actual Schur complement and coherent correction."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import numpy as np
import sympy as s
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[4]
RER=ROOT/'reverse-engineering-reality'
paths=[HERE/'PLAN.md',RER/'model/common_source/native_repair_mass.py',
       RER/'model/common_source/native_repair_mass_results.json',
       RER/'model/openwave_neutrino_oscillation/results.json',
       HERE.parent/'COMMON_CLOCK.json']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins={str(p.relative_to(ROOT)):sha(p) for p in paths}
receipt=json.loads(paths[2].read_text());nu=json.loads(paths[3].read_text())
assert receipt['architecture']['rotation']==[0,5,6,3,2,11,8,1,9,4,7,10]
assert receipt['architecture']['seed_port']==1
ports=receipt['native_feedback_source']['fixed_native_ports'];assert ports==[1,5,11,10,7]
p=s.Rational(1,12);pi=s.Matrix([(1-p)*p**j/(1-p**5) for j in range(5)])
S=s.zeros(5)
for j in range(5):S[(j+1)%5,j]=1
reset=s.zeros(5);reset[0,:]=s.ones(1,5)
T=p*S+(1-p)*reset
assert T*pi==pi and s.ones(1,5)*T==s.ones(1,5)
phi=(1+s.sqrt(5))/2
raw=s.Matrix([[-1,phi,0],[1,phi,0],[-1,-phi,0],[1,-phi,0],[0,-1,phi],[0,1,phi],
              [0,-1,-phi],[0,1,-phi],[phi,0,-1],[phi,0,1],[-phi,0,-1],[-phi,0,1]])
clean=lambda M:M.applyfunc(s.simplify)
P=[clean(raw[j,:].T*raw[j,:]/(raw[j,:]*raw[j,:].T)[0]) for j in ports]
Y=clean(sum((pi[j]*P[j] for j in range(5)),s.zeros(3)))
assert clean(Y-s.Matrix(nu['source_family_Y']))==s.zeros(3)
# Exact memory coordinates: marginal first, then zero-sum deviations.
memory_basis=s.Matrix.hstack(pi,*(s.eye(5)[:,j]-s.eye(5)[:,0] for j in range(1,5)))
mbi=memory_basis.inv();tb=mbi*T*memory_basis
assert tb[0,0]==1 and tb[0,1:]==s.zeros(1,4) and tb[1:,0]==s.zeros(4,1)
x=s.Symbol('x');pol=s.factor(tb[1:,1:].charpoly(x).as_expr())
assert s.expand(pol-(x**4+p*x**3+p*p*x*x+p**3*x+p**4))==0

tn=np.array(T,float);pin=np.array(pi,float).reshape(5);yn=np.array(Y.evalf(30),float)
pn=np.array([np.array(q.evalf(30),float) for q in P])
X=np.array([[0,1],[1,0]],complex);I6=np.eye(6);I36=np.eye(36)
H=np.array([np.kron(X,q) for q in pn])
L=np.array([-1j*(np.kron(q,I6)-np.kron(I6,q.T)) for q in H])
M0=np.kron(tn,I36)
Ldiag=np.zeros((180,180),complex)
for j in range(5):Ldiag[j*36:(j+1)*36,j*36:(j+1)*36]=L[j]
B=np.kron(np.array(memory_basis,float),I36);Bi=np.kron(np.array(mbi,float),I36)
M1=Bi@Ldiag@M0@B
M2=Bi@Ldiag@Ldiag@M0@B/2
D0=np.kron(np.array(tb[1:,1:],float),I36)
A=M1[:36,:36]
K2=M2[:36,:36]+M1[:36,36:]@np.linalg.solve(np.eye(144)-D0,M1[36:,:36])-A@A/2
HY=np.kron(X,yn)
expected=-1j*(np.kron(HY,I6)-np.kron(I6,HY.T))
assert np.max(abs(A-expected))<1e-14
anti=(K2-K2.conj().T)/2
# Partial trace of -i[H,.] yields -i*d*H for traceless H (row-major vec).
ha=1j*np.einsum('abcb->ac',anti.reshape(6,6,6,6))/6
assert np.max(abs(ha-ha.conj().T))<1e-13
recon=-1j*(np.kron(ha,I6)-np.kron(I6,ha.T))
hamerr=float(np.max(abs(anti-recon)))
assert hamerr<1e-13
assert np.max(abs(ha.imag))>1e-8 and np.max(abs(ha.real))<1e-13
# Each memory H_j=X tensor P_j; commutators live on I tensor imaginary skew-family.
ha3=(ha[:3,:3]+ha[3:,3:])/2
assert np.max(abs(ha-np.kron(np.eye(2),ha3)))<1e-13
mass,Fy=np.linalg.eigh(yn)
assert np.max(abs(np.diag(Fy.T@ha3@Fy)))<1e-13

records=[]
for step in [1/8,1/64,1/512]:
    U=np.array([I6+(np.cos(step)-1)*np.kron(np.eye(2),q)-1j*np.sin(step)*np.kron(X,q) for q in pn])
    M=np.block([[tn[new,old]*np.kron(U[new],U[new].conj()) for old in range(5)]for new in range(5)])
    assert np.max(abs(M-expm(step*Ldiag)@M0))<1e-14
    mb=Bi@M@B
    eig,V=np.linalg.eig(M);order=np.lexsort((np.angle(eig),-abs(eig)));eig=eig[order];V=V[:,order]
    invV=np.linalg.inv(V)
    marginal=np.kron(np.ones((1,5)),I36);injection=np.kron(pin[:,None],I36)
    right=marginal@V;left=invV@injection
    residues=[np.outer(right[:,j],left[j]) for j in range(180)]
    assert np.max(abs(sum(residues)-I36))<1e-10
    schurerr=[];spectralerr=[]
    for z in [1.1+.2j,.8+.3j,-.3+.2j]:
        full=marginal@np.linalg.solve(z*np.eye(180)-M,injection)
        sc=z*I36-mb[:36,:36]-mb[:36,36:]@np.linalg.solve(z*np.eye(144)-mb[36:,36:],mb[36:,:36])
        schurerr.append(float(np.max(abs(full-np.linalg.inv(sc)))))
        spectralerr.append(float(np.max(abs(full-sum(r/(z-e)for r,e in zip(residues,eig))))))
    assert max(schurerr)<1e-11 and max(spectralerr)<1e-10
    slow=abs(eig)>.5
    assert slow.sum()==36
    predictions=np.linalg.eigvalsh(HY+step*ha)
    rows=[]
    for j,e in enumerate(eig):
        lg=np.log(e)/step
        rows.append(dict(index=j,lambda_real=float(e.real),lambda_imag=float(e.imag),modulus=float(abs(e)),
             generator_decay_rate=float(-lg.real),generator_frequency=float(lg.imag),
             marginal_residue_Frobenius=float(np.linalg.norm(residues[j])),slow=bool(slow[j])))
    records.append(dict(step=step,complete_poles=rows,complete_count=len(rows),slow_poles=int(slow.sum()),
        peripheral_count_at_1e_minus_10=int(np.sum(abs(abs(eig)-1)<1e-10)),
        schur_identity_max_error=max(schurerr),spectral_residue_identity_max_error=max(spectralerr),
        residue_sum_identity_max_error=float(np.max(abs(sum(residues)-I36))),
        matrix_eigenbasis_condition=float(np.linalg.cond(V)),
        second_order_coherent_H_eigenvalues=predictions.tolist(),
        second_order_coherent_positive_masses=predictions[3:].tolist(),
        exact_channel_is_unitary=bool(np.max(abs(M.conj().T@M-np.eye(180)))<1e-12)))

out=dict(schema='oph.complete_five_memory_complement.v1',input_sha256=pins,
    domain='five recorded classical labels x six-dimensional Nambu/family density matrix, full180-dimensional CQ response',
    memory_p_exact=str(p),stationary_exact=[str(v)for v in pi],complement_polynomial=str(pol),
    exact_complement_radius='1/12',unit_circle_gap='11/12',
    gap_is_relaxation_not_sterile_mass=True,
    source_Y_eigenvalues=mass.tolist(),first_generator_source_Y_error=float(np.max(abs(A-expected))),
    second_generator_coherent_reconstruction_error=hamerr,
    second_generator_coherent_H_real=ha.real.tolist(),second_generator_coherent_H_imag=ha.imag.tolist(),
    second_generator_dissipative_HS_eigenvalues=np.linalg.eigvalsh((K2+K2.conj().T)/2).tolist(),
    first_order_in_step_coherent_mass_eigenvalue_correction_zero=True,
    source_generation_register_replaced=False,shared_charged_register_correction_derived=False,
    unitary_three_neutral_PMNS_matrix_inferred=False,physical_neutrino_masses_solved=False,
    records=records)
assert pins=={str(p.relative_to(ROOT)):sha(p)for p in paths}
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['complement_polynomial','source_Y_eigenvalues','second_generator_coherent_reconstruction_error']},indent=2))
for r in records:print({k:v for k,v in r.items() if k!='complete_poles'})
