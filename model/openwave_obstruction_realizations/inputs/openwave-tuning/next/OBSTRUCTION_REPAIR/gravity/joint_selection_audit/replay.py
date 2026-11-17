"""Independent replay of existing reciprocal thermal geometry selection.

Uses the source's stated finite matrices and complete transition, not its
implementation. It diagnoses admissible-state stationarity, never fits G.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
import hashlib,json,math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
X=np.array([[0.,1.],[1.,0.]])
I=np.eye(2)
D=np.diag([2*i.bit_count()+3 for i in range(8)])
coupling=np.kron(np.kron(X,X),I)+np.kron(np.kron(I,X),X)
rvals=[.9,1.1]
beta=1/math.sqrt(14)
matrices=[(r/4+1/(7*r)+3*r*r/28)*D-coupling/(7*r) for r in rvals]
levels=np.concatenate([np.linalg.eigvalsh(M) for M in matrices])
w=np.exp(-levels);pi=w/w.sum();g=np.repeat([0,1],8)
area=np.where(g==0,20/9,20/11)
T=np.zeros((16,16))
for i in range(16):
    for j in range(16):
        if g[i]!=g[j]:T[j,i]=w[j]/(16*(w[i]+w[j]))
    T[i,i]=1-T[:,i].sum()
pg=np.array([pi[g==a].sum() for a in (0,1)])
conditional_energy=np.array([np.sum(pi[g==a]*levels[g==a])/pg[a] for a in (0,1)])
eta=(conditional_energy[1]-conditional_energy[0])/(20/9-20/11)
f=g-pg[1];variance=np.sum(pi*f*f)
P=T.T
a=np.sum(pi*f*(P@f))/variance
b=P@f-a*f
k0=np.sum(pi*b*b)/variance
# This same-geometry matter tangent is admitted by the full source state
# space, preserves trace, volume and area, but changes matter entropy.
tangent=np.zeros(16);tangent[0]=-1;tangent[7]=1
entropy_slope=float(-np.dot(tangent,np.log(pi)))
area_slope=float(tangent@area)
energy_ex_slope=float(tangent@levels)
eps=min(pi[0],pi[7])/1000
entropy=lambda q:float(-np.sum(q[q>0]*np.log(q[q>0])))
fd=(entropy(pi+eps*tangent)-entropy(pi-eps*tangent))/(2*eps)
# Reproduce the published field randomization without discarding a branch.
q0=pg[g]/8
states=[q0]
records=[]
bar_energy=np.where(g==0,conditional_energy[0],conditional_energy[1])
for step in range(5):
    q=states[-1]
    qnext=T@q
    Sdiff=entropy(qnext)-entropy(q)
    Adiff=float(area@(qnext-q))
    Ediff=float(levels@(qnext-q))
    Exdiff=float((levels-bar_energy)@(qnext-q))
    records.append(dict(step=step+1,area_change=Adiff,entropy_change=Sdiff,
                        naive_generalized_entropy_change=Sdiff+eta*Adiff,
                        excess_energy_correction=Exdiff,
                        corrected_generalized_entropy_change=Sdiff+eta*Adiff-Exdiff,
                        heat_entropy_balance=Sdiff-Ediff))
    states.append(qnext)
# Every i has 16 proposals and two bath outcomes. Sum complete 512 labels.
# Opposite geometries have a thermal acceptance coin; same geometry proposals
# are retained no-ops. This enumeration is an ideal record law, not new native
# pulse execution or physical controller-work verification.
column_sums=[]
count=0
record_transition=np.zeros_like(T)
max_bath_balance=0.
for i in range(16):
    total=0
    for j in range(16):
        allowed=g[i]!=g[j]
        if allowed:
            low,high=sorted((i,j),key=lambda a:(levels[a],a))
            p1=w[high]/(w[low]+w[high]);gap=levels[high]-levels[low]
        else:
            low=high=i;p1=.5;gap=0.
        for before in (0,1):
            prob=(p1 if before else 1-p1)/16
            accepted=allowed and before==int(i==low)
            final=j if accepted else i
            after=1-before if accepted else before
            record_transition[final,i]+=prob
            max_bath_balance=max(max_bath_balance,abs(levels[final]-levels[i]+gap*(after-before)))
            count+=1;total+=prob
    column_sums.append(total)
out=dict(scope='Existing 16-state thermal geometry source, no spatial Einstein transfer',
         eta=eta,geometry_equilibrium_probability=pg[1],memory_a=a,memory_k0=k0,
         complete_ideal_records=count,record_normalization_error=max(abs(x-1) for x in column_sums),
         record_transition_reconstruction_error=float(np.max(abs(record_transition-T))),
         record_bath_energy_balance_error=float(max_bath_balance),
         detailed_balance_error=float(np.max(np.abs(T*pi[None,:]-(T*pi[None,:]).T))),
         fixed_geometry_tangent=dict(entropy_slope=entropy_slope,area_slope=area_slope,
                                    volume_slope=float(tangent.sum()),finite_difference=fd,
                                    corrected_slope=entropy_slope+eta*area_slope-energy_ex_slope),
         first_five_updates=records,
         source_derived_gravity_stationarity=False)
out['checks']={
    'eta_reproduced':.578<eta<.579,
    'memory_coefficient_reproduced':abs(k0-.02170170003)<1e-10,
    'positive_fixed_geometry_entropy_slope':entropy_slope>1,
    'area_and_volume_fixed':area_slope==0 and tangent.sum()==0,
    'thermal_correction_restores_identity':abs(out['fixed_geometry_tangent']['corrected_slope'])<1e-12,
    'finite_difference_reproduces_slope':abs(fd-entropy_slope)<1e-5,
    'complete_records_normalized':out['record_normalization_error']<1e-14,
    'complete_records_reconstruct_transition':out['record_transition_reconstruction_error']<1e-14,
    'all_record_bath_energy_balances':max_bath_balance<1e-14,
    'published_initial_area_response':abs(records[0]['area_change']+.0006270566)<1e-10,
    'exact_balance_all_updates':max(abs(x['corrected_generalized_entropy_change']-x['heat_entropy_balance']) for x in records)<1e-12,
}
out['checks']={key:bool(value) for key,value in out['checks'].items()}
out['pass']=all(out['checks'].values())
inputs={}
for rel in ['reverse-engineering-reality/model/reciprocal_area_dynamics/DERIVATION.md',
            'reverse-engineering-reality/model/reciprocal_area_dynamics/README.md',
            'reverse-engineering-reality/model/native_equilibration/reciprocal_vacuum.py',
            'reverse-engineering-reality/model/coherent_defect_wave/PROOF.md',
            'reverse-engineering-reality/model/particle_geometric_response/README.md',
            'reverse-engineering-reality/model/higgs_curvature_response/PROOF.md']:
    blob=(ROOT/rel).read_bytes();sha=hashlib.sha256(blob).hexdigest()
    (HERE/'inputs').mkdir(exist_ok=True);(HERE/'inputs'/sha).write_bytes(blob);inputs[rel]=sha
out['inputs']=inputs
(HERE/'RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
if not out['pass']:raise SystemExit(1)
