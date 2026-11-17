import itertools,json,time
import numpy as np
from source import prepare,region_matrix,literal_scaled_region,s
start=time.process_time()
cells=np.array([c for c in itertools.product(range(-2,3),repeat=3) if sum(x*x for x in c)<=4])

def first(C,dC):
    e,u=np.linalg.eigh(C);e=np.clip(e,1e-13,1-1e-13)
    return float(np.log((1-e)/e)@np.diag(u.conj().T@dC@u).real)

rows=[]
for L in [15,25]:
    for aliases in [0,1,2]:
        source=prepare(L,.5,1.,1e-4,aliases)
        C=region_matrix(source,cells);metric=region_matrix(source,cells,'metric');shape=region_matrix(source,cells,'shape')
        Kmetric=first(C,metric);Kshape=first(C,shape)
        row=dict(L=L,h=.5,R=1.,mass=1e-4,aliases=aliases,modes=4*len(cells),entropy=s.entropy(C),covariance_min=float(np.linalg.eigvalsh(C)[0]),covariance_max=float(np.linalg.eigvalsh(C)[-1]),metric_response=Kmetric,shape_response=Kshape,uniform_fixed_volume_residual=Kmetric+Kshape,relative_uniform_residual=abs(Kmetric+Kshape)/abs(Kmetric))
        rows.append(row);print(json.dumps(row),flush=True)

# Literal moving-CAR entropy derivative, separate from the analytic derivative.
L=15;aliases=1;eps=1e-4;source=prepare(L,.5,1.,1e-4,aliases)
C=region_matrix(source,cells);shape=region_matrix(source,cells,'shape')
plus=literal_scaled_region(L,.5,1.,1e-4,aliases,cells,np.exp(-eps))
minus=literal_scaled_region(L,.5,1.,1e-4,aliases,cells,np.exp(eps))
finite=(s.entropy(plus)-s.entropy(minus))/(2*eps)
deriv=first(C,shape)
control=dict(epsilon=eps,finite_entropy_derivative=finite,analytic_shape_derivative=deriv,error=abs(finite-deriv),matrix_derivative_error=float(np.linalg.norm((plus-minus)/(2*eps)-shape)))
out=dict(status='CONTINUUM_CELL_CAR_VOLUME_AND_ALIAS_CONTROLS',rows=rows,literal_shape_control=control,cpu_seconds=time.process_time()-start,exact_continuum_regional_algebra=False,observer_basis='Four spinor components of normalized disjoint cube indicators',metric_and_shape_calculated_independently=True)
with open('CONTROL.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(control,indent=2),flush=True)
assert control['error']<1e-6
