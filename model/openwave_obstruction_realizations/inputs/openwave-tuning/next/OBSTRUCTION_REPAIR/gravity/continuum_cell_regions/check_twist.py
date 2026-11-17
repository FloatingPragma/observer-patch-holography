import itertools,json,time
import numpy as np
from source import prepare,region_matrix,literal_scaled_region,s
start=time.process_time()
cells=np.array([c for c in itertools.product(range(-2,3),repeat=3) if sum(x*x for x in c)<=4])
def first(C,dC):
    e,u=np.linalg.eigh(C);e=np.clip(e,1e-13,1-1e-13)
    return float(np.log((1-e)/e)@np.diag(u.conj().T@dC@u).real)
rows=[]
for L in [16,24,40]:
    for aliases in [1,2] if L<40 else [1]:
        source=prepare(L,.5,1.,0.,aliases,.5)
        C=region_matrix(source,cells);metric=region_matrix(source,cells,'metric');shape=region_matrix(source,cells,'shape')
        Km=first(C,metric);Ks=first(C,shape)
        row=dict(L=L,h=.5,R=1.,mass=0.,twist=.5,aliases=aliases,entropy=s.entropy(C),metric_response=Km,shape_response=Ks,uniform_fixed_volume_residual=Km+Ks,relative_uniform_residual=abs(Km+Ks)/abs(Km))
        rows.append(row);print(json.dumps(row),flush=True)
L=16;aliases=1;eps=1e-4;source=prepare(L,.5,1.,0.,aliases,.5)
C=region_matrix(source,cells);shape=region_matrix(source,cells,'shape')
plus=literal_scaled_region(L,.5,1.,0.,aliases,cells,np.exp(-eps),.5)
minus=literal_scaled_region(L,.5,1.,0.,aliases,cells,np.exp(eps),.5)
control=dict(finite_entropy_derivative=(s.entropy(plus)-s.entropy(minus))/(2*eps),analytic_shape_derivative=first(C,shape),matrix_derivative_error=float(np.linalg.norm((plus-minus)/(2*eps)-shape)))
control['entropy_derivative_error']=abs(control['finite_entropy_derivative']-control['analytic_shape_derivative'])
with open('TWIST_CONTROL.json','w') as f:json.dump(dict(status='ANTIPERIODIC_MASSLESS_CONTINUUM_CELL_CONTROL',rows=rows,literal_control=control,cpu_seconds=time.process_time()-start),f,indent=2)
print(json.dumps(control,indent=2),flush=True)
assert control['entropy_derivative_error']<1e-6
