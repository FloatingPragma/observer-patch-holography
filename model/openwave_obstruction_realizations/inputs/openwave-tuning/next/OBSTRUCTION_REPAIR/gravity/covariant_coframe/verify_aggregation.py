"""Small independent full-pair trace versus separation-aggregated trace."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
import json,itertools
from pathlib import Path
import numpy as np
from run import compute,STENCILS,prior
HERE=Path(__file__).resolve().parent
prior.L=5;prior.R=1;prior.WIDTH=3
source=prior.Source()
rec=compute(source,1)
cells=np.array([x for x in itertools.product(range(-1,2),repeat=3) if sum(t*t for t in x)<=1])
N=len(cells);delta=cells[:,None,:]-cells[None,:,:]
def covariance(d):
    ix=d%5
    return source.kernel[ix[:,:,0],ix[:,:,1],ix[:,:,2]]*np.exp(1j*np.pi*d.sum(-1)/5)[:,:,None,None]
C=covariance(delta).transpose(0,2,1,3).reshape(4*N,4*N)
e,u=np.linalg.eigh(C);e=np.clip(e,1e-13,1-1e-13)
K=((u*np.log((1-e)/e))@u.conj().T).reshape(N,4,N,4)
checks=[]
for name,weights in STENCILS.items():
    dc=np.zeros((N,N,4,4),complex)
    for j in range(3):
        for m,w in enumerate(weights,1):
            shift=np.zeros(3,int);shift[j]=m
            dc-=w*((delta[:,:,j]+m)[:,:,None,None]*covariance(delta+shift)
                   -(delta[:,:,j]-m)[:,:,None,None]*covariance(delta-shift))
    full=np.einsum('iajb,jiba->',K,dc)
    agg=complex(*rec['pullbacks'][name]['entropy_inward_unit_dilation'])
    checks.append(dict(stencil=name,full_pair_trace=[full.real,full.imag],aggregated_trace=[agg.real,agg.imag],error=abs(full-agg)))
dense=json.loads((HERE/'DENSE_CONTROL.json').read_text())
mean=1.4*(5+2*np.cos(2*np.pi/5))/7
metric=1.4*rec['modes'][1]['responses']['covariant_Dsquare'][0]
shape=mean*rec['pullbacks']['central2']['entropy_inward_unit_dilation'][0]
out=dict(aggregation=checks,independent_dense_metric_error=abs(metric-dense['metric_entropy_tangent']),
         independent_dense_shape_error=abs(shape-dense['geometric_algebra_entropy_tangent']))
out['pass']=bool(max(x['error'] for x in checks)<1e-10 and out['independent_dense_metric_error']<1e-9 and out['independent_dense_shape_error']<1e-9)
(HERE/'AGGREGATION_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2),flush=True)
if not out['pass']:raise SystemExit(1)
