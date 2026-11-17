import json,time
import numpy as np
from kernels import cell,response_kernels,regional_response
from regions import build,aggregates,trace_coefficients,trace_response
start=time.process_time()
source=cell.prepare(16,.5,1.,0.,1,.5)
region=build(source,2)
shifts=[(0,0,0),(1,0,0)]
d,acc=aggregates(region,source,shifts)
C=cell.region_matrix(source,region['cells'])
e,u=np.linalg.eigh(C);e=np.clip(e,1e-13,1-1e-13);K=(u*np.log((1-e)/e))@u.conj().T
rows=[]
for key in ['metric','shape']:
    full=float(np.trace(K@cell.region_matrix(source,region['cells'],key)).real)
    reduced=trace_coefficients(acc[0],source,key,d)
    rows.append(dict(kind=key,full=full,reduced=[reduced.real,reduced.imag],error=abs(full-reduced)))
for i,shift in enumerate(shifts):
    km,kn=response_kernels(16,.5,1.,0.,1,.5,shift)
    for kind,kernel in [('metric_vertex',km),('lapse_vertex',kn)]:
        full=float(np.trace(K@regional_response(kernel,source,region['cells'],shift)).real)
        reduced=trace_response(acc[i],kernel,source,d)
        rows.append(dict(kind=kind,shift=shift,full=full,reduced=[reduced.real,reduced.imag],error=abs(full-reduced)))
        if i==0 and kind=='metric_vertex':
            direct=float(np.trace(K@cell.region_matrix(source,region['cells'],'metric')).real)
            rows.append(dict(kind='uniform_coframe_vertex_crosscheck',full=direct,reduced=[reduced.real,reduced.imag],error=abs(direct-reduced)))
out=dict(status='COMPLETE_SECTOR_KERNEL_AND_AGGREGATE_CONTROL',rows=rows,entropy_difference=abs(cell.s.entropy(C)-region['entropy']),all_irreps_retained=True,cpu_seconds=time.process_time()-start)
with open('KERNEL_CONTROL.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps(out,indent=2),flush=True)
assert max(r['error'] for r in rows)<1e-8
