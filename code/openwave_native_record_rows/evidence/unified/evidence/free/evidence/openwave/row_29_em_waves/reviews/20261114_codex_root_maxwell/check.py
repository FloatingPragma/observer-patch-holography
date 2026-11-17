#!/usr/bin/env python3
"""Independent record replay and Fourier covariance check, no producer imports."""
from pathlib import Path
import gzip,hashlib,importlib.util,json,sys,itertools
import numpy as np
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[4]
P=ROOT/'openwave/row_29_em_waves/work/2026-11-14_codex_record_maxwell/experiment'
raw=gzip.decompress((P/'RECORDS.json.gz').read_bytes()); data=json.loads(raw); receipt=json.loads((P/'RECEIPT.json').read_text())
assert hashlib.sha256(raw).hexdigest()==receipt['raw_uncompressed_sha256']
sp=importlib.util.spec_from_file_location('maxwell_review_source',P/'inputs/native.py')
mod=importlib.util.module_from_spec(sp);sys.modules[sp.name]=mod;sp.loader.exec_module(mod)
edges=data['edges']; W=len(edges); cells=list(itertools.product(range(3),range(3),range(3),range(8)))
C=np.zeros((12,216,216),dtype=np.int64);calls=0
for eid,slot,c,j,target in data['branches']:
    initial=[int(i==slot) for i in range(3)]
    out,record=mod.native(initial,(0,1),c);calls+=1
    assert record==(c,j) and out.index(1)==target
    a,b,axis=edges[eid]
    if axis>=0 and slot<2 and target!=slot:
        old=(a,b)[slot]; new=(a,b)[target]
        C[old%12,old//12,new//12]+=1
for port in range(12):
    C[port]+=np.diag(2*W-C[port].sum(axis=1))
assert all(np.array_equal(C[0],x) for x in C)
L=2*W*np.eye(216)-C[0]
for start,tape,expected in zip(data['channel_initial'],data['tapes'],data['occupation_counts']):
    position=start; counts=[0]*216
    for eid,c,j in tape:
        a,b,_=edges[eid]
        out,record=mod.native((int(position==a),int(position==b)),(0,1),c);calls+=1
        assert record==(c,j)
        if j:position=b if j==1 else a
        assert out==(int(position==a),int(position==b))
        counts[position//12]+=1
    assert counts==expected
errors=[]; polarizations=[]; native_lag2=[]
for item in receipt['spectra']:
    spatial=item['spatial']; n=item['block_events']; p=np.array(spatial)*2*np.pi/3
    k=2*np.sin(p/2); vals=[]; fourier_residual=0
    for m in range(8):
        w=2*np.pi*m/8
        v=np.array([np.exp(1j*(np.dot(p,z[:3])+w*z[3])) for z in cells])/np.sqrt(216)
        lam=float(np.vdot(v,L@v).real)
        fourier_residual=max(fourier_residual,float(np.max(np.abs(L@v-lam*v))))
        z=1-lam/(2*W)
        vals.append(1/lam-1/(4*W)-z*(-np.expm1(n*np.log(z)))/(2*W*n*(1-z)**2))
    assert fourier_residual<1e-10
    cov=[sum(vals[m]*np.cos(2*np.pi*m*t/8) for m in range(8))/8 for t in range(4)]
    errors.append(max(abs(a-b) for a,b in zip(cov,item['native_slice_covariance'])))
    K=np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    BB=cov[2]*(K@K.T)
    reflected_EE=-sum(vals[m]*np.cos(2*np.pi*m*2/8)*((2*np.sin(np.pi*m/8))**2*np.eye(3)+np.outer(k,k)) for m in range(8))/8
    assert np.max(np.abs(np.linalg.eigvalsh(BB)-item['magnetic_eigenvalues']))<1e-10
    assert np.max(np.abs(np.linalg.eigvalsh(reflected_EE)-item['reflected_electric_eigenvalues']))<1e-10
    if n==10**9:
        polarizations.append(float(np.max(np.abs(np.linalg.eigvalsh(reflected_EE)-np.linalg.eigvalsh(BB)))/np.linalg.eigvalsh(BB)[-1]))
        native_lag2.append(cov[2])
control=sum(np.cos(2*np.pi*m*2/8)/(4*W) for m in range(8))/8
checks={'all_stored_native_events_replayed':calls==166624,'source_fourier_covariance':max(errors)<1e-10,
        'two_electric_magnetic_modes':max(polarizations)<.01,'changed_law_separation':abs(control)<1e-12<min(native_lag2)}
checks={k:bool(v) for k,v in checks.items()}
out={'checks':checks,'native_calls':calls,'max_covariance_error':max(errors),'max_EE_BB_relative_error':max(polarizations),
     'changed_law_lag2':control,'scope':'Same-group component check; no native billion-event execution, finite Gaussianity, physical clock or signal intervention.'}
Path(__file__).with_name('RESULT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,sort_keys=True));sys.exit(0 if all(checks.values()) else 1)
