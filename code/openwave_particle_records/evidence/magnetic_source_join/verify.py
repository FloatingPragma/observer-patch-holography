#!/usr/bin/env python3
"""Short reviewer replay; exhaustive all-seam receipts are retained separately."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import headline as h

started=time.monotonic()
registered=json.loads((h.HERE/'PREREG.json').read_text())
assert registered['sha256']==h.dependencies()
custody=json.loads((h.HERE/'SOURCE_PINS.json').read_text())
for group in ('runtime_dependency_sha256','cited_proof_dependency_sha256'):
    assert all(h.sha(h.ROOT/p)==digest for p,digest in custody[group].items())
exact=json.loads((h.HERE/'exact_results.json').read_text())
source=json.loads((h.HERE/'source_results.json').read_text())
for r in (exact,source):
    assert r['preregistration_sha256']==h.sha(h.HERE/'PREREG.json')
    assert r['source_sha256']==h.dependencies()
    assert r['exact_integer_checks_pass']
    assert not r['physical_row5_complete'] and not r['physical_row6_complete']

# Relabel covariance is proved for arbitrary sites in PROOF.md; representative
# seams check the complete local/crossing law. The exhaustive132-case audit
# remains in both receipts and is reproduced by headline.py, not inferred here.
for sites,seam in ((2,(0,1)),(4,(0,3))):
    for law in ('C','E'):
        new=h.exact_audit(sites,law,seam)
        retained=next(r for r in exact['full_Fock_audits']
                      if r['sites']==sites and r['law']==law)
        assert new==retained
for law in ('C','E'):
    new=h.exact_audit(12,law,(0,11),[1<<i for i in range(24)])
    retained=next(r for r in exact['twelve_port_one_particle_audits']
                  if r['seam']==[0,11] and r['law']==law)
    assert new==retained
assert len(exact['twelve_port_one_particle_audits'])==132
assert source['full_Fock_audits']==exact['full_Fock_audits']
assert source['twelve_port_one_particle_audits']==exact['twelve_port_one_particle_audits']

execution=source['source_execution']
assert execution['source_arrays_sha256']==h.sha(h.HERE/'source_arrays.npz')
arrays=np.load(h.HERE/'source_arrays.npz')
spins=np.array([sum(((m>>(2*x))&1)-((m>>(2*x+1))&1)
                  for x in range(2))/2 for m in range(16)])
pointer=spins[np.arange(8192)>>9]
maximum=0.0
for case in execution['cases']:
    law=case['law']
    for name,sign in (('native',1),('spin_reversal',-1)):
        final=arrays[f'{law}_{name}_final_amplitudes']
        assert final.shape==(8192,16)
        probabilities=abs(final)**2
        assert np.max(abs(probabilities.sum(axis=0)-1))<2e-7
        wrong=[float(probabilities[pointer!=sign*spins[col],col].sum())
               for col in range(16)]
        maximum=max(maximum,max(wrong))
        assert max(wrong)<2e-12
        assert np.array_equal(wrong,case['pointer'][name]['wrong_spin_probability_by_input'])
assert execution['control_source_matrix_error']<2e-7
assert exact['changed_law']['max_doubled_Sz_intertwiner_squared_Frobenius']>0
print(json.dumps({'status':'PASS_NATIVE_FIELD_SPIN_CURRENT',
                  'retained_complete_source_laws':4,
                  'addresses_per_input':8192,'Fock_inputs_per_law':16,
                  'wrong_spin_probability_max':maximum,
                  'elapsed_seconds':time.monotonic()-started,
                  'physical_row5_complete':False,'physical_row6_complete':False}))
