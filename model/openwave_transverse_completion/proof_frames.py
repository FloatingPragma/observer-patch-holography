"""Complete dyadic proof coordinates; every Gram and residual is reverified."""
import hashlib
import json
from pathlib import Path
import numpy as np

WITNESSES=Path(__file__).with_name('complete_frames.json')
FRAME_KEYS={'charged_18_-3','vacuum_18','vacuum_24','vacuum_48'}

def fixed_eigen(A,key,record=None):
    """Supply exact binary64 dyadics to the unchanged complete-frame checker."""
    if record is None:record=json.loads(WITNESSES.read_text())
    assert record['schema']=='oph.transverse.complete_dyadic_proof_frames.v1'
    dv,de=record['vector_denominator'],record['eigenvalue_denominator']
    assert type(dv)is int and type(de)is int and (dv,de)==(1<<24,1<<20)
    assert set(record['frames'])==FRAME_KEYS and key in FRAME_KEYS
    w=record['frames'][key];n=len(A)
    assert A.shape==(n,n) and np.issubdtype(A.dtype,np.integer)
    assert np.array_equal(A,A.T)
    assert type(w['dimension'])is int and w['dimension']==n
    assert hashlib.sha256(json.dumps(A.tolist(),separators=(',',':')).encode()).hexdigest()==w['matrix_sha256']
    vectors,centres=w['integer_frame'],w['integer_eigenvalue_centres']
    assert isinstance(vectors,list) and isinstance(centres,list)
    assert len(vectors)==len(centres)==n
    assert all(isinstance(row,list)and len(row)==n for row in vectors)
    assert all(type(x)is int and abs(x)<=dv for row in vectors for x in row)
    assert all(type(x)is int and abs(x)<(1<<40)for x in centres)
    assert all(a<=b for a,b in zip(centres,centres[1:]))
    Z=np.asarray(vectors,dtype=np.int64);E=np.asarray(centres,dtype=np.int64)
    frame=Z.astype(float)/dv;energy=E.astype(float)/de
    assert np.array_equal(np.rint(frame*dv).astype(np.int64),Z)
    assert np.array_equal(np.rint(energy*de).astype(np.int64),E)
    return energy,frame
