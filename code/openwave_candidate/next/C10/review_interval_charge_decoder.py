"""Independent incidence/interval/PVM review, including idle-reference branches."""
from pathlib import Path
import json,hashlib
import numpy as np
import interval_charge_decoder as producer


def run():
    receipt=Path(producer.__file__).with_name('INTERVAL_CHARGE_DECODER.json');saved=json.loads(receipt.read_text())
    m=producer.FermionicSharedRecords(producer.Construction(3,((0,1),(1,2),(2,0)),.23),7,6)
    incidence=np.zeros((6,6),int)
    for e,(a,b) in enumerate(m.construction.seams):
        incidence[a,2*e]=1;incidence[b,2*e+1]=1;incidence[3+e,2*e:2*e+2]=-1
    occupations=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in m.states]);records=np.array([f for _,f in m.states])
    gauss=(occupations+records@incidence.T)%6
    assert np.all(gauss==np.array([1,1,1,0,0,0]))
    rng=np.random.default_rng(234);psi=rng.normal(size=(m.n,2))+1j*rng.normal(size=(m.n,2));psi/=np.linalg.norm(psi)
    branch_error=0.;comparisons=0
    for bits in range(64):
        region=np.array([(bits>>z)&1 for z in range(6)])
        true=occupations@region;boundary=records@(region@incidence);gr=gauss[0]@region
        residue=(gr-boundary)%6
        decoded=np.array([next(q for q in range(4) if q%6==r) for r in residue])
        assert np.array_equal(decoded,true)
        for q in range(4):
            direct=psi*(true==q)[:,None];indirect=psi*(decoded==q)[:,None]
            branch_error=max(branch_error,float(np.max(abs(direct-indirect))))
        comparisons+=m.n
    # General interval proof by difference, independently enumerate all widths.
    for total in range(7):
        for modulus in range(max(2,total+1),10):assert len({q%modulus for q in range(total+1)})==total+1
    assert 0%3==3%3 and producer.decode_charge(3,0,3,6)==3
    assert saved['projector_basis_comparisons']==comparisons
    files=[Path(__file__),Path(producer.__file__),receipt,producer.BASE/'next/C9/fermionic_shared_records.py']
    return {'schema':'oph.c10.interval_charge_independent_review.v1','outcome':'PASS_WITH_KNOWN_NONNEGATIVE_INTERVAL_SCOPE','independent_method':'Independent six-vertex oriented half-edge incidence; exhaustive all64masks/all120states; unique residue intersection with0..3; entangled two-dimensional idle-reference branch vectors','basis_comparisons':comparisons,'external_reference_branch_error':branch_error,'known_interval_condition':'N>L; arbitrary signed[-L,L] remains differentdomain','Z3_alias_independent':True,'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(producer.BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
if __name__=='__main__':
 r=run();Path(__file__).with_name('INTERVAL_CHARGE_REFERENCE_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
