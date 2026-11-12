"""Independent exact all-cut review from original Gauss columns and link matrices."""
from collections import Counter
from pathlib import Path
import hashlib,json,random
import numpy as np
import sympy as s
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.chiral_joint.oriented_links.model import endpoint_current
from model.coherent_defect_wave.chiral_joint.periodic_family.model import Configuration,kinetic_terms
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.model import gauss_column


def run():
    root=Path(__file__).resolve().parents[3]
    matter=representation()[2][11];scalar=higgs_memory()[1][11]
    ends=[endpoint_current(11,k) for k in (0,1)]
    assert ends[0]+ends[1]==s.zeros(119)
    edges=[(t.path[0],t.left_cell,t.right_cell) for t in kinetic_terms(Configuration())]
    rng=random.Random(20261112);checks=0
    for _ in range(24):
        occupation=rng.getrandbits(240);links=tuple(rng.randrange(119) for _ in range(24));scalars=tuple(rng.randrange(5) for _ in range(8))
        key=(occupation,links,scalars)
        density=[sum((matter[j%15,j%15] for j in range(30) if occupation>>(30*x+j)&1),s.Integer(0))+scalar[scalars[x],scalars[x]] for x in range(8)]
        source=[gauss_column(key,x,11).get(key,s.Integer(0)) for x in range(8)]
        density[0]+=1;source[0]+=1
        for mask in range(256):
            cut=sum((ends[0 if mask>>a&1 else 1][links[e],links[e]] for e,a,b in edges if ((mask>>a)^(mask>>b))&1),s.Integer(0))
            assert sum(source[x] for x in range(8) if mask>>x&1)==sum(density[x] for x in range(8) if mask>>x&1)+cut
            checks+=1
    # Direct full-basis endpoint constraints for both oriented ec routes.
    ranks=[]
    for origin in (0,1):
        indices=[i for i in range(119) if ends[origin][i,i]==-1]
        assert len(indices)==1
        i=indices[0];vec=s.eye(119)[:,i]
        for a,t in enumerate(representation()[2]):
            ec=next(j for j in range(15) if matter[j,j]==1)
            assert endpoint_current(a,origin)*vec==-t[ec,ec]*vec
            assert endpoint_current(a,1-origin)*vec==t[ec,ec]*vec
        ranks.append(i)
    total_Q=sum((endpoint_current(a,k) for a in (10,11) for k in (0,1)),s.zeros(119))
    assert total_Q!=s.zeros(119)
    # Independent polynomial convolution, with exact Python integers.
    weights=Counter(int(6*z) for z in ends[0].diagonal());offset=-min(weights)
    polynomial=np.array([weights.get(i-offset,0) for i in range(max(weights)+offset+1)],dtype=object)
    q=np.array([1],dtype=object);degrees=[]
    recorded=json.loads((root/'model/openwave_native_charge_flux/results.json').read_text())['actual_route_attachment']['full_boundary_flux_spectral_multiplicities_by_edge_count']
    assert weights==Counter({-k:v for k,v in weights.items()})
    for n in range(25):
        assert sum(q)==119**n
        assert all(x>=0 for x in q)
        if str(n) in recorded:
            expected={str(i-n*offset):str(int(x)) for i,x in enumerate(q) if x}
            assert expected==recorded[str(n)]
        degrees.append(len(q))
        q=np.convolve(q,polynomial)
    paths=['model/openwave_native_charge_flux/PROOF.md','model/openwave_native_charge_flux/algebra.py','model/openwave_native_charge_flux/results.json',str(Path(__file__).relative_to(root))]
    result={'verdict':'PASS_INDEPENDENT_NATIVE_GAUSS_FLUX_REVIEW','random_complete_basis_states':24,'exact_all_cut_operator_checks':checks,
        'oriented_unique_flux_indices':ranks,'full_boundary_tensor_dimensions_checked':[str(119**n) for n in range(25)],
        'omitted_weak_gauge_current_control_rejected':True,
        'scope':'Exact original source identities and actual endpoint constraints; coarse positive probabilities and unevaluated failure-flux distributions retain the producer scope.',
        'source_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}}
    (root/'model/openwave_native_charge_flux/joint_flux_review.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['verdict'])


if __name__=='__main__':run()
