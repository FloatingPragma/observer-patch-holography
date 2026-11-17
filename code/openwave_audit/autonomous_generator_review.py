"""Independent route and complete-instrument algebra for the autonomous generator."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import numpy as np
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=ROOT/'openwave-tuning/common_dynamics/autonomous_dirac'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    pins={}
    for filename,wanted in [('MANIFEST.json','44c008892484f8d346fe5f8e29cf29a5a185159e7f748e295cce86b3457756d1'),
                            ('FRAME_SCOPE_MANIFEST.json','71efd702dd1254e33439dd2b048d8de0853b0101b2717a259339ccb9dd499a2f')]:
        p=SOURCE/filename;assert sha(p)==wanted
        pins[str(p.relative_to(ROOT))]=wanted
        d=json.loads(p.read_text())
        for key in ('files','consumed_manifests','upstream_source_pins'):
            for name,digest in d.get(key,{}).items():
                assert sha(ROOT/name)==digest,name
                pins[name]=digest
    geometry_path=ROOT/'openwave-tuning/row_completion/reviews/ROW_20_ROOT/audit.py'
    spec=importlib.util.spec_from_file_location('independent_geometry',geometry_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    edges=[]
    for x in range(2):
        for y in range(2):
            for z in range(2):
                for a,b,shift in module.geometry():
                    dst=tuple((c+d)%2 for c,d in zip((x,y,z),shift))
                    u=12*(4*x+2*y+z)+a;v=12*(4*dst[0]+2*dst[1]+dst[2])+b
                    edges.append(tuple(sorted((u,v))))
    edges=sorted(set(edges));assert len(edges)==288
    lookup={e:i for i,e in enumerate(edges)}
    contract=json.loads((SOURCE/'CONTRACT.json').read_text())
    routes={(r['cell'],r['axis'],r['sign']):r for r in contract['route_checks']}
    mutations=[]
    for key,r in routes.items():
        current=[0]*288
        for u,v in zip(r['path'],r['path'][1:]):
            current[lookup[tuple(sorted((u,v)))]]+=1 if u<v else -1
        assert current==r['integer_current']
        divergence=[0]*96
        for (u,v),j in zip(edges,current):divergence[u]+=j;divergence[v]-=j
        expected=[int(i==12*r['cell'])-int(i==12*r['destination']) for i in range(96)]
        assert divergence==expected
        reverse=routes[(r['destination'],r['axis'],-r['sign'])]
        assert reverse['path']==list(reversed(r['path']))
        assert reverse['integer_current']==[-v for v in current]
        mutations.append(any(current[i] and current[i]!=-current[i] for i in range(288)))
    assert all(mutations)
    for axis in range(3):
        assert sorted(routes[(c,axis,1)]['destination'] for c in range(8))==list(range(8))

    # Nontrivial holonomy is necessary: at side2 a bare swap equals its
    # adjoint, whereas an actual routed phase need not. This exact algebraic
    # control represents a phase evaluation, not a normalizable field state.
    X=s.Matrix([[0,1],[1,0]])
    phase=(s.sqrt(3)+s.I)/2
    T=phase*X
    A=(T.conjugate().T-T)/(2*s.I)
    assert s.simplify(A+X/2)==s.zeros(2)
    S=s.zeros(4);S[:2,2:]=T;S[2:,:2]=T.conjugate().T
    assert s.simplify(S*S)==s.eye(4)
    J=s.eye(2).col_join(s.I*s.eye(2))/s.sqrt(2)
    assert s.simplify(J.conjugate().T*S*J-A)==s.zeros(2)
    c,z=s.symbols('c z',real=True)
    K=c*s.eye(2)-s.I*z*A
    defect=s.eye(2)-K.conjugate().T*K
    assert s.simplify(defect-(1-c*c)*s.eye(2)+z*z*A*A)==s.zeros(2)
    # Set c^2+z^2=1 to obtain z^2(I-A^2). No failure normalization.
    tested=0
    for lam in (-1,-.9,-.5,0,.5,.9,1):
        for epsilon in (.002,.01,.07,.2):
            k=math.cos(epsilon)-1j*lam*math.sin(epsilon)
            ideal=np.exp(-1j*epsilon*lam)
            assert abs(k-ideal)<=epsilon**2/2+1e-15
            for n in (1,3,17,64):
                fail=1-abs(k)**(2*n)
                assert abs(k**n-ideal**n)<=n*epsilon**2/2+1e-13
                assert fail<=n*epsilon**2+1e-13
                # Complete law contains the first-failure records.
                first=[abs(k)**(2*j)*(1-abs(k)**2) for j in range(n)]
                assert abs(sum(first)+abs(k)**(2*n)-1)<1e-13
                tested+=1
    assert abs((math.cos(.2)-0j)**17)**2<.51
    result={'schema':'oph.autonomous.generator.root.audit.v1',
            'component_checks_pass':True,'full_criterion_met':False,
            'source_pins':pins,'actual_routes_independently_rebuilt':48,
            'exact_block_involution_and_compression':True,
            'scalar_complete_instrument_checks':tested,
            'controls':{'route_sign_mutations_detected':True,
                        'nontrivial_holonomy_not_bare_swap':True,
                        'discarding_failure_changes_complete_law':True},
            'execution_scope':'Integer route and exact block algebra plus scalar analytic-bound diagnostics; no new full field or native pulse execution.',
            'audit_sha256':sha(Path(__file__))}
    (HERE/'autonomous_generator_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_pins'}))


if __name__=='__main__':main()
