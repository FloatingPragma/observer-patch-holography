#!/usr/bin/env python3
"""Independent finite-claim audit. Only writes beside this script with --write.

Uses endpoint floor/ceiling averaging and a solved absorbing Markov chain,
not producer flow/moment/energy routines. Reads graph and frozen arrays as data.
Does not award OpenWave cells or claim an independent producing group.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import sys
import time
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKET = ROOT / 'openwave/shared/static_candidate_join'
PROVIDER = ROOT / 'openwave/shared/localized_rest_join'
INPUTS = [
    PACKET / 'PROOF.md', PACKET / 'CONTRACT.json', PACKET / 'CLAUSES.md',
    PACKET / 'run.py', PACKET / 'source.py',
    ROOT / 'openwave/shared/compact_static_adjudication/INTERFACE_THEOREM.md',
    ROOT / 'openwave/shared/compact_static_adjudication/encoded_source.py',
    PROVIDER / 'PROOF.md', PROVIDER / 'FINITE_STACK_SCOPE.md',
    PROVIDER / 'evidence/graph.json',
    PROVIDER / 'evidence/raw_source_arrays.npz',
    PROVIDER / 'evidence/complete_outcome_laws.npz',
    PROVIDER / 'standalone/reverse-engineering-reality/model/observer/native.py',
    ROOT / 'reverse-engineering-reality/Lean/EventAlgebra/NativeBinarySourceEndpoints.lean',
    ROOT / 'reverse-engineering-reality/Lean/Screen/RegionalContinuity.lean',
    ROOT / 'reverse-engineering-reality/Lean/Screen/SeamCurrentEdge30Moment.lean',
]

INPUT_PATHS = [p.relative_to(ROOT) for p in INPUTS if not str(p.relative_to(ROOT)).startswith('reverse-engineering-reality/Lean/')]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def averaged(q, edge, coin):
    """Independent floor/ceiling endpoint construction, including signed loads."""
    a, b = edge
    out = list(q)
    total = q[a] + q[b]
    lo, hi = total // 2, -((-total) // 2)
    out[a], out[b] = (lo, hi) if coin == 0 else (hi, lo)
    return tuple(out), q[a] - out[a]

def source_pair(p, r, n):
    return tuple(int(i == p) - int(i == r) for i in range(n))

def extract_native():
    """Execute only the actual native endpoint function, without module imports."""
    p = PROVIDER / 'standalone/reverse-engineering-reality/model/observer/native.py'
    tree = ast.parse(p.read_text())
    definition = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'native')
    ns = {}
    exec(compile(ast.Module(body=[definition], type_ignores=[]), str(p), 'exec'), ns)
    return ns['native']

def compare_results(actual, expected, path="result"):
    """Exact structure and discrete evidence; finite floats within absolute 1e-10."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError("Dictionary fields differ: " + path)
        for key, value in expected.items():
            compare_results(actual[key], value, path + "." + str(key))
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError("List structure differs: " + path)
        for i, (got, wanted) in enumerate(zip(actual, expected)):
            compare_results(got, wanted, path + "[" + str(i) + "]")
    elif isinstance(expected, float):
        if (type(actual) not in (int, float) or not math.isfinite(actual)
                or not math.isfinite(expected) or abs(actual - expected) > 1e-10):
            raise ValueError("Finite numeric comparison failed: " + path)
    elif type(actual) is not type(expected) or actual != expected:
        raise ValueError("Exact comparison failed: " + path)


def main():
    global ROOT, PACKET, PROVIDER
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root',type=Path,default=ROOT,help='Workspace-shaped evidence root containing openwave/shared/')
    parser.add_argument('--lean-root',type=Path,help='Optional RER root; verifies the inspected Lean bytes, without building Lean')
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    ROOT=args.evidence_root.resolve()
    PACKET=ROOT/'openwave/shared/static_candidate_join'
    PROVIDER=ROOT/'openwave/shared/localized_rest_join'
    lean_sources=json.loads((HERE/'LEAN_SOURCES.json').read_text())
    if args.lean_root:
        for name,h in lean_sources.items():assert sha(args.lean_root/name)==h,name
    started = time.perf_counter()
    inputs_before = {str(p): sha(ROOT/p) for p in INPUT_PATHS}
    graph = json.loads((PROVIDER / 'evidence/graph.json').read_text())
    edges = [tuple(e) for e in graph['seams']]
    cell = [e for e in edges if e[1] < 12]
    assert len(edges) == 630 and len(cell) == 30
    raw = extract_native()
    endpoint_cases = 0
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for c in (0, 1):
                out, j = averaged((x, y), (0, 1), c)
                actual, record = raw((x, y), (0, 1), c)
                assert (actual, record) == (out, (c, j))
                endpoint_cases += 1

    # Rebuild faces from the actual graph, without importing the producer list.
    es = set(cell)
    faces = [t for t in combinations(range(12), 3)
             if all(e in es for e in combinations(t, 2))]
    assert len(faces) == 20
    B = sp.zeros(12, 30)
    for k, (a, b) in enumerate(cell):
        B[a, k], B[b, k] = 1, -1
    D = sp.zeros(30, 20)
    for k, (a, b, c) in enumerate(faces):
        D[cell.index((a,b)), k] = 1
        D[cell.index((b,c)), k] = 1
        D[cell.index((a,c)), k] = -1
    assert B * D == sp.zeros(12, 20)
    smith = smith_normal_form(D, domain=sp.ZZ)
    factors = [abs(int(smith[i,i])) for i in range(20)]
    assert factors == [1] * 19 + [0] and B.rank() == 11

    # Solve the full transient-chain linear system independently. The producer
    # substitutes the proposed flow into Bellman's equation instead.
    states = [source_pair(p,r,12) for p in range(12) for r in range(12) if p != r]
    index = {q:i for i,q in enumerate(states)}
    Qcounts, reward = sp.zeros(132), sp.zeros(132,30)
    killed = sp.zeros(132,1)
    partial = defaultdict(set)
    branches = 0
    for i,q in enumerate(states + [(0,)*12]):
        for e,edge in enumerate(cell):
            for c in (0,1):
                out,j = averaged(q,edge,c)
                key = (e,c,j)
                assert out not in partial[key]
                partial[key].add(out)
                assert sp.Matrix(q)-sp.Matrix(out) == B[:,e]*j
                if i < 132:
                    reward[i,e] += j
                    if any(out): Qcounts[i,index[out]] += 1
                    else: killed[i,0] += 1
                branches += 1
    A = 60*sp.eye(132)-Qcounts
    fundamental = A.inv(method='DM')
    absorption = fundamental*killed
    duration = fundamental*sp.ones(132,1)*60
    flow = fundamental*reward
    L = B*B.T
    plus = (L+sp.ones(12)/12).inv()-sp.ones(12)/12
    predicted = sp.Matrix.hstack(*(B.T*plus*sp.Matrix(q) for q in states)).T
    assert flow == predicted
    assert absorption == sp.ones(132,1)
    assert min(duration) > 0
    assert all(B*flow[i,:].T == sp.Matrix(q) for i,q in enumerate(states))

    # Source controls cannot absorb, checked exhaustively on every cell state.
    control_branches = 0
    for q in states:
        for a,b in cell:
            for c in (0,1):
                out=list(q)
                if c == 0: out[a],out[b] = out[b],out[a]
                assert any(out) and sorted(out) == sorted(q)
                control_branches += 1
    split = [Fraction(1),Fraction(0),Fraction(-1)]
    split[0]-=Fraction(1,2); split[1]+=Fraction(1,2)
    assert split[:2] == [Fraction(1,2)]*2

    # Full graph Laplacian, positive-transfer bound, and compact vectors.
    Lg = np.zeros((240,240),dtype=np.int64)
    for a,b in edges:
        Lg[a,a]+=1; Lg[b,b]+=1; Lg[a,b]-=1; Lg[b,a]-=1
    f=np.zeros(240,dtype=np.int64); f[[4,5,6,7]]=[1,-1,-1,1]
    g=np.zeros(240,dtype=np.int64); g[[16,17,18,19]]=[1,-1,-1,1]
    assert np.array_equal(Lg@f,6*f) and np.array_equal(Lg@g,6*g)
    assert np.max(np.diag(Lg)) == 6
    constant=np.ones(240,dtype=np.int64)
    assert not np.any(Lg@constant)
    contrast=np.zeros(240,dtype=np.int64); contrast[:12]=1;contrast[12:24]=-1
    low_rayleigh=Fraction(int(contrast@Lg@contrast),int(contrast@contrast))
    assert low_rayleigh == Fraction(1,3)
    # Sum actual independent endpoint maps on every support column, preserving
    # all off-support destinations and vacuum, rather than compressing first.
    transfer=Counter(); touched=set(); direct_calls=0
    for p in np.flatnonzero(f):
        for r in np.flatnonzero(g):
            q=source_pair(p,r,240)
            for edge in edges:
                for c in (0,1):
                    out,j=averaged(q,edge,c);direct_calls+=1
                    nz=tuple(i for i,v in enumerate(out) if v)
                    key=(-1,-1) if not nz else (out.index(1),out.index(-1))
                    touched.add(key)
                    transfer[key]+=int(f[p]*g[r])
    occupied={(p,r) for p in np.flatnonzero(f) for r in np.flatnonzero(g)}
    assert all(transfer[k] == 1248*int(f[k[0]]*g[k[1]]) if k in occupied else transfer[k] == 0 for k in touched)
    assert transfer[(-1,-1)] == 0
    exterior = touched - occupied - {(-1,-1)}
    assert len(exterior) == 64
    E6=-1260*np.log1p(-6/1260)
    E12=-1260*np.log1p(-12/1260)
    anchor=511/E6
    excess=1260*np.log1p(1/43680)
    assert abs(E12-2*E6-excess)<1e-12

    arrays=np.load(PROVIDER/'evidence/raw_source_arrays.npz')
    pointer_errors={}
    for name,vec in [('f',f),('g',g)]:
        target=np.zeros(256,dtype=complex);target[:240]=1j*vec/2
        U=arrays['U_'+name];V=arrays['inverse_'+name]
        error=max(float(np.max(abs(U[:,0]-target))),float(np.max(abs(V@U-np.eye(256)))))
        assert error<1e-10
        pointer_errors[name]=error
    # Derive ideal complete-certificate control probabilities directly.
    ray=f/2;plant=ray.copy();plant[4]*=-1
    ideal_certificate={'single_coherent':float((ray@ray)**2),
        'single_dephased':float(np.sum(ray**4)),
        'single_sign_plant':float((ray@plant)**2),
        'pair_coherent':float((ray@ray)**4),
        'pair_dephased':float(np.sum(ray**4)**2),
        'pair_sign_plant':float((ray@plant)**2)}
    laws=np.load(PROVIDER/'evidence/complete_outcome_laws.npz')
    saved_errors={}
    for mode in ('coherent','dephased','sign_plant'):
        for law in ('C','E'):
            arr=laws[f'{law}_energy_{mode}']
            err=max(abs(float(arr.sum())-1),abs(float(arr[0,0])-ideal_certificate['pair_'+mode]),
                abs(float(arr.sum(axis=1)[0])-ideal_certificate['single_'+mode]))
            assert err<1e-10;saved_errors[f'{law}_{mode}']=err

    # Recompute exact radical moment tensors from endpoint averaging. This
    # uses SymPy expressions, independently of the producer's Q(phi) algebra.
    phi=(1+sp.sqrt(5))/2
    pts=[sp.Matrix(v) for v in [(-1,phi,0),(1,phi,0),(-1,-phi,0),(1,-phi,0),
        (0,-1,phi),(0,1,phi),(0,-1,-phi),(0,1,-phi),(phi,0,-1),(phi,0,1),(-phi,0,-1),(-phi,0,1)]]
    expected_edges=[(a,b) for a in range(12) for b in range(a+1,12) if sp.expand((pts[a]-pts[b]).dot(pts[a]-pts[b])) == 4]
    assert expected_edges == cell
    anti=[next(j for j,x in enumerate(pts) if x == -v) for v in pts]
    R=sp.diag(-1,-1,1)
    permutation=[next(j for j,x in enumerate(pts) if x == R*v) for v in pts]
    assert all(tuple(sorted((permutation[a],permutation[b]))) in es for a,b in cell)
    moments=[];fixed_bad=0
    par=(3-phi)/9;perp=(phi+2)/15
    for p,v in enumerate(pts):
        q=source_pair(p,anti[p],12)
        vectors=[];fixed=[]
        for edge in cell:
            for c in (0,1):
                _,j=averaged(q,edge,c)
                d=j*(pts[edge[1]]-pts[edge[0]])
                vectors.append(d)
                if edge == (0,5): fixed.append(d)
        mean=sum(vectors,sp.zeros(3,1))/60
        cov=(sum((d*d.T for d in vectors),sp.zeros(3))/60-mean*mean.T).applyfunc(sp.simplify)
        assert (mean-(phi-3)*v/15).applyfunc(sp.simplify) == sp.zeros(3,1)
        fit=perp*sp.eye(3)+(par-perp)*(v*v.T)/v.dot(v)
        assert (cov-fit).applyfunc(sp.simplify) == sp.zeros(3)
        fm=sum(fixed,sp.zeros(3,1))/2
        fc=(sum((d*d.T for d in fixed),sp.zeros(3))/2-fm*fm.T).applyfunc(sp.simplify)
        fp=sp.simplify((v.T*fc*v)[0]/v.dot(v));ft=sp.simplify((sp.trace(fc)-fp)/2)
        ff=ft*sp.eye(3)+(fp-ft)*v*v.T/v.dot(v)
        aligned=(fm*v.dot(v)-v*fm.dot(v)).applyfunc(sp.simplify)==sp.zeros(3,1)
        uniaxial=(fc-ff).applyfunc(sp.simplify)==sp.zeros(3)
        fixed_bad+=int(not(aligned and uniaxial))
        moments.append((mean,cov))
    m,C=moments[0];mr,Cr=moments[permutation[0]]
    assert (mr+m).applyfunc(sp.simplify)==sp.zeros(3,1)
    assert (Cr-C).applyfunc(sp.simplify)==sp.zeros(3)
    assert (R*mr-m).applyfunc(sp.simplify)==sp.zeros(3,1)
    assert (moments[1][1]-C).applyfunc(sp.simplify)!=sp.zeros(3)
    assert fixed_bad==4

    # Reconstruct the retained classical ancestry and its failed quantum join.
    old_branches=0;different=0
    for p,r,s in ((p,r,s) for p in range(12) for r in range(12) for s in range(12) if len({p,r,s})==3):
        x=tuple(int(i in (p,s)) for i in range(12));y=tuple(int(i in (r,s)) for i in range(12))
        q=tuple(a-b for a,b in zip(x,y))
        for edge in cell:
            a,b=edge
            for c in (0,1):
                xo,jx=averaged(x,edge,c);yo,jy=averaged(y,edge,c)
                target=tuple(u-v for u,v in zip(xo,yo));k=jx-jy
                candidates=[c2 for c2 in (0,1) if averaged(q,edge,c2)==(target,k)]
                assert candidates
                if q[a]-q[b] in (-1,1):assert len(candidates)==1
                different+=int(averaged(q,edge,c)!=(target,k));old_branches+=1
    assert old_branches==79200 and different==12000
    old_records=[]
    old_joint=[]
    for p,r,c in ((0,2,0),(3,4,1)):
        x=tuple(int(i in (p,1)) for i in range(12));y=tuple(int(i in (r,1)) for i in range(12))
        xo,jx=averaged(x,(0,1),c);yo,jy=averaged(y,(0,1),c)
        old_records.append((c,jx,jy))
        environment=(c,jx,jy,tuple(i for i,(a,b) in enumerate(zip(xo,yo)) if a==b==1))
        old_joint.append((p,r,environment,sp.Rational(1,2)))
        assert averaged(source_pair(p,r,12),(0,1),1)==(source_pair(p,r,12),0)
    assert old_records[0]!=old_records[1]
    old_reduced=sp.Matrix(2,2,lambda i,j: old_joint[i][3]*old_joint[j][3] if old_joint[i][2]==old_joint[j][2] else 0)
    new_reduced=sp.ones(2)/4
    assert old_reduced == sp.eye(2)/4 and new_reduced[0,1] == sp.Rational(1,4)
    assert sp.trace(old_reduced) == sp.trace(new_reduced) == sp.Rational(1,2)
    for name,h in inputs_before.items():assert sha(ROOT/name)==h,name
    result={'schema':'oph.static-submission-independent-audit.v1','all_checks_pass':True,
      'audit_scope':'Independent implementation by a second Codex agent in the producing group; no G6 or platform award.',
      'native_endpoint_cases':endpoint_cases,
      'row_1':{'verdict':'CONFIRMED finite source; QUALIFIED electric-charge identification',
        'cell_branches':branches,'face_smith_factors':factors,'boundary_rank':11,
        'solved_transient_states':132,'solved_flow_entries':132*30,'absorption_probabilities':'all exactly 1',
        'expected_absorption_events_min':str(min(duration)),'expected_absorption_events_max':str(max(duration)),
        'complete_swap_control_branches':control_branches,'half_current_plant':['1/2','1/2']},
      'row_2':{'verdict':'CONFIRMED finite localized energy certificate; QUALIFIED rest-energy interpretation',
        'source_positions':240,'source_seams':630,'support_per_arm':4,'laplacian_eigenvalue':6,
        'full_pair_support_branch_calls':direct_calls,'cancelled_exterior_coordinates':len(exterior),
        'vacuum_output_amplitude':0,'transfer_single':'209/210','transfer_pair':'104/105',
        'lower_nonconstant_rayleigh':str(low_rayleigh),'E6':float(E6),'E12':float(E12),
        'pair_keV_at_511_anchor':float(E12*anchor),'finite_pair_excess_record_units':float(excess),
        'preparation_inverse_max_errors':pointer_errors,'ideal_certificate':ideal_certificate,
        'saved_complete_law_max_errors':saved_errors},
      'row_7':{'verdict':'CONFIRMED finite tensor/frame statement; QUALIFIED physical-rotation dictionary',
        'source_axes':12,'enumerated_alternatives':720,'mean_coefficient':str(sp.simplify((phi-3)/15)),
        'variance_parallel':str(sp.simplify(par)),'variance_transverse':str(sp.simplify(perp)),
        'half_turn_tensor_return':True,'half_turn_frame_reversal':True,'squared_frame_return':True,
        'covariance_nonconstant':True,'fixed_seam_failed_inputs':fixed_bad},
      'join':{'verdict':'CONFIRMED classical ancestry; REFUTED coherent equivalence (as disclosed)',
        'old_pair_branches':old_branches,'original_coin_disagreements':different,
        'old_distinguishing_records':old_records,'new_offdiagonal':str(new_reduced[0,1]),'old_reduced_offdiagonal':str(old_reduced[0,1]),
        'branch_probability':str(sp.trace(new_reduced))},
      'inputs_sha256':inputs_before,'inspected_lean_sources_sha256':lean_sources,
      'not_checked':['No new Lean build','No full rerun of the 82,575,360-column provider audit',
        'No laboratory precision validation','No independent producing group','No acceptance of the physical dictionary']}
    text=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:(HERE/'AUDIT.json').write_text(text)
    if args.check:compare_results(json.loads(text),json.loads((HERE/'AUDIT.json').read_text()))
    print(json.dumps({'all_checks_pass':True,'elapsed_seconds':time.perf_counter()-started,
        'lean_source_custody_checked':bool(args.lean_root),'row_1':result['row_1'],'row_2':result['row_2'],'row_7':result['row_7'],'join':result['join']},indent=2))
if __name__=='__main__':main()
