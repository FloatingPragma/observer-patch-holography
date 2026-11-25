"""Native material-record local channels on the registered Interval Echo board."""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np
import sympy as sp

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier
from ..signed_effects.sequential import compare

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECEIPT = HERE / 'results.json'
LOADS = tuple(product(range(3), repeat=2))
RECORDS = LOADS
DIM = 99


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


@lru_cache(None)
def exact_echo():
    """Independent exact assembly from the public native-event formula."""
    w = sp.MutableSparseMatrix(90, 9, {})
    for j, (a, b) in enumerate(LOADS):
        s = a + b
        for c in (0, 1):
            y = (s // 2 + c * (s % 2), s // 2 + (1-c) * (s % 2))
            label = law.LABELS.index((c, a-b))
            w[9*label+LOADS.index(y), j] = sp.sqrt(2)/2
    u = sp.MutableSparseMatrix(DIM, DIM, {})
    u[:9, 9:] = w.T
    u[9:, :9] = w
    u[9:, 9:] = sp.SparseMatrix.eye(90)-w*w.T
    return sp.ImmutableSparseMatrix(u)


@lru_cache(None)
def material_projectors():
    return tuple(sp.ImmutableSparseMatrix(DIM, DIM,
                 {(9*m+r, 9*m+r): 1 for m in range(11)}) for r in range(9))


@lru_cache(None)
def exact_instrument(adaptive=True):
    """All81 terminal two-read branches, without resetting either register."""
    u = exact_echo()
    p = material_projectors()
    identity = sp.SparseMatrix.eye(DIM)
    return tuple((r, s, sp.ImmutableSparseMatrix(p[s]*u*(u if a < b and adaptive else identity)*p[r]*u))
                 for r, (a, b) in enumerate(RECORDS) for s in range(9))


def numeric(matrix):
    out=np.zeros(matrix.shape,complex)
    for (row,col),v in matrix.todok().items():out[row,col]=complex(v)
    return out


@lru_cache(None)
def branch_maps(adaptive=True):
    """Return every public (first,second) material record and full99 postmap."""
    return {(r, s): numeric(k) for r, s, k in exact_instrument(adaptive)}


def zero(matrix):
    return all(sp.expand(v) == 0 for v in matrix.todok().values())


@lru_cache(None)
def exact_certificate():
    u = exact_echo()
    p = material_projectors()
    n = sp.diag(*[a+b for m in range(11) for a, b in LOADS])
    checks = {
        'echo_unitary_full99': zero(sp.ImmutableSparseMatrix(u.T*u-sp.eye(DIM))),
        'echo_preserves_local_charge': zero(sp.ImmutableSparseMatrix(n*u-u*n)),
        'nine_material_projectors_partition_full99': zero(sp.ImmutableSparseMatrix(sum(p, sp.SparseMatrix.zeros(DIM))-sp.eye(DIM))),
        'all_material_projectors_preserve_charge': all(zero(sp.ImmutableSparseMatrix(n*x-x*n)) for x in p),
    }
    for adaptive in (False, True):
        branches = exact_instrument(adaptive)
        gram = sp.SparseMatrix.zeros(DIM)
        for _, _, k in branches:
            gram += k.T*k
        checks[f'{"adaptive" if adaptive else "fixed"}_81_kraus_complete'] = zero(sp.ImmutableSparseMatrix(gram-sp.eye(DIM)))
        checks[f'{"adaptive" if adaptive else "fixed"}_all_branches_charge_preserving'] = all(
            zero(sp.ImmutableSparseMatrix(n*k-k*n)) for _, _, k in branches)
    require(all(checks.values()), 'exact native channel certificate failed')
    return checks


def geometry():
    board = carrier.federation(0)
    active, remote = (0, 1), (2, 3)
    require(active in board.seams and remote in board.seams, 'seams absent from pinned board')
    boundary = [list(e) for e in board.seams if len(set(e) & set(active)) == 1]
    degree_ports = int(sum(board.degree == 5))
    # The actual register index is a permutation, not an empirical factor fit.
    count = Counter()
    images = set()
    for a, b, c, d in product(range(3), repeat=4):
        for m, t in product(range(11), repeat=2):
            A, B = 9*m+3*a+b, 9*t+3*c+d
            q, Q = a+b, a+b+c+d
            images.add((A, B))
            count[Q] += 1
            require((A % 9)//3 == a and A % 3 == b and A//9 == m, 'split inverse failed')
            require(0 <= q <= 4, 'charge block invalid')
    require(len(images) == 99**2, 'register split is not bijective')
    return dict(ports=board.n, seams=len(board.seams), reference_qutrits=degree_ports,
                reference_collars=degree_ports, active_ports=list(active), active_collar=list(active),
                remote_ports=list(remote), remote_collar=list(remote),
                local_dimension=99, local_charge_dimensions=[11,22,33,22,11],
                boundary_seams=boundary, boundary_degree_references=[v for v in active if board.degree[v] == 5],
                selected_union_boundary_seams=[list(e) for e in board.seams if len(set(e)&{0,1,2,3})==1],
                selected_union_degree_references=[v for v in range(4) if board.degree[v]==5],
                unaffected_quantum_register_dimensions={'port_qutrits':238,'seam_collars':629,'reference_qutrit_collar_pairs':180},
                four_port_two_collar_sector_dimensions={str(k):count[k] for k in sorted(count)},
                all_register_split_words=len(images))


def balanced_reference():
    """Full Schmidt rank, fixed total four; reference is actual disjoint registers."""
    c = np.zeros((99,99), complex)
    for m in range(11):
        for j, (a,b) in enumerate(LOADS):
            c[9*m+j, 9*m+LOADS.index((2-a,2-b))] = 1/np.sqrt(99)
    return c


def remote_density(c):
    return c.T @ c.conj()


def reference_audit():
    c = balanced_reference()
    target = remote_density(c)
    maps = branch_maps(True)
    outputs = [k@c for k in maps.values()]
    after = sum((remote_density(x) for x in outputs), np.zeros((99,99),complex))
    # A complete Gram identity establishes this for arbitrary references and spectators.
    gram = sum((k.conj().T@k for k in maps.values()), np.zeros((99,99),complex))
    # Keep just first public01 and discard the other records: a trace-losing filter.
    accepted = [x for (r,s),x in zip(maps,outputs,strict=True) if r == LOADS.index((0,1))]
    sub = sum((remote_density(x) for x in accepted), np.zeros((99,99),complex))
    acceptance = float(np.trace(sub).real)
    conditioned = sub/acceptance
    # Second service on the disjoint material/collar factor commutes branchwise.
    b = numeric(material_projectors()[4]*exact_echo())
    order_error = max(np.linalg.norm(k@(c@b.T)-(k@c)@b.T) for k in maps.values())
    return dict(balanced_total_load=4, input_schmidt_rank=int(np.linalg.matrix_rank(c)),
                fixed_total_four_dimension=2299,
                normalized_remote_marginal_error=float(np.max(np.abs(after-target))),
                full_choi_trace_preservation_error=float(np.max(np.abs(gram-np.eye(99)))),
                reference_with_two_collar_states=99,
                disjoint_native_branch_order_error=float(order_error),
                postselected_first01_probability=acceptance,
                postselected_remote_trace_distance=float(np.linalg.eigvalsh(conditioned-target).__abs__().sum()/2))


def probabilities(vector, adaptive=True):
    result = []
    for (r,s), k in branch_maps(adaptive).items():
        x = k@vector
        result.append(float(np.vdot(x,x).real))
    require(abs(sum(result)-1)<1e-12, 'incomplete public outcomes')
    return np.asarray(result).reshape(9,9)


def apply_echo_sparse(state, edge):
    """Independent embedding on three physical ports and two distinct collars."""
    u = exact_echo()
    columns = {}
    for (row,col),v in u.todok().items():
        columns.setdefault(col,[]).append((row, complex(v)))
    out = {}
    ia, ib, im = edge
    for word, amplitude in state.items():
        col = 9*word[im]+3*word[ia]+word[ib]
        for row, v in columns[col]:
            new = list(word)
            new[im], j = divmod(row,9)
            new[ia], new[ib] = LOADS[j]
            key = tuple(new)
            out[key] = out.get(key,0)+amplitude*v
    return {k:v for k,v in out.items() if abs(v)>1e-15}


def controls():
    v = np.zeros(99,complex)
    v[LOADS.index((0,1))] = 1
    adaptive, fixed = probabilities(v,True), probabilities(v,False)
    # Distinct collars do not make overlapping endpoint services disjoint.
    start = {(0,2,1,0,0): 1.0}
    e, f = (0,1,3), (1,2,4) # physical seams(0,1) and(1,5)
    ef = apply_echo_sparse(apply_echo_sparse(start,e),f)
    fe = apply_echo_sparse(apply_echo_sparse(start,f),e)
    def loads(state):
        p = np.zeros(27)
        for w,a in state.items():p[9*w[0]+3*w[1]+w[2]] += abs(a)**2
        return p
    false_tv = float(np.abs(loads(ef)-loads(fe)).sum()/2)
    # Blank code boundary contact on physical(1,5), using the actual native law.
    from scipy.linalg import expm
    w = np.zeros(9,complex); w[LOADS.index((1,0))] = 1
    h0,h1=law.local_term_parts()
    physical=expm(-1j*np.pi/2*(h0+8*h1))@np.r_[w,np.zeros(90)]
    flowed=physical[:9]
    require(np.linalg.norm(physical[9:])<1e-12,'blank boundary contact leaked')
    outside_one = sum(abs(flowed[i])**2 for i,(_,b) in enumerate(LOADS) if b == 1)
    require(false_tv>.1 and outside_one>.49, 'boundary/locality controls failed')
    blank=np.zeros(99);blank[LOADS.index((1,1))]=1
    dirty=np.zeros(99);dirty[9*(1+law.LABELS.index((0,-2)))+LOADS.index((1,1))]=1
    u=numeric(exact_echo())
    def material(v):return (np.abs(v.reshape(11,9))**2).sum(axis=0)
    hidden_memory_tv=float(np.abs(material(u@blank)-material(u@dirty)).sum()/2)
    return dict(sharp_input=[0,1], adaptive_records=adaptive.tolist(),
                omit_feedback_records=fixed.tolist(), feedback_record_total_variation=float(np.abs(adaptive-fixed).sum()/2),
                distinct_collars_shared_port_order_record_tv=false_tv,
                open_boundary_seam=[1,5], open_boundary_dwell='pi/2',
                sealed_outside_load1_probability=0.0, open_outside_load1_probability=float(outside_one),
                same_material11_different_collar_output_tv=hidden_memory_tv)


@lru_cache(None)
def selected_support_branches():
    """Source-record-selected support, as exact tensor factors on9801 inputs."""
    u=exact_echo(); p=material_projectors(); identity=sp.ImmutableSparseMatrix.eye(99)
    result=[]
    for r,(a,b) in enumerate(RECORDS):
        first=p[r]*u
        for t in range(9):
            if a<b:
                result.append((r,t,(2,3),sp.ImmutableSparseMatrix(first),sp.ImmutableSparseMatrix(p[t]*u)))
            else:
                result.append((r,t,(0,1),sp.ImmutableSparseMatrix(p[t]*u*first),identity))
    return tuple(result)


def selected_support_certificate():
    # For each first record, second-read completeness reduces its Gram sum
    # to first^*first tensor I. Sum over first yields I tensor I.
    u=exact_echo(); p=material_projectors(); identity=sp.SparseMatrix.eye(99)
    checks={}
    for r,(a,b) in enumerate(RECORDS):
        group=[(left,right) for rr,_,_,left,right in selected_support_branches() if rr==r]
        first=p[r]*u
        if a<b:
            right=sum((right.T*right for _,right in group),sp.SparseMatrix.zeros(99))
            good=zero(sp.ImmutableSparseMatrix(right-identity)) and all(left==first for left,_ in group)
        else:
            left=sum((left.T*left for left,_ in group),sp.SparseMatrix.zeros(99))
            good=zero(sp.ImmutableSparseMatrix(left-first.T*first)) and all(right==identity for _,right in group)
        checks[str(r)]=good
    require(all(checks.values()),'record-selected tensor completeness failed')
    return checks


def selected_apply(state, left, right):
    from scipy.sparse import csr_matrix
    def sparse(matrix):
        items=matrix.todok()
        rows=[k[0] for k in items];cols=[k[1] for k in items]
        return csr_matrix(([complex(v) for v in items.values()],(rows,cols)),shape=(99,99))
    out=(sparse(left)@state.reshape(99,-1)).reshape(state.shape)
    return (sparse(right)@out.transpose(1,0,2).reshape(99,-1)).reshape(99,99,-1).transpose(1,0,2)


def selected_support_audit():
    # An arbitrary complex input within fixed total6, with arbitrary dirty
    # local collars and an actually accessible material reference on(4,5).
    q=np.array([a+b for a,b in LOADS]); localq=np.tile(q,11)
    rng=np.random.default_rng(440728)
    c=rng.normal(size=(99,99,9))+1j*rng.normal(size=(99,99,9))
    c*=localq[:,None,None]+localq[None,:,None]+q[None,None,:]==6
    c/=np.linalg.norm(c)
    initial=c.reshape(-1,9).T@c.reshape(-1,9).conj()
    after=np.zeros((9,9),complex); chosen=np.zeros((9,9),complex)
    for r,t,support,left,right in selected_support_branches():
        out=selected_apply(c,left,right).reshape(-1,9)
        marginal=out.T@out.conj()
        after+=marginal
        if support==(2,3):chosen+=marginal
    weight=float(np.trace(chosen).real)
    trace_distance=float(np.abs(np.linalg.eigvalsh(chosen/weight-initial)).sum()/2)
    # Simple public control: first01 consult chooses B; fixed-A cannot emit
    # second11 from A's total1. B's total2 balances to11 with certainty.
    sharp=np.zeros((99,99,1),complex)
    sharp[LOADS.index((0,1)),LOADS.index((2,0)),0]=1
    native=np.zeros((9,9)); forced=np.zeros((9,9))
    u=exact_echo(); p=material_projectors(); identity=sp.ImmutableSparseMatrix.eye(99)
    for r,t,_,left,right in selected_support_branches():
        native[r,t]=float(np.linalg.norm(selected_apply(sharp,left,right))**2)
        forced[r,t]=float(np.linalg.norm(selected_apply(sharp,p[t]*u*p[r]*u,identity))**2)
    require(np.max(np.abs(after-initial))<1e-12 and abs(native.sum()-1)<1e-12,'selected-support normalization failed')
    return dict(exact_per_first_record_completeness=selected_support_certificate(),
                union_ports=[0,1,2,3], union_collars=[[0,1],[2,3]], union_dimension=99**2,
                complement_material_reference=[4,5], reference_input_total_load=6,
                full_dirty_collar_fixed_charge_input=True, complement_marginal_max_error=float(np.max(np.abs(after-initial))),
                selected_B_probability=weight, selected_B_conditional_complement_trace_distance=trace_distance,
                sharp_input_A=[0,1],sharp_input_B=[2,0], source_selected_records=native.tolist(),
                forced_A_records=forced.tolist(), wrong_support_record_total_variation=float(np.abs(native-forced).sum()/2),
                consultation_rule='first material a<b selects seam(2,3); otherwise seam(0,1)',
                no_signalling_domain='complement of union(0,1,2,3) and their two enrolled collars; no conditioning')


def report():
    geo = geometry()
    native_error = float(np.max(np.abs(numeric(exact_echo()).real-law.echo())))
    refs = reference_audit()
    ctrl = controls()
    require(refs['normalized_remote_marginal_error']<1e-12 and native_error<1e-12,'native/reference audit failed')
    paths = [Path(__file__), Path(law.__file__), Path(carrier.__file__),
             ROOT/'model/candidate_laws/interval_echo/SPEC.md',
             ROOT/'model/candidate_laws/contrast_memory_spherical_v1/SPEC.md',
             carrier.HERE/'carriers/FEDERATION_LEVEL0.npz']
    return dict(source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                geometry=geo, exact_checks=exact_certificate(), native_echo_max_error=native_error,
                public_records={'first': [list(x) for x in RECORDS], 'second': [list(x) for x in RECORDS],
                                'all_pairs':81,'filtered_or_discarded_outcomes':0},
                scheduler={'first':'echo; material read', 'consultation':'first local material pair a<b',
                           'feedback':'native echo if a<b; identity otherwise',
                           'second':'echo; material read', 'reset':False,
                           'complete_native_echo_count':{'a<b':3,'a>=b':2}},
                reference=refs, controls=ctrl, selected_support=selected_support_audit(),
                scope='Protected addressed native echo services and supplied complete sharp material reads; no private collar read, free-boundary closure, source-selected architecture, clock, spacelike region or continuum inference.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    actual=report()
    if args.write:
        with RECEIPT.open('x') as f:json.dump(actual,f,indent=2,allow_nan=False);f.write('\n')
    else:
        compare(json.loads(RECEIPT.read_text()),actual)
        print('local interface: exact full99 native channels, all81 adaptive records, fixed-charge references and boundary controls replay')


if __name__=='__main__':main()
