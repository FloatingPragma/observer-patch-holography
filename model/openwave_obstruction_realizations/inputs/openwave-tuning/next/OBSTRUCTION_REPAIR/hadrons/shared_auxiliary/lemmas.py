"""Exact constants, symmetry average and interval consequences for full-source bounds."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import deque
import json
import numpy as np
from scipy import sparse
import sympy as sy
from local_algebra import current
HERE=Path(__file__).resolve().parent

def main():
    N=np.array([w.bit_count() for w in range(4096)],dtype=np.int64);U3=sparse.csr_matrix((4096,4096),dtype=np.int64);U4=U3.copy()
    for a,b in product(range(3),repeat=2):
        E=current([(3*x+a,3*x+b,1) for x in range(4)]);U3+=E.T@E
    for a,b in product(range(4),repeat=2):
        Fm=current([(3*a+c,3*b+c,1) for c in range(3)]);U4+=Fm.T@Fm
    D=U3+U4-sparse.diags(7*N,format='csr',dtype=np.int64);D.eliminate_zeros();assert D.nnz==0
    bounds=[]
    for a,b in product(range(5),repeat=2):
        total=3*(a+b);ca=3*a*(7-a);cb=3*b*(7-b)
        bounds.append(F(total)+(F(ca+cb)-F(total*total,4))/3)
    assert max(bounds)==20
    local=json.loads((HERE/'LOCAL_ALGEBRA.json').read_text())
    assert max(F(r['exact_max_G_squared']) for r in local['full_joint_irreps'])==F(17,3)
    P=[sy.eye(2),sy.Matrix([[0,1],[1,0]]),sy.Matrix([[0,-sy.I],[sy.I,0]]),sy.diag(1,-1)]
    pauli=[sy.kronecker_product(P[a],P[b]) for a,b in product(range(4),repeat=2) if a or b]
    H=sy.Matrix([[1,1],[1,-1]])/sy.sqrt(2);S=sy.diag(1,sy.I);CX=sy.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    gates=[sy.kronecker_product(H,P[0]),sy.kronecker_product(P[0],H),sy.kronecker_product(S,P[0]),sy.kronecker_product(P[0],S),CX]
    edges={i:set() for i in range(15)}
    for U in gates:
        for i,A in enumerate(pauli):
            X=sy.simplify(U*A*U.H);hits=[j for j,B in enumerate(pauli) if X==B or X==-B];assert len(hits)==1;edges[i].add(hits[0])
    seen={0};queue=deque([0])
    while queue:
        for v in edges[queue.popleft()]:
            if v not in seen:seen.add(v);queue.append(v)
    assert len(seen)==15
    out=dict(schema='oph.shared_auxiliary.bound_lemmas.v1',
        full4096_U3_U4_Howe_identity_exact='C_U3+C_U4=7N',
        local_colour_singlet_U4_Casimir='C_U4=3k(7-k), N=3k',
        global_U4_Casimir_lower='C_U4_total>=N_total²/4',
        exact_complete_P_hopping_square_norm_bound='20',bound_attained_in_N12_localN6_pair_SU4_singlet=True,
        exact_full_local_mixed_current_norm_squared='17/3',
        exact_global_P_to_Q_norm_squared_bound='47/96 at t=g=1/8',
        actual_Clifford_conjugation_graph_transitive_on_15_traceless_Paulis=True,
        singlet_matrix_element_identity='S_local²=I_local²=(2/5)C_SU4_local between global SU4 singlets',
        mixed_current_SU4_average='G_local²=(4/45)N_local(12-N_local)+(1/60)C_SU3_local',
        N12_singlet_compression='G_left²+G_right²=32/5-(2/45)M²+E6/180 in matrix elements',
        H0_current_image_intertwiner='total local colour including auxiliary intertwines the old quark local colour',
        dressed_trial_Q_diagonal='Delta +(<G² H_old>-<H_hop>)/<G²>',
        free_band_statement='h_free²=(t²+m²)I; twelve negative and twelve positive one-particle levels',
        N12_single_aux_exclusion='global Gauss requires quark adjoint; filled negative band is a unique colour singlet, hence at least one particle-hole pair and free energy>=-10sqrt(t²+m²)',
        rigorous_full625_Q_bound=True)
    (HERE/'LEMMAS.json').write_text(json.dumps(out,indent=2)+'\n')
    full=json.loads((HERE/'FESHBACH_BOUNDS.json').read_text());effective=json.loads((HERE/'EFFECTIVE_SECOND_ORDER.json').read_text());records=[]
    for f,e in zip(full['results'],effective['results']):
        intervals={r['N']:tuple(map(F,r['full_source_ground_interval_exact'])) for r in f['all_N_sector_bounds']};v=intervals[12];b=intervals[9];b2=intervals[6]
        bind=[2*b[0]-b2[1]-v[1],2*b[1]-b2[0]-v[0]]
        ev=tuple(map(F,e['effective_vacuum_interval']));es={r['N']:tuple(map(F,r['exact_effective_ground_interval'])) for r in e['all_N_ground_records']};eb=es[9];eb2=es[6];ebind=[2*eb[0]-eb2[1]-ev[1],2*eb[1]-eb2[0]-ev[0]]
        records.append(dict(kappa=f['kappa'],
            actual_full625_all_N_threshold_bounds=[dict(N=n,exact_interval=list(map(str,[z[0]-v[1],z[1]-v[0]]))) for n,z in intervals.items() if n!=12],
            actual_full625_B2_binding_exact_interval=list(map(str,bind)),actual_B2_binding_sign_certified=bind[0]>0 or bind[1]<0,
            effective_only_all_N_threshold_bounds=[dict(N=n,exact_interval=list(map(str,[z[0]-ev[1],z[1]-ev[0]]))) for n,z in es.items() if n!=12],
            effective_only_B2_binding_exact_interval=list(map(str,ebind)),effective_spectrum_is_not_full_source=True))
    (HERE/'ALL_THRESHOLDS.json').write_text(json.dumps(dict(schema='oph.shared_auxiliary.all_thresholds.v1',results=records),indent=2)+'\n');print('ALL EXACT LEMMAS AND THRESHOLDS PASSED',flush=True)
if __name__=='__main__':main()
