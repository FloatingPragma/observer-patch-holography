"""Exact N=1 coherence obstruction for the two-setting phase/parity pump.

The proof applies to uniform iid seam scheduling on any finite simple graph.
The exact carrier row norm is a sharper certificate, not a spectral-radius fit.
"""
import sympy as s
from .common import require
from model.common_source import slow_geometry


def certificate():
    # Verify the only nontrivial local coefficient from literal Kraus rows.
    plus=s.Matrix([1,1,0])/s.sqrt(2)
    minus=s.Matrix([1,-1,0])/s.sqrt(2)
    K0=plus*plus.T+s.diag(0,0,1)
    K1=plus*minus.T
    ops=[]
    for p in (0,1):
        V=s.eye(3);V[p,p]=s.I
        ops.extend(K*V for K in (K0,K1))
    gamma=(1+s.I)/4
    for p in range(3):
        for q in range(3):
            E=s.zeros(3);E[p,q]=1
            Y=s.simplify(sum((K*E*K.conjugate().T for K in ops),s.zeros(3))/2)
            if p!=q:
                require(all(Y[k,k]==0 for k in range(3)),'zero-diagonal invariant')
                if p<2 and q<2:require(Y==s.zeros(3),'active coherence erased')
                elif p<2:
                    target=s.zeros(3);target[0,2]=target[1,2]=gamma
                    require(Y==target,'endpoint/outside coherence coefficient')
                else:
                    target=s.zeros(3);target[2,0]=target[2,1]=s.conjugate(gamma)
                    require(Y==target,'outside/endpoint coherence coefficient')
    seams=list(slow_geometry.source().seams());M=len(seams);n=12
    degree=[sum(p in edge for edge in seams) for p in range(n)]
    row_norms={};bounds={}
    for p in range(n):
        for q in range(n):
            if p==q:continue
            coeff={}
            def add(key,z):coeff[key]=coeff.get(key,s.Integer(0))+z
            for a,b in seams:
                pin=p in (a,b);qin=q in (a,b)
                if pin and qin:continue
                if pin:
                    add((a,q),gamma);add((b,q),gamma)
                elif qin:
                    add((p,a),s.conjugate(gamma));add((p,b),s.conjugate(gamma))
                else:add((p,q),1)
            norm=s.simplify(sum(s.sqrt(s.expand_complex(z*s.conjugate(z))) for z in coeff.values())/M)
            A=int(tuple(sorted((p,q))) in [tuple(sorted(e)) for e in seams])
            bound=1-((degree[p]+degree[q]-2*A)*(1-1/s.sqrt(2))+A)/M
            require(s.simplify(bound-norm).is_nonnegative,'convex row norm bound')
            row_norms.setdefault(A,set()).add(norm);bounds.setdefault(A,set()).add(s.simplify(bound))
    exact=max((v for vals in row_norms.values() for v in vals),key=lambda x:float(x))
    general=max((v for vals in bounds.values() for v in vals),key=lambda x:float(x))
    require(exact==s.Rational(3,4)+s.sqrt(2)/12,'exact carrier infinity norm')
    require(general==1-(1-1/s.sqrt(2))/3,'carrier graph bound')
    require(exact<1,'strict coherence contraction')
    return {'status':'EXACT N=1 local identities and actual-carrier coefficient row sums',
      'local_matrix_units_checked':9,'carrier_off_diagonal_rows_checked':132,
      'gamma_endpoint_outside':str(gamma),
      'row_norm_by_adjacency':{str(k):list(map(str,v)) for k,v in row_norms.items()},
      'coherence_infinity_norm':str(exact),'coherence_norm_numeric':float(exact),
      'per_30_attempt_sweep_norm_upper_numeric':float(exact**M),
      'graph_only_bound':str(general),'graph_only_per_sweep_numeric':float(general**M),
      'general_bound':'1-((degree(p)+degree(q)-2*A_pq)*(1-1/sqrt(2))+A_pq)/M',
      'positive_min_degree_per_sweep_upper':'exp(-2*delta*(1-1/sqrt(2)))',
      'population_block':'I-L/(2*M), real nonnegative eigenvalues',
      'scope':'The unchanged N=1 pump with uniform iid seams cannot produce a long-lived infrared coherent wave in its sweep clock. No claim covers arbitrary controlled schedules, added coherent memory, or all many-load sectors.'}
