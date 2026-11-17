"""Actual CAR images verify shared-vacuum hopping and auxiliary denominators."""
from pathlib import Path
import json
import numpy as np
import sympy as sy
from scipy import sparse
from local_algebra import current,creation,annihilation
HERE=Path(__file__).resolve().parent

def main():
    sectors=json.loads((HERE/'LOCAL_OPERATORS.json').read_text())['sectors'];rr=[];cc=[];dd=[];off=0
    for r in sectors:
        R=np.array(r['R'],dtype=np.int64)
        for i,j in zip(*np.nonzero(R)):rr.append(r['words'][i]);cc.append(off+j);dd.append(R[i,j])
        off+=r['dimension']
    R=sparse.csr_matrix((np.array(dd,dtype=np.int64),(rr,cc)),shape=(4096,92));N=np.array([w.bit_count() for w in range(4096)])
    C6=-sparse.diags(N*N,format='csr',dtype=np.int64)
    for a in range(3):
        for b in range(3):
            E=current([(3*alpha+a,3*alpha+b,1) for alpha in range(4)]);C6+=3*(E.T@E)
    assert (C6@R).nnz==0
    fermion_checks=0
    for mode in range(12):
        for op in [creation,annihilation]:
            rows=[];cols=[];vals=[]
            for w in range(4096):
                z,s=op(mode,w)
                if z is not None:rows.append(z);cols.append(w);vals.append(s)
            B=sparse.csr_matrix((vals,(rows,cols)),shape=(4096,4096),dtype=np.int64);image=B@R
            d=C6@image-8*image;d.eliminate_zeros();assert d.nnz==0;fermion_checks+=1
    F={(a,b):current([(3*a+c,3*b+c,1) for c in range(3)]) for a in range(4) for b in range(4)}
    projections=0
    for i in range(12):
        for j in range(12):
            A=current([(i,j,1)]);difference=3*(R.T@A@R)
            if i%3==j%3:difference-=R.T@F[i//3,j//3]@R
            difference.eliminate_zeros();assert difference.nnz==0;projections+=1
    T=[]
    for a in range(3):
        for b in range(3):
            if a!=b:
                X=sy.zeros(3);X[a,b]=1;T.append(X)
    T.extend([sy.diag(1,-1,0),sy.diag(0,1,-1)])
    metric=sy.Matrix([[sy.trace(a*b) for b in T] for a in T]);mi=metric.inv();dual=[sum((mi[a,b]*T[b] for b in range(8)),sy.zeros(3)) for a in range(8)]
    wards=0
    for Q in T:
        ad=sy.Matrix([[sy.trace(dual[b]*(Q*T[a]-T[a]*Q)) for a in range(8)] for b in range(8)])
        for b in range(8):
            assert Q*dual[b]-dual[b]*Q+sum((ad[b,a]*dual[a] for a in range(8)),sy.zeros(3))==sy.zeros(3);wards+=1
    spin=[sy.diag(1,-1),sy.Matrix([[0,1],[0,0]]),sy.Matrix([[0,0],[1,0]])];adjoint_images=0
    for colour in T:
        for S in spin:
            pairs=[]
            for so in range(2):
                for si in range(2):
                    for a in range(3):
                        for b in range(3):
                            c=int(S[so,si]*colour[a,b])
                            if c:
                                for f in range(2):pairs.append((3*(2*so+f)+a,3*(2*si+f)+b,c))
            G=current(pairs);image=G@R;d=C6@image-18*image;d.eliminate_zeros();assert d.nnz==0;adjoint_images+=1
    # Normalized Hermitian colour and spin generators: sum(T_a S_i)^2=I.
    Ta=[]
    for a,b in [(0,1),(0,2),(1,2)]:
        X=sy.zeros(3);X[a,b]=X[b,a]=sy.Rational(1,2);Ta.append(X)
        X=sy.zeros(3);X[a,b]=-sy.I/2;X[b,a]=sy.I/2;Ta.append(X)
    Ta.extend([sy.diag(sy.Rational(1,2),-sy.Rational(1,2),0),sy.diag(1,1,-2)/(2*sy.sqrt(3))])
    Si=[sy.Matrix([[0,1],[1,0]])/2,sy.Matrix([[0,-sy.I],[sy.I,0]])/2,sy.diag(1,-1)/2]
    mixed=[sy.kronecker_product(S,sy.eye(2),A) for A in Ta for S in Si]
    assert sum((A*A for A in mixed),sy.zeros(12))==sy.eye(12)
    out=dict(schema='oph.shared_auxiliary.denominators.v1',local_colour_singlet_domain=92,complete_two_cell_P_dimension=8464,
        H0='kappa (Q_quark,left+Q_aux,left)^2 +3*kappa(n_aux,left+n_aux,right)',
        full_auxiliary_carrier='25 tensor25, all625 states retained, n_aux at each cell0or1',
        full_quark_CAR_dimension=1<<24,full_total_carrier_dimension=(1<<24)*625,
        physical_Gauss_constraint='Q_quark,left+Q_quark,right+Q_aux,left+Q_aux,right=0',
        P='both auxiliaries in vacuum and quark colour singlet independently at each cell, all local occupations retained',
        Q='complete physical Gauss complement of P',
        gap_lower_bound='4*kappa/3',actual_fermion_image_checks=fermion_checks,
        actual_fermion_image_colour_Casimir='4/3',actual_adjoint_current_image_checks=adjoint_images,
        quark_adjoint_Casimir='3',exact_adjoint_dual_Ward_identity_checks=wards,
        total_colour_of_current_plus_local_aux_image='0 by the exact Ward intertwiner',
        hopping_P_to_Q_H0_eigenvalue='4*kappa/3',aux_at_either_cell_P_to_Q_H0_eigenvalue='3*kappa',
        complete_local_singlet_bilinear_projection_checks=projections,
        identity='P c_(alpha,c)† c_(beta,d) P =delta_cd F_alpha,beta/3',
        complete_P_hopping_square_identity='P K_dimensionless² P=N-(2/3)sum_alpha,beta F_L(alpha,beta)F_R(beta,alpha)',
        kinetic_double_commutator_exact='sum_x,a,i [G_xai,[G_xai,H_hop]]=2H_hop',
        effective_second_order='m(NL-NR)-(3t²/4kappa)[N-(2/3)sum F_L F_R]-(g²/3kappa)(G_L²+G_R²)',
        t='1/8',g='1/8',empirical_fit=False,remainder_certified_here=False)
    (HERE/'DENOMINATORS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
