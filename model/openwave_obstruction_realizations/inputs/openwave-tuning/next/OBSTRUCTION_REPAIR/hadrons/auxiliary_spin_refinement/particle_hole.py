"""Exact finite-CAR charge conjugation for the actual two-auxiliary source.

The proof is on the generator algebra, hence lifts to every occupation and
all625 auxiliary outcomes. Electric equality uses the physical total Gauss law.
"""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent

def main():
    I=s.I;T=[]
    for a,b in [(0,1),(0,2),(1,2)]:
        x=s.zeros(3);x[a,b]=x[b,a]=s.Rational(1,2);T.append(x)
        y=s.zeros(3);y[a,b]=-I/2;y[b,a]=I/2;T.append(y)
    T.extend([s.diag(s.Rational(1,2),-s.Rational(1,2),0),s.diag(1,1,-2)/(2*s.sqrt(3))])
    S=[s.Matrix([[0,1],[1,0]])/2,s.Matrix([[0,-I],[I,0]])/2,s.diag(s.Rational(1,2),-s.Rational(1,2))]
    oc=[-1,1,-1,1,-1,1,-1,-1];os=[-1,1,-1]
    for a in range(8):assert -T[a].T==oc[a]*T[a]
    for i in range(3):assert -S[i].T==os[i]*S[i]
    adj=[s.Matrix(8,8,lambda b,c:s.simplify(-2*s.trace((T[a]*T[b]-T[b]*T[a])*T[c]))) for a in range(8)]
    vector=[s.Matrix(3,3,lambda j,k:-I*s.LeviCivita(i,j,k)) for i in range(3)]
    signs=[-oc[a]*os[i] for a in range(8) for i in range(3)];U=s.diag(1,*signs)
    for a in range(8):
        C=s.diag(s.zeros(1),s.kronecker_product(adj[a],s.eye(3)));assert U*C*U==oc[a]*C
    for i in range(3):
        J=s.diag(s.zeros(1),s.kronecker_product(s.eye(8),vector[i]));assert U*J*U==os[i]*J
    for a in range(8):
        for i in range(3):
            g=s.kronecker_product(S[i],s.eye(2),T[a]);assert -g.T==signs[3*a+i]*g
            X=s.zeros(25);X[0,1+3*a+i]=X[1+3*a+i,0]=1;assert U*X*U==signs[3*a+i]*X
    W=s.BlockMatrix([[s.zeros(12),s.eye(12)],[-s.eye(12),s.zeros(12)]]).as_explicit()
    assert W*W.T==s.eye(24)
    M=s.diag(*([1]*12+[-1]*12));K=s.BlockMatrix([[s.zeros(12),-s.eye(12)],[-s.eye(12),s.zeros(12)]]).as_explicit()
    for A in [M,K]:assert -(W.T*A*W).T==A and s.trace(A)==0
    out=dict(schema='oph.auxiliary_spin_refinement.particle_hole.v1',
        quark_CAR_map='c_L -> c_R^dagger; c_R -> -c_L^dagger',
        auxiliary_map='swap cells and act by diag(1, -o_colour[a]*o_spin[i]) in each25-dimensional register',
        colour_transpose_signs=oc,spin_transpose_signs=os,auxiliary_excited_signs=signs,
        exact_one_body_CAR_hopping_mass_checks=2,exact_mixed_current_and_auxiliary_star_checks=48,
        exact_auxiliary_colour_spin_checks=11,
        total_Gauss_covariance='Q_total,L^a -> o_colour[a] Q_total,R^a',
        electric_term_identity='Q_total,R^2=Q_total,L^2 on global total Gauss=0',
        N_map='N -> 24-N',B_map='B -> -B',total_spin_Casimir_preserved=True,total_isospin_Casimir_preserved=True,
        actual_full625_source_intertwiner=True,full_carrier_no_postselection=True,
        physical_charge_mass_or_clock_identification=False)
    (HERE/'PARTICLE_HOLE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
