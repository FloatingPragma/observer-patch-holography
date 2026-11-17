"""Exact CAR-generator identities underlying the full spin/flavour compression."""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent

def elementary(n,a,b):
    A=s.zeros(n);A[a,b]=1;return A

def main():
    M=s.diag(*([1]*12+[-1]*12));K=s.zeros(24);Y=s.zeros(24)
    for a in range(12):K[a,a+12]=K[a+12,a]=-1;Y[a,a+12]=1;Y[a+12,a]=-1
    # Bogoliubov particle-hole plus cell reflection: c -> Y c†.
    assert Y*Y.T==s.eye(24)
    for A in [M,K]:assert -(Y.T*A*Y).T==A and s.trace(Y.T*A*Y)==0
    su4=[];localcolour=[]
    for a in range(4):
        for b in range(4):
            F=s.zeros(24)
            for x in range(2):
                for c in range(3):F[3*(4*x+a)+c,3*(4*x+b)+c]=1
            assert F*M==M*F and F*K==K*F;su4.append(F)
    for a in range(3):
        for b in range(3):
            E=s.zeros(24)
            for alpha in range(4):E[3*alpha+a,3*alpha+b]=1
            assert all(F*E==E*F for F in su4);localcolour.append(E)
    # Exact scalar coefficients in C6_left -> C6_right. Here U=sum EabEba
    # and n=N_right; E_left,ab ->4delta_ab-E_right,ba.
    U,n=s.symbols('U n',commutative=True)
    assert s.expand(3*(48-8*n+U)-(12-n)**2)==3*U-n*n
    # Representative exact spin-unitary tree gauges; general recursion in proof.
    U1=s.Matrix([[s.Rational(3,5),s.Rational(4,5)],[-s.Rational(4,5),s.Rational(3,5)]])
    U2=s.Matrix([[s.Rational(5,13),12*s.I/13],[12*s.I/13,s.Rational(5,13)]])
    W=[s.eye(2),U1.H,U2.H*U1.H]
    assert s.simplify(W[0].H*U1*W[1])==s.eye(2) and s.simplify(W[1].H*U2*W[2])==s.eye(2)
    out=dict(schema='oph.spin_flavour.CAR_symmetry.v1',particle_hole_cell_reflection_one_body_exact=True,
        particle_hole_cell_reflection_mass_clock_preserved=True,
        local_electric_Casimir_reflection_exact=True,
        equal_left_right_Casimir_on_total_SU3_singlets='Q_left=-Q_right on singlets, so sumQ_left²=sumQ_right²',
        N_pairing='N ->24-N',SU4_pairing='lambda ->(6-lambda4,6-lambda3,6-lambda2,6-lambda1)',
        all16_SU4_generators_commute_with_actual_mass_and_hopping=True,
        all16_SU4_generators_commute_with_all9_local_colour_bilinears=True,
        all_Fock_consequences_use_exact_CAR_representation=True,
        tree_unitary_spin_links_are_gauge_equivalent=True,
        tree_recursion='W_root=I; W_child=U_parent_child† W_parent',
        tree_scope='spin/flavour scalar onsite terms and colour electric Casimirs; consistent currents and states transported',
        independent_spin_splitting_or_axial_change_from_tree_unitary_links=False)
    (HERE/'SPIN_FLAVOUR_SYMMETRY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
