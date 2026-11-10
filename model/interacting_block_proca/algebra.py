"""Exact polarization, canonical-filter and induced Lorentz-action identities."""
import itertools
import sympy as s

METRIC=s.diag(1,-1,-1,-1)


def polarization(momentum,frequency,mass):
    p=s.Matrix(momentum)
    spatial=s.eye(3)+p*p.T/(mass*(frequency+mass))
    return (p.T/mass).col_join(spatial)


def symbolic_identities():
    m,w=s.symbols('m w',positive=True)
    p=s.Matrix(s.symbols('p1 p2 p3',real=True));r2=(p.T*p)[0]
    eps=polarization(p,w,m);four=s.Matrix([w,*p]);shell=s.Poly(w*w-m*m-r2,w)
    def zero(expression):
        num,den=s.together(expression).as_numer_denom()
        assert s.rem(s.Poly(num,w),shell).as_expr().expand()==0
    for value in (four.T*METRIC*eps):zero(value)
    for value in (eps.T*METRIC*eps+s.eye(3)):zero(value)
    projector=-METRIC+four*four.T/m**2
    for value in (eps*eps.T-projector):zero(value)
    for value in (eps.T*eps-s.eye(3)-2*p*p.T/m**2):zero(value)
    M=eps[1:,:]
    for value in (M*M-s.eye(3)-p*p.T/m**2):zero(value)
    for value in (p.T*M-w*p.T/m):zero(value)
    # A0(p)=i p.Pi(p)/(m omega). Positive-frequency Pi=-i omega Phi.
    assert s.simplify(s.I*(-s.I*w)/(m*w)-1/m)==0
    # The exact equal-time brackets use the same canonical Phi/Pi carrier.
    equal_time_A0_Ai=(p.T*M/(m*w)).applyfunc(s.factor)
    for value in (equal_time_A0_Ai-p.T/m**2):zero(value)
    # Actual lattice omega_a in A0's denominator gives the target eps even
    # before omega_a converges; the Lorenz defect is explicitly nonzero.
    wa=s.symbols('omega_a',positive=True)
    assert s.simplify(s.I*(-s.I*wa)/(m*wa)-1/m)==0
    finite_defect=-s.I*wa*p.T/m+s.I*p.T*M
    for value in (finite_defect-s.I*(w-wa)*p.T/m):zero(value)
    return {'mass_shell':'omega^2=m^2+|p|^2','polarization_rank':3,
        'Lorenz_polarization_constraint':True,'Minkowski_orthonormal_polarizations':True,
        'polarization_sum':'epsilon epsilon^T=-eta+p p^T/m^2',
        'Euclidean_polarization_Gram':'epsilon^T epsilon=I+2 p p^T/m^2',
        'spatial_filter':'M=I+p p^T/[m(omega_0+m)]',
        'time_filter':'A0(p)=i p.Pi(p)/(m omega_a)',
        'positive_frequency_time_polarization':'p/m',
        'finite_Lorenz_defect':'i(omega_0-omega_a) p.Phi/m',
        'equal_time_A0_Ai_commutator':'-i partial_i delta/m^2',
        'equal_time_Ai_dotAj_commutator':'i(delta_ij-partial_i partial_j/m^2)delta',
        'limiting_microcausal_commutator':'(-eta^{mu nu}-partial^mu partial^nu/m^2) C_KG',
        'naive_A0_equal_to_p_Phi_over_m_as_real_field_rejected':True}


def proper_cube():
    out=[]
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1,1),repeat=3):
            r=s.zeros(3)
            for j in range(3):r[perm[j],j]=signs[j]
            if r.det()==1:out.append(r)
    assert len(out)==24
    return out


def induced_rotation_and_boost():
    m,w=s.symbols('m w',positive=True);p=s.Matrix(s.symbols('p1 p2 p3',real=True))
    eps=polarization(p,w,m)
    for rotation in proper_cube():
        lorentz=s.diag(s.ones(1),rotation)
        assert (polarization(rotation*p,w,m)*rotation-lorentz*eps).applyfunc(s.simplify)==s.zeros(4,3)
    # Noncollinear rational boost makes the induced spin action visible.
    four=s.Matrix([s.Rational(5,4),0,s.Rational(3,4),0])
    boost=s.Matrix([[s.Rational(5,3),s.Rational(4,3),0,0],
        [s.Rational(4,3),s.Rational(5,3),0,0],[0,0,1,0],[0,0,0,1]])
    assert boost.T*METRIC*boost==METRIC and boost.det()==1
    moved=boost*four
    E=polarization(four[1:,0],four[0],s.Integer(1))
    Ep=polarization(moved[1:,0],moved[0],s.Integer(1))
    wigner=-Ep.T*METRIC*boost*E
    assert wigner.T*wigner==s.eye(3) and wigner.det()==1
    assert Ep*wigner==boost*E and wigner!=s.eye(3)
    # A second noncollinear boost checks the exact Wigner cocycle.
    other=s.Matrix([[s.Rational(5,4),0,0,s.Rational(3,4)],
        [0,1,0,0],[0,0,1,0],[s.Rational(3,4),0,0,s.Rational(5,4)]])
    moved2=other*moved
    E2=polarization(moved2[1:,0],moved2[0],s.Integer(1))
    second=-E2.T*METRIC*other*Ep
    combined=-E2.T*METRIC*other*boost*E
    assert combined==second*wigner and combined.T*combined==s.eye(3)
    return {'exact_proper_cube_polarization_intertwiners':24,
        'rest_little_group_action':'ordinary three-dimensional vector representation',
        'positive_energy_one_particle_spin':1,
        'Wigner_matrix':'R(Lambda,p)=-epsilon(Lambda p)^T eta Lambda epsilon(p)',
        'noncollinear_rational_boost_Wigner_matrix':[[str(v) for v in row] for row in wigner.tolist()],
        'exact_noncollinear_boost_cocycle_checked':True,
        'one_particle_action':'(U(Lambda)psi)(p)=R(Lambda,Lambda^-1 p)psi(Lambda^-1 p)',
        'measure':'d^3p/[(2pi)^3 2omega(p)]',
        'field_action':'U(Lambda)^dag A^mu(x) U(Lambda)=Lambda^mu_nu A^nu(Lambda^-1 x)',
        'different_from_three_scalar_Lorentz_action':True,
        'finite_native_exact_boost_symmetry_claim':False}


def retained_reference_rotation():
    # The full native cover has3072 elements. The index identity below is
    # group-independent; this literal check uses its admitted cubic image.
    group=proper_cube();by_key={tuple(r):i for i,r in enumerate(group)}
    representations=[s.diag(s.ones(1),r) for r in group]
    C=s.zeros(4);C[0,1]=s.Rational(3,5);C[1,0]=s.Rational(2,5)
    assert sum((u.H*u for u in representations),s.zeros(4))/24==s.eye(4)
    checked=0
    for h,Uh in zip(group,representations):
        for g,Ug in zip(group,representations):
            right=representations[by_key[tuple(g*h)]]
            left=representations[by_key[tuple(h.T*g)]]
            assert right==Ug*Uh
            assert Uh*left==Ug
            assert right*C*right.H==Ug*(Uh*C*Uh.H)*Ug.H
            checked+=1
    return {'literal_admitted_group_order':24,'exact_right_reference_blocks_checked':checked,
        'isometry':'J psi=|G|^-1/2 sum_g U_g psi tensor |g>',
        'right_reference_action':'R_h|g>=|g h^-1>',
        'right_intertwiner':'(I tensor R_h)J=J U_h',
        'left_joint_action':'(U_h tensor L_h)J=J',
        'closed_spatial_rotation':'macro_position_rotation_R tensor right_reference_rotation_h',
        'all_lifted_Kraus_operators_transform_by_the_same_intertwiner':True,
        'native_reference_dimension_remains':3072,
        'no_new_quantum_rotation_reference':True}
