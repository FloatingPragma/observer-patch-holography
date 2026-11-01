"""Independent full-CAR and Clifford controls for the symmetry constraint."""
from itertools import product
import sympy as s


def annihilator(mode):
    if mode not in (0, 1):
        raise ValueError('two-mode control requires mode zero or one')
    a = s.zeros(4)
    for state in range(4):
        if state & (1 << mode):
            sign = (-1) ** ((state & ((1 << mode) - 1)).bit_count())
            a[state ^ (1 << mode), state] = sign
    return a


def full_fock_control(p, scalar):
    p = s.Rational(p); scalar = s.Rational(scalar)
    if not p > scalar >= 0:
        raise ValueError('control requires p>scalar>=0')
    a, b = annihilator(0), annihilator(1)
    n0, n1 = a.T*a, b.T*b
    H = (p+scalar)*n0 + (-p+scalar)*n1
    ground = s.eye(4)[:, 2]
    L = H - (-p+scalar)*s.eye(4)
    swap = s.Matrix([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,-1]])
    U = swap*(a+a.T)*(b+b.T)
    # All matrices here are real: Theta=U complex-conjugation is antiunitary.
    if U*ground == -ground:
        U = -U
    q = s.Matrix([[0,-1],[-1,0]])
    fields = (a,b)
    assert U*U.T == s.eye(4) and U*ground == ground
    assert all(U*field*U.T == sum((q[j,k]*fields[k].T for k in range(2)), s.zeros(4))
               for j,field in enumerate(fields))
    assert all(fields[j]*fields[k].T + fields[k].T*fields[j] == (s.eye(4) if j==k else s.zeros(4))
               for j in range(2) for k in range(2))
    return {
        'ground_energy': H[2,2], 'positive_GNS_generator': L,
        'CPT_generator_defect': U*L*U.T-L,
        'addition_cost': p+scalar, 'removal_cost': p-scalar,
        'Theta_matrix': U, 'charge_conjugation_matrix': q,
    }


def centering_failure():
    a, b = annihilator(0), annihilator(1)
    number = a.T*a+b.T*b
    H = a.T*a+3*b.T*b
    centered = H-2*number
    return {'original_generator_diagonal': list(H.diagonal()),
            'centered_generator_diagonal': list(centered.diagonal()),
            'original_reference': [1,0,0,0],
            'mean_before_centering': 2, 'mean_after_centering': 0,
            'same_reference_is_centered_ground': False}


def clifford_control():
    X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
    sigma=(X,Y,Z); identity=s.eye(2)
    alpha=[s.kronecker_product(Z,m) for m in sigma]
    beta=s.kronecker_product(X,identity)
    px,py,pz,m=s.symbols('px py pz m', real=True)
    H=sum((p*A for p,A in zip((px,py,pz),alpha)),s.zeros(4))+m*beta
    assert s.simplify(H*H-(px*px+py*py+pz*pz+m*m)*s.eye(4))==s.zeros(4)
    assert s.trace(H)==0 and all(A*beta+beta*A==s.zeros(4) for A in alpha)
    rows=[]
    for bits in product((0,1),repeat=3):
        chirality=-(-1)**sum(bits)
        orientation=[-(-1)**bit for bit in bits]
        proper_frame=[chirality*entry for entry in orientation]
        assert s.prod(proper_frame)==1
        rows.append({'valley':list(bits),'chirality':chirality,'proper_spin_frame_diagonal':proper_frame})
    return {'Dirac_mass_shell_polynomial':str(px*px+py*py+pz*pz+m*m),
            'Dirac_first_moment_trace':0, 'mass_term_is_scalar_energy_shift':False,
            'free_valleys':rows,'positive_chirality_count':4,'negative_chirality_count':4,
            'free_mass':0}


def certificate():
    rows=[]
    for p, scalar in ((2,0),(2,s.Rational(1,2))):
        control=full_fock_control(p,scalar)
        rows.append({'p':str(p),'scalar':str(scalar),
                     'positive_GNS_diagonal':[str(v) for v in control['positive_GNS_generator'].diagonal()],
                     'CPT_generator_defect_diagonal':[str(v) for v in control['CPT_generator_defect'].diagonal()],
                     'addition_cost':str(control['addition_cost']),
                     'removal_cost':str(control['removal_cost'])})
    return {'full_Fock_CAR_controls':rows,
            'centering_failure':{k:[str(x) for x in v] if isinstance(v,list) else v
                                 for k,v in centering_failure().items()},
            'Clifford_and_valley_controls':clifford_control()}
