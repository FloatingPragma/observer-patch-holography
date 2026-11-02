"""Exact full-Fock controls for exterior conditioning and entropy signs."""
import sympy as s


def annihilator(mode, modes=4):
    if type(mode) is not int or type(modes) is not int or not 0 <= mode < modes <= 8:
        raise ValueError('valid small CAR mode required')
    out = s.zeros(2**modes)
    for state in range(2**modes):
        if state & (1 << mode):
            out[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode)-1)).bit_count())
    return out


def partial_trace_exterior(matrix):
    if matrix.shape != (16, 16):
        raise ValueError('four-mode matrix required')
    return s.Matrix(4, 4, lambda i,j: sum(matrix[i+4*k,j+4*k] for k in range(4)))


def finite_example():
    R = s.Rational
    u = s.Matrix([R(3,5), 0, R(4,5), 0])
    v = s.Matrix([0, R(5,13), 0, R(12,13)])
    g = s.Matrix([0, 0, R(3,5), R(4,5)])
    P = u*u.T+v*v.T
    E = s.diag(1,1,0,0)
    A = P[:2,:2]
    F = (E*P*g)[:2,:]
    B = P[:2,2:]
    cs = [annihilator(i) for i in range(4)]
    cu = sum((u[i]*cs[i] for i in range(4)), s.zeros(16))
    cv = sum((v[i]*cs[i] for i in range(4)), s.zeros(16))
    cg = sum((g[i]*cs[i] for i in range(4)), s.zeros(16))
    vac = s.zeros(16,1);vac[0]=1
    psi = cu.T*cv.T*vac
    rho = psi*psi.T
    p = (g.T*P*g)[0]
    centered = cg.T*cg-p*s.eye(16)
    theta = s.Symbol('theta', real=True)
    local = partial_trace_exterior((s.eye(16)+theta*centered)*rho)
    expected = A-theta*F*F.T
    local_c = [annihilator(i,2) for i in range(2)]
    actual = s.Matrix(2,2,lambda i,j:s.trace(local*local_c[j].T*local_c[i]))
    four = s.trace(local*local_c[0].T*local_c[1].T*local_c[1]*local_c[0])
    return {'P':P,'A':A,'F':F,'B':B,'g':g,'occupation':p,'rho':local,
            'theta':theta,'actual_covariance':actual,'expected_covariance':expected,
            'four_point':four,'psi':psi,'global_density':rho}


def certificate():
    e=finite_example();t=e['theta']
    checks = {
        'pure_projection': e['P']**2==e['P'],
        'full_Fock_state_normalized': (e['psi'].T*e['psi'])[0]==1,
        'block_Douglas_identity': e['B']*e['B'].T==e['A']*(s.eye(2)-e['A']),
        'local_state_normalized': s.simplify(s.trace(e['rho'])-1)==0,
        'local_state_Hermitian': e['rho']==e['rho'].T,
        'rank_one_Wick_covariance': s.simplify(e['actual_covariance']-e['expected_covariance'])==s.zeros(2),
        'four_point_Wick_determinant': s.simplify(e['four_point']-e['expected_covariance'].det())==0,
        'no_quadratic_parameter_term': s.expand(e['expected_covariance'].det()).coeff(t,2)==0,
    }
    if not all(checks.values()):
        raise ArithmeticError('exact finite CAR conditioning control failed')
    # Nine existing directions, checked here independently rather than importing
    # any volatile predecessor receipt or its current validation booleans.
    R=s.Rational
    ns=[s.Matrix(x) for x in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),
        (R(3,5),R(4,5),0),(R(3,5),0,R(4,5)),(0,R(3,5),R(4,5))]]
    D=s.Matrix(9,9,lambda i,j:(R(5,4)+R(3,4)*ns[i].dot(ns[j]))**2)
    norm=max(sum(abs(x) for x in D.inv().row(i)) for i in range(9))
    if D.det()!=R(24407490807,2000000000) or norm!=R(68,5):
        raise ArithmeticError('observer response matrix differs')
    return {'exact_controls':checks,'regional_covariance_diagonal':[str(e['A'][i,i]) for i in range(2)],
            'exterior_occupation':str(e['occupation']),
            'cell_norm_and_occupation_are_distinct':True,
            'observer_response_determinant':str(D.det()),'observer_response_inverse_norm':str(norm),
            'global_affine_pure_vacuum_family_asserted':False}
