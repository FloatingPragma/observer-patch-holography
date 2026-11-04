"""Exact metric variation, nonminimal Higgs force and Wald controls."""
from functools import lru_cache

import sympy as s

from .geometry import Geometry, clean


def _zero(value):
    return all(s.simplify(entry) == 0 for entry in value)


def _pair(left, right):
    """Twice the real Hermitian pairing, valid for real coordinate symbols."""
    return clean((left.H*right+right.H*left)[0])


@lru_cache(maxsize=1)
def metric_variation_certificate():
    """Independent lapse/scale Euler variation includes second metric jets."""
    t = s.Symbol('t', real=True)
    a, N, F = [s.Function(name)(t) for name in ('a', 'N', 'F')]
    da, dd, dN, dF = s.diff(a, t), s.diff(a, t, 2), s.diff(N, t), s.diff(F, t)
    lagrangian = 6*F*(a*a*dd/N+a*da*da/N-a*a*da*dN/N**2)
    lapse_euler = s.diff(lagrangian, N)-s.diff(s.diff(lagrangian, dN), t)
    scale_euler = (s.diff(lagrangian, a)-s.diff(s.diff(lagrangian, da), t)
                   +s.diff(s.diff(lagrangian, dd), t, 2))
    C00 = 3*F*da**2/a**2+3*da*dF/a
    C11 = (-F*(2*a*dd+da**2-2*a*da*dN/N)-a*a*s.diff(F, t, 2)
           -2*a*da*dF+a*a*dN*dF/N)/N**2
    if clean(lapse_euler-2*a**3*C00/N**2) != 0 or clean(scale_euler+6*N*C11) != 0:
        raise ArithmeticError('lapse/scale metric variation failed')
    return {'lapse_Euler_identity_exact': True, 'scale_Euler_identity_exact': True,
            'metric_second_derivatives_and_boundary_integration_retained': True,
            'action': 'sqrt(-g)*F(X)*R',
            'metric_equation': '2F*G_ab+2(g_ab*box-nabla_a_nabla_b)F=T_min_ab',
            'improvement': 'T_xi_ab=2xi[X*G_ab+(g_ab*box-nabla_a_nabla_b)X]',
            'constant_vacuum_energy_selected': False}


@lru_cache(maxsize=1)
def wald_certificate():
    eta = s.diag(-1, 1, 1, 1)
    u = s.Matrix([s.Rational(5, 3), s.Rational(4, 3), 0, 0])
    n = s.Matrix([s.Rational(4, 3), s.Rational(5, 3), 0, 0])
    if (u.T*eta*u)[0] != -1 or (n.T*eta*n)[0] != 1 or (u.T*eta*n)[0] != 0:
        raise ArithmeticError('binormal frame is not Lorentz normalized')
    lower_u, lower_n = eta*u, eta*n
    epsilon = lower_u*lower_n.T-lower_n*lower_u.T
    F = s.Symbol('F', real=True)
    contraction = sum(F*(eta[a, c]*eta[b, d]-eta[a, d]*eta[b, c])
                      *epsilon[a, b]*epsilon[c, d]/2
                      for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    projector = eta+lower_u*lower_u.T
    volume = sum(F*(eta[a, c]*eta[b, d]-eta[a, d]*eta[b, c])
                 *lower_u[a]*projector[b, c]*lower_u[d]/2
                 for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    if s.expand(contraction+2*F) != 0 or s.expand(volume-3*F/2) != 0:
        raise ArithmeticError('Wald or generalized-volume curvature contraction failed')
    return {'binormal_contraction': str(s.expand(contraction)),
            'Wald_area_coefficient': '4*pi*F0-4*pi*xi*X',
            'slice_curvature_contraction': str(s.expand(volume)),
            'reference_curvature_coefficient_E0': 'Fstar/2',
            'generalized_volume_density': '(3F-Fstar)/(2Fstar)',
            'entropy_of_a_source_state_identified': False}


@lru_cache(maxsize=1)
def curved_ward_certificate():
    """Direct coordinate replay, with curved metric and nonzero source forces."""
    t, x, y, z = s.symbols('t x y z', real=True)
    q = (t, x, y, z)
    metric = (1+t*t)**2*s.diag(-1, 1, 1, 1)
    geometry = Geometry(q, metric)
    H = s.Matrix([x+s.I*t, 1+t+x*x])
    J = s.Matrix([s.Rational(1, 3)+s.I*x, t+s.Rational(1, 5)])
    sigma1, sigma3 = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    A = (s.I*x*sigma3/3, s.I*t*sigma1/5, s.zeros(2), s.zeros(2))
    if any(a.H != -a for a in A) or _zero(A[0]*A[1]-A[1]*A[0]):
        raise ArithmeticError('noncommuting anti-Hermitian weak control lost')
    DH = tuple((H.diff(coordinate)+a*H).applyfunc(clean) for coordinate, a in zip(q, A, strict=True))
    DJ = tuple((J.diff(coordinate)+a*J).applyfunc(clean) for coordinate, a in zip(q, A, strict=True))
    curvature = [[(A[b].diff(q[a])-A[a].diff(q[b])+A[a]*A[b]-A[b]*A[a]).applyfunc(clean)
                  for b in range(4)] for a in range(4)]
    covbox = s.zeros(2, 1)
    for a in range(4):
        for b in range(4):
            covbox += geometry.inverse[a, b]*(DH[b].diff(q[a])+A[a]*DH[b]
                         -sum((geometry.gamma[c][a][b]*DH[c] for c in range(4)), s.zeros(2, 1)))
    X = clean((H.H*H)[0])
    xi, lam, w = s.Rational(1, 6), s.Rational(1, 4), s.S.One
    potential = lam*(X-w*w)**2
    yukawa = _pair(H, J)
    kinetic = clean(sum(geometry.inverse[a, b]*(DH[a].H*DH[b])[0]
                       for a in range(4) for b in range(4)))
    minimal = s.Matrix(4, 4, lambda a, b: clean(_pair(DH[a], DH[b])-
                        metric[a, b]*(kinetic+potential+yukawa)))
    improvement = 2*xi*geometry.curvature_response(X)
    EH = (covbox-2*lam*(X-w*w)*H-J-xi*geometry.scalar*H).applyfunc(clean)
    eom = s.Matrix([_pair(EH, DH[b]) for b in range(4)])
    force = s.Matrix([clean(sum(geometry.inverse[a, c]*_pair(DH[c], curvature[a][b]*H)
                              for a in range(4) for c in range(4))) for b in range(4)])
    exchange = s.Matrix([-_pair(H, DJ[b]) for b in range(4)])
    divergence = geometry.divergence(minimal+improvement)
    residual = (divergence-eom-force-exchange).applyfunc(clean)
    improvement_expected = s.Matrix([-xi*geometry.scalar*s.diff(X, coordinate) for coordinate in q])
    improvement_divergence = geometry.divergence(improvement)
    if not _zero(residual) or not _zero(improvement_divergence-improvement_expected):
        raise ArithmeticError('curved nonminimal Higgs Ward identity failed')
    if not _zero(geometry.divergence(geometry.einstein)):
        raise ArithmeticError('direct curved metric violates contracted Bianchi')
    point = {t: s.Rational(1, 3), x: s.Rational(1, 5), y: 0, z: 0}
    evaluations = {name: [str(clean(v.subs(point))) for v in values]
                   for name, values in [('total_divergence', divergence), ('Higgs_EOM_term', eom),
                                        ('nonabelian_gauge_force', force), ('Yukawa_exchange', exchange),
                                        ('improvement_divergence', improvement_divergence)]}
    if any(all(s.sympify(v) == 0 for v in evaluations[name])
           for name in ('nonabelian_gauge_force', 'Yukawa_exchange', 'improvement_divergence')):
        raise ArithmeticError('a required force-deletion control became trivial')
    return {'metric': '(1+t^2)^2*diag(-1,1,1,1)',
            'curvature_at_point': str(clean(geometry.scalar.subs(point))),
            'point': {str(k): str(v) for k, v in point.items()},
            'xi': str(xi), 'exact_off_shell_residual': list(map(str, residual)),
            'noncommuting_weak_connections': True,
            'field_equations_imposed': False,
            'scalar_sector_conserved_in_isolation': False,
            'nonzero_force_terms': evaluations,
            'fermion_bilinear_role': 'formal covariant source; exchange retained',
            'full_Grassmann_or_quantum_Ward_replay_claimed': False}
