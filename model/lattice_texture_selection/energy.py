"""Source bond identities and finite energy controls for lattice textures."""
from fractions import Fraction as F
from itertools import product

import sympy as s

from .shells import integer


def positive(value):
    if type(value) not in (int, F) or value <= 0:
        raise ValueError('strictly positive exact rational input required')
    return F(value)


def shell_cost(degrees, exchange=F(1), kappa=F(1, 2000000)):
    J, k = positive(exchange), positive(kappa)
    if type(degrees) not in (list, tuple) or any(type(d) is not int for d in degrees):
        raise ValueError('complete integer shell degrees required')
    return J*k*sum(abs(d) for d in degrees)  # coefficient of pi


def mixed_shell_cost(nonacute_count, degrees, exchange=F(1), kappa=F(1, 2000000)):
    if type(nonacute_count) is not int or nonacute_count < 0:
        raise ValueError('a nonnegative exact nonacute-shell count required')
    J, k = positive(exchange), positive(kappa)
    return {'constant': J*k*nonacute_count/2,
            'coefficient_of_pi': shell_cost(degrees, J, k)}


def cut_sizes(side):
    side = integer(side, positive=True)
    vertices = set(product(range(side), repeat=3))
    done, cut, values = set(), 0, [0]
    for point in sorted(vertices):
        for axis in range(3):
            for sign in (-1, 1):
                other = list(point)
                other[axis] += sign
                other = tuple(other)
                if other in vertices:
                    cut += -1 if other in done else 1
        done.add(point)
        values.append(cut)
    if cut != 0 or min(values) < 0 or max(values) > side*side+side+1:
        raise ArithmeticError('lexicographic cut bound failed')
    return tuple(values)


def deformation_extra_bound(side, exchange=F(1)):
    side = integer(side, positive=True)
    return positive(exchange)*(side*side+side+7)


def rotation(quaternion):
    if len(quaternion) != 4:
        raise ValueError('four real exact quaternion coordinates required')
    q = [s.sympify(x) for x in quaternion]
    if any(x.has(s.Float) or x.is_real is not True for x in q):
        raise ValueError('four real exact quaternion coordinates required')
    a, b, c, d = q
    norm = sum(x*x for x in q)
    if norm == 0:
        raise ValueError('nonzero quaternion required')
    return s.Matrix([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])/norm


def bond_excess(left, right, a, b):
    return s.simplify(s.trace(s.diag(a, b, b)*(left.T*right-s.eye(3))))


def certificate():
    a, b = s.symbols('a b', real=True)
    q = s.symbols('q0:4', real=True)
    R = rotation(q)
    kappa = -(a+b)
    exact = s.cancel(s.trace(s.diag(a, b, b)*(R-s.eye(3))))
    lower = kappa*(1-R[0, 0])
    remainder = -b*(1+R[0, 0]-R[1, 1]-R[2, 2])
    if s.cancel(exact-lower-remainder) != 0:
        raise ArithmeticError('actual anisotropic covariance identity failed')
    if s.cancel(remainder+4*b*q[1]**2/sum(x*x for x in q)) != 0:
        raise ArithmeticError('nonnegative frame contribution changed')
    sigma = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    H = sum((s.kronecker_product(x, x)/4 for x in sigma), s.zeros(4))-s.eye(4)/4
    if H.eigenvals() != {-1: 1, 0: 3}:
        raise ArithmeticError('the complete source bond spectrum changed')
    alpha, beta, gamma = s.symbols('alpha beta gamma', nonnegative=True)
    u, v, w = [s.Matrix(s.symbols(name+'0:3', real=True)) for name in ('u', 'v', 'w')]
    zero = s.expand((alpha*u+beta*v+gamma*w).dot(v+w)
                    -alpha*u.dot(v+w)-beta*(v.dot(v)+v.dot(w))-gamma*(w.dot(w)+v.dot(w)))
    if zero != 0:
        raise ArithmeticError('convex-crossing identity changed')
    crossing = (s.Matrix([0, 1, 0]), s.Matrix([1, 0, 0]), s.Matrix([-1, 0, 0]))
    crossing_chords = sum((crossing[0]-x).dot(crossing[0]-x) for x in crossing[1:])
    if crossing_chords != 4 or crossing[1]+crossing[2] != s.zeros(3, 1):
        raise ArithmeticError('sharp local crossing chord threshold failed')
    cuts = [{'side': n, 'maximum_cut': max(cut_sizes(n)), 'cut_upper': n*n+n+1,
             'extra_energy_upper_at_J1': str(deformation_extra_bound(n))} for n in range(1, 9)]
    return {'uniform_reference_energy_per_bond': 'J*(a+2b-1/4)',
            'exact_covariance_identity': 'J*Tr[C*(Rx^T*Ry-I)]',
            'nonnegative_frame_remainder': '-4*J*b*q1^2/|q|^2',
            'bond_energy_spectrum_at_J1': {'-1': 1, '0': 3},
            'each_rotated_ground_bond_excess_between_zero_and_J': True,
            'actual_ground_kappa_lower': 'kappa>=m_order^2/2>1/2000000',
            'charged_shell_energy_lower': 'pi*J*kappa*abs(q)',
            'nonacute_shell_energy_strict_lower': 'J*kappa/2',
            'mixed_shell_lower': 'J*kappa*N_bad/2+pi*J*kappa*sum_good|q|',
            'finite_energy_forces_eventual_admissibility_and_zero_degree': True,
            'four_unit_charge_shells_lower_coefficient_of_pi_at_order_bound': str(shell_cost([1]*4)),
            'opposite_two_shell_families_lower_coefficient_of_pi_at_order_bound': str(shell_cost([1, 1, -1, -1])),
            'crossing_two_source_chord_squares_lower': 4,
            'crossing_absolute_energy_floor': '2*J*kappa>J/1000000',
            'crossing_floor_is_barrier_above_initial_energy': False,
            'finite_lexicographic_cut_controls': cuts,
            'deformation_peak_upper': 'initial_energy+J*(N^2+N+7)',
            'deformation_is_autonomous_source_evolution': False}
