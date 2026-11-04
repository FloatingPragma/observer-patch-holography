"""Finite support, exact shell degrees, and uniform source-energy bounds."""
from fractions import Fraction as F
from itertools import product

from model.lattice_texture_selection.shells import degree, shell


def separation(value):
    if type(value) is not int or value < 1:
        raise ValueError('positive integer half-separation required')
    return value


def direction(point, half_separation=1):
    a = separation(half_separation)
    if type(point) is not tuple or len(point) != 3 or any(type(x) is not int for x in point):
        raise ValueError('three integer lattice coordinates required')
    if max(map(abs, point)) > a+1:
        return (0, 0, 1)
    x, y, z = point
    return (x, y, z*z-a*a)


def topology(half_separation=1):
    a = separation(half_separation)
    rows = []
    for z in (-a, a):
        triangles, edges = shell(1, (0, 0, z))
        values = {p: direction(p, a) for t in triangles for p in t}
        value = degree(triangles, values)
        if value['degree'] != (1 if z > 0 else -1):
            raise ArithmeticError('opposite local degrees failed')
        rows.append(dict(center=[0, 0, z], radius=1, source_edges=len(edges), **value))
    triangles, _ = shell(a+2)
    exterior = degree(triangles, {p: direction(p, a) for t in triangles for p in t})
    if exterior['degree']:
        raise ArithmeticError('constant exterior must have degree zero')
    first = shell(1, (0, 0, -a))[1]
    second = shell(1, (0, 0, a))[1]
    side = 2*a+3
    support = set(product(range(-a-1, a+2), repeat=3))
    touched = set()
    for p in support:
        for axis in range(3):
            for sign in (-1, 1):
                q = tuple(x+sign*int(j == axis) for j, x in enumerate(p))
                touched.add(tuple(sorted((p, q))))
    if len(touched) != 3*side*side*(side+1):
        raise ArithmeticError('complete source bond count failed')
    return {'half_separation': a, 'support_radius': a+1, 'support_sites': side**3,
            'local_shells': rows, 'exterior_shell': exterior,
            'local_shell_edges_disjoint': not bool(first & second),
            'spin_energy_upper_coefficient_of_J': len(touched),
            'spin_energy_lower_coefficient_of_pi_J_kappa': 2 if not first & second else 1,
            'fermion_Fock_energy_difference_upper_coefficient_of_M': 2*(side**3-1),
            'bounds_are_uniform_in_exterior_volume': True,
            'bound_is_constituent_particle_mass': False}


def spectral_bounds(mass, half_separation=1, radius=3):
    a = separation(half_separation)
    if type(mass) not in (int, F) or mass < 8:
        raise ValueError('exact rational mass at least eight required')
    if type(radius) is not int or radius < a+2:
        raise ValueError('cube must include a constant exterior layer')
    m = F(mass)
    b2 = F(7, 4) if a == 1 else F(3, 2)
    r = b2/(m-4)
    outside = m-2
    leakage2 = b2/(outside-r)**2
    q = 2/(m-r)
    distance = radius+1-(a+1)
    error = 2*b2/(m-r-2)*q**(2*distance)
    return {'exact_isolated_band_dimension': 4,
            'band_abs_energy_upper': str(r), 'complement_abs_energy_lower': str(outside),
            'core_to_exterior_norm_squared': str(b2),
            'exterior_to_core_weight_ratio_upper': str(leakage2),
            'core_probability_lower': str(1/(1+leakage2)),
            'ordered_band_energy_infinite_volume_error_upper': str(error),
            'boundary_minimum_closed_path_length': 2*distance,
            'zero_eigenvalues_required_by_shell_degrees': False}

def leading_schur(half_separation=1):
    """Exact zero-energy one-hop exterior resolvent term, coefficient of 1/M."""
    import sympy as s
    a = separation(half_separation)
    sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    cores = ((0, 0, -a), (0, 0, a))
    rows = {}
    for index, core in enumerate(cores):
        for axis in range(3):
            for sign in (-1, 1):
                p = tuple(x+sign*int(j == axis) for j, x in enumerate(core))
                block = rows.setdefault(p, s.zeros(2, 4))
                block[:, 2*index:2*index+2] = -sign*s.I*sigma[axis]/2
    result = s.zeros(4)
    for point, block in rows.items():
        f = direction(point, a)
        norm = s.sqrt(sum(x*x for x in f))
        v = sum((s.Rational(f[j], 1)/norm*sigma[j] for j in range(3)), s.zeros(2))
        result -= block.H*v*block
    expected = s.zeros(4)
    if a == 1:
        expected[:2, 2:] = expected[2:, :2] = -sigma[2]/4
    if result != expected:
        raise ArithmeticError('exact leading Schur term or self-term cancellation failed')
    return {'half_separation': a, 'coefficient_of_inverse_M': [[str(x) for x in row] for row in result.tolist()],
            'eigenvalue_multiplicities': {str(x): count for x, count in result.eigenvals().items()},
            'all_single_core_diagonal_terms_cancel': True,
            'shared_neighbor_term_present': a == 1,
            'this_term_is_the_full_spectrum': False}
