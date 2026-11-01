"""Exact onsite mass algebra, bulk phases and finite texture-degree witnesses.

Coefficients describe a supplied quadratic perturbation of the Weyl hopping.
Neither contact-induced order nor a protected fermion defect is inferred.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd

import sympy as s


SITES = tuple(product((0, 1), repeat=3))
SIGMA = (s.ImmutableMatrix([[0, 1], [1, 0]]),
         s.ImmutableMatrix([[0, -s.I], [s.I, 0]]),
         s.ImmutableMatrix([[1, 0], [0, -1]]))
I2 = s.ImmutableMatrix.eye(2)
VERTICES = ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
            (0, -1, 0), (0, 0, 1), (0, 0, -1))


def rational(value):
    if type(value) is int:
        return F(value)
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or gcd(value.numerator, value.denominator) != 1):
        raise ValueError('an exact integer or canonical Fraction is required')
    return F(value.numerator, value.denominator)


def vector(values, size):
    if type(values) not in (tuple, list) or len(values) != size:
        raise ValueError(f'exactly {size} rational coordinates are required')
    return tuple(rational(x) for x in values)


def _sym(value):
    return s.Rational(value.numerator, value.denominator)


def _site(site):
    if (type(site) is not tuple or len(site) != 3
            or any(type(x) is not int for x in site)):
        raise ValueError('three integer lattice coordinates are required')
    return site


def onsite_basis(site):
    x = _site(site)
    return ((-1)**(sum(x) % 2)*I2,
            *((-1)**(x[j] % 2)*SIGMA[j] for j in range(3)))


def onsite(site, coefficients):
    m = vector(coefficients, 4)
    return s.ImmutableMatrix(sum((_sym(a)*b for a, b in zip(m, onsite_basis(site))), s.zeros(2)))


def edge_constraint_matrix():
    """Every real Pauli coefficient of every period-two directed edge."""
    rows = []
    for x in SITES:
        for j in range(3):
            y = tuple(v ^ (axis == j) for axis, v in enumerate(x))
            for a in range(4):
                row = [0]*32
                row[4*SITES.index(y)+a] += 1
                row[4*SITES.index(x)+a] += 1 if a in (0, j+1) else -1
                rows.append(row)
    return s.ImmutableMatrix(rows)


def exhaustion_certificate():
    constraints = edge_constraint_matrix()
    columns = []
    for a in range(4):
        columns.append([((-1)**(sum(x) % 2) if a == 0 else (-1)**x[a-1])
                        if b == a else 0 for x in SITES for b in range(4)])
    basis = s.ImmutableMatrix.hstack(*(s.ImmutableMatrix(c) for c in columns))
    rank = constraints.rank()
    if rank != 28 or basis.rank() != 4 or constraints*basis != s.zeros(96, 4):
        raise ArithmeticError('onsite edge-nullspace classification failed')
    return {'real_unknowns': 32, 'edge_equations': 96, 'constraint_rank': rank,
            'solution_dimension': 4, 'basis_columns': columns,
            'all_lattice_extension': 'B(x+e_j)=-sigma_j B(x) sigma_j',
            'extra_internal_states_added': False}


def source_anchor_certificate():
    """Reconstruct the actual open-cube one-body action, without eigentolerances."""
    from ..field_symmetry.symmetry import hidden_spin_matrix
    from ..full_fock_field.field import FieldGeometry

    def exact(matrix):
        return s.ImmutableMatrix([[s.Rational(float(z.real))+s.I*s.Rational(float(z.imag))
                                   for z in row] for row in matrix])

    geometry = FieldGeometry((2, 2, 2), speed=.125, interaction=1.)
    h = exact(geometry.single_h().toarray())
    if h != h.H or h.shape != (16, 16):
        raise ArithmeticError('source one-body domain changed')
    for a in range(4):
        mass = s.diag(*(onsite_basis(x)[a] for x in geometry.sites))
        if h*mass+mass*h != s.zeros(16):
            raise ArithmeticError('onsite family no longer anticommutes with the actual source')
    columns = []
    for site in range(8):
        for pauli in (I2, *SIGMA):
            local = s.zeros(16)
            local[2*site:2*site+2, 2*site:2*site+2] = pauli
            residual = h*local+local*h
            columns.append(s.Matrix([s.re(z) for z in residual]+[s.im(z) for z in residual]))
    rank = s.Matrix.hstack(*columns).rank()
    if rank != 28:
        raise ArithmeticError('actual source onsite nullspace dimension changed')
    for x in geometry.sites:
        for j in range(3):
            if onsite_basis(x)[j+1] != (-1)**sum(x)*exact(hidden_spin_matrix(x, j)):
                raise ArithmeticError('source hidden-frame Neel identity changed')
    return {'source': 'full_fock_field.field.FieldGeometry.single_h',
            'sides': [2, 2, 2], 'speed': '1/8', 'one_body_dimension': 16,
            'real_onsite_unknowns': 32, 'actual_source_constraint_rank': rank,
            'all_four_anticommutators_zero': True, 'source_hidden_frame_identity': True,
            'numeric_eigensolver_used': False}


def folded_generators():
    """Sixteen physical spin/valley rows; no valley projection."""
    kinetic, masses = [], []
    for j in range(3):
        matrix = s.zeros(16)
        for b in SITES:
            start = 2*SITES.index(b)
            matrix[start:start+2, start:start+2] = -(-1)**b[j]*SIGMA[j]
        kinetic.append(s.ImmutableMatrix(matrix))
    for a in range(4):
        matrix = s.zeros(16)
        flip = (1, 1, 1) if a == 0 else tuple(int(j == a-1) for j in range(3))
        for b in SITES:
            target = tuple(x ^ y for x, y in zip(b, flip))
            i, j = 2*SITES.index(target), 2*SITES.index(b)
            matrix[i:i+2, j:j+2] = I2 if a == 0 else SIGMA[a-1]
        masses.append(s.ImmutableMatrix(matrix))
    return tuple(kinetic), tuple(masses)


def folded_hamiltonian(kinetic_components, coefficients):
    """Inputs k_j=c sin(p_j), supplied as exact rationals for finite checks."""
    k, m = vector(kinetic_components, 3), vector(coefficients, 4)
    gamma, beta = folded_generators()
    return s.ImmutableMatrix(sum((_sym(a)*b for a, b in zip(k+m, gamma+beta)), s.zeros(16)))


def phase(coefficients):
    m0, *m = vector(coefficients, 4)
    r2 = sum((x*x for x in m), F())
    discriminant = m0*m0-r2
    if not discriminant:
        return 'gapless'
    if discriminant < 0:
        return 'spin_texture'
    return 'positive_scalar' if m0 > 0 else 'negative_scalar'


def spectrum(kinetic_components, coefficients):
    """Exact bulk levels, counting merged or zero levels with multiplicity."""
    k, m = vector(kinetic_components, 3), vector(coefficients, 4)
    r = s.sqrt(_sym(sum((x*x for x in m[1:]), F())))
    kinetic2 = _sym(sum((x*x for x in k), F()))
    levels = {}
    for sign in (-1, 1):
        energy = s.sqrt(kinetic2+(_sym(m[0])+sign*r)**2)
        for charge in (-1, 1):
            key = s.simplify(charge*energy)
            levels[key] = levels.get(key, 0)+4
    return levels


def bulk_gap_squared(coefficients):
    m = vector(coefficients, 4)
    r = s.sqrt(_sym(sum((x*x for x in m[1:]), F())))
    return s.simplify((abs(_sym(m[0]))-r)**2)


def component_retraction(coefficients, parameter):
    """Exact rational deformation retract within the declared mass family.

    The scalar cones retract to (+/-1,0,0,0). The central component retracts
    to m0=0, |m|_1=1, an octahedral sphere homeomorphic to S2.
    """
    m, t = vector(coefficients, 4), rational(parameter)
    label = phase(m)
    if not 0 <= t <= 1 or label == 'gapless':
        raise ValueError('a gapped coefficient and 0<=parameter<=1 are required')
    if label == 'spin_texture':
        scale = 1-t+t/sum(abs(x) for x in m[1:])
        return ((1-t)*m[0], *(scale*x for x in m[1:]))
    target = 1 if label == 'positive_scalar' else -1
    return ((1-t)*m[0]+t*target, *((1-t)*x for x in m[1:]))


def clifford_certificate():
    gamma, beta = folded_generators()
    identity, zero = s.eye(16), s.zeros(16)
    if any(a*a != identity or a != a.H for a in gamma+beta):
        raise ArithmeticError('a physical folded generator changed')
    if any(a*b+b*a != zero for a in gamma for b in beta):
        raise ArithmeticError('mass does not anticommute with original kinetic hopping')
    if any(a*b+b*a != zero for family in (gamma, beta[1:])
           for i, a in enumerate(family) for b in family[i+1:]):
        raise ArithmeticError('Clifford relation failed')
    if any(beta[0]*b != b*beta[0] for b in beta[1:]):
        raise ArithmeticError('scalar channel commutation changed')
    if (any(a.trace() != 0 for a in gamma+beta)
            or any((beta[0]*b).trace() != 0 for b in beta[1:])
            or any((a*beta[0]*b).trace() != 0 for a in gamma+beta for b in beta[1:])):
        raise ArithmeticError('signed-level multiplicity trace identities failed')
    return {'folded_dimension': 16, 'weyl_valleys': 8,
            'positive_chiralities': 4, 'negative_chiralities': 4,
            'mass_sector_dimension': 8, 'signed_level_trace_identities_exact': True,
            'generic_signed_level_multiplicities': [4, 4, 4, 4],
            'dispersion': '+/-sqrt(c^2 sum_j sin(p_j)^2+(m0 +/- |m|)^2)',
            'bulk_gap': 'abs(abs(m0)-|m|)',
            'mass_generated_by_contact': False}


def dirac_decomposition_certificate():
    """Four invariant irreducible four-row blocks for the third spin mass."""
    gamma, beta = folded_generators()
    generators = gamma+(beta[3],)
    blocks = []
    for bx, by in product((0, 1), repeat=2):
        indices = [2*SITES.index((bx, by, bz))+spin for bz in (0, 1) for spin in (0, 1)]
        outside = [i for i in range(16) if i not in indices]
        if any(g.extract(outside, indices) != s.zeros(12, 4) for g in generators):
            raise ArithmeticError('pure spin mass Dirac block is not invariant')
        small = [g.extract(indices, indices) for g in generators]
        words = []
        for mask in range(16):
            word = s.eye(4)
            for j, generator in enumerate(small):
                if mask & (1 << j):
                    word *= generator
            words.append(s.Matrix(list(word)))
        if s.Matrix.hstack(*words).rank() != 16:
            raise ArithmeticError('four-row Clifford block is not irreducible')
        blocks.append(indices)
    return {'reference_spin_axis': 2, 'dirac_block_count': 4, 'block_dimension': 4,
            'block_rows': blocks, 'each_block_matrix_algebra_dimension': 16,
            'other_spin_directions': 'unitary rotations generated by beta_i beta_j commute with every gamma',
            'physical_species_identified': False}


def wall_rotation(parameter, magnitude=1):
    """Rational constant-gap path from +x to -x through +y, t in [0,2].

    This joins uniform bulk coefficients. It is not a global spectral-gap
    theorem for a spatially varying interface Hamiltonian.
    """
    t, m = rational(parameter), rational(magnitude)
    if not 0 <= t <= 2 or m <= 0:
        raise ValueError('0<=parameter<=2 and positive magnitude are required')
    if t <= 1:
        x, y = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    else:
        u = t-1
        x, y = -2*u/(1+u*u), (1-u*u)/(1+u*u)
    return F(), m*x, m*y, F()


def octahedron_faces():
    faces = []
    for signs in product((0, 1), repeat=3):
        face = (signs[0], 2+signs[1], 4+signs[2])
        matrix = s.Matrix.hstack(*(s.Matrix(VERTICES[i]) for i in face))
        if matrix.det() < 0:
            face = face[0], face[2], face[1]
        faces.append(face)
    return tuple(faces)


def texture_degree(values, direction=(1, 2, 3)):
    """Degree of the radial PL octahedral map by exact oriented ray count.

    This bounded evaluator admits nonzero face determinants and a ray
    strictly off every image edge. It never rounds a numerical solid angle.
    """
    if type(values) not in (tuple, list) or len(values) != 6:
        raise ValueError('six ordered octahedral target vertices are required')
    vertices = tuple(vector(v, 3) for v in values)
    ray = vector(direction, 3)
    if not any(ray):
        raise ValueError('nonzero ray required')
    degree, intersections = 0, []
    for face in octahedron_faces():
        matrix = s.Matrix.hstack(*(s.Matrix([_sym(x) for x in vertices[i]]) for i in face))
        determinant = matrix.det()
        if not determinant:
            raise ValueError('the exact degree evaluator requires nondegenerate image faces')
        inverse_ray = matrix.inv()*s.Matrix([_sym(x) for x in ray])
        if any(x == 0 for x in inverse_ray):
            raise ValueError('choose a ray off image face edges')
        if all(x > 0 for x in inverse_ray):
            sign = 1 if determinant > 0 else -1
            degree += sign
            intersections.append({'face': list(face), 'sign': sign,
                                  'barycentric': [str(x/sum(inverse_ray)) for x in inverse_ray]})
    return {'degree': degree, 'intersections': intersections, 'faces_checked': 8,
            'fermionic_zero_mode_proved': False, 'repair_protection_proved': False}


def hedgehog_stability(values, orientation=1):
    """A rational separating-halfspace certificate for a degree +/-1 map.

    Each vertex is within delta<1/3 in coordinate max norm of orientation
    times its reference vertex. Every triangular linear homotopy stays
    separated from zero by the reference face normal: n.F >=1-3delta>0.
    """
    if type(orientation) is not int or orientation not in (-1, 1):
        raise ValueError('orientation must be exactly +1 or -1')
    if type(values) not in (tuple, list) or len(values) != 6:
        raise ValueError('six ordered octahedral target vertices are required')
    vertices = tuple(vector(v, 3) for v in values)
    delta = max(abs(a-orientation*b) for row, base in zip(vertices, VERTICES)
                for a, b in zip(row, base))
    margin = 1-3*delta
    if margin <= 0:
        raise ValueError('vertex perturbation exceeds the certified strict margin')
    return {'degree': orientation, 'maximum_coordinate_perturbation': str(delta),
            'face_separation_lower': str(margin), 'homotopy_avoids_zero': True,
            'scope': 'continuous interpolated mass-direction texture, not a fermion index'}


def certificate():
    return {'schema': 'oph.native_mass_channels.algebra.v1',
            'exhaustion': exhaustion_certificate(), 'clifford': clifford_certificate(),
            'source_anchor': source_anchor_certificate(),
            'pure_spin_dirac_decomposition': dirac_decomposition_certificate(),
            'gapped_components': {
                'positive_scalar': {'condition': 'm0>|m|', 'homotopy_type': 'point', 'pi1': '0', 'pi2': '0'},
                'negative_scalar': {'condition': 'm0<-|m|', 'homotopy_type': 'point', 'pi1': '0', 'pi2': '0'},
                'spin_texture': {'condition': '|m0|<|m|', 'homotopy_type': 'S2', 'pi1': '0', 'pi2': 'Z'}},
            'gapped_component_count': 3,
            'classification_domain': 'Hermitian onsite matrices anticommuting with the fixed Weyl hopping',
            'hedgehog': texture_degree(VERTICES),
            'antihedgehog': texture_degree(tuple(tuple(-x for x in v) for v in VERTICES)),
            'physical_species_identified': False, 'native_repair_fixed_algebra_extended': False,
            'spontaneous_order_proved': False, 'interacting_mass_shell_proved': False}
