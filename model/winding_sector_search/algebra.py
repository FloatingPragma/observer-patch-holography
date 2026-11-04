"""Exact full-valley comparison with the three-dimensional Callias operator."""
from functools import lru_cache
from itertools import product

import sympy as s

from ..native_mass_channels.algebra import SIGMA
from ..topological_defect_spectrum import algebra as source
from ..radial_defect_gap.bounds import spectral_window, encode


def winding(value):
    if type(value) is not int or value == 0:
        raise ValueError('a nonzero integer degree required')
    return value


def _fock_basis(gamma, beta):
    annihilators = [(b+s.I*g)/2 for g, b in zip(gamma, beta, strict=True)]
    vacuum = s.eye(16)
    for a in annihilators:
        vacuum *= a*a.H
    seeds = s.GramSchmidt(vacuum.columnspace(), orthonormal=True)
    if len(seeds) != 2:
        raise ArithmeticError('the full16-row Clifford multiplicity must be two')
    columns = []
    for occupation in product((0, 1), repeat=3):
        for seed in seeds:
            column = seed
            for j in reversed(range(3)):
                if occupation[j]:
                    column = annihilators[j].H*column
            columns.append(column)
    return s.ImmutableMatrix.hstack(*columns)


@lru_cache(maxsize=1)
def callias_intertwiner():
    eye = s.eye(2)
    gamma = tuple(s.kronecker_product(p, eye, SIGMA[0], eye) for p in SIGMA)
    beta = tuple(s.kronecker_product(eye, p, SIGMA[1], eye) for p in SIGMA)
    grading = s.kronecker_product(eye, eye, SIGMA[2], eye)
    sg, sb, _, _, sc = source.generators()
    change = source.internal_fock_basis()*_fock_basis(gamma, beta).H
    if change.H*change != s.eye(16) or change*change.H != s.eye(16):
        raise ArithmeticError('the Callias comparison must retain all sixteen rows')
    for actual, normal in zip(sg+sb, gamma+beta, strict=True):
        if actual*change != change*normal:
            raise ArithmeticError('a kinetic or mass generator fails to intertwine')
    if sc*change != change*grading:
        raise ArithmeticError('the chiral orientation fails to intertwine')
    return s.ImmutableMatrix(change)


def index_certificate():
    change = callias_intertwiner()
    n, u, v = (s.Matrix(s.symbols(prefix+'0:3', real=True)) for prefix in ('n', 'u', 'v'))
    pauli = lambda x: sum((x[j]*SIGMA[j] for j in range(3)), s.zeros(2))
    trace = s.expand(s.trace(pauli(n)*pauli(u)*pauli(v)))
    if s.expand(trace-2*s.I*n.dot(u.cross(v))) != 0:
        raise ArithmeticError('the oriented Callias boundary density changed')
    return {
        'source_dimension': 16, 'all_eight_valleys_retained': True,
        'Callias_4row_blocks_in_positive_to_negative_grading': 2,
        'intertwiner': [[str(x) for x in change.row(i)] for i in range(16)],
        'all_six_generators_and_grading_intertwined_exactly': True,
        'standard_gamma': 'sigma_j tensor I2 tensor sigma_x tensor I2',
        'standard_beta': 'I2 tensor sigma_j tensor sigma_y tensor I2',
        'standard_grading': 'I2 tensor I2 tensor sigma_z tensor I2',
        'positive_to_negative_block': 'i*(-v*sigma.dot(grad)+Phi) tensor I2',
        'Phi': 'm*N(x).sigma',
        'trace_boundary_density': 'tr[(n.sigma)(u.sigma)(v.sigma)]=2i*n.(u cross v)',
        'Callias_L_index': '-degree(N)',
        'source_chiral_index': '2*degree(N)',
        'zero_mode_dimension_lower': '2*abs(degree(N))',
        'zero_mode_dimension_exhausted_for_arbitrary_degree': False,
        'analytic_dependency': 'https://arxiv.org/abs/1506.05144',
        'theorem': 'Gesztesy-Waurick, Callias index formula, equation(0.1) and Theorem10.2',
        'hypotheses': ['bounded Hermitian C2 mass after a compact core smoothing',
                       'invertible exterior mass', 'first derivative O(1/r), second derivative O(1/r^2)',
                       'self-adjoint H on H1, chiral Fredholm off-diagonal domain'],
        'singular_bounded_core_extension': 'relative compactness and Fredholm-index homotopy',
        'native_repair_preserves_grading': False,
        'index_is_observed_family_count': False,
    }


def degree_certificate(degree):
    degree = winding(degree)
    k = abs(degree)
    u = s.symbols('u', nonnegative=True)
    primitive = -1/(1+u)
    if s.diff(primitive, u) != 1/(1+u)**2:
        raise ArithmeticError('rational-map area primitive changed')
    area = s.limit(primitive, u, s.oo)-primitive.subs(u, 0)
    if area != 1:
        raise ArithmeticError('rational-map degree integral failed')
    return {
        'degree': degree, 'map': f'stereographic w -> w^{k}, target orientation {1 if degree > 0 else -1}',
        'degree_integral_after_u_equals_r_to_2k': f'{degree}*integral_0^infinity du/(1+u)^2',
        'integral_exact': '1',
        'source_continuum_chiral_index': 2*degree,
        'source_continuum_kernel_dimension_lower': 2*k,
        'antipodal_equivariant': k % 2 == 1,
        'lattice_inversion_index_under_exterior_gap_hypothesis': 2 if k % 2 else None,
    }


def strong_mass_certificate(degree, **parameters):
    degree = winding(degree)
    if degree % 2 == 0:
        raise ValueError('the inversion-protected theorem requires an odd degree')
    return {'degree': degree, 'bounds': encode(spectral_window(**parameters)),
            'required_texture': 'unit-length and antipodal-equivariant at every nonzero lattice site; B0=0',
            'profile_specific_radial_equation_used': False,
            'kernel_count_is_degree_dependent': False,
            'continuum_limit_at_fixed_m_v_stays_in_strong_mass_domain': False}
