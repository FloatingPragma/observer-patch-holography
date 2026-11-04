"""All-valley sampling into the original bounded hedgehog lattice source."""
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s

from ..native_mass_channels import algebra as source


def parameters(spacing, velocity=1, mass=1):
    values = tuple(source.rational(x) for x in (spacing, velocity, mass))
    if any(x <= 0 for x in values):
        raise ValueError('positive exact spacing, velocity and mass required')
    return tuple(float(x) for x in values)


def positions(points):
    values = np.asarray(points)
    if values.ndim != 2 or values.shape[1] != 3 or values.dtype.kind not in ('i', 'u'):
        raise ValueError('an integer N-by-three lattice coordinate array is required')
    return values


def sample(points, spacing, field):
    """J_a f(n)=a^(3/2) sum_b (-1)^(b.n) f_b(a*n), including all rows.

    The diagnostic field callable returns shape (N,16,K), where K counts
    independent columns. Analytic inclusion is proved in PROOF.md, not by
    finite-grid numerical acceptance.
    """
    points = positions(points)
    a, _, _ = parameters(spacing)
    values = np.asarray(field(a*points), dtype=complex)
    if values.ndim != 3 or values.shape[:2] != (len(points), 16):
        raise ValueError('the field must retain all sixteen folded rows and its columns')
    signs = 1-2*((points@np.asarray(source.SITES).T) % 2)
    return a**1.5*np.einsum('nb,nbsk->nsk', signs, values.reshape(len(points), 8, 2, -1))


def lattice_action(points, spacing, field, velocity=1, mass=1):
    """Literal original Weyl hops and bounded staggered onsite hedgehog."""
    points = positions(points)
    a, v, m = parameters(spacing, velocity, mass)
    sigma = tuple(np.array(x, dtype=complex) for x in source.SIGMA)
    values = sample(points, spacing, field)
    answer = np.zeros_like(values)
    length = np.linalg.norm(points, axis=1)
    directions = np.divide(points, length[:, None], out=np.zeros(points.shape, float),
                           where=length[:, None] != 0)
    for axis in range(3):
        step = np.eye(3, dtype=int)[axis]
        difference = sample(points+step, spacing, field)-sample(points-step, spacing, field)
        answer += 1j*v/(2*a)*np.einsum('st,ntk->nsk', sigma[axis], difference)
        coefficient = m*directions[:, axis]*(1-2*(points[:, axis] % 2))
        answer += coefficient[:, None, None]*np.einsum('st,ntk->nsk', sigma[axis], values)
    return answer


@lru_cache(maxsize=1)
def _matrices():
    gamma, masses = source.folded_generators()
    return tuple(np.array(x, dtype=complex) for x in gamma), tuple(np.array(x, dtype=complex) for x in masses[1:])


def sampled_continuum_action(points, spacing, field, gradient, velocity=1, mass=1):
    points = positions(points)
    _, v, m = parameters(spacing, velocity, mass)
    gamma, beta = _matrices()

    def action(x):
        values = np.asarray(field(x), dtype=complex)
        derivatives = np.asarray(gradient(x), dtype=complex)
        if derivatives.shape != (len(x), 3, 16, values.shape[-1]):
            raise ValueError('all three derivatives and sixteen continuum rows are required')
        length = np.linalg.norm(x, axis=1)
        directions = np.divide(x, length[:, None], out=np.zeros(x.shape, float),
                               where=length[:, None] != 0)
        answer = np.zeros_like(values)
        for axis in range(3):
            answer += -1j*v*np.einsum('st,ntk->nsk', gamma[axis], derivatives[:, axis])
            answer += m*directions[:, axis, None, None]*np.einsum('st,ntk->nsk', beta[axis], values)
        return answer

    return sample(points, spacing, action)


def certificate():
    gamma, masses = source.folded_generators()
    for axis in range(3):
        kinetic, mass = s.zeros(16), s.zeros(16)
        for column, valley in enumerate(source.SITES):
            flipped = tuple(bit ^ int(j == axis) for j, bit in enumerate(valley))
            row = source.SITES.index(flipped)
            kinetic[2*column:2*column+2, 2*column:2*column+2] = -(-1)**valley[axis]*source.SIGMA[axis]
            mass[2*row:2*row+2, 2*column:2*column+2] = source.SIGMA[axis]
        if kinetic != gamma[axis] or mass != masses[axis+1]:
            raise ArithmeticError('the original source valley or mass-shift convention changed')
    u, following = s.symbols('u0:3', real=True), s.symbols('w0:3', real=True)
    for parity in product((-1, 1), repeat=3):
        present = sum((parity[k]*u[k]*source.SIGMA[k] for k in range(3)), s.zeros(2))
        for axis in range(3):
            neighbor = sum((parity[k]*(-1 if k == axis else 1)*following[k]*source.SIGMA[k]
                            for k in range(3)), s.zeros(2))
            actual = source.SIGMA[axis]*neighbor+present*source.SIGMA[axis]
            expected = -parity[axis]*(following[axis]-u[axis])*s.eye(2)
            expected += sum((parity[k]*(following[k]-u[k])*source.SIGMA[axis]*source.SIGMA[k]
                             for k in range(3) if k != axis), s.zeros(2))
            if (actual-expected).applyfunc(s.expand) != s.zeros(2):
                raise ArithmeticError('a nondecaying hopping-mass anticommutator was omitted')
    return {'sampling_map': 'J_a f(n)=a^(3/2)*sum_b (-1)^(b.n)*f_b(a*n)',
            'source_valleys_retained': 8, 'continuum_rows_retained': 16,
            'onsite_mass_intertwining_exact': True,
            'smooth_compact_test_field_kinetic_consistency_order': 2,
            'finite_test_field_Gram_converges_to_continuum_Gram': True,
            'bounded_source_mass': 'm*sum_j (n_j/|n|)*(-1)^n_j*sigma_j; B_0=0',
            'bounded_lattice_essential_spectrum': '[-sqrt(m^2+3*v^2/a^2),-m] union [m,sqrt(m^2+3*v^2/a^2)]',
            'all_continuum_subgap_levels_approached_by_lattice_discrete_spectrum': True,
            'near_zero_lattice_multiplicity_lower_bound': 2,
            'near_each_signed_Coulomb_level_lattice_multiplicity_lower_bound': '2*n^2 for integer n>=2',
            'spectral_inclusion_quantifiers': 'For each continuum gap level E and each epsilon>0 with (E-epsilon,E+epsilon) inside (-m,m), every sufficiently small a has at least the continuum multiplicity of lattice eigenvalues in that interval.',
            'global_spectral_inclusion_rate_proved': False,
            'additional_lattice_gap_levels_excluded': False,
            'exact_bounded_lattice_zero_modes_proved': False,
            'bounded_lattice_and_continuum_level_multiplicities_identified': False,
            'supplied_texture_required': True}
