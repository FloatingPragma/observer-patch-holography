"""Exact source moments and controlled finite spectral coarse-graining."""
from fractions import Fraction as F
from math import isqrt

from model.color_loop_readout.source import Controls, Sector, sqrt_upper


def matrix_add(a, b, scale=F(1)):
    return [[x+scale*y for x, y in zip(r, s, strict=True)]
            for r, s in zip(a, b, strict=True)]


def matrix_product(a, b):
    return [[sum((x*y for x, y in zip(r, c, strict=True)), F(0))
             for c in zip(*b)] for r in a]


def matrix_scale(a, scale):
    return [[scale*x for x in row] for row in a]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def difference(a, b):
    return matrix_add(a, b, -1)


def serialized(matrix):
    return [[str(x) for x in row] for row in matrix]


class Elimination:
    """The entire 432-dimensional complement is included in every product."""

    def __init__(self, controls=Controls()):
        self.source = Sector(0, controls)
        self.controls = controls
        self.p = [i for i, row in enumerate(self.source.labels) if row['N_H'] == 0]
        self.q = set(range(self.source.dimension))-set(self.p)
        if self.p != list(range(6)) or len(self.q) != 432:
            raise ArithmeticError('complete scalar-vacuum split required')
        self.scales = []
        for i in self.p:
            n = self.source.norms[i]
            root = isqrt(n.numerator)
            if n.denominator != 1 or root*root != n:
                raise ArithmeticError('exact normalized flux coordinates required')
            self.scales.append(root)
        self.vectors = [{i:F(1, self.scales[j])} for j, i in enumerate(self.p)]
        hp = [self.source.apply(v) for v in self.vectors]
        self.a = [[self.source.inner(v, w) for w in hp] for v in self.vectors]
        self.c = [self.keep_q(v) for v in hp]
        self.c_squared = 12*controls.higgs**2
        gram = [[self.source.inner(v, w) for w in self.c] for v in self.c]
        if gram != matrix_scale(identity(6), self.c_squared):
            raise ArithmeticError('exact complete Higgs Gram identity failed')
        bounds = {i:F(0) for i in self.q}
        for (i, j), value in self.source.numerators.items():
            if i in self.q and j in self.q:
                bounds[i] += sqrt_upper(value*value/(self.source.norms[i]*self.source.norms[j]))
        self.d_bound = max(bounds.values())

    def keep_q(self, vector):
        return {i:x for i, x in vector.items() if i in self.q}

    def apply_d(self, vector):
        return self.keep_q(self.source.apply(vector))

    def moments(self, order):
        if type(order) is not int or order < 0:
            raise ValueError('nonnegative integer moment order required')
        powers = self.c
        result = []
        for k in range(order+1):
            moment = [[self.source.inner(v, w) for w in powers] for v in self.c]
            if moment != list(map(list, zip(*moment))):
                raise ArithmeticError('Schur moment must be Hermitian')
            result.append(moment)
            if k < order:
                powers = [self.apply_d(v) for v in powers]
        return result

    def enclosure(self, mu, order=24, moments=None):
        if type(mu) not in (int, F):
            raise TypeError('exact rational spectral resolution required')
        mu = F(mu)
        if mu <= self.source.norm_bound:
            raise ValueError('spectral resolution must exceed full-source norm bound')
        moments = self.moments(order) if moments is None else moments
        if len(moments) != order+1:
            raise ValueError('complete requested moment list required')
        sigma = matrix_scale(identity(6), F(0))
        for k, moment in enumerate(moments):
            sigma = matrix_add(sigma, moment, F((-1)**k, 1)/mu**(k+1))
        remainder = self.c_squared*self.d_bound**(order+1)/(mu**(order+1)*(mu-self.d_bound))
        effective = matrix_add(matrix_add(self.a, identity(6), mu), sigma, -1)
        return effective, remainder


ELECTRIC = [[F(0) for _ in range(6)] for _ in range(6)]
for index, value in enumerate((0, 8, F(9, 4), F(25, 4), 6, 6)):
    ELECTRIC[index][index] = F(value)
MAGNETIC = [[F(0) for _ in range(6)] for _ in range(6)]
for index, value in ((1, 3), (3, 6), (4, 3)):
    MAGNETIC[0][index] = MAGNETIC[index][0] = F(value)
INDUCED = [[F(0) for _ in range(6)] for _ in range(6)]
for i, row in zip((0, 2, 5), ((36, 12, 12), (12, 6, 0), (12, 0, 12)), strict=True):
    for j, value in zip((0, 2, 5), row, strict=True):
        INDUCED[i][j] = F(value)


def component_certificate():
    """Check coefficient identities, giving arbitrary-control polynomial laws."""
    base = Elimination()
    fields = {
        'electric':Sector(0, Controls(1, 0, 0, 0, 0)),
        'magnetic':Sector(0, Controls(0, -1, 0, 0, 0)),
        'hopping':Sector(0, Controls(0, 0, 1, 0, 0)),
        'higgs':Sector(0, Controls(0, 0, 0, 1, 0)),
        'scalar':Sector(0, Controls(0, 0, 0, 0, 1)),
    }
    c = [base.keep_q(fields['higgs'].apply(v)) for v in base.vectors]
    for name, action in (('electric', ELECTRIC), ('magnetic', MAGNETIC),
                         ('scalar', identity(6))):
        for j, vector in enumerate(c):
            actual = fields[name].apply(vector)
            expected = {}
            for k in range(6):
                for i, value in c[k].items():
                    expected[i] = expected.get(i, F(0))+action[k][j]*value
            if actual != {i:x for i, x in expected.items() if x}:
                raise ArithmeticError('same-source electric/loop/scalar intertwiner failed')
    if any(base.keep_q(fields['higgs'].apply({i:F(1)})) for i in base.q):
        raise ArithmeticError('occupied-scalar Higgs compression must vanish')
    if any(fields['hopping'].apply(v) for v in base.vectors):
        raise ArithmeticError('empty-matter hopping must vanish')
    tc = [fields['hopping'].apply(v) for v in c]
    if any(i not in base.q for v in tc for i in v):
        raise ArithmeticError('hopping must preserve the scalar complement')
    first = [[base.source.inner(v, w) for w in tc] for v in c]
    second = [[base.source.inner(v, w) for w in tc] for v in tc]
    if first != matrix_scale(identity(6), F(0)) or second != INDUCED:
        raise ArithmeticError('complete hopping Gram identity failed')
    return {'intertwiner_columns':18, 'occupied_scalar_Higgs_columns':432,
            'hopping_Gram_R':serialized(second),
            'identities_hold_for_all_real_source_coefficients':True}


def exact_moment_certificate(controls=Controls()):
    model = Elimination(controls)
    e, m, h, g, s = (getattr(controls, key) for key in ('electric', 'magnetic', 'hopping', 'higgs', 'scalar'))
    expected_a = matrix_add(matrix_scale(ELECTRIC, e), MAGNETIC, -m)
    if model.a != expected_a:
        raise ArithmeticError('source flux action differs')
    b = matrix_add(model.a, identity(6), s)
    expected = [matrix_scale(identity(6), 12*g*g), matrix_scale(b, 12*g*g),
                matrix_add(matrix_scale(matrix_product(b, b), 12*g*g), INDUCED, g*g*h*h)]
    if model.moments(2) != expected:
        raise ArithmeticError('full-complement moment identity failed')
    return {'controls':controls.packet(), 'P_dimension':6, 'Q_dimension':432,
            'source_norm_bound':str(model.source.norm_bound),
            'complement_norm_bound':str(model.d_bound),
            'C_adjoint_C':serialized(expected[0]),
            'C_adjoint_D_C':serialized(expected[1]),
            'C_adjoint_D_squared_C':serialized(expected[2]),
            'induced_operator_R':serialized(INDUCED)}


def linear_fit_weights():
    """Hilbert--Schmidt projection, without asserting closure of the action."""
    tr_e = sum(ELECTRIC[i][i] for i in range(6))
    tr_e2 = sum(ELECTRIC[i][i]**2 for i in range(6))
    det = 6*tr_e2-tr_e*tr_e
    if (tr_e, tr_e2, det) != (F(57, 2), F(1441, 8), F(537, 2)):
        raise ArithmeticError('source normalization changed')
    offset = matrix_scale(identity(6), tr_e2/det)
    offset = matrix_add(offset, ELECTRIC, -tr_e/det)
    electric = matrix_add(matrix_scale(ELECTRIC, 6/det), identity(6), -tr_e/det)
    magnetic = matrix_scale(MAGNETIC, F(-1, 108))
    return {'offset':offset, 'electric':electric, 'magnetic':magnetic}


def interval_packet(center, radius, bits=96):
    scale = 1 << bits
    low, high = center-radius, center+radius
    lo = (low*scale).numerator//(low*scale).denominator
    hi = -((-high*scale).numerator//(-high*scale).denominator)
    return [str(F(lo, scale)), str(F(hi, scale))]


def inverse_readout_error(mu, source_bound, entry_error):
    """Convert complete six-by-six resolvent readout errors to action errors."""
    mu, source_bound, entry_error = map(F, (mu, source_bound, entry_error))
    if mu <= source_bound or entry_error < 0:
        raise ValueError('positive spectral margin and nonnegative readout error required')
    floor = 1/(mu+source_bound)
    eta = 6*entry_error
    if eta >= floor:
        raise ValueError('readout error exceeds the invertibility margin')
    return eta/(floor*(floor-eta))


def uniform_family_remainders(mu, complement_bound, higgs, normalization=F(1)):
    """Actual finite-source Schur family, with explicit dimensionless margins."""
    values = (mu, complement_bound, higgs, normalization)
    if any(type(x) not in (int, F) for x in values):
        raise TypeError('exact rational action scales required')
    mu, bound, g, scale = map(F, values)
    if not mu > bound >= 0 or scale <= 0:
        raise ValueError('positive action and spectral margins required')
    full = 12*g*g/(scale*(mu-bound))
    after_linear = 12*g*g*bound*bound/(scale*mu*mu*(mu-bound))
    return {'normalized_full_self_energy_upper':str(full),
            'normalized_remainder_after_first_two_moments':str(after_linear),
            'source_family_or_scale_selection_inferred':False}


def witness(mu=F(64), order=24, controls=Controls()):
    model = Elimination(controls)
    effective, remainder = model.enclosure(mu, order)
    coarse = matrix_add(effective, identity(6), -F(mu))
    fit = {}
    for name, weights in linear_fit_weights().items():
        center = sum((weights[i][j]*coarse[i][j] for i in range(6) for j in range(6)), F(0))
        error = remainder*sum((abs(x) for row in weights for x in row), F(0))
        fit[name] = interval_packet(center, error)
    intervals = [[interval_packet(x, remainder) for x in row] for row in effective]
    e_lo, e_hi = map(F, fit['electric'])
    m_lo, m_hi = map(F, fit['magnetic'])
    ratio = None
    if 0 < e_lo <= e_hi and 0 < m_lo <= m_hi:
        lo, hi = e_lo/m_hi, e_hi/m_lo
        ratio = interval_packet((lo+hi)/2, (hi-lo)/2)
    return {'spectral_resolution_mu':str(mu), 'neumann_degree':order,
            'source_norm_bound':str(model.source.norm_bound),
            'complement_norm_bound':str(model.d_bound),
            'operator_norm_remainder_upper':interval_packet(remainder, F(0))[1],
            'complete_effective_matrix_intervals':intervals,
            'Hilbert_Schmidt_action_coefficients':fit,
            'projected_electric_over_positive_magnetic':ratio,
            'induced_vacuum_to_weak_flux':intervals[0][2],
            'induced_vacuum_to_charged_singlet_flux':intervals[0][5],
            'scalar_complement_states_eliminated_exactly_in_the_formula':432,
            'single_Yang_Mills_coupling_ansatz_used':False,
            'physical_alpha_s':None}
