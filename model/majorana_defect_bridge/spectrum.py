"""Exact Takagi data, actual K45 and the complete momentum CAR symbol."""
from dataclasses import dataclass
from functools import lru_cache
from math import isqrt

import sympy as s

from ..common_source import calibrated_matter as calibrated
from ..common_source.chiral_candidate import action, kinetic
from ..conditional_standard_model import spectrum as charged


def clean(value):
    return s.ImmutableSparseMatrix(value)


@dataclass(frozen=True)
class MajoranaInputs:
    """Supplied Takagi masses, exact unitary factors and reference Higgs norm.

    The factor list fixes U in M=conjugate(U)*diag(m)*U dagger. It need
    not have a connected mixing graph. A zero mass is retained as a Weyl
    sector; no right-neutrino field is appended.
    """
    masses: tuple
    unitary_factors: tuple
    reference_higgs_norm: object
    provenance: tuple

    def __post_init__(self):
        if type(self.masses) is not tuple or len(self.masses) != 3:
            raise ValueError('three exact nonnegative masses required')
        masses = tuple(charged._exact(m, nonnegative=True) for m in self.masses)
        w = charged._exact(self.reference_higgs_norm, nonnegative=True)
        if w.is_positive is not True:
            raise ValueError('positive exact reference Higgs norm required')
        if type(self.unitary_factors) is not tuple or not self.unitary_factors:
            raise ValueError('a nonempty tuple of exact unitary factors required')
        factors = tuple(charged._matrix(f) for f in self.unitary_factors)
        for factor in factors:
            if not charged._zero(factor.H*factor-s.eye(3)):
                raise ValueError('each declared mixing factor must be exactly unitary')
        if (type(self.provenance) is not tuple or not self.provenance
                or any(type(row) is not str or not row for row in self.provenance)):
            raise ValueError('explicit nonempty provenance strings required')
        object.__setattr__(self, 'masses', masses)
        object.__setattr__(self, 'unitary_factors', factors)
        object.__setattr__(self, 'reference_higgs_norm', w)

    def mixing(self):
        value = s.eye(3)
        for factor in self.unitary_factors:
            value *= factor
        return clean(value)

    def family_kernel(self):
        u = self.mixing()
        return clean(u.conjugate()*s.diag(*self.masses)*u.H)


@lru_cache(maxsize=1)
def calibrated_inputs():
    choices = calibrated.neutrino_parameters()['unmeasured_architecture_choices']
    phases = s.diag(1, *(s.exp(s.I*s.Rational(choices[key])/2)
                       for key in ('Majorana_phase21_radians', 'Majorana_phase31_radians')))
    return MajoranaInputs(tuple(calibrated.neutrino_mass_coordinates()),
                          (*calibrated._pmns_rotations(), phases), calibrated.neutral_higgs_norm(),
                          ('existing optional low-energy NuFIT normal-ordering calibration',
                           'lightest mass and Majorana phases are supplied architecture choices',
                           'phenomenological tree matching at the common Higgs background; no certified EFT running'))


def _inputs(value):
    value = calibrated_inputs() if value is None else value
    if type(value) is not MajoranaInputs:
        raise ValueError('MajoranaInputs required')
    return value


def _higgs(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError('two exact finite complex Higgs components required')
    return tuple(charged._exact(x) for x in value)


def majorana_kernel(higgs, inputs=None):
    inputs, higgs = _inputs(inputs), _higgs(higgs)
    ell = calibrated._LH_tensor()*s.Matrix(higgs)
    return clean(s.kronecker_product(ell*ell.T, inputs.family_kernel()/inputs.reference_higgs_norm**2))


def full_kernel(higgs=None, *, matter=None, neutrinos=None):
    matter = charged.calibrated_inputs() if matter is None else matter
    if type(matter) is not charged.ActionInputs:
        raise ValueError('ActionInputs required')
    higgs = _higgs((0, matter.vEW/s.sqrt(2)) if higgs is None else higgs)
    value = s.zeros(45)
    for name, y in matter.yukawas().items():
        components = calibrated._oriented_components(name)
        conjugate = action._representation()[4][name]['conjugate_Higgs']
        for h, component in zip(higgs, components, strict=True):
            coefficient = s.conjugate(h) if conjugate else h
            directed = coefficient*s.kronecker_product(component, y.conjugate())
            value += directed+directed.T
    return clean(value+majorana_kernel(higgs, neutrinos))


def pairing_kernel(kernel):
    if kernel.shape != (45, 45) or not charged._zero(kernel-kernel.T):
        raise ValueError('the full symmetric K45 required')
    return clean(s.kronecker_product(action.SPIN_EPSILON, kernel.conjugate()))


def nambu_symbol(momentum, *, matter=None, neutrinos=None):
    """The existing complete180-row finite symbol, with this actual K45."""
    return kinetic.bdg_symbol(momentum, kernel=full_kernel(matter=matter, neutrinos=neutrinos), speed=1.)


def universal_symbol_certificate():
    # A generic symmetric family block has all six complex entries. The
    # proof uses these dimension-independent spin identities for all45 rows.
    x, y, z, v = s.symbols('px py pz v', real=True)
    entries = s.symbols('k00 k01 k02 k11 k12 k22')
    a, b, c, d, e, f = entries
    k = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    epsilon = action.SPIN_EPSILON
    if any(g*epsilon+epsilon*g.T != s.zeros(2) for g in sigma):
        raise ArithmeticError('the Weyl spin contraction changed')
    h = -v*s.kronecker_product(x*sigma[0]+y*sigma[1]+z*sigma[2], s.eye(3))
    delta = s.kronecker_product(epsilon, k.conjugate())
    symbol = s.Matrix.vstack(s.Matrix.hstack(h, delta), s.Matrix.hstack(delta.H, h.T))
    expected = s.diag(s.kronecker_product(s.eye(2), k.conjugate()*k),
                      s.kronecker_product(s.eye(2), k*k.conjugate()))
    expected += v*v*(x*x+y*y+z*z)*s.eye(12)
    if not charged._zero(symbol*symbol-expected) or delta.T != -delta:
        raise ArithmeticError('the generic complex symmetric CAR shell identity failed')
    return {'generic_complex_symmetric_family_entries': 6,
            'spin_contraction_all_three_directions_exact': True,
            'opposite_momentum_hole_block': '-h(-p)^T=h(p)^T',
            'creation_pairing': 'epsilon_spin tensor conjugate(K45)',
            'squared_upper_block': '|p|^2*I90+I2 tensor (conjugate(K45)*K45)',
            'squared_lower_block': '|p|^2*I90+I2 tensor (K45*conjugate(K45))',
            'same_identity_for_every_symmetric_internal_dimension': True,
            'full_internal_dimension': 45, 'physical_CAR_annihilation_modes': 90,
            'full_Nambu_dimension': 180, 'Nambu_rows_are_independent_species': False,
            'additional_sterile_right_neutrinos': 0,
            'native_continuous_time_execution_or_interacting_pole_solver': False}


def _sqrt_interval(square, digits=14):
    square = s.Rational(square)
    if square < 0:
        raise ValueError('nonnegative rational square required')
    scale = 10**digits
    lower = isqrt(int(square.p)*scale*scale//int(square.q))
    if s.Rational(lower*lower, scale*scale) == square:
        upper = lower
    else:
        upper = lower+1
    lo, hi = s.Rational(lower, scale), s.Rational(upper, scale)
    if not lo*lo <= square <= hi*hi:
        raise ArithmeticError('exact mass interval failed')
    return [str(lo), str(hi)]


def mixing_support_certificate():
    """Phase-independent positive lower bounds for all calibrated PMNS entries."""
    packet = calibrated.neutrino_parameters()['oscillation_inputs']
    a, b, c = (s.Rational(packet['sin2_theta'+ij]) for ij in ('12', '23', '13'))
    # Mixed entries have |A+B exp(i delta)| >= ||A|-|B||.
    # With |A|,|B|<=1 this gives squared modulus >=(A²-B²)²/4.
    lower = [[(1-a)*(1-c), a*(1-c), c],
             [(a*(1-b)-(1-a)*b*c)**2/4,
              ((1-a)*(1-b)-a*b*c)**2/4, b*(1-c)],
             [(a*b-(1-a)*(1-b)*c)**2/4,
              ((1-a)*b-a*(1-b)*c)**2/4, (1-b)*(1-c)]]
    if any(value <= 0 for row in lower for value in row):
        raise ArithmeticError('the calibrated full mixing graph has no positive support bound')
    data = calibrated_inputs()
    r23, r13, r12 = calibrated._pmns_rotations()
    phase = s.exp(s.I*s.Rational(packet['delta_degrees'])*s.pi/180)
    sa, sb, sc, ca, cb, cc = map(s.sqrt, (a, b, c, 1-a, 1-b, 1-c))
    expected = s.Matrix([[ca*cc, sa*cc, sc/phase],
                        [-sa*cb-ca*sb*sc*phase, ca*cb-sa*sb*sc*phase, sb*cc],
                        [sa*sb-ca*cb*sc*phase, -ca*sb-sa*cb*sc*phase, cb*cc]])
    if not charged._zero(r23*r13*r12-expected):
        raise ArithmeticError('the exact PMNS entry formulas changed')
    return {'squared_modulus_lower': [[str(x) for x in row] for row in lower],
            'nonzero_entries': 9, 'charged_lepton_neutrino_graph_connected': True,
            'bounds_hold_for_every_Dirac_and_Majorana_phase_with_these_angles': True,
            'arbitrary_unitary_mixing_has_this_support': False,
            'supplied_mixing_factor_count': len(data.unitary_factors)}


def takagi_block_certificate(kernel, data, matter):
    """Check the actual45-row congruence and its90-row CAR basis map."""
    u = data.mixing()
    if (not charged._zero(u.H*u-s.eye(3))
            or not charged._zero(u.T*data.family_kernel()*u-s.diag(*data.masses))):
        raise ArithmeticError('the actual supplied Takagi congruence failed')
    indices = action.ordering()['massless_neutrino_internal_family_indices']
    w = s.eye(45)
    for i, source_i in enumerate(indices):
        for j, source_j in enumerate(indices):
            w[source_i, source_j] = u[i, j]
    w = clean(w)
    expected = s.MutableSparseMatrix(charged.mass_kernel(matter))
    scale = (matter.vEW/s.sqrt(2)/data.reference_higgs_norm)**2
    for i, mass in zip(indices, data.masses, strict=True):
        expected[i, i] = scale*mass
    transformed = w.T*kernel*w
    if not charged._zero(transformed-expected):
        raise ArithmeticError('the complete K45 Takagi block differs from the actual action')
    physical = clean(s.kronecker_product(s.eye(2), w))
    actual_pairing = physical.H*pairing_kernel(kernel)*physical.conjugate()
    expected_pairing = s.kronecker_product(action.SPIN_EPSILON, expected.conjugate())
    if not charged._zero(actual_pairing-expected_pairing):
        raise ArithmeticError('the complete physical CAR transformation changed conjugation')
    sigma = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    for g in sigma:
        principal = s.kronecker_product(g, s.eye(45))
        if not charged._zero(physical*principal-principal*physical):
            raise ArithmeticError('a Takagi factor changed the common kinetic principal symbol')
    return {'actual_family_congruence_U_transpose_M_U_diagonal': True,
            'complete_K45_block_congruence_exact': True,
            'complete_Delta90_creation_conjugation_exact': True,
            'all_three_kinetic_principal_matrices_intertwined': True,
            'neutrino_internal_indices': list(indices),
            'charged_K45_action_retained': True,
            'unitary_factors_are_inputs_not_a_source_selected_mixing_matrix': True}


def catalogue(inputs=None, *, matter=None):
    default = inputs is None and matter is None
    data = _inputs(inputs)
    matter = charged.calibrated_inputs() if matter is None else matter
    if type(matter) is not charged.ActionInputs:
        raise ValueError('ActionInputs required')
    w = matter.vEW/s.sqrt(2)
    scale = s.simplify((w/data.reference_higgs_norm)**2)
    kernel = full_kernel(matter=matter, neutrinos=data)
    delta = pairing_kernel(kernel)
    if (not charged._zero(delta.T+delta)
            or not charged._zero(action.electric_charge().T*kernel+kernel*action.electric_charge())):
        raise ArithmeticError('the complete Majorana CAR/charge identity failed')
    rank = sum(m != 0 for m in data.masses)
    charged_rank = sum((2 if key == 'charged_lepton' else 6)*y.rank()
                       for key, y in matter.yukawas().items())
    rows = []
    for j, mass in enumerate(data.masses):
        value = s.simplify(scale*mass)
        rows.append({'label': 'nu_'+str(j+1), 'supplied_Takagi_column': j,
                     'action_mass_exact': str(value), 'action_mass_decimal_display': str(s.N(value, 15)),
                     'mass_squared_exact': str(s.simplify(value**2)),
                     'massive_Majorana': mass != 0, 'massless_Weyl': mass == 0,
                     'physical_on_shell_degrees': 2,
                     'signed_Nambu_multiplicity_at_nonzero_energy': 2,
                     'quadratic_shell': 'E^2=|p|^2+('+str(value)+')^2',
                     'physical_pole_mass': None, 'assigned_defect_class': None})
    result = {'input_provenance': list(data.provenance), 'reference_Higgs_norm': str(data.reference_higgs_norm),
              'actual_Higgs_norm': str(w), 'neutrino_mass_scaling': str(scale),
              'coefficient': 'kappa=Mnu(reference)/w_reference^2',
              'quadratic_Higgs_scaling_exact': True,
              'mixing': action._serialize(data.mixing()),
              'Takagi_congruence_from_exact_unitary_factors': True,
              'actual_action_intertwiner': takagi_block_certificate(kernel, data, matter),
              'neutrino_rows': rows, 'neutrino_rank': rank, 'full_K45_rank': int(charged_rank+rank),
              'full_K45_nullity': int(45-charged_rank-rank),
              'full_K45_symmetric_and_electric_neutral': True,
              'full_Delta90_antisymmetric': True,
              'neutrino_mass_phase_or_ordering_selected_by_source': False,
              'low_energy_to_common_scale_EFT_running_certified': False,
              'interacting_relativistic_continuum_or_pole_masses_proved': False,
              'physical_defect_to_field_intertwiner_constructed': False}
    if default:
        if (not charged._zero(kernel-calibrated.mass_kernel(include_neutrinos=True))
                or not charged._zero(data.mixing()-calibrated.pmns_matrix())):
            raise ArithmeticError('the bridge changed the retained calibrated action')
        result['same_existing_calibrated_K45_exact'] = True
        result['calibration'] = calibrated.neutrino_parameters()
        result['masses_eV_exact_intervals'] = [_sqrt_interval(s.simplify((m*10**9)**2)) for m in data.masses]
        result['mass_squared_splittings_eV2'] = [str(s.simplify((data.masses[j]**2-data.masses[0]**2)*10**18)) for j in (1, 2)]
        result['calibrated_mixing_support'] = mixing_support_certificate()
        # Retain the actual T+ action; no external PMNS table replaces it.
        labels = action.ordering()['all_left_components']
        n = next(row['index'] for row in labels if row['internal_type'] == 'L' and row['electric_Q'] == '0')
        e = next(row['index'] for row in labels if row['internal_type'] == 'L' and row['electric_Q'] == '-1')
        raw = calibrated.weak_raising_current().extract(range(3*n, 3*n+3), range(3*e, 3*e+3))
        actual = data.mixing().H*raw
        if raw != s.eye(3) or actual != data.mixing().H:
            raise ArithmeticError('the source charged-current PMNS dictionary changed')
        result['actual_neutrino_to_charged_lepton_Tplus_block'] = action._serialize(actual)
        result['Tplus_block_is_PMNS_dagger_exact'] = True
    return result
