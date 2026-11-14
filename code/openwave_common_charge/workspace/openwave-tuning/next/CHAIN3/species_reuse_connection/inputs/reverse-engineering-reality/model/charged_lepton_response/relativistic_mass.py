"""Same-coefficient free Dirac dispersion, including all lattice valleys."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import sympy as s

from model.common_source import native_matter_attachment as matter
from model.common_source import native_record_yukawa_field as field
from model.common_source.chiral_candidate import kinetic
from model.standard_model_source_map import koide_source
from model.native_flavor_orbit import intervals as iv


def _rational(value):
    if type(value) not in (str, int, F):
        raise ValueError('an exact rational is required')
    return F(value)


def continuum_bound(spacing, momentum_box, speed, time=0, mass_floor=None):
    """Uniform bounds for |p_j|<=P, with unchanged source mass coefficients."""
    a, p, c, t = map(_rational, (spacing, momentum_box, speed, time))
    if a <= 0 or p < 0 or c <= 0 or t < 0:
        raise ValueError('positive spacing/speed and nonnegative momentum/time required')
    delta = c*a*a*p**3/2
    out = {'symbol_norm_error_upper': str(delta),
           'each_energy_error_upper': str(delta),
           'unitary_norm_error_upper': str(min(F(2), t*delta)),
           'all_eight_valleys_retained': True}
    if mass_floor is not None:
        m = _rational(mass_floor)
        if m <= 0:
            raise ValueError('a proved positive charged mass floor is required')
        out['negative_energy_projector_error_upper'] = str(min(F(1), delta/m))
        out['free_two_point_operator_error_upper'] = str(min(F(2), delta*(t+1/m)))
    return out


def source_kernel(phase=1):
    """The actual native neutral charged-Higgs line, with all other rows zero."""
    rows = [r for r in matter.invariant_higgs_lines() if r['neutral_rank'] == 2]
    if len(rows) != 1:
        raise ArithmeticError('the unique native lepton-type Higgs line is required')
    row = rows[0]
    representation = matter.representation()
    selected = list(representation['selected'])
    weak = representation['operators']['weak_C2'].extract(selected, selected)
    left = s.diag(*(int(weak[i, i] != 0) for i in range(15)))
    d = left*row['components'][1]*(s.eye(15)-left)
    if s.trace(d.H*d) != 1 or d*d != s.zeros(15):
        raise ArithmeticError('the charged Higgs orientation differs')
    y = koide_source.exact_source(phase)['action_coefficient']
    directed = s.kronecker_product(d, y.conjugate())
    k = s.ImmutableSparseMatrix(directed+directed.T)
    if k != k.T:
        raise ArithmeticError('all-left pairing must be symmetric before spin epsilon')
    return k, y


def clifford_certificate():
    """Exact spin identities proving the square for every symmetric kernel."""
    sigmas = (s.Matrix([[0, 1], [1, 0]]),
              s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1))
    epsilon = s.Matrix([[0, 1], [-1, 0]])
    for i, j in product(range(3), repeat=2):
        if sigmas[i]*sigmas[j]+sigmas[j]*sigmas[i] != 2*int(i == j)*s.eye(2):
            raise ArithmeticError('spatial Clifford identity failed')
    if any(a*epsilon+epsilon*a.T != s.zeros(2) for a in sigmas):
        raise ArithmeticError('actual all-left kinetic/pairing cross term is nonzero')
    beta = s.kronecker_product(sigmas[0], s.eye(2))
    alpha = tuple(s.kronecker_product(sigmas[2], a) for a in sigmas)
    gamma = (beta, *(beta*a for a in alpha))
    signature = (1, -1, -1, -1)
    for i, j in product(range(4), repeat=2):
        if gamma[i]*gamma[j]+gamma[j]*gamma[i] != 2*signature[i]*int(i == j)*s.eye(4):
            raise ArithmeticError('the massive Dirac Clifford convention differs')
    return {'spatial_Clifford_checks': 9, 'spin_pairing_cross_term_checks': 3,
            'massive_Dirac_spacetime_Clifford_checks': 16,
            'massive_Dirac_signature': '+---',
            'dimension': 'three spatial directions and one common evolution parameter',
            'exact_square': 'H(p)^2=diag(c^2|p|^2 I+Delta Delta^dagger,c^2|p|^2 I+Delta^dagger Delta)',
            'creation_pairing': 'Delta=epsilon_spin tensor conjugate(K45)',
            'lower_Nambu_block': '-h(-p)^T',
            'free_energies': 'plus/minus sqrt(c^2|p|^2+w_H^2 sigma_j(Y)^2) on the charged block'}


def phase_zero_gap():
    value = (15-10*iv.sqrt_rational(F(2)))/216
    if value.lo <= F(1, 1000):
        raise ArithmeticError('the zero-phase example lacks its stated positive gap')
    return {'program_phase': '0', 'neutral_Higgs_norm': '1',
            'smallest_charged_action_mass_interval': value.strings(),
            'certified_mass_floor': '1/1000',
            'applies_to_other_phases_without_a_separate_check': False}


def certificate():
    controls = []
    for phase in (s.S.One, (s.Integer(99)+20*s.I)/101):
        k, y = source_kernel(phase)
        kn = np.array(k.evalf(17), complex)
        yn = np.array(y.evalf(17), complex)
        # This calls the owning full 180-dimensional source symbol, not a fitted spectrum.
        report = field.dispersion_certificate(kn, momentum=(.23, -.31, .17), speed=1)
        if report['Hermiticity_error'] > 1e-12 or report['dispersion_spectral_error'] > 1e-12:
            raise ArithmeticError('the full native neutral-Higgs symbol failed its numerical control')
        masses = np.linalg.eigvalsh(yn)
        expected = np.sort(np.r_[np.zeros(39), np.repeat(masses, 2)])
        if np.max(abs(np.sort(np.linalg.svd(kn, compute_uv=False))-expected)) > 1e-12:
            raise ArithmeticError('the native Higgs lift changed its family coefficients')
        for bits in product((0, 1), repeat=3):
            p = tuple((-1)**b*x for b, x in zip(bits, (.23, -.31, .17)))
            h = kinetic.bdg_symbol(p, kn, speed=1)
            if np.max(abs(np.linalg.eigvalsh(h)-np.linalg.eigvalsh(kinetic.bdg_symbol((.23, -.31, .17), kn, speed=1)))) > 1e-12:
                raise ArithmeticError('a retained valley changed the free mass spectrum')
        controls.append({'phase': str(phase), 'full_native_symbol_control': report,
                         'positive_charged_action_coefficients_approx': masses.tolist(),
                         'positive_charged_Nambu_energy_rows': 12,
                         'spectator_internal_zero_coefficients': 39,
                         'all_eight_valley_energy_controls': 8})
    return {'exact_spin_proof': clifford_certificate(), 'finite_controls': controls,
            'compact_momentum_example': continuum_bound('1/32', 2, 1, 3, '1/1000'),
            'compact_momentum_example_source_gap': phase_zero_gap(),
            'continuum_source_mass_rule': 'm_j=w_H sigma_j(Y); the lattice spacing introduces no additional mass factor',
            'reference': 'declared filled negative-energy Gaussian vacuum on the charged free block',
            'scope': {'same_coherent_family_coefficient_and_native_Higgs_line': True,
                'controlled_free_relativistic_valley_limit': True,
                'mirror_valleys_deleted': False,
                'interacting_gauge_field_continuum_proved': False,
                'Higgs_vacuum_selected_by_source': False,
                'native_global_attempt_clock_identified': False,
                'physical_charged_leptons_identified': False,
                'physical_measured_mass_ratios': [],
                'floating_controls_are_roundoff_certificates': False}}
