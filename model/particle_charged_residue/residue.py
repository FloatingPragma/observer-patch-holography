"""Exact source-root counts, ground trial and complete rational window detector."""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.native_flavor_orbit.intervals import sqrt_rational
from model.particle_ground_spectral_records import records as inherited
from model.spectral_readout import kernels


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
SOURCE_FILES = (
    'model/particle_ground_spectral_records/records.py',
    'model/particle_ground_spectral_records/PROOF.md',
    'model/particle_field_state_response/state_response.py',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_candidate_atlas/interacting_chiral.py',
    'model/native_flavor_orbit/intervals.py',
    'model/spectral_readout/kernels.py',
)
SOURCE_INPUTS = SOURCE_FILES+('model/particle_charged_residue/source_pins.json',)
BITS = 48
GROUND_CUT = F(-2)
CHARGED_CUT = F(-3, 5)
TRIAL_SHIFT = F(-9, 4)
FILTER_WIDTH = F(1, 16)
FILTER_POWER = 4


def _q(value):
    value = F(value)
    return s.Rational(value.numerator, value.denominator)


def authenticate_sources():
    expected = json.loads((PACKAGE/'source_pins.json').read_text())
    actual = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCE_FILES}
    if actual != expected:
        raise ArithmeticError('an inherited spectral-source input differs from its frozen pin')
    return actual


@lru_cache(maxsize=256)
def _square_free_factors(matrix):
    x = s.Symbol('x')
    polynomial = s.Matrix(matrix).charpoly(x).as_poly()
    return tuple(polynomial.sqf_list()[1])


def root_count(matrix, lower=None, upper=None):
    """Exact eigenvalue multiplicity in an interval with nonspectral endpoints.

    Sympy counts distinct roots. Square-free decomposition restores each
    algebraic multiplicity, equal to eigenspace dimension for Hermitian H.
    """
    lo = -s.oo if lower is None else _q(lower)
    hi = s.oo if upper is None else _q(upper)
    if lo >= hi:
        raise ValueError('an increasing root-count interval is required')
    total = 0
    for polynomial, multiplicity in _square_free_factors(matrix):
        if (lo != -s.oo and polynomial.eval(lo) == 0) or (hi != s.oo and polynomial.eval(hi) == 0):
            raise ValueError('root-count endpoints must lie off the spectrum')
        total += multiplicity*int(polynomial.count_roots(lo, hi))
    return total


@lru_cache(maxsize=1)
def spectral_windows():
    minima = inherited.sector_minima(BITS)
    ground = inherited.global_ground_interval(BITS)
    blocks = inherited.rational_blocks()
    ground_counts = [root_count(matrix, None, GROUND_CUT) for _, _, matrix in blocks]
    if sum(ground_counts) != 1:
        raise ArithmeticError('the complete source has no unique ground below the declared cut')
    winner = next(j for j, count in enumerate(ground_counts) if count == 1)
    if blocks[winner][1] != 2:
        raise ArithmeticError('the unique source ground changes its charge sector')
    if sum(root_count(matrix, *ground) for _, _, matrix in blocks) != 1:
        raise ArithmeticError('the inherited ground interval does not isolate one state')
    result = {}
    for occupation in (1, 3):
        selected = [matrix for _, n, matrix in blocks if n == occupation]
        window = minima[occupation]
        in_window = sum(root_count(matrix, *window) for matrix in selected)
        below_cut = sum(root_count(matrix, None, CHARGED_CUT) for matrix in selected)
        if in_window != 2 or below_cut != in_window or window[1] >= CHARGED_CUT:
            raise ArithmeticError('the lowest charged window is not isolated from its complement')
        result[occupation] = {'energy_interval': window, 'projector_rank': in_window,
                              'complement_energy_lower': CHARGED_CUT}
    x = s.Symbol('x')
    ground_poly = s.Poly(x**5-54*x**4+958*x**3-5328*x**2-5816*x+23184, x)
    charged_poly = s.Poly(x**3-35*x**2+286*x+341, x)
    if ground_poly.count_roots(*map(_q, ground)) != 1:
        raise ArithmeticError('the source ground polynomial root differs')
    for row in result.values():
        if charged_poly.count_roots(*map(_q, row['energy_interval'])) != 1:
            raise ArithmeticError('the source charged polynomial root differs')
    # Divisibility attaches these displayed polynomials to the literal source.
    winner_poly = s.Matrix(blocks[winner][2]).charpoly(x).as_poly()
    if not winner_poly.rem(ground_poly).is_zero:
        raise ArithmeticError('the displayed ground polynomial is not a source factor')
    for occupation in (1, 3):
        carrying = [(matrix, root_count(matrix, *result[occupation]['energy_interval']))
                    for _, n, matrix in blocks if n == occupation]
        if any(not s.Matrix(matrix).charpoly(x).as_poly().rem(charged_poly).is_zero
               for matrix, count in carrying if count):
            raise ArithmeticError('a lowest charged source block has a different polynomial')
    return {'ground_interval': ground, 'ground_rank': 1,
            'ground_complement_energy_lower': GROUND_CUT,
            'ground_block_indices': blocks[winner][0],
            'ground_block_matrix': blocks[winner][2],
            'ground_polynomial_descending': tuple(int(z) for z in ground_poly.all_coeffs()),
            'charged_polynomial_descending': tuple(int(z) for z in charged_poly.all_coeffs()),
            'charged': result, 'source_dimension': 208,
            'root_counts_include_all_algebraic_multiplicities': True}


def _primitive_vector(vector):
    denominator = s.ilcm(*[value.q for value in vector])
    integers = [int(value*denominator) for value in vector]
    divisor = s.igcd(*integers)
    if divisor == 0:
        raise ArithmeticError('zero inverse-iteration source trial')
    return s.ImmutableMatrix([value//divisor for value in integers])


@lru_cache(maxsize=8)
def ground_trial(iterations=8):
    """Deterministic rational inverse iteration, admitted by its energy bound."""
    if type(iterations) is not int or not 1 <= iterations <= 32:
        raise ValueError('inverse iterations must be an integer between1 and32')
    data = spectral_windows()
    lower, upper = data['ground_interval']
    if not TRIAL_SHIFT < lower < upper < GROUND_CUT:
        raise ArithmeticError('the source trial shift or ground complement cut is invalid')
    h = s.Matrix(data['ground_block_matrix'])
    inverse = (h-_q(TRIAL_SHIFT)*s.eye(h.rows)).inv()
    vector = s.eye(h.rows)[:, 0]
    for _ in range(iterations):
        vector = _primitive_vector(inverse*vector)
    norm = int((vector.T*vector)[0])
    energy = F((vector.T*h*vector)[0])/norm
    leakage = max(F(0), min(F(1), (energy-lower)/(GROUND_CUT-upper)))
    distance = sqrt_rational(leakage, 32).hi
    return {'indices': data['ground_block_indices'], 'vector': vector,
            'norm_squared': norm, 'mean_energy': energy,
            'ground_leakage_upper': leakage, 'half_trace_distance_upper': distance,
            'inverse_iterations': iterations, 'inverse_shift': TRIAL_SHIFT,
            'occupation': 2, 'source_density_exactly_stationary': False}


def _source_chart_probe(k, output, input_):
    phases = inherited.source.real_computation_chart()[0]
    return s.ImmutableMatrix(len(output), len(input_), lambda i, j:
        s.conjugate(phases[output[i]])*k[output[i], input_[j]]*phases[input_[j]])


@lru_cache(maxsize=256)
def collision_block(matrix, energy_lower, width=FILTER_WIDTH):
    """Rational source-derived collision; its pointer dilation is unitary."""
    h = s.Matrix(matrix)
    eye = s.eye(h.rows)
    d = (h-_q(energy_lower)*eye)/_q(width)
    effect = s.ImmutableMatrix((eye+d*d).inv())
    a = s.ImmutableMatrix(effect-s.I*d*effect)
    b = s.ImmutableMatrix(eye-effect+s.I*d*effect)
    return {'D': s.ImmutableMatrix(d), 'success_effect': effect, 'A': a, 'B': b}


@lru_cache(maxsize=1)
def complete_window_detector():
    """Forty original probabilities, retaining each first-failure time."""
    data = spectral_windows(); trial = ground_trial()
    vector, norm = trial['vector'], trial['norm_squared']
    output = []
    for row in inherited.instrument():
        target = 2+row['sign']
        lower = data['charged'][target]['energy_interval'][0]
        total = s.S.Zero; passed = s.S.Zero
        failed = [s.S.Zero for _ in range(FILTER_POWER)]
        for indices, occupation, matrix in inherited.rational_blocks():
            probe = _source_chart_probe(row['K'], indices, trial['indices'])
            if not any(probe):
                continue
            value = probe*vector
            weight = sum(s.conjugate(x)*x for x in value)
            total += weight
            if occupation == target:
                effect = collision_block(matrix, lower)['success_effect']
                retained = s.eye(len(indices))
                for j in range(FILTER_POWER):
                    failure_effect = retained*(s.eye(len(indices))-effect)
                    failed[j] += (value.H*failure_effect*value)[0]
                    retained = retained*effect
                passed += (value.H*retained*value)[0]
            else:
                # This is a declared charge-and-energy apparatus. Other
                # output charges are retained in the first-failure record.
                failed[0] += weight
        total, passed = F(total)/norm, F(passed)/norm
        failed = [F(value)/norm for value in failed]
        if not 0 <= passed <= total or any(value < 0 for value in failed) or passed+sum(failed) != total:
            raise ArithmeticError('the source window filter is not contractive')
        outcomes = [('window', passed)]+[(f'failure_{j}', value) for j, value in enumerate(failed)]
        for label, probability in outcomes:
            output.append({'sign': row['sign'], 'mode': row['mode'], 'outcome': label,
                           'original_probability': probability,
                           'compensating_reference_charge': -6*row['sign'],
                           'compensating_reference_parity': 1})
    if sum(row['original_probability'] for row in output) != 1:
        raise ArithmeticError('the complete original source detector is not normalized')
    for sign in (1, -1):
        if sum(row['original_probability'] for row in output if row['sign'] == sign) != F(1, 2):
            raise ArithmeticError('the complete original charged branch differs from one half')
    return {'records': tuple(output), 'original_probability_sum': F(1),
            'all_40_detector_records_retained': True, 'failures_retained': True,
            'reference_retained': True, 'conditioned_on_success': False,
            'filter_width': FILTER_WIDTH, 'filter_power': FILTER_POWER}


@lru_cache(maxsize=1)
def certified_residues():
    data = spectral_windows(); trial = ground_trial(); detector = complete_window_detector()
    result = {}
    for sign in (1, -1):
        target = 2+sign
        lower, upper = data['charged'][target]['energy_interval']
        outside_success = (FILTER_WIDTH**2/((CHARGED_CUT-lower)**2+FILTER_WIDTH**2))**FILTER_POWER
        inside_success_lower = (FILTER_WIDTH**2/((upper-lower)**2+FILTER_WIDTH**2))**FILTER_POWER
        passed = sum(row['original_probability'] for row in detector['records']
                     if row['sign'] == sign and row['outcome'] == 'window')
        # Both states lie in N=2. Each sign's effect is bounded by I/2 there.
        state_error = trial['half_trace_distance_upper']/2
        weight_lower = max(F(0), passed-state_error-outside_success/2)
        weight_upper = min(F(1, 2), (passed+state_error)/inside_success_lower)
        if not 0 < weight_lower <= weight_upper <= F(1, 2):
            raise ArithmeticError('the original source residue is not certified nonzero')
        result[sign] = {'target_occupation': target,
                       'target_source_energy_interval': (lower, upper),
                       'energy_transfer_interval': (lower-data['ground_interval'][1],
                                                    upper-data['ground_interval'][0]),
                       'projector_rank': 2,
                       'original_ground_residue_interval': (weight_lower, weight_upper),
                       'trial_original_window_probability': passed,
                       'outside_window_success_upper': outside_success,
                       'inside_window_failure_upper': 1-inside_success_lower,
                       'state_transfer_probability_error_upper': state_error,
                       'original_charge_branch_weight': F(1, 2)}
    # The exact antiunitary maps the unique ground to itself up to phase.
    if not all(inherited.charge_conjugation()['identities'].values()):
        raise ArithmeticError('source charge conjugation failed')
    common = (max(result[sign]['original_ground_residue_interval'][0] for sign in (1, -1)),
              min(result[sign]['original_ground_residue_interval'][1] for sign in (1, -1)))
    if common[0] > common[1]:
        raise ArithmeticError('exact conjugate residue intervals are inconsistent')
    return {'channels': result, 'exact_charge_conjugate_ground_weights_equal': True,
            'common_original_residue_interval': common,
            'physical_particle_identification': False}


def finite_clock_window_plan(count=131072, duration=F(1, 16), steps=32768,
                             radius=F(1, 32), resolution=F(1, 64)):
    """Same-H finite Fourier energy detector with every clock label retained."""
    if type(count) is not int or count < 2:
        raise ValueError('a complete Fourier register has at least two labels')
    radius, resolution, duration = F(radius), F(resolution), F(duration)
    if not 0 < resolution < radius:
        raise ValueError('positive energy resolution must be below the reported radius')
    residues = certified_residues(); data = spectral_windows()
    window = data['charged'][3]['energy_interval']
    center = sum(window)/2
    if (window[1]-window[0])/2+resolution >= radius:
        raise ValueError('the certified source window is not inside the eroded readout interval')
    if center+radius+resolution >= CHARGED_CUT:
        raise ValueError('the enlarged readout interval reaches the spectral complement')
    leakage = kernels.real_energy_leakage(count, duration, resolution, F(22**2), support_upper=22)
    clock_error = min(F(1), (count-1)*inherited.source.cayley_bound(duration, steps))
    preparation_error = ground_trial()['half_trace_distance_upper']/2
    # Each sign has original weight1/2 on both N=2 input states.
    event_error = (leakage['total']+clock_error)/2+preparation_error
    low, high = residues['common_original_residue_interval']
    ground = data['ground_interval']
    return {'Fourier_bins': count, 'clock_step': duration,
            'Cayley_substeps_per_clock_step': steps,
            'maximum_source_time': (count-1)*duration,
            'final_energy_window_center': center,
            'reported_energy_radius': radius, 'energy_resolution': resolution,
            'source_ground_reference_interval': ground,
            'reported_transfer_center_interval': (center-ground[1], center-ground[0]),
            'alias_probability_upper': leakage['alias'],
            'filter_leakage_upper': leakage['total'],
            'complete_clock_error_upper': clock_error,
            'original_event_error_per_sign_upper': event_error,
            'original_readout_weight_per_sign_interval': (max(F(0), low-event_error),
                                                         min(F(1, 2), high+event_error)),
            'all_clock_bins_and_compensating_references_retained': True,
            'full_clock_histogram_executed': False}


def _outward_pair(pair):
    return (inherited.dyadic_interval(pair[0])[0], inherited.dyadic_interval(pair[1])[1])


def certificate():
    pins = authenticate_sources()
    windows = spectral_windows(); trial = ground_trial(); detector = complete_window_detector()
    residues = certified_residues(); plan = finite_clock_window_plan()
    rows = [{key: value for key, value in row.items() if key != 'original_probability'} |
            {'original_probability_interval': inherited.dyadic_interval(row['original_probability'])}
            for row in detector['records']]
    channels = {}
    for sign, row in residues['channels'].items():
        channels[sign] = {key: (inherited.dyadic_interval(value) if isinstance(value, F)
                               and key not in ('original_charge_branch_weight',) else value)
                         for key, value in row.items() if key != 'original_ground_residue_interval'}
        channels[sign]['original_ground_residue_interval'] = _outward_pair(row['original_ground_residue_interval'])
    # The exact large rational values remain executable; receipts enclose them
    # outward to avoid printing arbitrary-size inverse-iteration denominators.
    plan['original_event_error_per_sign_upper'] = inherited.dyadic_interval(plan['original_event_error_per_sign_upper'])[1]
    plan['original_readout_weight_per_sign_interval'] = _outward_pair(plan['original_readout_weight_per_sign_interval'])
    return inherited.serial({
        'schema': 'oph.native_charged_residue.v1', 'source_sha256': pins,
        'source_dimension': 208, 'all_fock_states': 16, 'all_flux_states': 13,
        'ground_rank': 1, 'ground_energy_interval': windows['ground_interval'],
        'ground_complement_lower': GROUND_CUT,
        'ground_polynomial_descending': windows['ground_polynomial_descending'],
        'charged_polynomial_descending': windows['charged_polynomial_descending'],
        'ground_trial': {'inverse_iterations': trial['inverse_iterations'], 'shift': TRIAL_SHIFT,
                         'mean_energy_interval': inherited.dyadic_interval(trial['mean_energy']),
                         'half_trace_distance_upper': inherited.dyadic_interval(trial['half_trace_distance_upper'])[1],
                         'source_density_exactly_stationary': False,
                         'prepared_by_postselection': False},
        'charged_channels': channels,
        'exact_charge_conjugate_ground_weights_equal': True,
        'common_original_residue_interval': _outward_pair(residues['common_original_residue_interval']),
        'complete_rational_window_detector': {'records': rows,
            'original_probability_sum': 1, 'original_probability_per_sign': '1/2',
            'filter_width': FILTER_WIDTH, 'filter_power': FILTER_POWER,
            'failures_and_references_retained': True, 'conditioned_on_success': False},
        'finite_clock_window_plan': plan,
        'same_source_H_and_clock': True,
        'physical_particle_identification': False,
        'relativistic_mass_shell': False, 'GeV_mass': None,
    })
