"""Existing positive clock bridges and an exact reference-substitution audit.

The pending clock continuation and the coherent wave use the same binary
carrier. Common operational time does not require equal clock and field
generators. Reuse the native-return wave reference and finite-record phase
inverse before testing the stronger, optional equality of full inner flows.
Compiled apparatus certificates remain those of their original producers.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import sympy as sp

from model.native_modular_reference import history_clock as clock
from model.native_modular_reference import history_transport as transport
from model.observer_time_fields import spatial_dynamics as fields
from model.observer_time_fields import wave_reference
from model.native_modular_reference import phase_inference, reference_tomography
from tools import strict_json
from .contract import binding, check_receipt as check_receipt_file

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('clock_wave_results.json')
PENDING_PR_HEAD = '3d282b3d59763fa1a5171ad9a4fbc09fc0abf35a'
CORES = ((0, 1), (0, 8), (1, 5), (5, 8))


def _fraction(value):
    """Own canonical exact input values; bools, floats and subclasses fail."""
    if type(value) is int:
        return F(value)
    if type(value) is not F:
        raise ValueError('an original exact integer or Fraction is required')
    n, d = value.numerator, value.denominator
    if (type(n) is not int or type(d) is not int or d <= 0
            or math.gcd(abs(n), d) != 1):
        raise ValueError('canonical exact rational required')
    return F(n, d)


def _sqrt_bounds(value, bits=80):
    value = _fraction(value)
    if value < 0 or type(bits) is not int or not 8 <= bits <= 256:
        raise ValueError('nonnegative rational and 8..256 precision bits required')
    scale = 1 << bits
    n = math.isqrt(value.numerator*scale*scale//value.denominator)
    lo, hi = F(n, scale), F(n+1, scale)
    if lo*lo == value:
        hi = lo
    if not lo*lo <= value <= hi*hi:
        raise ArithmeticError('rational square-root enclosure failed')
    return lo, hi


def _log_bounds(value, terms=48):
    value = _fraction(value)
    if not 1 <= value <= 5 or type(terms) is not int or not 1 <= terms <= 128:
        raise ValueError('log certificate requires a rational in [1,5] and 1..128 terms')
    z = (value-1)/(value+1)
    lo = 2*sum((z**(2*k+1)/(2*k+1) for k in range(terms)), F())
    tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return lo, lo+tail


def _wave_frequency_bounds(laplacian_bounds):
    """arccos(1-lambda/60)=2 asin(sqrt(lambda/120)), exactly enclosed."""
    if type(laplacian_bounds) is not tuple or len(laplacian_bounds) != 2:
        raise ValueError('two ordered exact Laplacian endpoints required')
    lo, hi = map(_fraction, laplacian_bounds)
    if not 0 < lo <= hi <= 12:
        raise ValueError('the elementary source band must lie in (0,12]')
    xlo, _ = _sqrt_bounds(lo/120)
    _, xhi = _sqrt_bounds(hi/120)
    terms = 24
    def partial(x):
        return sum((F(math.comb(2*k, k), 4**k*(2*k+1))*x**(2*k+1)
                    for k in range(terms)), F())
    # Every remaining positive Taylor coefficient is at most one.
    return 2*partial(xlo), 2*(partial(xhi)+xhi**(2*terms+1)/(1-xhi*xhi))


def _ratio(numerator, denominator):
    a, b = numerator
    c, d = denominator
    if not 0 <= a <= b or not 0 < c <= d:
        raise ValueError('positive ordered interval denominator required')
    return a/d, b/c


def _difference(a, b):
    return a[0]-b[1], a[1]-b[0]


def _outward(interval, places=12):
    scale = 10**places
    lo, hi = interval
    return [str(F(lo.numerator*scale//lo.denominator, scale)),
            str(F(-((-hi.numerator*scale)//hi.denominator), scale))]


def spectral_shape():
    """Affine-invariant level-gap shapes, even after reversing clock orientation.

    This is an additional comparison of one vector from each of the three
    density bands with one vector from each reference eigenspace. It does
    not assert that the actual density subspace reduces the reference clock.
    """
    root = _sqrt_bounds(F(5))
    bands = ((5-root[1], 5-root[0]), (F(6), F(6)),
             (5+root[0], 5+root[1]))
    omega = tuple(_wave_frequency_bounds(b) for b in bands)
    shape = _ratio(_difference(omega[2], omega[0]), _difference(omega[1], omega[0]))
    # Reference weights, ordered by increasing -log(weight): 2/45,1/36,2/225.
    # Total gap log(5); lower adjacent gap log(8/5), upper log(25/8).
    positive = _ratio(_log_bounds(F(5)), _log_bounds(F(8, 5)))
    reversed_shape = _ratio(_log_bounds(F(5)), _log_bounds(F(25, 8)))
    margins = [min(abs(shape[0]-r[1]), abs(r[0]-shape[1])) for r in (positive, reversed_shape)]
    if any(not (shape[1] < r[0] or r[1] < shape[0]) for r in (positive, reversed_shape)):
        raise ArithmeticError('clock and wave gap-shape intervals overlap')
    return {'wave_laplacian_bands': ['5-sqrt(5)', '6', '5+sqrt(5)'],
            'wave_band_multiplicities': [3, 5, 3],
            'wave_phase_per_coherent_step': 'arccos(1-lambda/60)',
            'wave_frequency_intervals': [_outward(x) for x in omega],
            'ratio_definition': '(highest_level-lowest_level)/(middle_level-lowest_level)',
            'wave_gap_shape_interval': _outward(shape),
            'modular_positive_orientation_shape_interval': _outward(positive),
            'modular_reversed_orientation_shape_interval': _outward(reversed_shape),
            'positive_separation_lower': _outward((min(margins), min(margins)))[0],
            'either_affine_orientation_excluded': True,
            'common_units_or_positive_incident_clock_factors_change_shapes': False}


def _source():
    data = fields.native_density()
    counts = sp.ImmutableMatrix(data['counts'])
    z = sp.ImmutableMatrix(data['field'])
    laplacian = sp.ImmutableMatrix(data['laplacian'])
    if counts.shape != (66, 66) or counts != counts.T:
        raise ArithmeticError('native binary field kernel changed')
    if z*counts != (60*sp.eye(12)-laplacian)*z:
        raise ArithmeticError('native density intertwiner changed')
    # The inherited producer proves these actual source bands; recheck the
    # characteristic polynomial here before using the analytic gap intervals.
    x = sp.Symbol('lambda')
    if laplacian.charpoly(x).as_expr() != sp.expand(x*(x-6)**5*(x*x-10*x+20)**3):
        raise ArithmeticError('source density spectrum changed')
    return counts, z, laplacian


def commutator_certificate(weights, counts):
    """A necessary functional-calculus test, deliberately not a sufficiency test."""
    weights = transport._reference(weights)
    if not isinstance(counts, sp.ImmutableMatrix) or counts.shape != (66, 66):
        raise ValueError('the immutable exact 66-state count matrix is required')
    if counts != counts.T or any(not x.is_Integer or x < 0 for x in counts):
        raise ValueError('symmetric nonnegative integer counts required')
    if any(sum(counts[:, j]) != 60 for j in range(66)):
        raise ValueError('complete sixty-branch columns required')
    sigma = tuple(sp.Rational(x.numerator, x.denominator) for x in weights)
    entries = [(i, j, (sigma[i]-sigma[j])*counts[i, j]/60)
               for i in range(66) for j in range(66)
               if (sigma[i]-sigma[j])*counts[i, j] != 0]
    return {'nonzero_entries': len(entries),
            'squared_frobenius_norm': str(sum((value**2 for _, _, value in entries), sp.Rational(0))),
            'first_matrix_unit_witness': list(entries[0][:2])+[str(entries[0][2])] if entries else None,
            'commutation_necessary_condition_passes': not entries}


def transported_comparison():
    counts, z, _ = _source()
    sigma = transport.reference()
    states = tuple(transport.native.states(binary=True))
    source_order = tuple(tuple(int(x) for x in state) for state in states)
    if tuple(fields.coherent.STATES) != source_order:
        raise ArithmeticError('clock and wave use different native basis orders')
    projector = z.T*z/10+sp.ones(66)/66
    if projector*projector != projector or sp.trace(projector) != 12:
        raise ArithmeticError('constant plus native density projection failed')
    seen, law_rows = {}, []
    for law in ('C', 'E'):
        branches = clock.branches(law)
        cores = set()
        for branch in branches:
            permutation = transport.port_permutation(law, branch)
            core = tuple(sorted(permutation[p] for p in (0, 1)))
            current = transport.extended_reference(law, branch, sigma)
            expected = tuple((F(2, 225), F(1, 36), F(2, 45))[sum(x[p] for p in core)] for x in states)
            if current != expected or (core in seen and current != seen[core]):
                raise ArithmeticError('transported reference disagrees with its native current tape')
            seen[core] = current
            cores.add(core)
        if tuple(sorted(cores)) != CORES:
            raise ArithmeticError('the pending clock no longer uses the declared four cores')
        law_rows.append({'law': law, 'complete_two_attempt_branches': len(branches),
                         'transported_cores': [list(x) for x in sorted(cores)]})
    rows = []
    for core in CORES:
        current = seen[core]
        witness = commutator_certificate(current, counts)
        if witness['commutation_necessary_condition_passes']:
            raise ArithmeticError('the native reference unexpectedly commutes with the wave kernel')
        diagonal = sp.diag(*(sp.Rational(x.numerator, x.denominator) for x in current))
        leakage = (sp.eye(66)-projector)*diagonal*z.T
        leak = sum((x*x for x in leakage), sp.Rational(0))
        if leak <= 0:
            raise ArithmeticError('constant-plus-density subspace unexpectedly reduces the reference')
        rows.append({'core': list(core), **witness,
                     'constant_plus_density_reference_leakage_squared_frobenius': str(leak),
                     'reference_weights_and_multiplicities': {str(k): v for k, v in sorted(Counter(current).items())}})
    # Four independent Krylov vectors imply >=4 distinct eigenvalues of the
    # self-adjoint K. arccos is injective on its real spectrum. No numerical rank.
    vector, vectors = sp.eye(66)[:, 0], []
    for _ in range(4):
        vectors.append(vector)
        vector = counts*vector
    columns = sp.Matrix.hstack(*vectors)
    gram = columns.T*columns
    determinant = gram.det()
    if determinant <= 0:
        raise ArithmeticError('four-level exact wave spectral witness failed')
    return {'native_dimension': 66, 'native_event_count_per_column': 60,
            'kernel_definition': 'K=counts/60; H_wave=arccos(K)',
            'reference_generator_definition': 'H_sigma=-log(sigma); global shifts cancel from conjugation',
            'laws': law_rows, 'core_certificates': rows,
            'full_shell_spectrum': {'krylov_seed_native_index': 0,
                'krylov_powers_of_counts': [0, 1, 2, 3],
                'exact_gram_matrix': [[str(x) for x in row] for row in gram.tolist()],
                'exact_gram_determinant': str(determinant),
                'wave_distinct_generator_levels_lower': 4,
                'transported_reference_distinct_generator_levels': 3,
                'full_shell_affine_identification_even_after_unitary_reencoding': False},
            'tracial_control': commutator_certificate((F(1, 66),)*66, counts)}


def compatible_two_mode_control():
    """Actual source modes, explicitly different reference and encoding contract."""
    counts, z, laplacian = _source()
    port_mode = sp.Matrix([1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    mode = z.T*port_mode
    constant = sp.ones(66, 1)
    if (laplacian*port_mode != 6*port_mode or counts*mode != 54*mode
            or counts*constant != 60*constant or (constant.T*mode)[0] != 0):
        raise ArithmeticError('the actual source two-mode control failed')
    return {'source_mode_port_coefficients': list(map(int, port_mode)),
            'constant_mode_squared_norm': str((constant.T*constant)[0]),
            'density_mode_squared_norm': str((mode.T*mode)[0]),
            'kernel_on_normalized_two_mode_code': ['1', '9/10'],
            'wave_generator_on_code': ['0', 'theta=arccos(9/10)'],
            'supplied_reference_on_code': ['2/3', '1/3'],
            'gibbs_inverse_parameter': 'b=log(2)/theta',
            'generator_identity': '-log(sigma_control)=b*H_wave+log(3/2)*I',
            'predual_clock_conversion': 's=-t/b=-t*theta/log(2) for U_sigma(s)=sigma^(-is) and U_wave(t)=exp(-it H_wave)',
            'one_wave_step_balanced_probe_plus_probability_exact': '19/20',
            'same_as_transported_native_reference': False,
            'new_native_preparation_or_compiled_service_executed': False}


def native_return_bridge():
    """Reuse the existing positive same-K reference, without preparing it."""
    counts, z, _ = _source()
    kernel = counts/60
    projector = z.T*z/10+sp.ones(66)/66
    rows = []
    for cap in (2, 4):
        reference = wave_reference.reference(cap)
        sigma, transfer = reference['density'], reference['transfer']
        if (sigma*kernel != kernel*sigma or sp.trace(sigma) != 1
                or (sp.eye(66)-projector)*sigma*z.T != sp.zeros(66, 12)):
            raise ArithmeticError('existing native-return reference lost its wave dictionary')
        probabilities = tuple(reference[key] for key in
                              ('return_probability', 'different_probability', 'cap_probability'))
        if sum(probabilities, F()) != 1:
            raise ArithmeticError('native-return reference discarded a terminal outcome')
        difference = sigma-reference['native_return_state']
        if not any(difference):
            raise ArithmeticError('reconstruction was silently called a prepared source state')
        # Positivity follows analytically from positive polynomial coefficients
        # and 0<=K<=I; this exact scalar floor makes faithfulness quantitative.
        floor = 1/sp.trace(sum((sp.Rational(b)*kernel**n for n, b in
                    enumerate(wave_reference.coefficients(cap))), sp.zeros(66)))
        rows.append({'cap': cap, 'return_probability': str(probabilities[0]),
                     'different_probability': str(probabilities[1]),
                     'cap_probability': str(probabilities[2]),
                     'reference_minimum_eigenvalue_lower': str(floor),
                     'reference_kernel_commutator_exact_zero': True,
                     'constant_plus_density_subspace_preserved': True,
                     'reconstructed_minus_actual_return_density_squared_frobenius': str(sum(x*x for x in difference)),
                     'reconstructed_density_sha256': fields._digest(sigma),
                     'finite_reference_transfer_sha256': fields._digest(transfer)})
    bounds = wave_reference.cap_bounds(512)
    return {'reused_theorem': 'observer_time_fields/README.md:459-540; wave_reference.py',
            'same_native_configuration_kernel': True,
            'ideal_reference': 'sigma_wave=exp(-Theta)/Tr(exp(-Theta)), Theta=arccos(K)',
            'positive_length_law': 'exp(arcsin(z)-pi/2)',
            'finite_cap_reference': 'sigma_h=P_h(K)/Tr(P_h(K))',
            'finite_cap_generator': 'Theta_h=-log(P_h(K)/P_h(1))',
            'generator_bound': '0<=Theta-Theta_h<=-log(1-delta_h) I',
            'same_positive_wave_phase_magnitudes_and_cosine_amplitudes': True,
            'exact_finite_cap_checks': rows,
            'cap_512_reference_halftrace_upper': str(bounds['reconstructed_reference_trace_distance_upper']),
            'cap_512_generator_norm_error_upper': str(bounds['wave_modular_generator_norm_error_upper']),
            'ideal_uncapped_mean_length': 'infinite',
            'coherent_reference_prepared_or_vertex_modular_unitary_executed': False,
            'new_theorem_claimed': False}


def operational_parameter_bound(interval, conversion, frequency_bound):
    """Propagate an identified clock interval to a separately specified field.

    Conditional on t=conversion*s and ||H_field||<=frequency_bound, Duhamel
    gives the unitary/half-diamond and cosine bounds. No equality with the
    clock's Hamiltonian is assumed. Hardware, coupling and field errors add
    their own allowances; the clock interval retains its own failure event.
    """
    if type(interval) is not tuple or len(interval) != 2:
        raise ValueError('an ordered exact parameter interval is required')
    lo, hi = map(_fraction, interval)
    rate, omega = map(_fraction, (conversion, frequency_bound))
    if lo > hi or rate == 0 or omega < 0:
        raise ValueError('ordered endpoints, nonzero conversion and nonnegative band bound required')
    middle, radius = (lo+hi)/2, (hi-lo)/2
    error = min(F(2), abs(rate)*omega*radius)
    return {'parameter_interval': [str(lo), str(hi)],
            'midpoint': str(middle), 'radius': str(radius),
            'declared_wave_parameter_per_clock_parameter': str(rate),
            'field_generator_norm_upper': str(omega),
            'wave_parameter_interval': list(map(str, sorted((rate*lo, rate*hi)))),
            'unitary_and_cosine_operator_error_upper': str(error),
            'quantum_channel_halfdiamond_and_detector_probability_error_upper': str(min(F(1), error)),
            'clock_and_field_hamiltonians_required_equal': False}


def operational_readout_bridge():
    """Replay the existing blind inverse, then propagate its complete interval.

    This reuses retained counts and the producer's admitted systematic budget;
    it neither resamples nor reruns compiled apparatus. The programmed target
    is not consumed. The field conversion is an explicit conditional choice.
    """
    path = ROOT/'model/native_modular_reference/phase_clock_results.json'
    retained = strict_json.loads(path.read_text())
    for name, digest in retained['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('retained phase-clock source binding changed')
    if retained['schema'] != 'oph.native_modular_reference.phase_clock.v1':
        raise ValueError('the existing complete-record phase-clock receipt is required')
    rows = []
    # An independent certified density-band bound, not a property of the clock.
    root = _sqrt_bounds(F(5))
    if _wave_frequency_bounds((5+root[0], 5+root[1]))[1] >= F(1, 2):
        raise ArithmeticError('actual elementary wave band exceeds the declared bound')
    for alternative in retained['alternatives']:
        intervals, failures = [], []
        for row in alternative['readouts']:
            counts = {quadrature: {tuple(record): count for record, count in
                                  zip(sample['records'], sample['counts'], strict=True)}
                      for quadrature, sample in row['sampling'].items()}
            previous = row['inference']
            inferred = phase_inference.infer_counts(
                counts, alternative['training']['execution']['counts'],
                law=alternative['law'], index=row['index'],
                systematic_halftrace_upper=F(previous['budget']['systematic_halftrace_upper']),
                reference_sampling_bias=reference_tomography.SAMPLER_GATE)
            if (inferred['status'] != 'identified'
                    or list(map(str, inferred['parameter_interval'])) != previous['parameter_interval']):
                raise ArithmeticError('clock parameter interval does not follow its complete records')
            intervals.append(inferred['parameter_interval'])
            failures.append(inferred['confidence_failure_upper'])
        interval = max(x[0] for x in intervals), min(x[1] for x in intervals)
        if list(map(str, interval)) != alternative['common_cap_parameter_intersection']:
            raise ArithmeticError('overlapping-cap clock intervals disagree')
        # One training event is shared by the two caps, not charged twice.
        training_failure = inferred['reference_inference']['confidence_failure_upper']
        joint_failure = sum(failures, F())-training_failure
        if str(joint_failure) != alternative['both_caps_joint_confidence_failure_upper']:
            raise ArithmeticError('shared-training confidence accounting changed')
        rows.append({'law': alternative['law'],
                     'joint_clock_confidence_failure_upper': str(joint_failure),
                     'conditional_density_band_prediction': operational_parameter_bound(interval, F(1), F(1, 2))})
    return {'reused_clock': 'native_modular_reference.phase_clock: learned enacted reference gap, complete native quadrature scores',
            'replayed_from_retained_counts_without_programmed_parameter': True,
            'clock_spectrum_matched_to_field_spectrum': False,
            'declared_conversion_for_this_conditional_control': 't=s; not inferred from the two separate experiments',
            'field_band': 'constant plus the three elementary native density bands of arccos(K), norm<1/2',
            'readouts': rows,
            'new_joint_clock_field_instrument_executed': False,
            'original_source_or_sampler_execution_reaudited': False,
            'confidence_is_conditional_on_existing_source_and_sampling_contract': True,
            'calibration_or_coupling_is_selected_by_SI_c': False}


def source_bindings():
    paths = set(fields.SOURCE_PATHS) | set(transport.source_bindings()) | set(wave_reference.SOURCE_PATHS)
    phase = strict_json.loads((ROOT/'model/native_modular_reference/phase_clock_results.json').read_text())
    paths.update(phase['source_sha256'])
    paths.update(('model/native_modular_reference/history_clock.py',
                  'model/native_modular_reference/history_clock_results.json.gz',
                  'model/native_modular_reference/phase_clock_results.json',
                  'model/observer_time_fields/wave_reference_results.json',
                  'model/clock_defect_energy/contract.py',
                  'model/clock_defect_energy/clock_wave.py',
                  'model/clock_defect_energy/test_clock_wave.py',
                  'model/clock_defect_energy/CLOCK_WAVE.md'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def report():
    result = {'schema': 'oph.clock_defect_energy.clock_wave.v2',
              'contract': binding(),
              'pending_input': {'pull_request': 1149, 'head': PENDING_PR_HEAD,
                                'merge_status_asserted': False, 'full_pr_audit_asserted': False},
              'source_bindings': source_bindings(),
              'existing_positive_native_return_bridge': native_return_bridge(),
              'operational_clock_readout_to_field_bound': operational_readout_bridge(),
              'transported_clock_vs_wave': transported_comparison(),
              'three_band_affine_shape': spectral_shape(),
              'changed_reference_positive_control': compatible_two_mode_control(),
              'scope': {'new_result': 'reuse existing positive wave-reference and operational-clock bridges; guard against substituting the pending diagonal reference for the different wave reference',
                        'equal_full_inner_flows_required_for_common_operational_time': False,
                        'negative_result_excludes_a_common_operational_clock': False,
                        'earlier_common_clock_test_framing_withdrawn': True,
                        'physical_c_or_lab_duration_selected': False,
                        'new_empirical_comparisons': 0,
                        'finite_coherent_wave_is_passive_density_dynamics': False,
                        'compiled_clock_word_receipt_reexecuted': False,
                        'all_restricted_encodings_or_future_clock_field_joins_excluded': False,
                        'completion_gate_promoted': False}}
    return result


def check_receipt(retained, fresh):
    # JSON-type-sensitive comparison: Python equality alone equates 0 and False.
    if json.dumps(retained, sort_keys=True) != json.dumps(fresh, sort_keys=True):
        raise ValueError('retained clock/wave certificate differs from its exact producer')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    group.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = report()
    payload = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.write:
        args.write.write_text(payload)
    elif args.check:
        check_receipt_file(args.check, result)
        print('Existing positive clock bridges and scoped reference-substitution guard verified.')
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
