"""Fixed corrected-repair response and the first surviving correlation probe.

The existing record-selected corrector protects code6 exactly in its ideal
target contract. Here its actual C/E maps and fixed completion selector are
used to calculate the full logical transfer, rather than inserting another
logical Hamiltonian. That transfer is identity: its static response seed is
zero. Native one-port coordinates and every native-edge diagonal agreement
cost are scalar on this code.

A two-hop diagonal correlation cost, fixed by the same accepted seam
weights, has a nonscalar compression. Its full syndrome leakage is retained:
neither that compression nor the availability of source controls selects a
physical Yukawa or a closed coherent evolution. No measured parameters enter.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as s

from tools import strict_json
from . import record_code as code, slow_geometry
from ..native_modular_reference import local_preparation as preparation
from ..observer import native
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..repair_scheduler import coins


ROOT = Path(__file__).resolve().parents[2]
SOURCE_GATE = 1e-8
SOURCE_REPLAY_TOLERANCE = 2e-9
PERTURBATION = F(1, 64)


def _serialize(matrix):
    return [[str(value) for value in row] for row in s.Matrix(matrix).tolist()]


def completion_weights():
    """One existing finite proposal; identity timeout is never restarted."""
    packet = preparation.selectors(1)
    value = dict(packet['completion_accepted_seam_weights'])
    if any(type(w) is not F or w <= 0 for w in value.values()):
        raise ArithmeticError('exact positive existing completion weights required')
    if sum(value.values(), F())+packet['completion_timeout_probability'] != 1:
        raise ArithmeticError('complete selector lost its terminal probability')
    return value, packet['completion_timeout_probability']


def code_geometry():
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    support = tuple(tuple(index[state] for state in row) for row in code.SUPPORTS)
    if len(set(j for row in support for j in row)) != 36 or any(len(row) != 6 for row in support):
        raise ArithmeticError('the fixed six matching code changed')
    v = np.zeros((66, 6), complex)
    for logical, row in enumerate(support):
        v[list(row), logical] = 1/math.sqrt(6)
    frame = s.ImmutableMatrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5)
                               for a, b in row] for row in slow_geometry.raw_frame()])
    if frame.T*s.ones(12, 1) != s.zeros(3, 1):
        raise ArithmeticError('the existing native slow frame is not centered')
    port_compressions = []
    for port in range(12):
        diagonal = tuple(F(sum(states[i][port] for i in row), 6) for row in support)
        if diagonal != (F(1, 6),)*6:
            raise ArithmeticError('native one-port compression ceased to be scalar')
        port_compressions.append(diagonal)
    means = s.ones(12, 1)/6
    gram = (frame.T*frame)[0, 0]
    feature = s.simplify(frame.T*s.diag(*means)*frame/gram)
    if feature != s.eye(3)/6:
        raise ArithmeticError('corrected-code native mean feature is anisotropic')
    core = set(preparation.CORE)
    histograms = [Counter(sum(states[i][port] for port in core) for i in row) for row in support]
    if any(histogram != Counter({0: 4, 1: 2}) for histogram in histograms):
        raise ArithmeticError('core-biased modular-reference compression changed')
    sigma = preparation.reference()
    reference_levels = tuple({sigma[i] for i, state in enumerate(states)
                              if sum(state[port] for port in core) == q} for q in (0, 1, 2))
    if reference_levels != ({F(2, 225)}, {F(1, 36)}, {F(2, 45)}):
        raise ArithmeticError('the existing faithful reference is not the declared core-count function')
    original_means = s.Matrix([sum(s.Rational(sigma[i])*state[port] for i, state in enumerate(states))
                               for port in range(12)])
    original_feature = s.simplify(frame.T*s.diag(*original_means)*frame/gram)
    return {'states': states, 'index': index, 'supports': support, 'V': v,
            'frame': frame, 'feature': feature, 'original_feature': original_feature,
            'core_histograms': histograms, 'sigma': tuple(sigma)}


@lru_cache(maxsize=64)
def _reflection(mask):
    """Exact matching-block recovery; an independent six-row calculation."""
    selected = tuple((mask >> i) & 1 for i in range(6))
    count = sum(selected)
    if not count:
        return s.eye(6)
    target = s.ones(6, 1)/s.sqrt(6)
    incoming = s.Matrix(selected)/s.sqrt(count)
    difference = incoming-target
    norm = s.simplify((difference.T*difference)[0])
    result = s.eye(6) if norm == 0 else s.eye(6)-2*difference*difference.T/norm
    if (result*incoming-target).applyfunc(s.simplify) != s.zeros(6, 1):
        raise ArithmeticError('exact matching-block recovery does not return its code direction')
    return s.ImmutableMatrix(result)


def _raw_maps(law, seam, geometry):
    """Actual native maps cross-checked against floor arithmetic on all66 rows."""
    records = 1 << circuit.RECORD_BITS[law]
    maps = {(coin, raw): {} for coin in (0, 1) for raw in range(records)}
    for coin, label, mapping in coins.coin_maps(law, seam):
        reverse = {value: raw for raw, value in circuit.record_values(law, coin).items()}
        if label not in reverse:
            raise ArithmeticError('source emitted an unrecognized native record')
        maps[coin, reverse[label]] = dict(mapping)
    oracle = {(coin, raw): {} for coin in (0, 1) for raw in range(records)}
    a, b = seam
    for coin in (0, 1):
        for before, state in enumerate(geometry['states']):
            total = state[a]+state[b]
            final = list(state)
            final[a], final[b] = (total+coin)//2, (total+1-coin)//2
            current = state[a]-final[a]
            raw = 0 if current == 0 else 1 if current == 1 else 2
            if law == 'E':
                raw += 4*(total % 2)
            oracle[coin, raw][before] = geometry['index'][tuple(final)]
    if maps != oracle:
        raise ArithmeticError('native recorded maps differ from independent floor arithmetic')
    return maps


def diagonal_cost(power=1, *, uniform=False):
    """H_k=sum_(a<b)(W**k)_ab(n_a-n_b)^2, with W the fixed seam matrix.

    H2 is an additional correlation observable, NOT an already present term
    in corrected repair or a selected coherent interaction. Its normalization
    is exactly the existing one-proposal weight convention, with no mass fit.
    """
    if type(power) is not int or power not in (1, 2) or type(uniform) is not bool:
        raise ValueError('fixed one-hop/two-hop integer and Boolean control required')
    weights = ({edge: F(1) for edge in code.SEAMS} if uniform else completion_weights()[0])
    adjacency = s.zeros(12)
    for (a, b), weight in weights.items():
        adjacency[a, b] = adjacency[b, a] = s.Rational(weight)
    coefficients = adjacency**power
    geometry = code_geometry()
    values = tuple(sum(coefficients[a, b]*(state[a]-state[b])**2
                       for a in range(12) for b in range(a+1, 12))
                   for state in geometry['states'])
    means = tuple(sum(values[i] for i in row)/6 for row in geometry['supports'])
    variance = tuple(s.factor(sum((values[i]-mean)**2 for i in row)/6)
                     for row, mean in zip(geometry['supports'], means, strict=True))
    gaps = tuple(s.factor(means[2*j+1]-means[2*j]) for j in range(3))
    return {'W': s.ImmutableMatrix(adjacency), 'coefficients': s.ImmutableMatrix(coefficients),
            'values': values, 'compression': s.diag(*means), 'variance': variance,
            'occupation_differences': gaps, 'code_preserving': all(value == 0 for value in variance)}


def corrected_cycle():
    """Every native branch/recovery, complete selector, and36 logical operators."""
    geometry = code_geometry()
    weights, timeout = completion_weights()
    perturb = diagonal_cost(2)
    identity = np.eye(6)
    initial_reference_code_population = sum((sum(geometry['sigma'][i] for i in row)/6
                                            for row in geometry['supports']), F())
    reports = {}
    for law in ('C', 'E'):
        logical_transfer = float(timeout)*np.eye(36, dtype=complex)
        record_count = zero_count = checked = 0
        max_recovery_error = max_gram_error = 0.
        first_order = [s.Integer(0)]*6
        post_leakage = [s.Rational(timeout)*value for value in perturb['variance']]
        reference_code_population = timeout*initial_reference_code_population
        common = None
        digest = hashlib.sha256()
        for seam in code.SEAMS:
            maps = _raw_maps(law, seam, geometry)
            seam_rows = []
            for (coin, raw), mapping in maps.items():
                record_count += 1
                selected = [[j for j, state in enumerate(row) if state in mapping]
                            for row in geometry['supports']]
                counts = tuple(map(len, selected))
                if len(set(counts)) != 1:
                    raise ArithmeticError('corrected native record leaks logical information')
                count = counts[0]
                probability = F(count, 12)
                if not probability and mapping:
                    raise ArithmeticError('a zero code branch unexpectedly has an original binary-reference input')
                label = circuit.record_values(law, coin).get(raw)
                current = (label if law == 'C' else
                           (label-(2*coin-1)*(abs(label) % 2))//2) if label is not None else 0
                for before, after in mapping.items():
                    restored = list(geometry['states'][after])
                    if current:
                        a, b = seam
                        restored[a], restored[b] = restored[b], restored[a]
                    if tuple(restored) != geometry['states'][before]:
                        raise ArithmeticError('record-selected inverse native movement is not the declared endpoint swap')
                if probability != code.branch_weight(law, seam, coin, raw):
                    raise ArithmeticError('source code branch probability differs')
                zero_count += probability == 0
                if probability:
                    for positions in selected:
                        _reflection(sum(1 << j for j in positions))
                    actual = np.zeros((66, 6), complex)
                    for before, after in mapping.items():
                        actual[after] += geometry['V'][before]/math.sqrt(2)
                    recovery = code.declared_recovery_matrix66(law, seam, coin, raw)
                    returned = recovery@actual
                    target = math.sqrt(float(probability))*geometry['V']
                    error = float(np.linalg.norm(returned-target, 2))
                    gram_error = float(np.linalg.norm(returned.conj().T@returned-float(probability)*identity, 2))
                    max_recovery_error = max(max_recovery_error, error)
                    max_gram_error = max(max_gram_error, gram_error)
                    if max(error, gram_error) > 5e-13:
                        raise ArithmeticError('actual same-carrier recovery differs from exact identity branch')
                    checked += 6*66
                    actual_derivative = np.zeros((66, 6), complex)
                    for before, after in mapping.items():
                        actual_derivative[after] += float(perturb['values'][before])*geometry['V'][before]/math.sqrt(2)
                    logical_derivative = geometry['V'].conj().T@recovery@actual_derivative
                    derivative_diagonal = []
                    for logical, positions in enumerate(selected):
                        energies = [perturb['values'][geometry['supports'][logical][j]] for j in positions]
                        mean = sum(energies)/count
                        derivative_diagonal.append(float(mean)*math.sqrt(float(probability)))
                        if seam in weights:
                            post_leakage[logical] += s.Rational(weights[seam])*s.Rational(probability)*(sum(e*e for e in energies)/count-mean*mean)
                    if np.linalg.norm(logical_derivative-np.diag(derivative_diagonal), 2) > 5e-13:
                        raise ArithmeticError('same-carrier first derivative differs from exact selected energy moments')
                    local = geometry['V'].conj().T@returned
                    if seam in weights:
                        logical_transfer += float(weights[seam])*np.kron(local.conj(), local)
                        sigma_selected = sum((geometry['sigma'][geometry['supports'][logical][j]]
                                              for logical, positions in enumerate(selected) for j in positions), F())
                        reference_code_population += weights[seam]*sigma_selected/(2*count)
                        projection = geometry['V'].conj().T@recovery
                        actual_reference_population = sum(float(geometry['sigma'][before])*float(np.vdot(projection[:, after], projection[:, after]).real)/2
                                                          for before, after in mapping.items())
                        if abs(actual_reference_population-float(sigma_selected/(2*count))) > 5e-13:
                            raise ArithmeticError('full-reference projector witness differs from exact native subset calculation')
                # Logical first derivative is calculated from all native
                # records; no output is selected or conditioned on success.
                if probability and seam in weights:
                    for logical, positions in enumerate(selected):
                        selected_mean = sum(perturb['values'][geometry['supports'][logical][j]] for j in positions)/count
                        first_order[logical] += s.Rational(weights[seam])*s.Rational(probability)*selected_mean
                seam_rows.append({'coin': coin, 'raw': raw, 'weight': str(probability),
                                  'abort': raw not in circuit.record_values(law, coin)})
                digest.update(json.dumps([law, seam, coin, raw, str(probability), selected], separators=(',', ':')).encode())
            if sum((F(row['weight']) for row in seam_rows), F()) != 1:
                raise ArithmeticError('complete native record menu lost probability')
            if common is None:
                common = seam_rows
            elif common != seam_rows:
                raise ArithmeticError('corrected native records depend on selected seam')
        means = tuple(perturb['compression'][j, j] for j in range(6))
        first_order = tuple(s.factor(value+s.Rational(timeout)*mean)
                            for value, mean in zip(first_order, means, strict=True))
        if first_order != means or np.linalg.norm(logical_transfer-np.eye(36), 2) > 2e-12:
            raise ArithmeticError('full corrected response or first-order all-record identity failed')
        reference_defect = reference_code_population-initial_reference_code_population
        if reference_defect <= 0:
            raise ArithmeticError('active same-carrier faithful-reference obstruction disappeared')
        reports[law] = {'native_raw_records_checked': record_count, 'raw_zero_records_retained': zero_count,
            'same_carrier_output_entries_checked': checked, 'all_six_unknown_logical_columns': True,
            'maximum_full66_recovery_operator_error': max_recovery_error,
            'maximum_branch_gram_error': max_gram_error,
            'full_logical_transfer_dimension': 36,
            'logical_transfer_identity_error': float(np.linalg.norm(logical_transfer-np.eye(36), 2)),
            'exact_branch_and_support_sha256': digest.hexdigest(), 'common_raw_record_factor': common,
            'all_record_first_order_diagonal': list(map(str, first_order)),
            'perturbed_complete_record_syndrome_leakage_quadratic_coefficients': list(map(str, map(s.factor, post_leakage))),
            'original_faithful_reference_code_projector_probability_before': str(initial_reference_code_population),
            'original_faithful_reference_code_projector_probability_after': str(reference_code_population),
            'original_reference_half_trace_distance_lower_bound': str(reference_defect),
            'reference_projector_witness_independent_of_unitary_recovery_extension': True,
            'reference_witness_identity': 'R_h*K_h*V=sqrt(p_h)*V implies R_h.H*V=K_h*V/sqrt(p_h); zero-weight raw maps are empty on binary66',
            'original_same_carrier_faithful_reference_stationary': False}
    return reports


def graph_invariant_control():
    """Classify a*sum_i n_i+sum_(i<j)b_orbit(i,j)*n_i*n_j on totalN=2.

    These a,b coefficients name a normal-ordered invariant observable. They
    are not the edge coefficients of the separate (n_i-n_j)^2 cost probe.
    """
    transformations = slow_geometry.source().incidence_automorphisms()
    pairs = {(a, b) for a in range(12) for b in range(a+1, 12)}
    orbits = []
    while pairs:
        a, b = min(pairs)
        orbit = {tuple(sorted((permutation[a], permutation[b]))) for permutation in transformations}
        if not orbit <= pairs:
            raise ArithmeticError('pair action is not a complete orbit partition')
        orbits.append(orbit); pairs -= orbit
    if sorted(map(len, orbits)) != [6, 30, 30]:
        raise ArithmeticError('native graph pair orbits changed')
    antipodal = next(orbit for orbit in orbits if len(orbit) == 6)
    if antipodal != set(code.MATCHINGS[1]):
        raise ArithmeticError('antipodal matching no longer equals logical1')
    control = diagonal_cost(2, uniform=True)
    if control['compression'] != 36*s.eye(6)+4*s.diag(0, 1, 0, 0, 0, 0):
        raise ArithmeticError('symmetric graph correlation control changed')
    return {'actual_incidence_automorphisms': len(transformations), 'unordered_pair_orbit_sizes': sorted(map(len, orbits)),
            'antipodal_pairs': [list(pair) for pair in sorted(antipodal)],
            'invariant_observable_definition': 'a*sum_i n_i+sum_(i<j)b_graph_orbit(i,j)*n_i*n_j, restricted to total native occupation2',
            'arbitrary_invariant_diagonal_pair_compression': '(2*a+b_distance2)*I6+(b_antipodal-b_distance2)*|1><1|',
            'occupation_difference_rank_upper_for_fixed_2s_plus_n_factorization': 1,
            'uniform_two_hop_compression': _serialize(control['compression']),
            'uniform_two_hop_preserves_code_exactly': control['code_preserving'],
            'uniform_two_hop_occupation_differences': list(map(str, control['occupation_differences'])),
            'non_diagonal_or_context_breaking_operators_excluded': False}


def dictionary_control():
    """Hold H2 fixed while testing two CAR subsystem observable assignments.

    These are alternative observer dictionaries for ONE CAR mode with three
    multiplicity sectors. They are not three independent fermion fields and
    not countermodels after a complete physical dictionary has been fixed.
    """
    h = diagonal_cost(2)['compression']
    output = {}
    for name, pairs in (('existing_2s_plus_n', ((0, 1), (2, 3), (4, 5))),
                        ('alternative_same_operator', ((0, 3), (1, 4), (2, 5)))):
        annihilator = s.zeros(6)
        differences = []
        for empty, occupied in pairs:
            annihilator[empty, occupied] = 1
            differences.append(s.factor(h[occupied, occupied]-h[empty, empty]))
        charge = annihilator.T*annihilator
        if (annihilator**2 != s.zeros(6) or annihilator*annihilator.T+charge != s.eye(6)
                or charge*annihilator-annihilator*charge != -annihilator
                or h*charge != charge*h):
            raise ArithmeticError('same-operator observer dictionary fails CAR or conserved charge')
        output[name] = {'empty_occupied_pairs': [list(pair) for pair in pairs],
            'annihilator': _serialize(annihilator), 'occupation_charge': _serialize(charge),
            'occupation_differences': list(map(str, differences)),
            'one_CAR_mode_with_multiplicity_three': True}
    return {'fixed_native_instrument_and_H2_identical': True, 'dictionaries': output,
            'native_total_occupation_in_both_dictionaries': '2*I6',
            'native_raw_record_laws_distinguish_these_dictionaries': False,
            'physical_charge_or_family_dictionary_selected_here': False,
            'M15_gauge_currents_and_existing_observer_apparatus_preserved_or_replayed': False,
            'same_complete_Standard_Model_or_observer_premises_claimed': False,
            'scope': 'ambiguity before an observable attachment is supplied; does not alter a previously fixed field program'}


def _probe_energies():
    probe = diagonal_cost(2)
    exact_energies = [s.Integer(0)]*512
    for state in native.states():
        energy = sum(probe['coefficients'][a, b]*(state[a]-state[b])**2
                     for a in range(12) for b in range(a+1, 12))
        exact_energies[circuit.encode(state)] = energy
    return exact_energies


def verify_source_arrays(arrays):
    """Independent operator gates on retained binary arrays, without repair.

    Fresh platform output comparison is separate: an archived source result
    must itself meet the fixed all-column target gate, including absolute
    phase, raw zero outcomes and every syndrome/padding row.
    """
    shapes = {'source_encoded_six_columns': (512, 6), 'two_hop_phase_outputs': (512, 6),
              'native_C_corrected_outputs': (8, 512, 6), 'native_E_corrected_outputs': (16, 512, 6)}
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete quantum source output keys differ')
    for name, shape in shapes.items():
        value = arrays[name]
        if (type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype('complex128')
                or not np.isfinite(value).all()):
            raise ValueError('complete quantum source schema differs: '+name)
    ideal = code.ideal_code_matrix()
    phases = -float(PERTURBATION)*np.asarray(list(map(float, _probe_energies())))
    phase_target = np.exp(1j*phases)[:, None]*ideal
    metrics = {'encoder_operator_error': float(np.linalg.norm(arrays['source_encoded_six_columns']-ideal, 2)),
               'two_hop_source_operator_error': float(np.linalg.norm(arrays['two_hop_phase_outputs']-phase_target, 2))}
    for law, seam in (('C', (0, 1)), ('E', (2, 10))):
        value = arrays['native_'+law+'_corrected_outputs']
        targets = [math.sqrt(float(code.branch_weight(law, seam, coin, raw)))*ideal
                   for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law])]
        metrics[law] = {'maximum_branch_operator_error': max(float(np.linalg.norm(actual-target, 2))
                          for actual, target in zip(value, targets, strict=True)),
                       'aggregate_gram_error': float(np.linalg.norm(sum((b.conj().T@b for b in value), np.zeros((6, 6), complex))-np.eye(6), 2))}
    errors = [metrics['encoder_operator_error'], metrics['two_hop_source_operator_error']]
    errors += [value for law in ('C', 'E') for value in metrics[law].values()]
    if any(not math.isfinite(value) or value > SOURCE_GATE for value in errors):
        raise ValueError('retained source outputs exceed unchanged all-column operator/Gram gate')
    return metrics


def source_execution():
    """Existing emitted encoders/native instruments/correctors on full512 rows."""
    actual_code = code.code_matrix()
    ideal_code = code.ideal_code_matrix()
    exact_energies = _probe_energies()
    energies = np.asarray(list(map(float, exact_energies)))
    phases = -float(PERTURBATION)*energies
    # The existing convenience compiler is capped at8bits; the identical
    # complete commuting Walsh construction is written explicitly for this
    # existing9bit carrier. Every zero coefficient/scalar phase is retained.
    word = []
    for mask in range(512):
        coefficient = -s.Rational(PERTURBATION)*sum(value*(-1)**((address & mask).bit_count())
                                                  for address, value in enumerate(exact_energies))/512
        label = ['Z' if mask >> bit & 1 else 'I' for bit in range(9)]
        word.extend(compiler.pauli_rotation(label, float(coefficient)))
    perturbed = compiler.execute(word, 9, actual_code)
    independent = np.exp(1j*phases)[:, None]*ideal_code
    phase_error = float(np.linalg.norm(perturbed-independent, 2))
    if phase_error > SOURCE_GATE:
        raise ArithmeticError('actual two-hop phase source differs from independent full512 action')
    arrays = {'source_encoded_six_columns': actual_code, 'two_hop_phase_outputs': perturbed}
    rows = {}
    for law, seam in (('C', (0, 1)), ('E', (2, 10))):
        branches = code.repair_branches(law, seam, actual_code)
        outputs = np.stack([branch['amplitude'] for branch in branches])
        errors = [float(np.linalg.norm(branch['amplitude']-math.sqrt(float(branch['code_branch_weight']))*ideal_code, 2))
                  for branch in branches]
        if max(errors) > SOURCE_GATE:
            raise ArithmeticError('actual corrected source cycle differs from identity branch contract')
        arrays['native_'+law+'_corrected_outputs'] = outputs
        rows[law] = {'source_seam': list(seam), 'raw_branch_count': len(branches),
            'all512_rows_and_six_reference_columns_retained': True,
            'maximum_branch_operator_error': max(errors),
            'aggregate_gram_error': float(np.linalg.norm(sum((b.conj().T@b for b in outputs), np.zeros((6, 6), complex))-np.eye(6), 2)),
            'source_repair_resources': code.repair_resources(law, seam),
            'maximum_record_recovery_native_pulses': max(branch['recovery_resources']['native_pulses'] for branch in branches)}
    verify_source_arrays(arrays)
    return {'executed_source_native_seams': rows, 'encoder_resources': code.encoding_resources(),
            'two_hop_control_resources': compiler.resources(word), 'two_hop_source_operator_error': phase_error,
            'probe_parameter': str(PERTURBATION), 'parameter_is_programmed_coherent_probe_not_physical_time': True,
            'correlation_phase_extension': 'same graph cost on all78 total-two configurations; zero phase on all other padding rows',
            'exact_all30_ideal_certificate_is_separate_from_two_seam_floating_execution': True,
            'projection_or_padding_removal_or_probability_normalization_applied': False,
            'numerical_source_admission': SOURCE_GATE,
            'fresh_source_output_operator_replay_tolerance': SOURCE_REPLAY_TOLERANCE,
            'proper_services_and_physical_errors_authenticated_here': False}, arrays


def report():
    geometry = code_geometry()
    first, second = diagonal_cost(1), diagonal_cost(2)
    if first['compression'] != s.Rational(9, 32)*s.eye(6) or first['occupation_differences'] != (0, 0, 0):
        raise ArithmeticError('native weighted-edge agreement ceased to be scalar')
    if second['occupation_differences'] != (s.Rational(1, 256), -s.Rational(1, 1024), -s.Rational(1, 3072)):
        raise ArithmeticError('fixed completion two-hop candidate changed')
    cycles = corrected_cycle()
    source, arrays = source_execution()
    weights, timeout = completion_weights()
    return {'schema': 'oph.parameter_selection.v1',
        'classification': 'fixed corrected-native response obstruction and an unselected correlation-coupling candidate',
        'premises': {'same_native_C_E_primitives': True, 'fixed_code': 'record_code.MATCHINGS; six-dimensional sharp recorded code',
            'fixed_recovery': 'existing inverse movement plus matching Householder; existing independent512-source compilation',
            'fixed_selector_rounds': 1, 'accepted_probability': str(sum(weights.values(), F())),
            'terminal_identity_timeout_probability': str(timeout),
            'record_boundary': 'all raw native coin/value records and selector timeout; fine auxiliary selector transcripts marginalized under existing exact contract',
            'family_test_factorization': 'logical_label=2*spectator+occupation; this is one CAR mode with multiplicity3, not three independently occupiable species',
            'new_measured_inputs': [], 'new_logical_mass_hamiltonian_inserted_in_baseline': False},
        'fixed_cycle': {'C_E': cycles, 'exact_all_code_operator_transfer': 'I36',
            'complete_CQ_instrument': 'rho -> direct_sum_h p_h*rho, including all raw zeros and terminal identity timeout',
            'all_finite_adaptive_histories': 'record law depends on chosen seams and earlier records, never on logical state; corrected conditional work is unchanged',
            'inverse_response_for_exact_positive_Laplace_z': 'S(z)=z*I36', 'static_seed': '0_36x36',
            'resolvent_boundary': 'algebraic inverse of I-T+zI for one complete proposal; not continuation after a terminal timeout',
            'positive_source_family_gaps_generated': 0,
            'time_parameter': 'native completion proposal; no physical or modular frequency identified'},
        'native_observable_compressions': {'one_port': 'V.H*n_a*V=I6/6 for every a',
            'centered_native_slow_coordinates': 'V.H*sum_a(B_ai*n_a)*V=0 for each native slow axis i',
            'code_reference_mean_feature': _serialize(geometry['feature']),
            'original_faithful_reference_mean_feature': _serialize(geometry['original_feature']),
            'independent_reference_clock_copy_would_be_an_explicit_additional_tensor_factor': True,
            'core_count_histogram_in_each_matching': {'0': 4, '1': 2, '2': 0},
            'all_functions_of_existing_core_biased_reference_are_scalar_after_compression': True,
            'arbitrary_weighted_native_edge_diagonal_cost_is_scalar': True,
            'fixed_weighted_edge_cost': _serialize(first['compression'])},
        'first_surviving_correlation': {'definition': 'H2=sum_(a<b)(W**2)_ab*(n_a-n_b)**2; W_ab=existing accepted seam weight',
            'normalization_is_fixed_one_proposal_convention': True,
            'compressed_six_level_operator': _serialize(second['compression']),
            'occupation_differences_exact': list(map(str, second['occupation_differences'])),
            'absolute_difference_ratio': [12, 3, 1],
            'all_record_logical_block_first_derivative': '-i*[V.H*H2*V,rho]; verified by exact native record moments, including timeout',
            'syndrome_leakage_variance_by_logical_column': list(map(str, second['variance'])),
            'preserves_code_exactly': second['code_preserving'],
            'logical_block_compression_is_not_a_closed_finite_step_evolution': True,
            'already_present_as_coherent_term_in_fixed_native_cycle': False,
            'physical_coherent_coefficient_or_mass_scale_selected': False,
            'rapid_correction_limit_or_attractive_syndrome_elimination_proved': False},
        'symmetric_graph_control': graph_invariant_control(), 'same_operator_dictionary_control': dictionary_control(),
        'source_execution': source,
        'nearest_existing_results': {'passive_stationary_rank_bound': 'full_repair_response.stationary_selection; rank<=2 for same biased reference and nonnegative passive seam weights',
            'recorded_recovery_identity': 'record_code exact_certificate/execution_report; Knill-Laflamme and classical-feedback recovery',
            'addition_here': 'fixed active cycle logical inverse response and source-observable scalarization; explicit two-hop correlation compression and its unremoved syndrome leakage'},
        'primary_theory_sources': ['https://arxiv.org/abs/quant-ph/9604034', 'https://arxiv.org/abs/quant-ph/0403092'],
        'scope': {'Standard_Model_mass_derivation_completed': False, 'calibrated_packets_used': False,
            'active_correction_alone_selects_Yukawa_matrix': False,
            'other_codes_references_coherent_interactions_and_non_diagonal_native_correlations_excluded': False}}, arrays


def source_bindings():
    paths = {name for name in preparation.SOURCE_PATHS if not name.endswith('.md')}
    paths.update({'model/common_source/'+name for name in (
        'parameter_selection.py', 'test_parameter_selection.py', 'record_code.py', 'native_routing.py',
        'experiment.py', 'slow_geometry.py')})
    paths.update({'model/observer/native.py', 'model/observer/quantum_mixing.py',
        'model/observer/source_scalar/compiler.py', 'model/observer/source_discrimination/circuit.py',
        'model/observer/source_discrimination/experiment.py', 'model/repair_scheduler/coins.py',
        'tools/strict_json.py', 'requirements.txt', 'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'})
    paths.update(compiler.source_summary()['source_sha256'])
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if set(expected) != set(actual):
            raise ValueError(path+' complete keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' length differs')
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _compare(left, right, path+f'[{i}]')
    elif type(expected) is float:
        if path.endswith(('numerical_source_admission', 'fresh_source_output_operator_replay_tolerance')) and expected != actual:
            raise ValueError(path+' fixed admission differs')
        if not math.isfinite(actual) or not math.isclose(expected, actual, rel_tol=2e-10, abs_tol=2e-12):
            raise ValueError(path+' fresh numerical diagnostic differs')
    elif expected != actual:
        raise ValueError(path+' differs')


def check(path):
    retained = strict_json.load(path)
    if type(retained) is not dict or retained.get('source_bindings') != source_bindings():
        raise ValueError('retained closed source bindings differ')
    artifact = path.with_suffix('.npz')
    if retained.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('retained complete source artifact custody differs')
    with np.load(artifact, allow_pickle=False) as archive:
        if len(archive.files) != len(set(archive.files)):
            raise ValueError('duplicate source artifact keys')
        arrays = {key: archive[key] for key in archive.files}
    fresh, outputs = report()
    # The retained binary values define their own reported diagnostics. This
    # avoids requiring a fresh BLAS run to reproduce floating residuals while
    # the independent all-column gate remains unchanged on BOTH executions.
    retained_metrics = verify_source_arrays(arrays)
    fresh['source_execution']['two_hop_source_operator_error'] = retained_metrics['two_hop_source_operator_error']
    for law in ('C', 'E'):
        fresh['source_execution']['executed_source_native_seams'][law].update(retained_metrics[law])
    _compare(fresh, {k: v for k, v in retained.items() if k not in ('source_bindings', 'artifact')})
    if set(arrays) != set(outputs):
        raise ValueError('complete quantum source output keys differ')
    for name, expected in outputs.items():
        value = arrays[name]
        blocks = zip(value, expected, strict=True) if value.ndim == 3 else [(value, expected)]
        if (value.dtype != expected.dtype or value.shape != expected.shape or not np.isfinite(value).all()
                or any(np.linalg.norm(left-right, 2) > SOURCE_REPLAY_TOLERANCE for left, right in blocks)):
            raise ValueError('complete quantum source output differs: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('fixed corrected-native parameter-selection replay: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('choose fresh receipt and source artifact paths')
    result, arrays = report()
    result['source_bindings'] = source_bindings()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with artifact.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    result['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    print('fixed cycle: zero logical seed; two-hop candidate:12:3:1, with retained syndrome leakage')


if __name__ == '__main__':
    main()
