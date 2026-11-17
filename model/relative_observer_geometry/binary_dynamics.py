"""Recover spatial response and propagation inside the binary quantum sector.

Calibration records contain current, not the hidden placement coin. Fresh
known preparations and addressed seams are declared services. The native fair
coin, finite collection cap and every failed collection remain explicit.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path
import random

from ..observer import native
from ..observer_reconstruction import native_geometry as ng
from ..primitive_observer import transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from . import intrinsic_geometry as intrinsic
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
PREPARATIONS = tuple((0, i) for i in range(1, 12))+((1, 2),)
CAP = 24
SEED = 1729
SOURCE_PATHS = tuple(sorted(set((*transfer.SOURCE_PATHS, *intrinsic.SOURCE_PATHS,
    'model/relative_observer_geometry/binary_dynamics.py',
    'model/observer/quantum_mixing.py',
    'model/native_quantum_history/construction.py',
))))


def _loads(pair):
    return tuple(int(i in pair) for i in range(12))


def _matmul(left, right):
    return [[sum((x*y for x, y in zip(row, column, strict=True)), F(0))
             for column in zip(*right, strict=True)] for row in left]


def _encode(matrix):
    return [[str(x) for x in row] for row in matrix]


def execute_calibration(seams, *, seed=SEED, cap=CAP, coin_source=None):
    """Collect the deterministic current support with finite native attempts.

    Coin zero/one is internal. Only (preparation, oriented seam, current) is
    retained. Each attempt starts from the same disclosed binary preparation.
    A custom coin source is only for adversarial timeout tests.
    """
    if type(cap) is not int or cap < 1:
        raise ValueError('positive integer calibration cap required')
    rng = random.Random(seed)
    draw = (lambda: rng.randrange(2)) if coin_source is None else coin_source
    records, timeouts, attempts_by_setting = [], [], []
    for p, pair in enumerate(PREPARATIONS):
        initial = _loads(pair)
        for a, b in seams:
            delta = initial[a]-initial[b]
            expected = {0, delta}
            seen = set()
            for used in range(1, cap+1):
                coin = draw()
                if type(coin) is not int or coin not in (0, 1):
                    raise ValueError('the native placement coin must be binary')
                _, (_, current) = native.native(initial, (a, b), coin)
                records.append([p, a, b, current])
                seen.add(current)
                if seen == expected:
                    break
            attempts_by_setting.append([p, a, b, used])
            if seen != expected:
                timeouts.append([p, a, b])
    return {'preparations': [list(pair) for pair in PREPARATIONS],
            'current_records': records, 'timeout_settings': timeouts,
            'attempts_by_setting': attempts_by_setting,
            'cap_per_setting': cap, 'success': not timeouts}


def recover_response(calibration):
    """Infer the response from complete current-support records in M_66.

    The declared fair-coin rule assigns weights 1/2 to the two current
    outcomes when endpoint occupations differ. Their relative sample counts
    are not probability estimates, and fairness is not inferred from data.
    This is inference from a current table plus its reported completion flag;
    it does not authenticate cap metadata, preparation services or entropy.
    The execution driver enforces the reported finite collection protocol.
    """
    if type(calibration) is not dict or not calibration.get('success') or calibration.get('timeout_settings'):
        raise ValueError('a timed-out calibration supplies no recovered geometry')
    preparations = calibration.get('preparations')
    if preparations != [list(pair) for pair in PREPARATIONS]:
        raise ValueError('the disclosed spanning binary preparation family is required')
    records = calibration.get('current_records')
    if type(records) is not list or not records:
        raise ValueError('nonempty current records required')
    support = {}
    for row in records:
        if (type(row) is not list or len(row) != 4 or any(type(x) is not int for x in row)):
            raise ValueError('records contain only integer preparation, seam and current fields')
        p, a, b, current = row
        if p not in range(12) or not 0 <= a < b < 12:
            raise ValueError('valid preparation and oriented seam labels required')
        initial = _loads(PREPARATIONS[p])
        if current not in {0, initial[a]-initial[b]}:
            raise ValueError('current contradicts its known binary preparation')
        support.setdefault((p, a, b), set()).add(current)
    seam_sets = [{(a, b) for q, a, b in support if q == p} for p in range(12)]
    if not seam_sets[0] or any(seams != seam_sets[0] for seams in seam_sets):
        raise ValueError('all preparations require the same complete seam support')
    seams = sorted(seam_sets[0])
    S = [[F(int(v in pair)) for pair in PREPARATIONS] for v in range(12)]
    inverse = intrinsic._inverse(S)
    if _matmul(S, inverse) != [[F(i == j) for j in range(12)] for i in range(12)]:
        raise ArithmeticError('the binary preparation family does not span port space')
    output = [row[:] for row in S]
    for p, pair in enumerate(PREPARATIONS):
        initial = _loads(pair)
        for a, b in seams:
            expected = {0, initial[a]-initial[b]}
            if support[p, a, b] != expected:
                raise ValueError('both possible current outcomes must be retained before averaging')
            mean = sum(support[p, a, b], F(0))/len(support[p, a, b])
            output[a][p] -= mean/len(seams)
            output[b][p] += mean/len(seams)
    response = _matmul(output, inverse)
    laplacian = [[2*len(seams)*(F(i == j)-response[i][j])
                  for j in range(12)] for i in range(12)]
    checks, projector = ng._band_certificate(laplacian)
    if len(seams) != 30 or not all(checks.values()):
        raise ValueError('current response does not certify the native carrier geometry')
    gram = ng._encode_q5([[ng._qscale(x, 4) for x in row] for row in projector])
    return {'seams': [list(e) for e in seams], 'preparation_inverse_exact': _encode(inverse),
            'response_exact': _encode(response), 'laplacian_exact': _encode(laplacian),
            'gram_qsqrt5': gram, 'band_certificate': checks,
            'probability_or_coin_records_used': False,
            'declared_fair_coin_law_used': True,
            'record_estimator_authenticates_physical_services_or_cap_metadata': False,
            'preparations_all_in_binary_66_state_sector': True}


def stopping_certificate(seams, cap):
    crossing = sum(_loads(pair)[a] != _loads(pair)[b]
                   for pair in PREPARATIONS for a, b in seams)
    fixed = len(PREPARATIONS)*len(seams)-crossing
    # A failed crossing setting saw only zero or only nonzero current.
    per_setting = F(2, 2**cap)
    return {'crossing_settings': crossing, 'constant_current_settings': fixed,
            'timeout_probability_per_crossing_setting': str(per_setting),
            'probability_of_any_timeout_upper_bound': str(min(F(1), crossing*per_setting)),
            'maximum_native_attempts': fixed+crossing*cap,
            'unbounded_mean_attempts': fixed+3*crossing,
            'coin_access_required': False,
            'native_fairness_fresh_preparation_and_record_access_are_premises': True,
            'each_native_coin_is_fair_conditioned_on_the_full_past': True,
            'seeded_pseudorandom_replay_authenticates_independent_entropy': False,
            'timeout_is_retained_not_discarded': True,
            'successful_reconstruction_exact_within_declared_model': True,
            'laboratory_error_budget_supplied': False}


def heisenberg_and_fluctuation(kernel, response):
    """Exact quantum first-moment closure and stationary lag covariance.

    Each native Kraus branch is a partial injection, so its pullback of any
    diagonal occupation is diagonal. Explicit branch counts then prove the
    Heisenberg identity for every initial density, including coherence.
    """
    state_index = {state: i for i, state in enumerate(kernel.states)}
    branch_counts = {}
    for law in ('C', 'E'):
        branches = {}
        for a, b in kernel.seams:
            for coin in (0, 1):
                for i, state in enumerate(kernel.states):
                    output, (_, current) = native.native(state, (a, b), coin)
                    label = current if law == 'C' else state[a]-state[b]
                    branches.setdefault((a, b, coin, label), {})[i] = state_index[output]
        moment = [[F(0) for _ in range(12)] for _ in kernel.states]
        for branch in branches.values():
            if len(set(branch.values())) != len(branch):
                raise ArithmeticError('a native binary Kraus branch is not a partial injection')
            for i, j in branch.items():
                for v in range(12):
                    moment[i][v] += F(kernel.states[j][v], 60)
        expected = [[sum(response[v][u]*state[u] for u in range(12))
                     for v in range(12)] for state in kernel.states]
        if moment != expected:
            raise ArithmeticError('native quantum occupations fail the inferred response identity')
        branch_counts[law] = len(branches)
    Q = [[F(i == j)-F(1, 12) for j in range(12)] for i in range(12)]
    fields = [[F(value)-F(1, 6) for value in state] for state in kernel.states]
    propagated = [row[:] for row in fields]
    expected = [[F(5, 33)*x for x in row] for row in Q]
    lag_checks = []
    for lag in range(9):
        if lag in (0, 1, 2, 4, 8):
            covariance = [[sum(fields[s][u]*propagated[s][v] for s in range(66))/66
                           for v in range(12)] for u in range(12)]
            if covariance != expected:
                raise ArithmeticError('occupation fluctuations disagree with response propagation')
            lag_checks.append({'native_attempt_lag': lag,
                               'covariance_equals_5_over_33_T_power_Q_exact': True,
                               'covariance_sha256': hashlib.sha256(canonical(_encode(covariance)).encode()).hexdigest()})
        propagated = [[sum(F(count, 60)*propagated[j][v]
                           for j, count in enumerate(row) if count)
                       for v in range(12)] for row in kernel.counts]
        expected = _matmul(expected, response)
    contrast = [F(int(i == 0)-int(i == 1)) for i in range(12)]
    after = [sum(response[i][j]*contrast[j] for j in range(12)) for i in range(12)]
    before_norm = F(5, 33)*sum(x*x for x in contrast)
    after_norm = F(5, 33)*sum(x*x for x in after)
    if not 0 < after_norm < before_norm:
        raise ArithmeticError('native occupation contraction witness changed')
    return {'quantum_kraus_branch_counts': branch_counts,
            'all_C_and_E_branches_are_partial_injections': True,
            'Heisenberg_occupation_identity_exact_for_every_binary_density': True,
            'stationary_centered_occupation_gram': '(5/33)*(I-11^T/12)',
            'constant_free_GNS_dimension': 11,
            'slow_response_GNS_subspace_dimension': 3,
            'lag_covariance_formula': 'C_k=(5/33)*T^k*(I-11^T/12)',
            'lag_checks': lag_checks,
            'stationary_law_used_for_covariance_not_for_calibration': True,
            'explicit_occupation_contraction': {
                'observable': 'N_0-N_1', 'initial_GNS_norm_squared': str(before_norm),
                'one_attempt_GNS_norm_squared': str(after_norm),
                'one_native_attempt_is_not_a_state_preserving_automorphism_on_these_observables': True},
            'global_no_wave_or_no_quantum_field_theorem': False}


def propagation_certificate(response, seams):
    adjacency = {i: set() for i in range(12)}
    for a, b in seams:
        adjacency[a].add(b); adjacency[b].add(a)
    distance = []
    for start in range(12):
        reached, queue = {start: 0}, [start]
        for i in queue:
            for j in sorted(adjacency[i]):
                if j not in reached:
                    reached[j] = reached[i]+1; queue.append(j)
        if len(reached) != 12:
            raise ArithmeticError('the reconstructed graph is disconnected')
        distance.append(reached)
    powers = [[[F(i == j) for j in range(12)] for i in range(12)]]
    for _ in range(3):
        powers.append(_matmul(powers[-1], response))
    for i in range(12):
        for j in range(12):
            d = distance[i][j]
            if any(powers[k][i][j] for k in range(d)) or powers[d][i][j] <= 0:
                raise ArithmeticError('native first-moment response violates its exact graph cone')
    # This zero-total perturbation is the difference of two allowed binary
    # preparations related by one seam move, not an inadmissible added load.
    witness = None
    for a, b in seams:
        for receiver in range(12):
            if distance[receiver][a] == 2 and distance[receiver][b] == 3:
                spectator = next(v for v in range(12) if v not in (a, b, receiver))
                values = [power[receiver][a]-power[receiver][b] for power in powers]
                if values[0] != 0 or values[1] != 0 or values[2] <= 0:
                    raise ArithmeticError('physical binary preparation difference failed front readout')
                witness = {'preparation_A': sorted((a, spectator)),
                           'preparation_B': sorted((b, spectator)),
                           'receiver_port': receiver,
                           'occupation_mean_difference_by_attempt_0_to_3': list(map(str, values))}
                break
        if witness:
            break
    if witness is None:
        raise ArithmeticError('no native binary two-hop response witness')
    return {'all_port_pairs_and_native_attempts_before_graph_distance_checked': True,
            'first_nonzero_response_kernel_entry_occurs_at_graph_distance': True,
            'individual_kernel_columns_are_not_fixed_total_preparation_interventions': True,
            'binary_sector_preparation_contrast': witness,
            'time_parameter': 'native attempts with independent uniform seam choices',
            'metric_parameter': 'native seam hops, not calibrated length',
            'front_speed_in_hops_per_native_attempt': 1,
            'physical_Lorentz_cone_or_lab_speed_identified': False,
            'nonconstant_mean_multipliers': [
                {'value': '(55+sqrt(5))/60', 'multiplicity': 3},
                {'value': '9/10', 'multiplicity': 5},
                {'value': '(55-sqrt(5))/60', 'multiplicity': 3}],
            'all_nonconstant_occupation_modes_contract_without_phase_rotation': True,
            'scope_of_contraction_obstruction': 'direct unconditioned occupation/current means under positive scalar reparameterization; not Euclidean reconstruction or a global no-QFT claim'}


def run():
    kernel = transfer.native_counts()
    transfer.exact_audit(kernel)
    calibration = execute_calibration(kernel.seams)
    response = recover_response(calibration)
    T = [[F(x) for x in row] for row in response['response_exact']]
    if response['gram_qsqrt5'] != ng.recover(ng.execute_calibration())['gram_qsqrt5']:
        raise ArithmeticError('binary and integer calibration geometries disagree')
    return {'schema': 'oph.relative_observer_geometry.binary_dynamics.v1',
            'sources': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}),
            'calibration_seed': SEED, 'sampled_calibration': calibration,
            'calibration_stopping_certificate': stopping_certificate(kernel.seams, CAP),
            'binary_record_response': response,
            'same_gram_as_existing_integer_calibration_exact': True,
            'quantum_occupation_and_fluctuation_bridge': heisenberg_and_fluctuation(kernel, T),
            'propagation': propagation_certificate(T, kernel.seams),
            'scope': {'main_calibration_and_dynamic_states_in_binary_66_state_sector': True,
                      'integer_calibration_used_only_as_external_comparison_control': True,
                      'native_C_current_readout_used_for_calibration': True,
                      'hidden_placement_coin_is_not_an_observer_input': True,
                      'one_quantum_occupation_response_for_both_C_and_E': True,
                      'single_ongoing_unreset_calibration_history': False,
                      'modular_parameter_identified_with_native_attempts': False,
                      'spatial_continuum_or_physical_clock_units_recovered': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Binary spatial dynamics receipt reproduced')
    elif args.write:
        with args.write.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
