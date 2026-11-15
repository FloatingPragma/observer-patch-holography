"""Finite native force readout on the frozen reciprocal classical source.

Only the existing pair-table compiler and literal native rule are consumed.
No full pair generator or pair semigroup is constructed. Every first-event
mark, holding mark and cap output is retained in a compressed outcome archive.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BRIDGE = ROOT / 'openwave-tuning/auditor/native_field_bridge'
FIELD = HERE.parent / 'field_limit'
FROZEN = '4597dc9c226d203d890bc3d3a929bbe1be80e2dc57ec03898c77da4fe5b954e1'
FIELD_FROZEN = '612aedd40f13faea2b80b5c31d4afbb19cb2332ff22db5b459ac2ae936879008'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def dependencies():
    out = {}
    for directory, expected in ((BRIDGE, FROZEN), (FIELD, FIELD_FROZEN)):
        name = 'FINAL_MANIFEST.json' if directory == BRIDGE else 'MANIFEST.json'
        manifest_path = directory / name
        assert sha(manifest_path) == expected
        manifest = json.loads(manifest_path.read_text())
        out[str(manifest_path.relative_to(ROOT))] = expected
        for rel, pin in manifest['files'].items():
            path = directory / rel
            assert sha(path) == pin, str(path)
            out[str(path.relative_to(ROOT))] = pin
    # Resolve transitive pins from receipts whose bytes are themselves frozen.
    # Merely recording the current native file would admit a different source.
    for receipt in (BRIDGE / 'RECIPROCAL_FEEDBACK.json', FIELD / 'RESULT.json'):
        for rel, expected in json.loads(receipt.read_text())['source_pins'].items():
            assert sha(ROOT / rel) == expected, rel
            out[rel] = expected
    load('field_feedback', BRIDGE / 'field_feedback.py')
    integer = HERE.parent / 'coordination/INTEGER_FIELD_CONTROL.json'
    out[str(integer.relative_to(ROOT))] = sha(integer)
    return out


def counters(state, n=96):
    return [[int(i == p) for i in range(n)] for p in state]


def first_event(state, edges, rates, native, horizon=.001, target_seam=0):
    """Exact cap-one marked law, evaluated through the native endpoint rule.

    Columns are probability, copy, seam, coin, native current, measured
    current, terminal x, terminal y, cap flag. Row zero is no event.
    All events, including holds, carry the cap flag and their terminal state.
    """
    total = float(rates.sum())
    no_event = math.exp(-total * horizon)
    exposure = -math.expm1(-total * horizon) / total
    rows = [[no_event, -1, -1, -1, 0, 0, *state, 0]]
    initial = [tuple(c) for c in counters(state)]
    for side in range(2):
        for seam, edge in enumerate(edges):
            for coin in (0, 1):
                after, (record_coin, current) = native(initial[side], edge, coin)
                assert record_coin == coin
                final = list(state)
                final[side] = int(np.argmax(after))
                # Independent endpoint rule on the one-load domain.
                expected = (edge[0] if coin else edge[1]) if state[side] in edge else state[side]
                assert final[side] == expected
                assert current == int(state[side] == edge[0]) - int(expected == edge[0])
                measured = current if side == 0 and seam == target_seam else 0
                rows.append([exposure * rates[side, seam] / 2,
                             side, seam, coin, current, measured, *final, 1])
    law = np.array(rows, dtype=float)
    assert abs(law[:, 0].sum() - 1) < 3e-14
    return law, exposure, total


def measure(partner, edges, table, native, reciprocal, kappa=1., nu=1.,
            beta=4., horizon=.001, held_partner=None):
    """Equiprobable endpoint preparation, with an unconditional score.

    Each preparation has its own recorded exposure g. Score=current/g,
    including zero on no-event and all unrelated marks. No cap branch is
    removed, and no success-conditioned normalization is performed.
    """
    seam = edges.index((0, 1))
    laws = []
    exposures = []
    rates_total = []
    for endpoint in (0, 1):
        state = (endpoint, partner)
        rates = reciprocal.pair_rates(counters(state), edges, table, kappa, beta, nu)
        if held_partner is not None:
            held = reciprocal.pair_rates(counters((endpoint, held_partner)),
                                         edges, table, kappa, beta, nu)
            rates[0] = held[0]
        law, exposure, total = first_event(state, edges, rates, native, horizon, seam)
        laws.append(law)
        exposures.append(exposure)
        rates_total.append(total)
    drift = sum(float(law[:, 0] @ law[:, 5]) / g / 2
                for law, g in zip(laws, exposures))
    second_rate = sum(float(law[:, 0] @ (law[:, 5] ** 2)) / g / 2
                      for law, g in zip(laws, exposures))
    field_partner = partner if held_partner is None else held_partner
    force = float(kappa * (table[0, field_partner] - table[1, field_partner]))
    predicted = nu * math.sinh(beta * force / 2) / 2
    assert abs(drift - predicted) < 3e-14
    g_min = min(exposures)
    return {'partner': partner, 'held_partner': held_partner, 'kappa': kappa,
            'beta': beta, 'nu': nu, 'horizon': horizon, 'cap': 1,
            'seam': [0, 1], 'force_minus_energy_difference': force,
            'measured_drift_per_source_time': drift, 'predicted_drift': predicted,
            'measured_second_current_rate': second_rate,
            'exposure_by_preparation': exposures, 'total_attempt_rates': rates_total,
            'complete_probability': [float(v[:, 0].sum()) for v in laws],
            'cap_probabilities': [float(v[v[:, 8] == 1, 0].sum()) for v in laws],
            'hold_probabilities': [float(v[(v[:, 8] == 1) & (v[:, 4] == 0), 0].sum())
                                   for v in laws],
            'outcome_count': sum(len(v) for v in laws),
            'score_absolute_bound': 1 / g_min}, np.stack(laws)


def arithmetic_transfer(edges, neighbours, table, field, reciprocal, native):
    """Named comparison architecture, separate from the frozen scheduler.

    Each candidate pair uses its own zero-start 512-round 40-bit solve.
    This evaluation compiles only pairs required by the cap-one instrument;
    the analytic bound is uniform on all 9216 pairs.
    """
    rounds, bits = 512, 40
    alpha = F(23, 24) ** rounds
    eta = 2 * alpha + F(408, 1 << bits) * (1 - alpha)
    frozen_energy_error = 2 * alpha
    delta = eta + frozen_energy_error
    needed = set()
    for endpoint in (0, 1):
        state = (endpoint, 9)
        needed.add(state)
        needed.update((int(dest), 9) for dest in neighbours[endpoint])
        needed.update((endpoint, int(dest)) for dest in neighbours[9])
    energy = {}
    residuals = []
    for state in sorted(needed):
        charge = np.array(counters(state)[0]) - np.array(counters(state)[1])
        z = field.local_fixed_point([0] * 96, charge, neighbours, rounds, bits)
        value = F(sum(int(q) * zz for q, zz in zip(charge, z)), 2 * (1 << bits))
        residual = field.residual_squared_certificate(z, charge, neighbours, bits)
        energy[state] = float(value)
        assert abs(energy[state] - table[state]) <= float(delta) + 2e-14
        residuals.append({'state': state, 'exact_energy': str(value),
                          'residual_field_error_squared_upper': str(residual),
                          'frozen_energy_difference': energy[state] - float(table[state])})
    laws, frozen_laws, exposures = [], [], []
    for endpoint in (0, 1):
        state = (endpoint, 9)
        rr = np.ones((2, len(edges)))
        for side, position in enumerate(state):
            for seam, (a, b) in enumerate(edges):
                if position not in (a, b):
                    continue
                dest = b if position == a else a
                after = (dest, state[1]) if side == 0 else (state[0], dest)
                rr[side, seam] = math.exp(-2 * (energy[after] - energy[state]))
        law, _, _ = first_event(state, edges, rr, native)
        frozen_rates = reciprocal.pair_rates(counters(state), edges, table)
        frozen_law, exposure, _ = first_event(state, edges, frozen_rates, native)
        laws.append(law)
        frozen_laws.append(frozen_law)
        exposures.append(exposure)
    lbound = 4 * delta  # beta=4, kappa=1: beta*kappa*delta.
    # exp(l)-1 <= l/(1-l), an exact rational outward bound.
    assert 0 < lbound < 1
    rate_difference_bound = 96 * lbound / (1 - lbound)
    time = F(1, 1000)
    tv_bound = time * rate_difference_bound
    tv = [float(abs(a[:, 0] - b[:, 0]).sum() / 2)
          for a, b in zip(laws, frozen_laws)]
    assert max(tv) <= float(tv_bound)
    # The same frozen exposure score is used in both laws for this transfer.
    observed_error = abs(sum(float((a[:, 0] - b[:, 0]) @ a[:, 5]) / g / 2
                             for a, b, g in zip(laws, frozen_laws, exposures)))
    score_error_bound = 2 * float(tv_bound) / min(exposures)
    assert observed_error <= score_error_bound
    return {'variant': 'direct-pair-zero-start512-fixedpoint40-electric-scheduler',
            'same_native_rule_preparation_and_clock_parameter': True,
            'frozen_scheduler_replaced': False,
            'comparison_domain': 'beta=4, kappa=1, nu=1, horizon=1/1000; any finite native attempt cap',
            'zero_start_field_error_upper': field.rounded_upper(eta),
            'unit_energy_difference_to_ideal_frozen_upper': field.rounded_upper(delta),
            'log_seam_rate_ratio_absolute_upper': field.rounded_upper(lbound),
            'complete_marked_capped_path_TV_upper': field.rounded_upper(tv_bound),
            'measured_first_event_TV': tv,
            'measured_frozen_score_error': observed_error,
            'frozen_score_error_upper': score_error_bound,
            'compiled_pair_count': len(needed), 'compiled_pairs': residuals,
            'rounding_scope': 'Integer controller arithmetic is exact. The coupling certificate compares exact-real exponentials and the exact-arithmetic frozen local stencil. Floating-library evaluation error is separately unproved.'}, np.stack(laws)


def run():
    before = dependencies()
    fixed = sys.modules['field_feedback']
    reciprocal = load('reciprocal_feedback', BRIDGE / 'reciprocal_feedback.py')
    field = load('force_controlled_field', FIELD / 'controlled_field.py')
    edges, neighbours = fixed.graph()
    table = reciprocal.pair_table(neighbours)
    native = fixed.literal_native()
    archive = {}
    measurements = []
    for kappa in (0., .01, .1, 1., 2.):
        item, law = measure(9, edges, table, native, reciprocal, kappa=kappa)
        item['archive_key'] = f'kappa_{kappa:g}'
        archive[item['archive_key']] = law
        force = item['force_minus_energy_difference']
        z = 2 * force
        linear = force  # beta D0 = 4*(nu/4)=1 at nu=1.
        bound = abs(z) ** 3 * math.cosh(abs(z)) / 12
        error = abs(item['measured_drift_per_source_time'] - linear)
        assert error <= bound + 3e-14
        item['Einstein_linear_drift'] = linear
        item['linear_response_remainder'] = error
        item['Taylor_remainder_bound'] = bound
        measurements.append(item)
    zero, primary = measurements[0], measurements[3]
    assert abs(zero['measured_drift_per_source_time']) < 1e-14
    assert abs(zero['measured_second_current_rate'] / 2 - .25) < 1e-14
    assert primary['measured_drift_per_source_time'] > .03
    # The finite hop changes its own charge field. The midpoint/self term
    # is required when turning a present field into an energy increment.
    rho = np.array(counters((0, 9))[0]) - np.array(counters((0, 9))[1])
    phi = fixed.local_field(rho, neighbours, 512)
    seam_field = float(phi[0] - phi[1])
    self_increment = float(table[0, 1])
    assert abs(primary['force_minus_energy_difference'] - (seam_field - self_increment)) < 3e-14
    field_to_energy = {'present_oriented_seam_field': seam_field,
                       'half_edge_resistance_self_increment': self_increment,
                       'force_with_self_response': seam_field - self_increment,
                       'bare_present_field_force_error': self_increment,
                       'identity': 'Delta U/kappa = -sigma*(B_e^T phi) + (B_e^T G_K B_e)/2'}
    # Genuine graph reflection: x-coordinate reflection of the local frame.
    port_reflection = [1, 0, 3, 2, 4, 5, 6, 7, 10, 11, 8, 9]
    perm = [12 * (i // 12) + port_reflection[i % 12] for i in range(96)]
    assert {tuple(sorted((perm[a], perm[b]))) for a, b in edges} == set(edges)
    reflection_error = float(np.max(abs(table - table[np.ix_(perm, perm)])))
    assert reflection_error < 3e-14
    reflected, law = measure(11, edges, table, native, reciprocal)
    archive['reflection'] = law
    assert abs(reflected['measured_drift_per_source_time'] + primary['measured_drift_per_source_time']) < 1e-14
    # Real native move of the other copy; retain its current as intervention.
    candidate = max((int(p) for p in neighbours[9] if p not in (0, 1)),
                    key=lambda p: abs((table[0, p] - table[1, p]) - (table[0, 9] - table[1, 9])))
    edge = tuple(sorted((9, candidate)))
    coin = int(candidate == edge[0])
    moved, (_, current) = native(tuple(counters((0, 9))[1]), edge, coin)
    assert int(np.argmax(moved)) == candidate and abs(current) == 1
    changed, law = measure(candidate, edges, table, native, reciprocal)
    archive['other_copy_moved'] = law
    held, law = measure(candidate, edges, table, native, reciprocal, held_partner=9)
    archive['held_field'] = law
    cross_response = changed['measured_drift_per_source_time'] - primary['measured_drift_per_source_time']
    held_response = held['measured_drift_per_source_time'] - primary['measured_drift_per_source_time']
    assert abs(cross_response) > .001 and abs(held_response) < 1e-14
    doubled, law = measure(9, edges, table, native, reciprocal, nu=2.)
    archive['double_nu'] = law
    assert abs(doubled['measured_drift_per_source_time'] - 2 * primary['measured_drift_per_source_time']) < 1e-14
    # Actual short executions exercise the frozen generator's native/cap path.
    paths = []
    for cap in (1, 2):
        path = reciprocal.execute((0, 9), edges, table, 1., 89, horizon=.5, cap=cap)
        path['independent_replay'] = reciprocal.replay(path, edges, table)
        assert path['terminal']['kind'] == 'cap_exhausted'
        paths.append(path)
    variant, law = arithmetic_transfer(edges, neighbours, table, field, reciprocal, native)
    archive['fixedpoint40_variant'] = law
    archive_path = HERE / 'FIRST_EVENT_OUTCOMES.npz'
    np.savez_compressed(archive_path, **archive)
    assert dependencies() == before
    result = {'schema': 'oph.native_force_measurement.v1', 'source_pins': before,
              'producer_sha256': sha(Path(__file__)),
              'checks': {'actual_literal_native_marked_outcomes': True,
                         'full_first_event_law_and_cap_retained': True,
                         'finite_drift_force_identity': True,
                         'Einstein_linear_response_with_bound': True,
                         'zero_coupling_and_reflection': True,
                         'actual_partner_move_and_held_field': True,
                         'common_nu_scaling': True,
                         'bounded_fixedpoint_variant_transfer': True,
                         'source_pins_unchanged': True},
              'measurements': measurements,
              'field_to_energy': field_to_energy,
              'reflection': {'port_permutation': port_reflection, 'graph_automorphism': True,
                             'energy_invariance_error': reflection_error, 'measurement': reflected},
              'held_field_control': {'intervention': {'copy': 1, 'before': 9, 'after': candidate,
                                                      'seam': edge, 'coin': coin, 'current': current},
                                     'reciprocal': changed, 'held_field': held,
                                     'measured_cross_response': cross_response,
                                     'held_field_cross_response': held_response},
              'clock_control': doubled, 'actual_cap_paths': paths,
              'fixedpoint_variant_transfer': variant,
              'outcome_archive': {'path': archive_path.name, 'sha256': sha(archive_path),
                                   'keys': list(archive), 'shape_each': [2, 1153, 9],
                                   'preparation_probabilities': [.5, .5],
                                   'columns': ['probability', 'copy', 'seam', 'coin', 'native_current',
                                               'measured_current', 'terminal_x', 'terminal_y', 'cap_flag']},
              'official_rows_changed': []}
    (HERE / 'RESULT.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'drift': primary['measured_drift_per_source_time'],
                      'force': primary['force_minus_energy_difference'], 'cross_response': cross_response,
                      'variant_TV_bound': variant['complete_marked_capped_path_TV_upper']}, indent=2))


if __name__ == '__main__':
    run()
