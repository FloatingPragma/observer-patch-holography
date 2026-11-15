"""Adversarial semantic verification of the native first-event instrument."""
from __future__ import annotations

import copy
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np

import force_measurement as producer

HERE = Path(__file__).resolve().parent


def validate(laws, row, edges, table):
    assert laws.shape == (2, 1153, 9)
    assert row['held_partner'] is None
    y, kappa, beta, nu, time = [row[k] for k in ('partner', 'kappa', 'beta', 'nu', 'horizon')]
    drift = 0.
    for x, law in enumerate(laws):
        # Independent reconstruction of seam rates, native target and coin.
        rates = np.full((2, len(edges)), nu)
        for side, p in enumerate((x, y)):
            for e, (a, b) in enumerate(edges):
                if p in (a, b):
                    dest = a + b - p
                    target = (dest, y) if side == 0 else (x, dest)
                    rates[side, e] = nu * math.exp(-beta * kappa * (table[target] - table[x, y]) / 2)
        total = float(rates.sum())
        exposure = -math.expm1(-total * time) / total
        assert np.max(abs(law[0] - [math.exp(-total*time), -1, -1, -1, 0, 0, x, y, 0])) < 3e-14
        assert abs(law[:, 0].sum() - 1) < 3e-14
        for index, actual in enumerate(law[1:]):
            side, rem = divmod(index, 2 * len(edges))
            seam, coin = divmod(rem, 2)
            a, b = edges[seam]
            state = [x, y]
            old = state[side]
            target = (a if coin else b) if old in (a, b) else old
            current = int(old == a) - int(target == a)
            measured = current if side == 0 and seam == 0 else 0
            state[side] = target
            expected = [exposure * rates[side, seam]/2, side, seam, coin,
                        current, measured, *state, 1]
            assert np.max(abs(actual - expected)) < 3e-14
        drift += float(law[:, 0] @ law[:, 5]) / exposure / 2
    force = kappa * (table[0, y] - table[1, y])
    assert abs(drift - nu * math.sinh(beta*force/2)/2) < 3e-14
    assert abs(drift - row['measured_drift_per_source_time']) < 3e-14
    return drift


def rejected(fn):
    try:
        fn()
    except (AssertionError, ValueError):
        return True
    return False


def run():
    result = json.loads((HERE/'RESULT.json').read_text())
    assert producer.dependencies() == result['source_pins']
    assert producer.sha(HERE/'force_measurement.py') == result['producer_sha256']
    assert producer.sha(HERE/'FIRST_EVENT_OUTCOMES.npz') == result['outcome_archive']['sha256']
    fixed = producer.sys.modules['field_feedback']
    reciprocal = producer.load('reciprocal_feedback', producer.BRIDGE/'reciprocal_feedback.py')
    edges, neighbours = fixed.graph()
    table = reciprocal.pair_table(neighbours)
    archive = np.load(HERE/'FIRST_EVENT_OUTCOMES.npz')
    for row in result['measurements']:
        validate(archive[row['archive_key']], row, edges, table)
    primary = result['measurements'][3]
    original = archive['kappa_1']
    mutations = {}
    cap_removed = original.copy()
    cap_removed[:, 1:, 0] = 0
    mutations['cap_output_deleted'] = rejected(lambda: validate(cap_removed, primary, edges, table))
    holds_removed = original.copy()
    for law in holds_removed:
        law[(law[:, 8] == 1) & (law[:, 4] == 0), 0] = 0
    mutations['native_holds_deleted'] = rejected(lambda: validate(holds_removed, primary, edges, table))
    selected = cap_removed.copy()
    selected[:, 0, 0] = 1
    mutations['conclusive_only_renormalization'] = rejected(lambda: validate(selected, primary, edges, table))
    reversed_current = original.copy()
    reversed_current[:, :, 5] *= -1
    mutations['current_sign_reversed'] = rejected(lambda: validate(reversed_current, primary, edges, table))
    doubled_current = original.copy()
    doubled_current[:, :, 5] *= 2
    mutations['fair_coin_factor_omitted'] = rejected(lambda: validate(doubled_current, primary, edges, table))
    wrong_clock = copy.deepcopy(primary)
    wrong_clock['nu'] = 2
    mutations['native_clock_rate_mislabeled'] = rejected(lambda: validate(original, wrong_clock, edges, table))
    mutations['held_field_claimed_reciprocal'] = rejected(
        lambda: validate(archive['held_field'], result['held_field_control']['reciprocal'], edges, table))
    bare = result['field_to_energy']['present_oriented_seam_field']
    mutations['bare_field_omits_self_response'] = abs(bare - primary['force_minus_energy_difference']) > .1
    zero = result['measurements'][0]
    mutations['nonzero_Green_field_claimed_force_at_zero_coupling'] = (
        abs(bare) > .1 and abs(zero['measured_drift_per_source_time']) < 1e-14)
    reflected = result['reflection']['measurement']
    mutations['reflection_with_unchanged_drift'] = abs(reflected['measured_drift_per_source_time'] - primary['measured_drift_per_source_time']) > .06
    assert all(mutations.values())
    # The source envelope Lambda<=660 makes the score bound rational:
    # g=integral exp(-Lambda*t)dt >= T-Lambda*T²/2, at T=1/1000.
    time = F(1, 1000)
    exposure_lower = time - 660*time*time/2
    tv_upper = F(result['fixedpoint_variant_transfer']['complete_marked_capped_path_TV_upper'])
    score_upper = 2*tv_upper/exposure_lower
    assert result['fixedpoint_variant_transfer']['measured_frozen_score_error'] <= float(score_upper)
    review = {'schema': 'oph.native_force_measurement.adversarial_review.v1',
              'result_sha256': producer.sha(HERE/'RESULT.json'),
              'checks': {'all_pins': True, 'independent_probability_and_endpoint_reconstruction': True,
                         'all_five_couplings': True, 'full_cap_outcome_archive': True},
              'mutation_rejections': mutations,
              'rational_force_transfer': {'exposure_lower': str(exposure_lower),
                                           'TV_upper': str(tv_upper),
                                           'frozen_score_error_upper': str(score_upper)},
              'scope': 'Separate implementation by the same author; no external scientific review claimed.'}
    (HERE/'REVIEW.json').write_text(json.dumps(review, indent=2)+'\n')
    print(json.dumps(review, indent=2))


if __name__ == '__main__':
    run()
