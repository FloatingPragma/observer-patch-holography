"""Exhaustive sector proof and new complete-tape native executions."""
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import random

import neutral_record as nr

HERE = Path(__file__).resolve().parent


def coin_stream(rng):
    while True:
        yield rng.randrange(2)


def independent_step(n, edge, c):
    a, b = edge
    total = n[a]+n[b]
    na = total//2 + c*(total % 2)
    out = list(n)
    out[a], out[b] = na, total-na
    return tuple(out), n[a]-na


def main():
    branches = 0
    for q in (1, 7, 100000):
        configurations = [tuple([q]*12)]
        for a in range(12):
            for b in range(12):
                if a != b:
                    n = [q]*12
                    n[a] += 1
                    n[b] -= 1
                    configurations.append(tuple(n))
        assert len(configurations) == 133
        for n in configurations:
            for edge in nr.EDGES:
                for c in (0, 1):
                    out, j = independent_step(n, edge, c)
                    actual, (_, actualj) = nr.original.native(n, edge, c)
                    assert (out, j) == (actual, actualj)
                    assert nr.pair_sector(out, q) and abs(j) <= 1
                    shifted, js = independent_step(tuple(x-q for x in n), edge, c)
                    assert tuple(x+q for x in shifted) == out and js == j
                    branches += 1
    rng = random.Random(273511)  # diagnostic source, not a reserved seed
    executed = []
    controls = {'bad_contrast_rejected': True, 'missing_preparation_rejected': True,
                'fractional_current_rejected': True, 'hidden_source_ignored': True,
                'hidden_pointer_ignored': True, 'wrong_pointer_rejected': True}
    wrong_pointer_exercised = 0
    observed = set()
    outcomes = set()
    # Each of these source histories carries every possible regional value.
    # Both final source states and the preparation are used only in verification.
    for case in range(8):
        q = 1 + case % 3
        source, audit = nr.execute_source(q, case % 30, (-1)**case,
                         [(rng.randrange(30), rng.randrange(2)) for _ in range(case)],
                         [rng.randrange(2) for _ in range(30)])
        for region in ((case % 12,), tuple(v for v in range(12) if v != case % 12), tuple(range(12))):
            obs = nr.observe_source(q, source, region)
            charge = sum(audit['final_loads'][v]-q for v in region)
            assert obs['charge'] == charge
            for initial in range(12):
                pointer = nr.original.meter.execute(initial, (obs['degree'],), coin_stream(rng), 2)
                verdict = nr.admit(q, source, region, pointer)
                assert verdict['admitted'] and initial in verdict['compatible_pointer_initials']
                assert verdict['readout'] == (charge if pointer['success'] else None)
                assert Fraction(verdict['fine_record_effect'][initial]) == Fraction(1, 2**verdict['native_attempts'])
                outcomes.add(pointer['reason'])
                if pointer['success']:
                    observed.add(charge)
                bad = deepcopy(source)
                bad['records'][-1]['contrast'] += 1
                controls['bad_contrast_rejected'] &= not nr.admit(q, bad, region, pointer)['admitted']
                bad = deepcopy(source)
                bad['records'].pop(0)
                controls['missing_preparation_rejected'] &= not nr.admit(q, bad, region, pointer)['admitted']
                bad = deepcopy(source)
                bad['records'][0]['current'] = .5
                controls['fractional_current_rejected'] &= not nr.admit(q, bad, region, pointer)['admitted']
                hidden = deepcopy(source)
                hidden['final_loads'] = ['invalid diagnostic data']
                controls['hidden_source_ignored'] &= nr.admit(q, hidden, region, pointer) == verdict
                hiddenp = deepcopy(pointer)
                hiddenp['final_loads'] = ['invalid diagnostic data']
                controls['hidden_pointer_ignored'] &= nr.admit(q, source, region, hiddenp) == verdict
                wrong = nr.original.meter.execute(initial, (1 if obs['degree'] != 1 else -1,), coin_stream(rng), 2)
                if wrong['success']:
                    wrong_pointer_exercised += 1
                    controls['wrong_pointer_rejected'] &= not nr.admit(q, source, region, wrong)['admitted']
                executed.append({'case': case, 'filling': q, 'region': region,
                                 'source_records': source, 'pointer': pointer,
                                 'audit_only': audit, 'initial_pointer_for_audit_only': initial,
                                 'verdict': verdict})
    assert all(controls.values()), controls
    assert wrong_pointer_exercised > 0
    assert {-1, 0, 1} <= observed, observed
    assert len(outcomes) >= 2, outcomes
    files = [Path(nr.__file__), Path(__file__), nr.PROVIDER/'native_kernel.py',
             nr.PROVIDER/'native_cycle_meter.py', nr.PROVIDER/'strict_tape_validator.py']
    pins = {str(p.relative_to(nr.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    result = {'schema': 'oph.obstruction_repair.neutral_record.v1',
              'candidate_id': 'native-neutral-pair-record-only-v1',
              'new_preparation_domain': True, 'new_source_update_law': False,
              'exhaustive_native_branches': branches, 'fresh_complete_tapes': len(executed),
              'charge_values': sorted(observed), 'pointer_outcomes': sorted(outcomes),
              'controls': controls, 'source_pins': pins, 'official_rows_changed': [],
              'wrong_pointer_control_executions': wrong_pointer_exercised,
              'solved_component': 'neutral-start single-pair sector joined to record-only finite counters and native pointer',
              'retained_inputs': ['true scalar neutral filling', 'recorded pair intervention',
                                  'integer native alphabet', 'faithful emitted E channel',
                                  'declared source/pointer scheduler and spatial charge dictionary'],
              'not_solved': ['electric charge identification', 'independent physical coupling',
                             'spontaneous pair creation', 'full common-source physics', 'platform validation']}
    (HERE/'NEUTRAL_RECORD_EXAMPLES.json').write_text(json.dumps(executed, indent=2)+'\n')
    (HERE/'NEUTRAL_RECORD_RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('exhaustive_native_branches', 'fresh_complete_tapes', 'charge_values', 'pointer_outcomes', 'controls')}))


if __name__ == '__main__':
    main()
