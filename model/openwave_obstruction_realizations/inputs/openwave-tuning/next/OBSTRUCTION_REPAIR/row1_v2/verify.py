"""Independent-field mutants, unknown-start sweeps and exact instrument laws."""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import cmath
import hashlib
import json
import math
import random
import sys

sys.dont_write_bytecode = True
import phase_reader as phase
import closing_reader as closing

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
GEOMETRY = ROOT/'openwave-tuning/row_completion/ROW_01/GEOMETRY.json'


def exact_degree(geom, currents, region):
    # Incidence of the primal seams, independent of dual contour assembly.
    total = 0
    for (a, b), j in zip(geom['seams'], currents):
        total += j*((a in region)-(b in region))
    return total


def phase_controls():
    geom = json.loads(GEOMETRY.read_text())
    rng = random.Random(411719)  # new diagnostics, never held-out acceptance seeds
    counts = defaultdict(int)
    max_error = 0.0
    for case in range(72):
        amplitude = (0, 1, 3, 4, 16, 64)[case % 6]
        recorded = [rng.randint(-amplitude, amplitude) for _ in range(30)]
        region = [v for v in range(12) if rng.randrange(2)]
        if len(region) in (0, 12):
            region = [case % 12]
        cut_edge = next(e for e, (a, b) in enumerate(geom['seams']) if (a in region) != (b in region))
        expected = exact_degree(geom, recorded, region)
        for deformed in (False, True):
            got = phase.read_phase(geom, recorded, region, field=phase.declared_field(geom, recorded, deformed=deformed))
            assert got['status'] == 'read' and got['phase_degree'] == expected, got
            assert got['max_principal_increment'] < math.pi/4+1e-12
            assert not got['regularity_verified_from_samples']
            max_error = max(max_error, abs(got['float_degree']-expected))
            counts['declared_reads'] += 1
        # The source field can disagree with the record which fixed its sampling.
        reversed_field = [-j for j in recorded]
        got = phase.read_phase(geom, recorded, region, field=phase.declared_field(geom, reversed_field))
        assert got['status'] == 'read' and got['phase_degree'] == -expected
        if expected:
            assert not phase.compare_degree(got, expected)
            counts['reversed_nonzero_comparison_rejections'] += 1
        counts['reversed_reads'] += 1
        changed = recorded.copy()
        changed[cut_edge] += 1
        got = phase.read_phase(geom, recorded, region, field=phase.declared_field(geom, changed))
        assert got['status'] == 'read' and got['phase_degree'] == exact_degree(geom, changed, region)
        assert not phase.compare_degree(got, expected)
        counts['extra_turn_comparison_rejections'] += 1
        changed[cut_edge] = 3*abs(recorded[cut_edge])+3
        got = phase.read_phase(geom, recorded, region, field=phase.declared_field(geom, changed))
        assert got['status'] == 'increment_guard' and got['phase_degree'] is None, got
        counts['large_increment_rejections'] += 1
        changed = recorded.copy()
        changed[cut_edge] += 8*abs(recorded[cut_edge])+8
        got = phase.read_phase(geom, recorded, region, field=phase.declared_field(geom, changed))
        assert got['status'] == 'read' and got['phase_degree'] == expected
        assert got['phase_degree'] != exact_degree(geom, changed, region)
        # This is an EXPECTED failure of arbitrary-field authentication.
        counts['out_of_class_aliases_retained'] += 1
    zero = [0]*30
    basefield = phase.declared_field(geom, zero)
    badfield = lambda e, t: basefield(e, t)*cmath.exp(.17j*t) if e == 0 else basefield(e, t)
    assert phase.read_phase(geom, zero, [0], field=badfield)['status'] == 'open_curve'
    assert phase.read_phase(geom, zero, [0], field=lambda e, t: 0)['status'] == 'invalid_field'
    badgeom = deepcopy(geom)
    badgeom['dual_edges'][0][1] = -1
    assert phase.read_phase(badgeom, zero, [0], field=basefield)['status'] == 'rejected_geometry'
    for bad in (.5, F(1, 2), 4.0, True, '1', None):
        record = zero.copy()
        record[0] = bad
        assert phase.read_phase(geom, record, [0], field=basefield)['status'] == 'rejected_input'
    calls = []
    unused = lambda e, t: calls.append((e, t))
    huge = [10**100]+[0]*29
    assert phase.read_phase(geom, huge, [0], field=unused)['status'] == 'resource_stop' and not calls
    assert phase.read_phase(geom, zero, [0], field=unused, sample_budget=0)['status'] == 'resource_stop' and not calls
    assert phase.read_phase(geom, zero, [0], field=basefield, sample_budget=100001)['status'] == 'rejected_input'
    for m, status in ((12494, 'read'), (12495, 'resource_stop')):
        record = [m]+[0]*29
        got = phase.read_phase(geom, record, [0], field=phase.declared_field(geom, record))
        assert got['status'] == status, got
    return {'counts': dict(counts), 'max_float_degree_error': max_error,
            'closure_geometry_unit_circle_and_input_controls': True,
            'registered_cap_controls': True, 'arbitrary_field_authentication': False}


def step(n, seam, coin):
    """Own endpoint averaging rule; no producer native implementation."""
    a, b = seam
    result = list(n)
    total = n[a]+n[b]
    result[a] = total//2+coin*(total % 2)
    result[b] = total-result[a]
    return tuple(result), n[a]-result[a]


def configurations():
    result = [(0,)*12]
    for a in range(12):
        for b in range(12):
            if a != b:
                n = [0]*12
                n[a], n[b] = 1, -1
                result.append(tuple(n))
    return result


def make_sweep(initial, order, coins):
    state = initial
    records = []
    for e, coin in zip(order, coins):
        a, b = closing.EDGES[e]
        d = state[a]-state[b]
        state, j = step(state, (a, b), coin)
        records.append({'kind': 'E', 'seam': e, 'coin': coin, 'current': j, 'contrast': d})
    return records, state


def coin_stream(rng):
    while True:
        yield rng.randrange(2)


def closing_controls():
    rng = random.Random(411720)
    counts = defaultdict(int)
    outcomes = set()
    successful_values = set()
    examples = []
    for number, initial in enumerate(configurations()):
        for repetition in range(2):
            order = list(range(30))
            rng.shuffle(order)
            coins = [rng.randrange(2) for _ in order]
            records, final = make_sweep(initial, order, coins)
            shifted_records, shifted_final = make_sweep(tuple(v+100000 for v in initial), order, coins)
            assert records == shifted_records and tuple(v-100000 for v in shifted_final) == final
            counts['arbitrary_start_sweeps'] += 1
            for region in ([number % 12], [v for v in range(12) if v != number % 12], list(range(12))):
                observation = closing.observe_closing(records, region)
                assert observation['initial_relative_state'] == list(initial)
                assert observation['final_relative_state'] == list(final)
                charge = sum(final[v] for v in region)
                assert observation['charge'] == charge
                counts['closing_only_reads'] += 1
                for initial_pointer in (0, 5, 11):
                    pointer = closing.original.meter.execute(initial_pointer, (-charge,), coin_stream(rng), 2)
                    verdict = closing.admit_closing(records, region, pointer)
                    assert verdict['admitted'] and initial_pointer in verdict['compatible_pointer_initials'], verdict
                    assert verdict['source_initial_relative_column'] == list(initial)
                    assert F(verdict['fine_pointer_effect'][initial_pointer]) == F(1, 2**verdict['native_attempts'])
                    assert verdict['readout'] == (charge if pointer['success'] else None)
                    outcomes.add(verdict['reason'])
                    if pointer['success']:
                        successful_values.add(charge)
                    hidden = deepcopy(pointer)
                    hidden['final_loads'] = ['not an observer input']
                    assert closing.admit_closing(records, region, hidden) == verdict
                    counts['fine_effect_reads'] += 1
                    if len(examples) < 18 or (pointer['success'] and charge and len(examples) < 24):
                        examples.append({'records': records, 'region': region, 'pointer': pointer, 'verdict': verdict})
                # Corruption is tested at the actual observer interface.
                bad = deepcopy(records)
                bad[-1]['contrast'] += 1
                assert not closing.admit_closing(bad, region, pointer)['admitted']
                bad = deepcopy(records)
                bad[0]['current'] = F(1, 2)
                assert not closing.admit_closing(bad, region, pointer)['admitted']
                bad = deepcopy(records)
                bad[1] = deepcopy(bad[0])
                assert not closing.admit_closing(bad, region, pointer)['admitted']
                assert not closing.admit_closing(records[:-1], region, pointer)['admitted']
                counts['four_corruption_controls'] += 1
    assert successful_values == {-1, 0, 1}, successful_values
    assert len(outcomes) >= 3, outcomes
    # Two genuine pairs and a lone excess have valid native sweeps, but their
    # observed relative states are outside the declared instrument's domain.
    for initial in ((1, 1, -1, -1)+(0,)*8, (1,)+(0,)*11):
        records, _ = make_sweep(initial, list(range(30)), [0]*30)
        try:
            closing.observe_closing(records, [0])
        except ValueError:
            pass
        else:
            raise AssertionError('outside-sector source admitted')
    # No claim of detecting a transfer after the final recorded event.
    records, final = make_sweep((0,)*12, list(range(30)), [0]*30)
    result = closing.observe_closing(records, [0])
    altered_endpoint = list(final)
    altered_endpoint[0] += 1
    altered_endpoint[1] -= 1
    assert result['charge'] != altered_endpoint[0]
    (HERE/'EXAMPLES.json').write_text(json.dumps(examples, indent=2)+'\n')
    return {'counts': dict(counts), 'successfully_read_charges': sorted(successful_values),
            'retained_pointer_outcomes': sorted(outcomes), 'outside_sector_controls': True,
            'post_record_motion_detected': False, 'preparation_history_is_observer_input': False}


def exact_laws():
    states = configurations()
    state_set = set(states)
    branches = 0
    for state in states:
        for seam in closing.EDGES:
            for coin in (0, 1):
                output, _ = step(state, seam, coin)
                assert output in state_set
                branches += 1
    # Retain both coin outcomes even when they produce the same state.
    all_sweeps = []
    order = list(range(7, 30))+list(range(7))
    for initial in states:
        distribution = {initial: 1}
        for depth, e in enumerate(order, 1):
            nxt = defaultdict(int)
            for n, multiplicity in distribution.items():
                for c in (0, 1):
                    n2, _ = step(n, closing.EDGES[e], c)
                    nxt[n2] += multiplicity
            distribution = dict(nxt)
            assert sum(distribution.values()) == 2**depth
        all_sweeps.append({'initial_relative_state': initial,
                           'final_counts': [[list(n), count] for n, count in sorted(distribution.items())],
                           'denominator': 2**30})
    # Exact pointer aggregate law derived from discovery/route/measurement
    # stages, rather than replaying a favorable or postselected tape.
    pointer_laws = []
    h, r = F(1, 4), F(3, 4)
    for x in range(12):
        for degree in (-1, 0, 1):
            probability = {'discovery_timeout': h**5, 'routing_timeout': F(),
                           'measurement_timeout': F(), 'complete_closed_pointer_word': F()}
            k = 0
            for a, b in closing.EDGES:
                if x not in (a, b):
                    continue
                destination = b if a == x else a
                length = len(closing.original.meter.route(destination, 0))
                first_discovery = h**k*r
                probability['routing_timeout'] += first_discovery*(1-r**length)
                probability['measurement_timeout'] += first_discovery*r**length*(1-r**(6*abs(degree)))
                probability['complete_closed_pointer_word'] += first_discovery*r**(length+6*abs(degree))
                k += 1
            assert k == 5 and sum(probability.values()) == 1
            pointer_laws.append({'initial_pointer': x, 'degree': degree,
                                 'probabilities': {key: str(p) for key, p in probability.items()}})
    payload = {'native_branches': branches, 'complete_source_sweep_laws': all_sweeps,
               'pointer_laws': pointer_laws, 'all_joint_effects_sum_to_identity': True,
               'scope': '133 initial relative states times 12 pointer columns, conditional on sweep start'}
    (HERE/'EXACT_LAWS.json').write_text(json.dumps(payload, indent=2)+'\n')
    return {'exhaustive_native_branches': branches, 'complete_source_sweep_laws': len(all_sweeps),
            'source_coin_histories_per_initial_state': 2**30, 'pointer_laws': len(pointer_laws),
            'complete_joint_instrument': True}


def main():
    result = {'schema': 'oph.row1.interfaces.v2', 'phase': phase_controls(),
              'closing': closing_controls(), 'laws': exact_laws(),
              'official_rows_changed': [], 'physical_charge_identification_solved': False,
              'external_review': 'openwave-push/reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B/REVIEW.md'}
    (HERE/'RESULT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
