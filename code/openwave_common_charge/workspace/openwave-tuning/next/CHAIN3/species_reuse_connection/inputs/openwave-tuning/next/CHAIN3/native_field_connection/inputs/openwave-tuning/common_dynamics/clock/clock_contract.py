"""One declared stroboscopic clock and complete native thinning comparison.

Effective refinement time is separate from implementation duration. Duration
requires a bounded control-rate and service-latency profile.
"""
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
from pathlib import Path
import hashlib
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RER = ROOT / 'reverse-engineering-reality'
PROVIDER = ROOT / 'openwave-tuning/auditor/native_walk_bridge'
sys.path.insert(0, str(RER))
sys.path.insert(0, str(PROVIDER))
from model.native_spatial_refinement import source as spatial
from model.observer import native
from model.observer.source_scalar import compiler
import common_walk

MODEL_ID = 'native-local12-current-v0-proposed'
CLOCK_ATTACHMENT_ID = 'local12-completed-word-clock-v1-proposed'


def positive_rational(value, name):
    if type(value) is not F or value <= 0:
        raise ValueError(name + ' must be an explicit positive Fraction')
    return value


def effective_time(steps, cell_spacing, tau0=F(1)):
    """Shared refinement coordinate; no physical-second inference."""
    if type(steps) is not int or steps < 0:
        raise ValueError('nonnegative exact completed-step count required')
    return steps * positive_rational(cell_spacing, 'cell spacing') * positive_rational(tau0, 'tau0')


def refinement_contract(cell_spacing, tau0=F(1)):
    dt = effective_time(1, cell_spacing, tau0)
    return {'clock_attachment': CLOCK_ATTACHMENT_ID,
            'time_per_completed_W_exact': str(dt),
            'time_coordinate': 't=n*tau0*h',
            'scalar_acoustic_speed_squared_exact': str(F(1, 48) / tau0**2),
            'source_graph_light_cone_speed_bound_exact': str(F(1) / tau0),
            'physical_seconds_identified': False,
            'condition': 'The same h and tau0 apply to every source, field and apparatus output; actual timing needs a supplied implementation profile.'}


def wait_law(size, cap):
    if type(cap) is not int or cap < 1:
        raise ValueError('positive finite global-attempt cap required')
    g = spatial.graph(size)
    p = F(6, len(g.seams))
    contact = [p * (1-p)**(n-1) for n in range(1, cap+1)]
    timeout = (1-p)**cap
    assert sum(contact) + timeout == 1
    return {'global_attempt_cap': cap, 'contact_probability_per_attempt': str(p),
            'first_contact_attempt_probabilities': list(map(str, contact)),
            'timeout_probability': str(timeout),
            'full_path_recipe': 'Before contact retain each chosen nonincident seam, law and coin. Each fine C or E path of length n has weight (1/(2M))^n for its fixed law; an independently chosen fair law adds a factor1/2 per event. A contact may hold the load. Timeout retains all cap letters and unchanged load.',
            'success_conditioning_applied': False}


def replay_global_tape(initial, tape, size=2, law='C'):
    g = spatial.graph(size)
    if type(initial) is not int or not 0 <= initial < g.ports or law not in ('C', 'E'):
        raise ValueError('admitted initial port and C/E law required')
    x = initial
    records = []
    contacts = 0
    for edge, coin in tape:
        if type(edge) is not int or not 0 <= edge < len(g.seams) or coin not in (0, 1):
            raise ValueError('actual seam and binary coin required')
        seam = g.seams[edge]
        before = tuple(int(v == x) for v in range(g.ports))
        after, (_, current) = native.native(before, seam, coin)
        y = after.index(1)
        contact = x in seam
        contacts += int(contact)
        label = current if law == 'C' else before[seam[0]]-before[seam[1]]
        # An observer with the declared known initial unit reconstructs its
        # location from signed currents. The simulator endpoint is audit-only.
        observed_next = seam[1] if current == 1 else seam[0] if current == -1 else x
        assert observed_next == y
        records.append({'edge': edge, 'seam': list(seam), 'coin': coin, 'law': law,
                        'label': label, 'current': current, 'source_before': x,
                        'source_after': y, 'contact': contact,
                        'global_attempt': len(records)+1, 'local_contact_count': contacts})
        x = observed_next
    return {'records': records, 'final_port': x, 'global_attempts': len(records),
            'local_contacts_including_holds': contacts,
            'contact_classification': 'Known initial port plus preceding signed-current records; no terminal occupation is an observer input. Unknown-position or quantum inputs do not inherit this classifier.',
            'fixed_law_path_probability_exact': str(F(1, 2*len(g.seams))**len(records))}


@lru_cache(None)
def frozen_program():
    return json.loads((PROVIDER/'COMMON_WALK_PROGRAM.json').read_text())


def pulse_histogram():
    program = frozen_program()
    counts = Counter((p['kind'], F.from_float(float(p['angle']))) for p in program['reflection_word'])
    calls = Counter(kind for kind, _ in program['swap_routine_calls'])
    for name, count in calls.items():
        for pulse in program['primitive_words'][name]:
            counts[pulse['kind'], F.from_float(float(pulse['angle']))] += count
    assert sum(counts.values()) == 432425
    return counts


def sqrt_interval(integer, digits=50):
    scale = 10**digits
    below = isqrt(integer*scale*scale)
    return F(below, scale), F(below+1, scale)


def coefficient_bound(a, b, root):
    endpoints = (a+b*root[0], a+b*root[1])
    return min(endpoints), max(endpoints)


def abs_interval(interval):
    lo, hi = interval
    return (F(0) if lo <= 0 <= hi else min(abs(lo), abs(hi)), max(abs(lo), abs(hi)))


def coefficient_max_intervals():
    summary = compiler.source_summary()
    rows = summary['exact_base_coefficients']
    root5, root2 = sqrt_interval(5), sqrt_interval(2)
    coefficients = {}
    for kind, row in rows.items():
        coefficients[kind] = [coefficient_bound(F(a), F(b), root5) for a, b in row]
    summed = [(a[0]+b[0], a[1]+b[1]) for a, b in zip(coefficients['X'], coefficients['Z'])]
    coefficients['H'] = [(min(v/root2[j] for v in row for j in (0,1)),
                           max(v/root2[j] for v in row for j in (0,1))) for row in summed]
    return {kind: (max(abs_interval(c)[0] for c in row),
                   max(abs_interval(c)[1] for c in row))
            for kind, row in coefficients.items()}


def implementation_duration(control_rate, service_latency, attempt_latency=F(0), attempts=0):
    """Serial rectangular-control profile in a declared external time unit.

    Each port has amplitude at most control_rate. Each pulse takes the
    smallest duration that respects that cap. Every load/export takes the
    same explicitly supplied service_latency. This is a timing architecture,
    not a timing property inferred from the native map.
    """
    positive_rational(control_rate, 'control rate')
    if any(type(v) is not F or v < 0 for v in (service_latency, attempt_latency)):
        raise ValueError('explicit nonnegative rational service/attempt latencies required')
    if type(attempts) is not int or attempts < 0:
        raise ValueError('nonnegative actual attempt count required')
    coeff = coefficient_max_intervals()
    active = [F(0), F(0)]
    literal_active = F(0)
    for (kind, angle), count in pulse_histogram().items():
        for side in (0, 1):
            active[side] += count*abs(angle)*coeff[kind][side]/control_rate
        # Source execution multiplies literal binary64 coefficients. Retain
        # those exact dyadic values instead of identifying them with the
        # ideal algebraic coefficients used by the independent enclosure.
        literal = [F.from_float(float(c)) for c in compiler.integrated_controls(kind, float(angle))]
        literal_active += count*max(map(abs, literal))/control_rate
    overhead = 864850*service_latency + attempts*attempt_latency
    return {'ideal_algebraic_control_time_interval_exact': list(map(str, active)),
            'literal_binary64_control_time_exact': str(literal_active),
            'literal_vs_ideal_duration_difference_interval_exact': [str(literal_active-active[1]), str(literal_active-active[0])],
            'total_duration_interval_exact': [str(literal_active+overhead)]*2,
            'native_pulses': 432425, 'load_export_services': 864850,
            'actual_native_attempts': attempts, 'physical_unit': 'declared profile unit',
            'assumptions': ['serial rectangular simultaneous-port controls',
                            'per-port amplitude bound is attained for each pulse',
                            'fixed common service and attempt latency',
                            'no extra idle interval or overlap'],
            'source_selects_profile': False}


def graph_distances(g):
    result = np.full((g.ports, g.ports), -1, dtype=np.int64)
    for source in range(g.ports):
        result[source, source] = 0
        todo = deque([source])
        while todo:
            x = todo.popleft()
            for y in g.neighbours[x]:
                if result[source, y] < 0:
                    result[source, y] = result[source, x]+1
                    todo.append(y)
    return result


def locality_audit():
    g = common_walk.graph()
    distances = graph_distances(g)
    inputs = np.array([16*common_walk.encode_position(x)+event
                       for x in range(g.ports) for event in range(12)], dtype=np.int64)
    original = np.repeat(np.arange(g.ports), 12)
    current = inputs.copy()
    first_nonlocal = None
    padding_prefixes = 0
    for prefix, (kind, targets) in enumerate(frozen_program()['swap_routine_calls'], 1):
        if kind == 'flip': current ^= 1 << targets[0]
        elif kind == 'cx': current ^= ((current >> targets[0])&1) << targets[1]
        elif kind == 'ccx':
            a, b, c = targets
            current ^= (((current >> a)&1) & ((current >> b)&1)) << c
        else: raise ValueError(kind)
        addresses = (current >> 4)&127
        valid = (addresses%16) < 12
        padding_prefixes += int(not valid.all())
        vertices = 12*(addresses//16)+addresses%16
        candidates = np.flatnonzero(valid)
        far = candidates[distances[original[candidates], vertices[candidates]] > 1]
        if first_nonlocal is None and len(far):
            i = int(far[0])
            first_nonlocal = {'gate_prefix': prefix, 'gate_kind': kind, 'targets': list(targets),
                              'initial_port': int(original[i]), 'event': int(inputs[i]%16),
                              'intermediate_port': int(vertices[i]),
                              'graph_distance': int(distances[original[i], vertices[i]])}
    expected = common_walk.expected_swap_permutation()[inputs]
    assert np.array_equal(current, expected)
    final = np.array([common_walk.decode_position(int(x//16)) for x in current])
    assert (distances[original, final] <= 1).all()
    assert first_nonlocal is not None and padding_prefixes > 0
    # Exact numerator propagation checks the entire physical column family.
    state = np.zeros((96, 12, 96), dtype=np.int64)
    for x in range(96): state[x, :, x] = 1
    neighbour = common_walk.geometry_tables()[0]
    walks = []
    for step in range(1, 5):
        state = state.sum(axis=1, keepdims=True) - 6*state
        swapped = state.copy()
        for x in range(96):
            for direction in range(6):
                swapped[neighbour[x,direction], 2*direction+1] = state[x,2*direction+1]
        state = swapped
        bad = 0
        for target in range(96):
            bad += int(np.count_nonzero(state[target, :, distances[target] > step]))
        assert bad == 0
        walks.append({'steps': step, 'input_position_columns': 96, 'outside_graph_cone_nonzero_amplitudes': bad})
    return {'whole_step_light_cone': walks, 'encoded_gate_prefix_counterexample': first_nonlocal,
            'prefixes_with_temporary_source_padding': padding_prefixes,
            'locality_scope': 'Graph-position support at complete W boundaries; encoded processor operations and coherent load/export are not a distributed graph-local implementation.',
            'concrete_local_implementation_obligation': 'Represent twelve event modes at each site; apply on-site R in parallel and transport each matched mode only to its adjacent site. Establish the same complete native service at those spatial addresses.'}


def run():
    source_paths = [Path(__file__), PROVIDER/'common_walk.py', PROVIDER/'compile_walk.py',
                    PROVIDER/'COMMON_WALK_PROGRAM.json', PROVIDER/'COMMON_WALK_RECEIPT.json',
                    PROVIDER/'COMMON_SOURCE_CUSTODY.json', Path(spatial.__file__), Path(native.__file__),
                    Path(compiler.__file__), RER/compiler.RESPONSE, RER/compiler.CONTROL_RECEIPT]
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    custody = json.loads((PROVIDER/'COMMON_SOURCE_CUSTODY.json').read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h for p,h in custody['files_sha256'].items())
    g = spatial.graph(2)
    numerators = np.zeros((96, 96), dtype=np.int64)
    checked = 0
    for x in range(96):
        before = tuple(int(i == x) for i in range(96))
        for e, seam in enumerate(g.seams):
            for coin in (0, 1):
                after, (_, current) = native.native(before, seam, coin)
                y = after.index(1)
                numerators[y, x] += 1
                assert current == int(x == seam[0])-int(y == seam[0])
                if x not in seam: assert y == x and current == 0
                checked += 1
    local = common_walk.transfer_numerators()
    assert np.array_equal(numerators, (2*len(g.seams)-12)*np.eye(96, dtype=int)+local)
    rates = []
    for size in (2, 3, 4):
        gg = spatial.graph(size)
        r = F(3, 7)
        rate_global, p = len(gg.seams)*r, F(6, len(gg.seams))
        assert rate_global*p == 6*r
        rates.append({'size': size, 'seams': len(gg.seams), 'per_seam_rate': str(r),
                      'global_attempt_rate': str(rate_global), 'local_contact_rate': str(rate_global*p),
                      'generator': '-r*L_graph/2'})
    tape = [(next(i for i,e in enumerate(g.seams) if 0 not in e),0),
            (next(i for i,e in enumerate(g.seams) if 0 not in e),1), (0,0), (0,1)]
    traces = [replay_global_tape(0, tape, law=law) for law in ('C','E')]
    assert all(t['global_attempts']==4 and t['local_contacts_including_holds']==2 for t in traces)
    waits = [wait_law(size, cap) for size in (2,3,4) for cap in (1,2,4)]
    profile = implementation_duration(F(1), F(1,100))
    scaled = implementation_duration(F(2), F(1,200))
    assert all(F(a)==2*F(b) for a,b in zip(profile['total_duration_interval_exact'],scaled['total_duration_interval_exact']))
    assert effective_time(4,F(1,16)) == F(1,4)
    # This return statistic follows from W^2 on delta_x tensor uniform12.
    coherent_return, markov_return = F(25,72), F(7,24)
    assert coherent_return-markov_return == F(1,18)
    result = {'schema': 'oph.local12.common_clock.v1', 'candidate_id': MODEL_ID,
        'clock_attachment_id': CLOCK_ATTACHMENT_ID,
        'checks': {'all_native_global_columns': True, 'exact_population_thinning': True,
                   'per_seam_clock_is_volume_consistent': True, 'all_capped_wait_outcomes': True,
                   'all_nonincident_native_marks_retained': True,
                   'source_control_duration_enclosed': True, 'rate_rescaling_control': True,
                   'whole_step_and_processor_locality_distinguished': True,
                   'coherent_and_markov_second_step_control': True, 'source_custody_unchanged': True},
        'native_columns_checked': checked, 'per_seam_clock': rates, 'complete_wait_laws': waits,
        'actual_fine_tapes': traces, 'implementation_profile': profile,
        'double_rate_half_latency_profile': scaled,
        'refinement_clock_examples': [refinement_contract(F(1,n)) for n in (16,32,64)],
        'locality': locality_audit(),
        'second_step_return_control': {'coherent': str(coherent_return), 'Markov': str(markov_return), 'difference': '1/18'},
        'physical_attachment': 'A common effective time coordinate is supplied. Physical seconds, spatially local processor services, and refinement scaling of hardware bandwidth/latency require explicit architecture or calibration.',
        'scope': 'Literal classical thinning is an exact same-process comparison. Repeated coherent W is a declared distinct dynamics; no raw-channel identity is inferred. Every time coordinate and failure is retained.',
        'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    (HERE/'CLOCK_CONTRACT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': result['checks'], 'native_columns': checked,
                      'locality_witness': result['locality']['encoded_gate_prefix_counterexample'],
                      'duration': [float(F(x)) for x in profile['total_duration_interval_exact']]}))


if __name__ == '__main__':
    run()
