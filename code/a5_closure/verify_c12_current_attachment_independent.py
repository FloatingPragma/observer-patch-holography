"""Independent endpoint, spectral and partial-injection proof of C12 scope.

No producer, native executor, sampler or compiled-work module is imported.
The prior proper-code verifier supplies a separate Q(phi) arithmetic route.
"""
import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
from itertools import combinations, combinations_with_replacement
import json
import math
from pathlib import Path
import subprocess

import verify_proper_code_attachment_independent as prior
from verify_galois_port_frame_independent import ZERO, ONE, add, sub, mul, div, scalar, total, dot

ROOT = Path(__file__).resolve().parents[2]
BASELINE = '343ea42c64fb783c84c9fb4e461aa3c8ce054107'
# Current prerequisite receipt; historical work baseline remains above.
PRIOR_COMMIT = '2a96056fc2b3f3afc05bb4e16ba30b23d0fddfc5'
SOURCE_COMMIT = prior.SOURCE_COMMIT
INPUTS = sorted(set(prior.INPUTS + [
    'model/observer/native.py', 'model/observer/quantum_mixing.py',
    'model/observer/coherent_continuation.py',
    'model/identifiability_checks/rule_observability.py',
    'model/predictions/observer_instrument.py',
]))
PRIOR_FILES = ('code/a5_closure/proper_code_attachment_probe.py',
               'code/a5_closure/verify_proper_code_attachment_independent.py',
               'code/a5_closure/receipts/proper_code_attachment_probe.json')
IMPLEMENTATIONS = ('code/a5_closure/c12_current_attachment_probe.py',
                   'code/a5_closure/verify_c12_current_attachment_independent.py')
REFERENCE = ROOT/'code/a5_closure/receipts/c12_current_attachment_probe.json'
m = prior.m


def digest(value):
    return hashlib.sha256(json.dumps(prior.encoded(value), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def custody(paths, commit, repository, root):
    result = {}
    for path in paths:
        data = (root/path).read_bytes()
        m.need(data == subprocess.check_output(['git', 'show', commit+':'+path], cwd=ROOT), 'source custody: '+path)
        result[path] = {'repository': repository, 'commit': commit,
            'blob': hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(),
            'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
    return result


def configurations(full=False):
    pairs = combinations_with_replacement(range(12), 2) if full else combinations(range(12), 2)
    return sorted(tuple(int(p == a)+int(p == b) for p in range(12)) for a, b in pairs)


def endpoint(x, edge, coin):
    a, b = edge
    y = list(x)
    y[a] = (x[a]+x[b])//2+coin*((x[a]+x[b]) % 2)
    y[b] = x[a]+x[b]-y[a]
    return tuple(y), x[a]-y[a]


def validate_partial_injections(groups, dimension):
    counts = [0]*dimension
    for entries in groups.values():
        m.need(len({a for a, b in entries}) == len(entries) and
               len({b for a, b in entries}) == len(entries), 'native partial-injection condition')
        for a, b in entries:
            m.need(0 <= a < dimension and 0 <= b < dimension, 'native branch dimension')
            counts[a] += 1
    m.need(counts == [2]*dimension, 'two fair-coin branches per native column, including zero currents')


def native_profile(seams):
    states = configurations(full=True)
    lookup = {x: i for i, x in enumerate(states)}
    counts = Counter()
    for edge in seams:
        for law in ('C', 'E'):
            groups = defaultdict(list)
            for i, x in enumerate(states):
                for coin in (0, 1):
                    y, current = endpoint(x, edge, coin)
                    label = current if law == 'C' else x[edge[0]]-x[edge[1]]
                    groups[coin, label].append((i, lookup[y]))
            validate_partial_injections(groups, len(states))
            counts[law] += len(groups)
    return {'native_dimension': len(states), 'seams': len(seams),
        'matrix_units_checked': len(states)**2*len(seams),
        'C_recorded_branch_count': counts['C'], 'E_recorded_branch_count': counts['E'],
        'channel_identity': 'Phi=Phi after native-basis dephasing when terminal native coherence is inaccessible',
        'work_controls': 'arbitrary complete work-only CP instruments conditioned on classical native records',
        'additional_coherent_native_controls_admitted': False,
        'coherent_native_work_coupling_admitted': False,
        'coherent_uncommitted_native_record_access_admitted': False,
        'terminal_native_access': 'trace or diagonal read only',
        'finite_adaptive_histories_and_all_failure_flags': True}


def verify(packet, root=ROOT):
    pins = custody(INPUTS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography', root)
    previous = custody(PRIOR_FILES, PRIOR_COMMIT, 'FloatingPragma/observer-patch-holography', root)
    frozen = m.load(root/'code/a5_closure/manifests/marked_triplet_attachment_reference.json')
    frame = m.source_frame(root)
    vertices = m.response_vertices()
    mark = m.source_marked_action(root, vertices, [vertices.index(v) for v in frame])
    m.need(frozen['actions'].index(list(mark)) == 4, 'source marker fixed before characters')
    prior.verify(m.load(root/PRIOR_FILES[-1]), root=root)
    _, ds = m.incidence(root)
    seams = [(i, j) for i in range(12) for j in range(i+1, 12) if ds[i][j] == 1]
    states = configurations()
    lookup = {x: i for i, x in enumerate(states)}
    g = [[6*x[p]-1 for p in range(12)] for x in states]
    counts = [[0]*66 for _ in range(66)]
    for edge in seams:
        for i, x in enumerate(states):
            for coin in (0, 1):
                y, _ = endpoint(x, edge, coin)
                counts[lookup[y]][i] += 1
    gram = [[sum(r[i]*r[j] for r in g) for j in range(12)] for i in range(12)]
    m.need(gram == [[360*int(i == j)-30 for j in range(12)] for i in range(12)], 'G isometry on mean-free subspace')
    lap = [[5*int(i == j)-int(ds[i][j] == 1) for j in range(12)] for i in range(12)]
    for i in range(66):
        for p in range(12):
            lhs = sum(counts[i][j]*g[j][p] for j in range(66))
            rhs = 60*g[i][p]-sum(g[i][j]*lap[j][p] for j in range(12))
            m.need(lhs == rhs, 'native mean intertwining')
    join = {'native_binary_dimension': 66, 'work_site_dimension': 12, 'mean_zero_rank': 11,
        'G_formula': 'G[x,p]=6*x_p-1', 'G_star_G': '360 I12 - 30 J12',
        'scheduler': 'uniform thirty seams and fair placement coins; declared instance',
        'intertwining': 'K66 G = G (I12-L/60)', 'native_events_checked': 3960,
        'kernel_count_sha256': digest(counts), 'G_sha256': digest(g),
        'source_density_audit_replayed': True, 'state_encoding_operation_implemented': False}
    b = [[total(mul(scalar(x), frame[j][p]) for j, x in enumerate(row)) for p in range(3)] for row in g]
    for action in frozen['actions']:
        inverse = [action.index(i) for i in range(12)]
        for i, x in enumerate(states):
            y = tuple(x[inverse[p]] for p in range(12))
            m.need(g[lookup[y]] == [g[i][inverse[p]] for p in range(12)], 'native pair action transports the same G columns')
    scale = mul(scalar(360), (Q(8), Q(4)))  # 360*(10+2sqrt(5)).
    m.need(m.mm(prior.trans(b), b) == [[scale if i == j else ZERO for j in range(3)] for i in range(3)], 'source slow-code Gram')
    a, c = [r[0] for r in b], [r[1] for r in b]
    # For v_+=a+i b, v_-=a-i b: their inner product is
    # ||a||²-||b||²-2i<a,b>=0; each normalized outer product is pure.
    m.need(dot(a, a) == dot(c, c) == scale and dot(a, c) == ZERO, 'normalized orthogonal conjugate code states')
    diagonal = [div(add(mul(x, x), mul(y, y)), mul(scalar(2), scale)) for x, y in zip(a, c)]
    m.need(total(diagonal) == ONE, 'shared native diagonal')
    i, j = next((i, j) for i in range(66) for j in range(i+1, 66) if mul(c[i], a[j]) != mul(a[i], c[j]))
    coefficient = div(sub(mul(c[i], a[j]), mul(a[i], c[j])), mul(scalar(2), scale))
    control = {'current_marking_index': 4, 'historical_marking_index': 1,
        'current_native_action': list(mark), 'coordinate_choice': [0, 1],
        'native_code_Gram_scale': scale, 'actions_checked': 60,
        'witness_definition': 'rho_+-=(aa^T+bb^T +/- i(ba^T-ab^T))/(2k), [a,b]=first two columns of G raw_frame',
        'witness_conjugation': 'complex conjugation at fixed real embedding; not Galois conjugation',
        'vectors_sha256': digest([a, c]), 'shared_diagonal_sha256': digest(diagonal),
        'trace_each': 1, 'purity_each': 1, 'mutually_orthogonal': True,
        'nonzero_imaginary_entry': {'row': i, 'column': j, 'plus_coefficient': coefficient},
        'native_physical_preparation_of_witness_claimed': False, 'ideal_C12_images_orthogonal': True,
        'classical_record_only_outputs_equal': True,
        'complete_deterministic_transfer_worst_case_halftrace_floor': '1/2',
        'uniform_herald_success_p_floor': 'p/2; p>0 prescribed for the target channel',
        'arbitrary_untouched_references_covered': True,
        'coherent_native_output_retained_by_receiver': False,
        'work_input_independent_of_unknown_native_payload': True,
        'additional_coherent_native_work_access_admitted': False}
    # Spectral route to all filter traces; producer instead powers rational
    # matrices. The existing independent projector replay binds this spectrum.
    exponential = sum((Q(-1, 2)**j/Q(math.factorial(j)) for j in range(14)), Q())
    grid = 1 << 64
    weights = [Q((exponential*r*grid).__floor__(), grid) for r in (Q(1), Q(1, 4), Q(3, 32), Q(1, 64), Q(1, 512))]
    fs, ff = div((Q(114), Q(2)), scalar(120)), div((Q(116), Q(-2)), scalar(120))
    ps, pf = ONE, ONE
    traces, slow_mass = [], ZERO
    for n, w in enumerate(weights):
        t = add(scalar(1+5*Q(19, 20)**(2*n)), mul(scalar(3), add(ps, pf)))
        m.need(t[1] == 0, 'rational full trace')
        traces.append(t[0])
        slow_mass = add(slow_mass, mul(scalar(w), ps))
        ps, pf = mul(ps, mul(fs, fs)), mul(pf, mul(ff, ff))
    q = sum((w*t/12 for w, t in zip(weights, traces)), Q())
    successes = [w*t/12 for w, t in zip(weights, traces)]
    failures = [sum(weights[k:])*(traces[k-1]-traces[k])/12 for k in range(1, 5)]
    rejection = 1-sum(weights)
    m.need(rejection+sum(successes)+sum(failures) == 1 and weights[0] > 0, 'complete preparation probability law and full-rank floor')
    fraction = div(mul(scalar(3), slow_mass), scalar(12*q))
    initial = (0,)*10+(1, 1)
    m.need({p[0] for p in frozen['actions']} == set(range(12)), 'source group is transitive on work sites')
    native_centroid = [div(total(mul(scalar(n), row[p]) for n, row in zip(initial, frame)), scalar(2)) for p in range(3)]
    m.need(native_centroid != [ZERO]*3 and all(total(row[p] for row in frame) == ZERO for p in range(3)), 'source-centered centroid countercontrol')
    census = Counter()
    for e1 in seams:
        for c1 in (0, 1):
            middle, j1 = endpoint(initial, e1, c1)
            feedback = int(j1*(int(e1[1] == 10)-int(e1[0] == 10)) < 0)
            for e2 in seams:
                census[int(10 in e1), feedback, int(10 in e2)] += 2
    prep = {'weights': list(map(str, weights)), 'filter': 'F=I12-L/120',
        'preparation_services_are_declared': True,
        'reference': 'M=sum_n w_n F^(2n); rho=M/Tr(M)',
        'request_rejection': str(rejection), 'success_by_request': list(map(str, successes)),
        'first_failure_by_attempt': list(map(str, failures)), 'success_probability': str(q),
        'complete_probability_mass': '1', 'trace_F_even_powers': list(map(str, traces)),
        'rank_of_ideal_success_work_state': 12, 'minimum_unnormalized_eigenvalue_lower': str(weights[0]),
        'slow_population': fraction, 'outside_slow_population': sub(ONE, fraction),
        'normalized_slow_block': 'P_s/3; same for every retained request',
        'current_protocol_band': 2, 'current_active_pair': [0, 2],
        'slow_detector_addresses': [1, 3, 5], 'selected_slow_Hermitian_observable_rank': 3,
        'full_slow_Hermitian_dimension': 9, 'native_initial_configuration': list(initial),
        'work_initial_density': 'supplied I12/12 with work/reference purification',
        'centroid_dictionary_control': {'native_initial_centroid': native_centroid,
            'ideal_success_work_centroid': [ZERO]*3, 'automatic_identification_preserved': False,
            'scope': 'same native port frame, total-two centroid versus one-site work position, before continuing reads'},
        'native_tapes': 3600, 'controller_census': [{'control': list(k), 'count': v} for k, v in sorted(census.items())],
        'ideal_slow_terminal_each_address_given_request_and_tape': 'w_n*f_s^(2n)/(12*3600)',
        'ideal_slow_block_encodes_native_feedback': False, 'supplied_work_reference_is_native_reference': False,
        'failure_poststates_reconstructed_by_this_probe': False,
        'raw_filter_or_clock_words_reexecuted': False, 'raw_numerical_guards_promoted_to_rigorous_error_bounds': False}
    expected = prior.encoded({'schema': 'oph.c12-current-attachment-probe.v1', 'work_baseline': BASELINE,
        'source_pins': pins, 'prior_probe_pins': previous,
        'implementation_sha256': {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        'mean_observable_join': join, 'slow_code_control': control, 'C12_preparation': prep,
        'native_record_grammar': native_profile(seams),
        'verdict': 'SOURCE_INTERFACE_UNDERDETERMINED',
        'certified_subverdict': 'RECORD_ONLY_C12_COHERENT_ATTACHMENT_OBSTRUCTED',
        'classical_moment_attachment_present': True, 'universal_OPH_or_compiler_impossibility': False,
        'physical_golden_embedding_selected': False, 'phase_0_stop_preserved': True})
    m.need(json.dumps(packet, sort_keys=True, allow_nan=False) == json.dumps(expected, sort_keys=True, allow_nan=False),
           'independent C12-current probe mismatch')
    return {'verified': True, 'verdict': expected['verdict'], 'subverdict': expected['certified_subverdict'],
        'native_mean_events': 3960, 'native_matrix_units_covered': 182520,
        'method': 'independent endpoint averaging, spectral traces and partial-injection support proof'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=REFERENCE)
    args = parser.parse_args()
    print(json.dumps(verify(m.load(args.receipt)), indent=2, sort_keys=True))
