"""Exact native-current/C12 interface audit; no coherent native access added."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import proper_code_attachment_probe as prior
from galois_port_frame_certificate import F, Z, O, encode, mm, transpose, require
from model.observer import native, quantum_mixing, coherent_continuation
from model.identifiability_checks import rule_observability
from model.predictions import inverse_clock, prepared_clock

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


def digest(value):
    return hashlib.sha256(json.dumps(encode(value), sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def custody(paths, commit, repository):
    result = {}
    for path in paths:
        data = (ROOT/path).read_bytes()
        require(data == subprocess.check_output(['git', 'show', commit+':'+path], cwd=ROOT), 'source custody: '+path)
        result[path] = {'repository': repository, 'commit': commit,
            'blob': hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest(),
            'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
    return result


def native_join():
    states, seams = native.states(binary=True), quantum_mixing.seams()
    lookup = {x: i for i, x in enumerate(states)}
    counts = [[0]*66 for _ in states]
    for i, x in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                y, _ = native.native(x, seam, coin)
                counts[lookup[y]][i] += 1
    g = [[6*x[p]-1 for p in range(12)] for x in states]
    gram = [[sum(row[i]*row[j] for row in g) for j in range(12)] for i in range(12)]
    require(gram == [[360*int(i == j)-30 for j in range(12)] for i in range(12)], 'native centered-load Gram')
    lap = [[0]*12 for _ in range(12)]
    for a, b in seams:
        lap[a][a] += 1; lap[b][b] += 1
        lap[a][b] -= 1; lap[b][a] -= 1
    require(counts == list(map(list, zip(*counts))), 'native kernel self-adjointness')
    k12 = [[60*int(i == j)-lap[i][j] for j in range(12)] for i in range(12)]
    require([[sum(counts[i][j]*g[j][p] for j in range(66)) for p in range(12)] for i in range(66)]
            == [[sum(g[i][j]*k12[j][p] for j in range(12)) for p in range(12)] for i in range(66)],
            'same-source mean transfer intertwining')
    source = coherent_continuation.density_sector_audit()
    require(source['all_density_compressions_intertwined'], 'existing density source audit')
    return states, g, lap, {'native_binary_dimension': 66, 'work_site_dimension': 12,
        'mean_zero_rank': 11, 'G_formula': 'G[x,p]=6*x_p-1',
        'scheduler': 'uniform thirty seams and fair placement coins; declared instance',
        'G_star_G': '360 I12 - 30 J12', 'intertwining': 'K66 G = G (I12-L/60)',
        'native_events_checked': 3960, 'kernel_count_sha256': digest(counts), 'G_sha256': digest(g),
        'source_density_audit_replayed': True, 'state_encoding_operation_implemented': False}


def slow_pair(g):
    frame = prior.marked.raw_frame()
    vertices = prior.marked.response_vertices()
    relabel = [vertices.index(v) for v in frame]
    mark = prior.marked.source_marked_action(frame, vertices, relabel)
    frozen = prior.marked.build()
    index = frozen['actions'].index(list(mark))
    require(index == 4, 'source marking before character use')
    b = [[sum((F(gij)*frame[j][p] for j, gij in enumerate(row)), Z) for p in range(3)] for row in g]
    scale = F(360)*F(10, 2)
    require(mm(transpose(b), b) == [[scale*int(i == j) for j in range(3)] for i in range(3)], 'source slow-code Gram')
    states = native.states(binary=True)
    lookup = {x: i for i, x in enumerate(states)}
    rows = json.loads((ROOT/'code/source_selection_model/response.json').read_text())['closed_paths']
    for action, row in zip(frozen['actions'], rows, strict=True):
        u = [[F(Q(x[0]), Q(x[1])) for x in r] for r in row['observed_response'][0]]
        for i, x in enumerate(states):
            y = [0]*12
            for p, bit in enumerate(x):
                y[action[p]] = bit
            require(b[lookup[tuple(y)]] == [sum((u[p][q]*b[i][q] for q in range(3)), Z) for p in range(3)],
                    'native slow code has the same source marked upper action')
    a, c = [row[0] for row in b], [row[1] for row in b]
    diagonal = [(x*x+y*y)/(2*scale) for x, y in zip(a, c)]
    require(sum(diagonal, Z) == O and all(x.sign() >= 0 for x in diagonal), 'normalized equal diagonals')
    i, j = next((i, j) for i in range(66) for j in range(i+1, 66) if c[i]*a[j] != a[i]*c[j])
    coefficient = (c[i]*a[j]-a[i]*c[j])/(2*scale)
    return {'current_marking_index': index, 'historical_marking_index': 1,
        'current_native_action': list(mark), 'coordinate_choice': [0, 1],
        'native_code_Gram_scale': scale, 'actions_checked': 60,
        'witness_definition': 'rho_+-=(aa^T+bb^T +/- i(ba^T-ab^T))/(2k), [a,b]=first two columns of G raw_frame',
        'witness_conjugation': 'complex conjugation at fixed real embedding; not Galois conjugation',
        'vectors_sha256': digest([a, c]), 'shared_diagonal_sha256': digest(diagonal),
        'trace_each': 1, 'purity_each': 1, 'mutually_orthogonal': True,
        'nonzero_imaginary_entry': {'row': i, 'column': j, 'plus_coefficient': coefficient},
        'native_physical_preparation_of_witness_claimed': False,
        'ideal_C12_images_orthogonal': True, 'classical_record_only_outputs_equal': True,
        'complete_deterministic_transfer_worst_case_halftrace_floor': '1/2',
        'uniform_herald_success_p_floor': 'p/2; p>0 prescribed for the target channel',
        'arbitrary_untouched_references_covered': True,
        'coherent_native_output_retained_by_receiver': False,
        'work_input_independent_of_unknown_native_payload': True,
        'additional_coherent_native_work_access_admitted': False}


def preparation(lap):
    weights = inverse_clock.subordinated_filter_weights()
    f = [[Q(int(i == j))-Q(lap[i][j], 120) for j in range(12)] for i in range(12)]
    power = [[Q(int(i == j)) for j in range(12)] for i in range(12)]
    traces = []
    diagonal = [Q()]*12
    for n in range(len(weights)):
        traces.append(sum(power[i][i] for i in range(12)))
        for i in range(12):
            diagonal[i] += weights[n]*power[i][i]
        power = [[sum(power[i][k]*f[k][j] for k in range(12)) for j in range(12)] for i in range(12)]
        power = [[sum(power[i][k]*f[k][j] for k in range(12)) for j in range(12)] for i in range(12)]
    q = sum((w*t/12 for w, t in zip(weights, traces)), Q())
    successes = [w*t/12 for w, t in zip(weights, traces)]
    failures = [sum(weights[k:])*(traces[k-1]-traces[k])/12 for k in range(1, len(weights))]
    rejection = 1-sum(weights)
    require(rejection+sum(successes)+sum(failures) == 1, 'all request/failure/success mass')
    lam = F(5, -1)
    factor = O-lam/120
    slow = Z
    power_s = O
    for w in weights:
        slow += F(w)*power_s
        power_s = power_s*factor*factor
    fraction = 3*slow/F(12*q)
    require(fraction.sign() > 0 and (O-fraction).sign() > 0 and weights[0] > 0, 'all-band support, not a slow herald')
    census = Counter({key: len(rows) for key, rows in prepared_clock.native_histories()})
    require(sum(census.values()) == 3600 and prepared_clock.BAND == 2, 'actual continuing protocol scope')
    addresses = list(prepared_clock.source.BAND_ADDRESSES[1])
    active = prepared_clock.source.PAIR_ADDRESSES[prepared_clock.BAND-1]
    require(active not in addresses and 0 not in addresses and addresses == [1, 3, 5], 'selected menu leaves the slow block unmixed')
    require(diagonal == [q]*12, 'the prepared ideal work density is transitive on site diagonals')
    frame = prior.marked.raw_frame()
    native_centroid = [sum((F(n)*row[p]/2 for n, row in zip(prepared_clock.INITIAL, frame)), Z) for p in range(3)]
    work_centroid = [sum((F(v/(12*q))*row[p] for v, row in zip(diagonal, frame)), Z) for p in range(3)]
    require(work_centroid == [Z]*3 and native_centroid != work_centroid, 'natural centroid dictionary fails at this preparation checkpoint')
    return {'weights': list(map(str, weights)), 'filter': 'F=I12-L/120',
        'preparation_services_are_declared': True,
        'reference': 'M=sum_n w_n F^(2n); rho=M/Tr(M)',
        'request_rejection': str(rejection), 'success_by_request': list(map(str, successes)),
        'first_failure_by_attempt': list(map(str, failures)), 'success_probability': str(q),
        'complete_probability_mass': '1', 'trace_F_even_powers': list(map(str, traces)),
        'rank_of_ideal_success_work_state': 12, 'minimum_unnormalized_eigenvalue_lower': str(weights[0]),
        'slow_population': fraction, 'outside_slow_population': O-fraction,
        'normalized_slow_block': 'P_s/3; same for every retained request',
        'current_protocol_band': prepared_clock.BAND, 'current_active_pair': [0, active],
        'slow_detector_addresses': addresses, 'selected_slow_Hermitian_observable_rank': 3,
        'full_slow_Hermitian_dimension': 9, 'native_initial_configuration': list(prepared_clock.INITIAL),
        'work_initial_density': 'supplied I12/12 with work/reference purification',
        'centroid_dictionary_control': {'native_initial_centroid': native_centroid,
            'ideal_success_work_centroid': work_centroid, 'automatic_identification_preserved': False,
            'scope': 'same native port frame, total-two centroid versus one-site work position, before continuing reads'},
        'native_tapes': 3600, 'controller_census': [{'control': list(k), 'count': v} for k, v in sorted(census.items())],
        'ideal_slow_terminal_each_address_given_request_and_tape': 'w_n*f_s^(2n)/(12*3600)',
        'ideal_slow_block_encodes_native_feedback': False,
        'supplied_work_reference_is_native_reference': False,
        'failure_poststates_reconstructed_by_this_probe': False,
        'raw_filter_or_clock_words_reexecuted': False, 'raw_numerical_guards_promoted_to_rigorous_error_bounds': False}


def build():
    source_pins = custody(INPUTS, SOURCE_COMMIT, 'FloatingPragma/observer-patch-holography')
    prior_pins = custody(PRIOR_FILES, PRIOR_COMMIT, 'FloatingPragma/observer-patch-holography')
    _, g, lap, join = native_join()
    branch = rule_observability.intertwining_certificate()
    require(branch['matrix_units_checked'] == 182520, 'complete native diagonal-block control')
    return encode({'schema': 'oph.c12-current-attachment-probe.v1', 'work_baseline': BASELINE,
        'source_pins': source_pins, 'prior_probe_pins': prior_pins,
        'implementation_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in IMPLEMENTATIONS},
        'mean_observable_join': join, 'slow_code_control': slow_pair(g), 'C12_preparation': preparation(lap),
        'native_record_grammar': {'native_dimension': 78, 'seams': 30,
            'matrix_units_checked': branch['matrix_units_checked'],
            'C_recorded_branch_count': sum(r['recorded_partial_injections']['C'] for r in branch['seams']),
            'E_recorded_branch_count': sum(r['recorded_partial_injections']['E'] for r in branch['seams']),
            'channel_identity': 'Phi=Phi after native-basis dephasing when terminal native coherence is inaccessible',
            'work_controls': 'arbitrary complete work-only CP instruments conditioned on classical native records',
            'additional_coherent_native_controls_admitted': False,
            'coherent_native_work_coupling_admitted': False,
            'coherent_uncommitted_native_record_access_admitted': False,
            'terminal_native_access': 'trace or diagonal read only',
            'finite_adaptive_histories_and_all_failure_flags': True},
        'verdict': 'SOURCE_INTERFACE_UNDERDETERMINED',
        'certified_subverdict': 'RECORD_ONLY_C12_COHERENT_ATTACHMENT_OBSTRUCTED',
        'classical_moment_attachment_present': True,
        'universal_OPH_or_compiler_impossibility': False, 'physical_golden_embedding_selected': False,
        'phase_0_stop_preserved': True})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--write', type=Path)
    args = parser.parse_args()
    text = json.dumps(build(), indent=2, sort_keys=True)+'\n'
    if args.check:
        require(REFERENCE.read_text() == text, 'C12 current probe receipt differs')
        print('C12_CURRENT_ATTACHMENT_VALID: RECORD_ONLY_C12_COHERENT_ATTACHMENT_OBSTRUCTED')
    elif args.write:
        with args.write.open('x') as stream:
            stream.write(text)
    else:
        print(text, end='')
