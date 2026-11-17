"""Whole-sector identities, countercontrols and hostile receipt inputs."""

import copy
import fnmatch
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pytest

from . import build, circuits, experiments, geometry_check, model, preparation, spatial, verify, walk
from .check import annihilator, at, exact, exterior, verify_small, X, Z
from .experiment_check import verify_detector, verify_interaction, verify_noise, verify_resolved_reads
from .pauli import Word


@pytest.fixture(scope='module')
def evidence():
    return build.candidate()


@pytest.fixture(scope='module')
def small_cache(evidence):
    return verify_small(evidence['small_graphs'])


def test_receipt_and_parent_replay():
    verify.verify()


def test_rebuilt_physics_matches_receipt(evidence):
    stored = verify.parent.load(verify.HERE/'receipt.json')
    # BLAS cancellation can change a matrix component of order 1e-17.
    # Replay both payloads first; compare complex arrays in the same norm
    # tolerance as their independent matrix identities. Scalar resource
    # bounds, including the tiny positive noise cap, keep strict relative
    # comparison and never get this absolute floor.
    verify.verify_evidence(evidence)
    verify.verify_evidence(stored['evidence'])
    arrays = {'state', 'isometry', 'scalar', 'parameter', 'unitary', 'matrix', 'kraus',
              'kernel', 'code_effect', 'initial', 'encoded_output', 'output', 'banks'}
    def compare(actual, expected, key=None):
        assert type(actual) is type(expected)
        if key in arrays and isinstance(actual, list):
            np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-13)
        elif isinstance(actual, dict):
            assert actual.keys() == expected.keys()
            for name in actual:
                compare(actual[name], expected[name], name)
        elif isinstance(actual, list):
            assert len(actual) == len(expected)
            for a, e in zip(actual, expected):
                compare(a, e)
        else:
            exact(actual, expected)
    compare(evidence, stored['evidence'])


@pytest.mark.parametrize('n,edges', [(True, [(0, 1)]), (1, []), (3, [(0, 1)]),
    (2, [(0, 0)]), (2, [(1, 0)]), (2, [(0, 1), (0, 1)]), (2, [(False, 1)]), (2, [(0, 2)])])
def test_invalid_graphs(n, edges):
    with pytest.raises(ValueError):
        model.Encoding(n, edges)


@pytest.mark.parametrize('parity', [0, 1])
def test_nondefault_edge_order_and_every_gate(parity):
    enc = model.Encoding(4, [(2, 3), (0, 3), (1, 3), (0, 1), (1, 2)])
    basis, iso = enc.isometry(parity)
    a = [annihilator(4, i) for i in range(4)]
    for i, j in list(enc.edges)+[(j, i) for i, j in enc.edges]:
        gamma_i, gamma_j = a[i]+a[i].conj().T, a[j]+a[j].conj().T
        assert np.allclose(enc.a(i, j).matrix(enc.m)@iso,
                           iso@(-1j*gamma_i@gamma_j)[np.ix_(basis, basis)])
        for g in (np.eye(2), 1j*X, np.array([[0, np.exp(.4j)], [-np.exp(-.4j), 0]])):
            single = np.eye(4, dtype=complex)
            single[np.ix_([i, j], [i, j])] = g
            assert np.allclose(enc.execute(enc.su2(i, j, g, parity))@iso,
                               iso@exterior(single)[np.ix_(basis, basis)])


def test_local_order_change_is_paid_conjugation_not_state_fault():
    enc = model.Encoding(4, list(itertools.combinations(range(4), 2)))
    orders = copy.deepcopy(enc.incident)
    e, f = orders[0][:2]
    orders[0][:2] = [f, e]
    other = model.Encoding(4, enc.edges, orders)
    cz = np.diag([(-1.)**(((s >> e) & 1)*((s >> f) & 1)) for s in range(1 << enc.m)])
    for i, j in enc.edges:
        assert np.allclose(cz@enc.a(i, j).matrix(enc.m)@cz, other.a(i, j).matrix(enc.m))
    for s, t in zip(enc.loops, other.loops):
        assert np.allclose(cz@s.matrix(enc.m)@cz, t.matrix(enc.m))
    state, _ = enc.prepare()
    fault = Word(z=1 << enc.chords[0]).matrix(enc.m)@state
    assert np.linalg.norm(fault-state) > 1
    assert any(np.vdot(fault, s.matrix(enc.m)@fault).real < -.99 for s in enc.loops)
    z = Word(z=1 << e).matrix(enc.m)
    for k, (i, j) in enumerate(enc.edges):
        assert np.allclose(z@enc.a(i, j).matrix(enc.m)@z,
                           (-1 if k == e else 1)*enc.a(i, j).matrix(enc.m))


def test_ordinary_swap_fails_inside_fixed_two_particle_sector():
    # Superposition of two configurations, BOTH with N=2. This avoids a
    # counterexample that only compares an unobservable inter-parity phase.
    n = 3
    permutation = np.eye(n)
    permutation[[0, 1]] = permutation[[1, 0]]
    fermionic = exterior(permutation)
    ordinary = np.zeros((8, 8))
    for s in range(8):
        t = (s & ~3) | ((s & 1) << 1) | ((s & 2) >> 1)
        ordinary[t, s] = 1
    psi = np.zeros(8)
    psi[[3, 5]] = 1/np.sqrt(2)
    assert abs(np.vdot(fermionic@psi, ordinary@psi)) < 1e-12
    enc = model.Encoding(n, [(0, 1), (0, 2), (1, 2)])
    basis, iso = enc.isometry(0)
    assert np.allclose(enc.execute(enc.fswap(0, 1))@iso@psi[basis], iso@(fermionic@psi)[basis])


@pytest.mark.parametrize('raw,theta', [(Word(3, 5, 1), .77), (Word(5, 7, 0), -.38), (Word(0, 7, 2), .9)])
def test_native_parity_rotation_and_both_instrument_outcomes(raw, theta):
    n = 3
    p = raw.matrix(n)
    u = circuits.execute(circuits.pauli_tape(raw, theta, n), n+1)
    expected = np.cos(theta)*np.eye(16)+1j*np.sin(theta)*np.kron(Z, p)
    assert np.allclose(u, expected)
    v = circuits.execute(circuits.measurement_tape(raw, n), n+1)
    assert np.allclose(v[:8, :8], (np.eye(8)+p)/2)
    assert np.allclose(v[8:, :8], (np.eye(8)-p)/2)
    assert np.allclose(v.conj().T@v, np.eye(16))


@pytest.mark.parametrize('weight', [-.01, 1.01, True, float('nan'), float('inf')])
def test_invalid_weighted_acceptance(weight):
    with pytest.raises(ValueError):
        circuits.weighted_program(Word(0, 7), weight, 3)


def test_same_click_effect_cannot_hide_changed_postmeasurement_state(evidence):
    from .check import complex_array
    from .experiment_check import verify_weighted_measurements
    rows = copy.deepcopy(evidence['weighted_measurements'])
    k = complex_array(rows[0]['kraus'][3], (8, 8))
    # A data unitary after the click leaves its POVM effect identical but
    # changes the retained quantum output. Comparing effects alone misses it.
    wrong = at(Z, 1, 3)@k
    assert np.allclose(wrong.conj().T@wrong, k.conj().T@k)
    assert np.linalg.norm(wrong-k) > .5
    rows[0]['kraus'][3] = model.encode(wrong)
    with pytest.raises(ValueError):
        verify_weighted_measurements(rows)


def test_occupied_helper_must_be_reset_before_acceptance():
    program = circuits.weighted_program(Word(0, 3), .37, 2)
    proper = circuits.execute_weighted(program)
    program['reset_before_acceptance'] = False
    wrong = circuits.execute_weighted(program)
    assert np.linalg.norm(proper[3].conj().T@proper[3]-wrong[3].conj().T@wrong[3]) > .2


@pytest.mark.parametrize('q', [1, 2, 3, 4])
def test_local_decoder_all_syndromes_and_exact_counts(q):
    vertices, edges, cycles, _, correction, _ = preparation.plan(q)
    assert len(edges) == 63*q**3-3*q*q
    assert len(cycles) == 31*q**3-3*q*q+1
    lookup = {e: k for k, e in enumerate(edges)}
    # A linear identity covers every syndrome, including the all-negative one.
    for k, path in enumerate(cycles):
        decoded = 0
        for i, j in zip(path, path[1:]):
            decoded ^= correction[lookup[tuple(sorted((i, j)))]]
        assert decoded == 1 << k
    assert len(vertices) == 32*q**3


def test_color_schedule_has_no_clock_processor_collisions():
    q = 5
    vertices, edges, flights = spatial.graph(q)
    incident = [set() for _ in vertices]
    for e, (u, v) in enumerate(edges):
        incident[u].add(e)
        incident[v].add(e)
    ids = {v: k for k, v in enumerate(vertices)}
    for kind, modes, _ in walk.layers(.1):
        occupied = {}
        if kind == 'ferry':
            pairs = [(u, v) for u, v in flights if vertices[u][2] == modes[0]]
        else:
            pairs = [(ids[(x, 0, modes[0])], ids[(x, 1 if kind == 'restore' else 0, modes[-1])])
                     for x in itertools.product(range(q), repeat=3)]
        for u, v in pairs:
            anchor = vertices[u][0]
            color = tuple(t % 4 for t in anchor)
            support = incident[u] | incident[v]
            assert not (occupied.setdefault(color, set()) & support)
            occupied[color].update(support)
            for e in support:
                i, j = edges[e]
                midpoint = [(a+b)/2 for a, b in zip(vertices[i][0], vertices[j][0])]
                assert max(abs(s-t) for s, t in zip(midpoint, anchor)) <= 1.5


@pytest.mark.parametrize('q,a', [(1, .07), (2, .13), (3, .03)])
def test_entire_actual_spatial_matrix_and_blank_bank(q, a):
    matrix, bank = walk.execute(q, a)
    assert np.allclose(matrix, geometry_check.closed_walk(q, a))
    assert np.allclose(matrix.conj().T@matrix, np.eye(len(matrix)))
    assert np.max(np.abs(bank)) < 1e-12


def test_preparation_time_is_diameter_not_volume():
    times = [spatial.budget(q, 1/q)['preparation_time_upper'] for q in (4, 8, 16, 32)]
    assert all(a > b for a, b in zip(times, times[1:]))
    assert abs(spatial.budget(4, .1)['wall_factor']-spatial.budget(100, .003)['wall_factor']) < 1e-7
    for q in (1, 2, 8):
        b = spatial.budget(q, .1)
        assert b['event_time'] > 0 and b['native_drive_bound'] > 0
        assert b['phase_rate_for_point_zero_zero_one'] > 0


def test_complete_detector_and_additive_countercontrol(evidence, small_cache):
    verify_detector(evidence['detector'], small_cache)
    assert .37+.81 > 1  # Fully occupied: naive summed effect is not an effect.
    assert 0 < 1-(1-.37)*(1-.81) < 1


def test_complex_read_conjugation_matters(evidence, small_cache):
    verify_resolved_reads(evidence['resolved_reads'], small_cache)
    f = np.array([1., .3j, 0., -.4+.2j, .1])
    f /= np.linalg.norm(f)
    right = sum(f[i].conjugate()*annihilator(5, i) for i in range(5))
    wrong = sum(f[i]*annihilator(5, i) for i in range(5))
    assert np.linalg.norm(right.conj().T@right-wrong.conj().T@wrong) > .5


def test_actual_interaction_and_vacuum_noise(evidence, small_cache):
    verify_interaction(evidence['interaction'], small_cache)
    verify_noise(evidence['noise'], small_cache)
    assert all(r['failure_lower'] > 0 for r in evidence['noise'])


def test_interaction_input_is_prepared_from_even_vacuum():
    enc = model.Encoding(4, list(itertools.combinations(range(4), 2)))
    basis, iso = enc.isometry(0)
    vacuum, _ = enc.prepare()
    create_pair = enc.execute((1., [(enc.a(0, 2), np.pi/2)]))
    g = np.array([[1., -1.], [1., 1.]])/np.sqrt(2)
    actual = enc.execute(enc.su2(2, 3, g))@enc.execute(enc.su2(0, 1, g))@create_pair@vacuum
    target = np.zeros(16)
    target[[5, 6, 9, 10]] = .5
    assert np.allclose(actual, iso@target[basis])


def test_conversion_preserves_entangled_input(evidence):
    from .check import complex_array
    row = evidence['conversions'][1]
    maps = [complex_array(k, (8, 4)) for k in row['corrected_branches']]
    bell = np.eye(4).reshape(-1)/2
    states = [np.kron(k, np.eye(4))@bell for k in maps]
    rho = sum(np.outer(s, s.conj()) for s in states)
    assert np.isclose(np.trace(rho@rho), 1)
    # Tracing the new code qubit is deliberately not claimed as inverse conversion.
    restored = sum(rho[k*16:(k+1)*16, k*16:(k+1)*16] for k in range(2))
    assert np.linalg.norm(restored-np.outer(bell, bell)) > .1


MUTATIONS = [
    (('small_graphs',), []),
    (('small_graphs', 0, 'blocks'), []),
    (('small_graphs', 0, 'blocks', 1, 'parity'), 0),
    (('small_graphs', 0, 'branches'), []),
    (('small_graphs', 0, 'branches', 1, 'probability'), 1.),
    (('small_graphs', 0, 'branches', 0, 'state', 0, 0), .123),
    (('small_graphs', 0, 'loops', 0, 2), 0),
    (('small_graphs', 0, 'edge_operators', 0, 2), 2),
    (('small_graphs', 0, 'blocks', 0, 'isometry', 0, 0, 0), .321),
    (('small_graphs', 0, 'blocks', 0, 'gates', 2, 'scalar'), [1., 0.]),
    (('native_rotations', 0, 'tape'), []),
    (('native_rotations', 0, 'pulse_code_events'), 0),
    (('native_measurements', 0, 'unitary', 0, 0, 0), 0.),
    (('native_measurements', 2, 'tape'), [['h', 0]]),
    (('weighted_measurements', 0, 'program', 'outcomes'), [[1, 1]]),
    (('weighted_measurements', 0, 'program', 'reset_before_acceptance'), False),
    (('weighted_measurements', 0, 'program', 'acceptance_angles'), [0., 0.]),
    (('weighted_measurements', 0, 'program', 'code_events_per_path'), 0),
    (('weighted_measurements', 0, 'kraus'), []),
    (('geometry', 1, 'graph_sha256'), '0'*64),
    (('geometry', 1, 'flight_sha256'), '0'*64),
    (('preparation', 1, 'decoder_sha256'), '0'*64),
    (('preparation', 1, 'signed_loops_sha256'), '0'*64),
    (('preparation', 1, 'prefix_rounds'), 0),
    (('preparation', 1, 'maximum_slots_per_anchor'), 1),
    (('budgets', 0, 'event_time'), 0.),
    (('budgets', 0, 'preparation_time_upper'), 0.),
    (('budgets', 0, 'phase_rate_for_point_zero_zero_one'), 0.),
    (('budgets', 0, 'noisy_interfaces_upper'), 1),
    (('walks', 0, 'banks', 0, 0), .1),
    (('walks', 0, 'program', 0, 'parameter'), .1),
    (('walks', 0, 'program', 24, 'modes'), [0]),
    (('conversions', 1, 'corrected_branches'), []),
    (('resolved_reads', 0, 'blocks'), []),
    (('resolved_reads', 0, 'kernel', 1, 1), -.123),
    (('detector', 0, 'weights'), [.37, 1.81]),
    (('detector', 0, 'instruments'), []),
    (('detector', 1, 'instruments', 0, 'outcomes'), [[0, 0], [1, 0], [1, 1]]),
    (('interaction', 'encoded_output', 0, 0), .7),
    (('interaction', 'circuit', 'rotations', 0, 'theta'), 0.),
    (('noise', 0, 'vacuum_fidelity'), 1.),
    (('noise', 0, 'failure_lower'), 0.),
    (('clock_witness', 'phase_rate_cap'), 0.),
    (('clock_witness', 'report_time_upper'), 0.),
    (('clock_witness', 'primitive_events_upper'), 1),
    (('clock_witness', 'accounting_range_upper'), 1.),
    (('clock_witness', 'swing_lower'), float('nan')),
    (('clock_witness', 'ticks'), True),
]


@pytest.mark.parametrize('path,value', MUTATIONS, ids=[str(p) for p, _ in MUTATIONS])
def test_physics_mutations_rejected_even_with_custody_bypassed(evidence, path, value):
    packet = copy.deepcopy(evidence)
    target = packet
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        verify.verify_evidence(packet)


@pytest.mark.parametrize('field,factor', [('wall_tick', 1-1e-11), ('total_exposure_upper', 1-1e-11),
    ('phase_rate_cap', 1+1e-11), ('trace_error_upper', 1-1e-11),
    ('accounting_error_upper', 1-1e-11), ('swing_lower', 1+1e-11)])
def test_nearby_false_bounds_fail_even_inside_replay_tolerance(evidence, field, factor):
    row = copy.deepcopy(evidence['clock_witness'])
    row[field] *= factor
    with pytest.raises(ValueError):
        geometry_check.verify_clock_witness(row)


@pytest.mark.parametrize('text', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}'])
def test_strict_json(tmp_path, text):
    path = tmp_path/'bad.json'
    path.write_text(text, encoding='utf-8')
    with pytest.raises(ValueError):
        verify.parent.load(path)


@pytest.mark.parametrize('kind', ['missing', 'extra', 'integer_matrix', 'boolean_angle'])
def test_schema_and_numeric_types_fail_closed(evidence, kind):
    bad = copy.deepcopy(evidence)
    if kind == 'missing':
        del bad['noise']
    elif kind == 'extra':
        bad['assume_success'] = True
    elif kind == 'integer_matrix':
        bad['small_graphs'][0]['blocks'][0]['isometry'][0][0][0] = 0
    else:
        bad['small_graphs'][0]['blocks'][0]['gates'][0]['rotations'][0]['theta'] = True
    with pytest.raises(ValueError):
        verify.verify_evidence(bad)


def test_cli_optimized_rejects_custody_and_physics(tmp_path):
    packet = verify.parent.load(verify.HERE/'receipt.json')
    env = dict(os.environ, PYTHONPATH=str(verify.ROOT/'code'))
    for field in ('sources', 'evidence'):
        bad = copy.deepcopy(packet)
        if field == 'sources':
            bad[field][next(iter(bad[field]))] = '0'*64
        else:
            bad[field]['clock_witness']['phase_rate_cap'] = 0.
        path = tmp_path/(field+'.json')
        path.write_text(json.dumps(bad), encoding='utf-8')
        result = subprocess.run([sys.executable, '-O', '-m', 'm1_fermionic_source.verify', str(path)],
                                env=env, cwd=verify.ROOT, capture_output=True, text=True)
        assert result.returncode != 0 and 'ValueError' in result.stderr


def test_checker_runs_with_producers_disabled():
    script = "import sys\n"
    for name in ('build', 'circuits', 'experiments', 'model', 'pauli', 'preparation', 'spatial', 'walk'):
        script += f"sys.modules['m1_fermionic_source.{name}'] = None\n"
    script += "from m1_fermionic_source import verify\nverify.verify()\n"
    result = subprocess.run([sys.executable, '-O', '-c', script], cwd=verify.ROOT,
                            env=dict(os.environ, PYTHONPATH=str(verify.ROOT/'code')),
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


