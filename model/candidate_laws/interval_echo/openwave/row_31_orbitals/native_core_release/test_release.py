"""Native local interventions, finite source instrument and held orbital reads."""

import gzip
import json
from itertools import product

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import cell
from . import measure as m


def receipt():
    return json.loads(gzip.decompress((m.HERE / 'results.json.gz').read_bytes()))


def test_native_echo_on_full_memory_and_load_space():
    law = m.candidate.native_law
    # Independent sparse integer construction of the native write isometry.
    w = np.zeros((90, 9))
    for j, (a, b) in enumerate(law.SRC):
        for coin in (0, 1):
            amount = a + b
            x = (amount // 2, amount - amount // 2)
            if coin:
                x = x[::-1]
            w[9 * law.LABELS.index((coin, a - b)) + law.SRC.index(x), j] += 1 / np.sqrt(2)
    u = np.block([[np.zeros((9, 9)), w.T], [w, np.eye(90) - w @ w.T]])
    np.testing.assert_allclose(u, law.echo(), atol=1e-15)
    np.testing.assert_allclose(u @ u, np.eye(99), atol=6e-16)
    r = (u @ np.eye(99)[:, law.SRC.index((0, 0))]).reshape(11, 9)[:, 0]
    for word in ((0, 0), (1, 0), (0, 1)):
        j = law.SRC.index(word)
        initial = np.kron(r, np.eye(9)[:, j])
        expected = np.eye(99)[:, j] if word == (0, 0) else initial
        np.testing.assert_allclose(u @ initial, expected, atol=5e-16)


def test_finite_source_filter_from_explicit_clock_and_refusal_energy(monkeypatch):
    board = cell()
    seam = tuple(board.seams[0])
    h = m.transition.generator(board, seam, 8)
    e, v = np.linalg.eigh(h)
    source = np.eye(board.n)[:, seam[0]]
    bins, step = 17, .031
    # Sum a separate short coherent clock, without the spectral kernel formula.
    direct = sum(np.exp(1j * j * step * m.TARGET) * expm(-1j * j * step * h)
                 @ source / bins for j in range(bins))
    phase = step * (e - m.TARGET)
    kernel = np.exp(-.5j * (bins - 1) * phase) * np.sinc(bins * phase / (2 * np.pi)) / np.sinc(phase / (2 * np.pi))
    np.testing.assert_allclose(direct, v @ (kernel * v[seam[0]]), atol=2e-14)
    monkeypatch.setattr(m.transition, 'BINS', bins)
    monkeypatch.setattr(m.transition, 'STEP', step)
    prepared, _, report = m.prepare(h, seam[0])
    np.testing.assert_allclose(prepared, direct, atol=2e-14)
    assert abs(report['acceptance'] + report['refusal'] - 1) < 1e-12
    assert abs(report['admitted_relative_energy'] + report['refused_relative_energy']
               - report['initial_relative_energy']) < 1e-12


def test_release_chart_matches_native_pulse_not_a_potential_switch():
    board = cell()
    seam = tuple(board.seams[0])
    rng = np.random.default_rng(3131)
    psi = rng.normal(size=board.n) + 1j * rng.normal(size=board.n)
    psi /= np.linalg.norm(psi)
    blank, written = m.branch_states(psi, seam, 'one_native_echo')
    assert np.linalg.norm(blank[list(seam)]) == 0
    assert np.linalg.norm(written[[i for i in range(board.n) if i not in seam]]) == 0
    np.testing.assert_allclose(blank + written, psi)
    law = m.candidate.native_law
    r = np.sqrt(m.native_certificate()['formation_probabilities_by_memory'])
    actual = np.zeros((11, board.n), complex)
    for port in range(board.n):
        word = (int(port == seam[0]), int(port == seam[1]))
        local = (law.echo() @ np.kron(r, np.eye(9)[:, law.SRC.index(word)])).reshape(11, 9)
        actual[:, port] = psi[port] * local[:, law.SRC.index(word)]
    expected = r[:, None] * written[None, :]
    expected[0] += blank
    np.testing.assert_allclose(actual, expected, atol=3e-16)
    # A separate supplied unitary can erase the memory at every load position.
    direction = r - np.eye(11)[0]
    direction /= np.linalg.norm(direction)
    swap = np.eye(11) - 2 * np.outer(direction, direction)
    np.testing.assert_allclose(swap @ r, np.eye(11)[0], atol=3e-16)


def test_complete_records_persistent_native_shape_and_changed_law():
    data = receipt()
    assert sum(data['independent_formation_memory_diagnostics']['counts']) == m.SHOTS
    for panel in data['panels']:
        prep = panel['preparation']
        assert prep['acceptance'] > .96 and prep['conditional_ground_fidelity'] > .9999998
        native, changed = panel['arms']
        for arm in panel['arms']:
            for case in arm['cases']:
                assert len(case['readings']) == len(m.TIMES)
                energies = []
                for record in case['readings']:
                    assert len(record['probabilities']) == 11 * panel['ports'] + 1
                    assert sum(record['counts']) == m.SHOTS
                    assert abs(sum(record['probabilities']) - 1) < 2e-10
                    assert record['sampling_grid_l1'] < 11 * panel['ports'] / data['sampling_weight_denominator'] + 2e-10
                    assert min(record['probabilities']) >= 0
                    public_counts = np.array(record['counts'][:-1]).reshape(11, panel['ports']).sum(axis=0)
                    assert record['public_counts_by_port_then_source_refusal'] == public_counts.tolist() + [record['counts'][-1]]
                    assert sum(record['public_counts_by_port_then_source_refusal']) == m.SHOTS
                    core = int(public_counts[panel['seam']].sum())
                    tail = .01 / (2 * 2 * 2 * len(m.OPERATIONS) * len(m.TIMES))
                    lower, upper = m.count_intervals.clopper_pearson(core, int(public_counts.sum()), tail)
                    # The count is exact. scipy's inverse beta and forward
                    # binomial tails differ across CPU math libraries by a few
                    # units in the last place, and the tail bisection can move
                    # an endpoint by a similar relative amount. Each endpoint
                    # is widened outward by 1e-8 of its distance from the point
                    # estimate; 1e-12 relative stays far inside that margin.
                    np.testing.assert_allclose(
                        record['simultaneous_99_core_given_admission_interval'],
                        [lower, upper], rtol=1e-12, atol=0)
                    assert abs(record['probabilities'][-1] - prep['refusal']) < 1e-14
                    energies.append(record['released_relative_energy'])
                assert max(energies) - min(energies) < 2e-10
        retain, echo, double, supplied = native['cases']
        for a, b in zip(retain['readings'], double['readings'], strict=True):
            np.testing.assert_allclose(a['probabilities'], b['probabilities'], atol=1e-14)
            assert a['core_probability_per_start'] / prep['acceptance'] > .981
            assert a['simultaneous_99_core_given_admission_interval'][0] > .98
        assert .018 < echo['blank_memory_probability_before_flight'] / prep['acceptance'] < .019
        assert abs(supplied['blank_memory_probability_before_flight'] - prep['acceptance']) < 1e-13
        assert supplied['readings'][-1]['core_probability_per_start'] / prep['acceptance'] < .02
        assert changed['cases'][0]['readings'][-1]['core_probability_per_start'] / prep['acceptance'] < .03
        assert changed['cases'][0]['readings'][-1]['simultaneous_99_core_given_admission_interval'][1] < .05
        local, common, hydrogen = panel['transition_records']
        assert local['acceptance_per_start'] > .079
        assert common['acceptance_per_start'] < .0004
        assert hydrogen['acceptance_per_start'] < 3e-8
        assert sum(hydrogen['counts'][:-2]) == 0
        for record in panel['transition_records']:
            assert sum(record['counts']) == m.SHOTS
            assert len(record['counts']) == panel['ports'] + 2


def test_full_native_memory_and_port_embedding_energy_and_pulse():
    board = cell()
    seam = tuple(board.seams[0])
    law = m.candidate.native_law
    other = np.diag(8. * (6 - board.degree))
    for a, b in board.seams:
        if {a, b} == set(seam):
            continue
        other[a, a] += 8.5
        other[b, b] += 8.5
        other[a, b] -= .5
        other[b, a] -= .5

    def embed(local):
        answer = np.zeros((11 * board.n, 11 * board.n))
        for memory in range(11):
            for port in range(board.n):
                word = (int(port == seam[0]), int(port == seam[1]))
                column = 9 * memory + law.SRC.index(word)
                for row in np.flatnonzero(abs(local[:, column]) > 1e-14):
                    out_memory, out_word = divmod(int(row), 9)
                    loads = law.SRC[out_word]
                    if loads == (0, 0):
                        assert port not in seam
                        out_port = port
                    else:
                        assert port in seam and loads in ((1, 0), (0, 1))
                        out_port = seam[loads == (0, 1)]
                    answer[out_memory * board.n + out_port, memory * board.n + port] += local[row, column]
        return answer

    hx, hc = law.local_term_parts()
    full = np.kron(np.eye(11), other) + embed(hx + 8 * hc)
    pulse = embed(law.echo())
    r = np.sqrt(m.native_certificate()['formation_probabilities_by_memory'])
    blank_chart = np.kron(np.eye(11)[:, :1], np.eye(board.n))
    written_chart = np.kron(r[:, None], np.eye(board.n))
    h0 = board.laplacian().toarray() / 2
    hw = m.transition.generator(board, seam, 8)
    np.testing.assert_allclose(full @ blank_chart, blank_chart @ (48 * np.eye(board.n) + h0), atol=2e-14)
    np.testing.assert_allclose(full @ written_chart, written_chart @ (48 * np.eye(board.n) + hw), atol=2e-14)
    psi, _, _ = m.prepare(hw, seam[0])
    a, b = m.branch_states(psi, seam, 'one_native_echo')
    before = written_chart @ psi
    after = blank_chart @ a + written_chart @ b
    np.testing.assert_allclose(pulse @ before, after, atol=5e-16)
    work = np.vdot(after, full @ after).real - np.vdot(before, full @ before).real
    reduced_work = (np.vdot(a, h0 @ a) + np.vdot(b, hw @ b) - np.vdot(psi, hw @ psi)).real
    assert abs(work - reduced_work) < 3e-14


def test_seeded_draw_ignores_last_bit_roundoff():
    data = receipt()
    assert data['sampling_weight_denominator'] == m.GRID == 2**40
    record = data['panels'][1]['arms'][0]['cases'][1]['readings'][2]
    p = np.array(record['probabilities'])
    # Symmetric ports give exactly tied probabilities. Raise the first member
    # of every tie by one unit in the last place, as a different math library
    # could. A sampler ordered by probability would reorder its draws.
    values, first, multiplicity = np.unique(p, return_index=True, return_counts=True)
    tied = first[(multiplicity > 1) & (values > 0)]
    assert len(tied) > 100
    nudged = p.copy()
    nudged[tied] = np.nextafter(p[tied], 1.)
    fresh = m.sample(p, np.random.default_rng(31))
    assert m.sample(nudged, np.random.default_rng(31))['counts'] == fresh['counts']
    assert fresh['sampling_grid_l1'] < len(p) / m.GRID + 2e-10


def test_invalid_probabilities_and_receipt_mismatch_fail():
    with pytest.raises(ArithmeticError):
        m.sample([-1e-16, 1.], np.random.default_rng(0))
    with pytest.raises(ArithmeticError):
        m.sample([.1, .2], np.random.default_rng(0))
    with pytest.raises(ArithmeticError):
        m.transition.compare({'count': 1}, {'count': 2})
