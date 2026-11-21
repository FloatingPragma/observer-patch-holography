"""Tests of the cage: the construction in integers, the record laws by three routes, the readers, the register
words against dense fermions on a small graph, and the replay with counts that count against the claims."""
import copy
import json
import math

import numpy as np
import pytest
from scipy.linalg import block_diag

from model.candidate_laws.scalar_channel import compiler
from . import construction as cn
from . import instruments as ins
from . import registers, replay


@pytest.fixture(scope='module')
def cell():
    car = cn.cell()
    Phi, h, residual = cn.modes(car, cn.ring_profile(car, cn.RING))
    return dict(car=car, D=cn.law(car), Phi=Phi, h=h, residual=residual, chi=cn.circulation(car, 0))


@pytest.fixture(scope='module')
def run():
    sim = replay.simulate()
    counts = replay.load_counts()
    read, gates = replay.observe(counts, sim)
    return dict(sim=sim, counts=counts, read=read, gates=gates)


def test_the_incidence_identity_and_the_flat_branch(cell):
    car = cell['car']
    assert cn.line_graph_identity(car) == 0
    assert cn.flat_branch(car) == dict(rank_of_Q=12, profiles=18, rings=30, rank_of_the_rings=18)
    supports = cn.small_supports(car)
    assert [len(supports[k]) for k in (1, 2, 3, 4)] == [0, 0, 0, 30]
    e = cn.ring_profile(car, cn.RING)
    assert not np.any(car['Q'] @ e) and np.array_equal(car['C'] @ car['C'].T @ e, 10 * e)
    assert cn.kept_profiles(car) == [18]


@pytest.mark.parametrize('rest, hop', [(1.0, 1 / math.sqrt(2)), (0.3, 1 / math.sqrt(2)), (1.0, 1.0), (2.5, 0.2)])
def test_the_pair_reduces_the_law_for_every_rest_level_and_hop(cell, rest, hop):
    car = cell['car']
    Phi, h, residual = cn.modes(car, cn.ring_profile(car, cn.RING), rest, hop)
    assert residual < 1e-12
    assert h[0, 1] == pytest.approx(math.sqrt(10) * hop, abs=1e-12)
    D = cn.law(car, rest, hop)
    assert np.abs(D @ Phi - Phi @ h).max() < 1e-12
    assert np.count_nonzero(Phi[:30, 0]) == 4 and np.count_nonzero(np.abs(Phi[30:, 1]) > 1e-12) == 28


def test_profiles_without_the_alternating_signs_are_not_kept(cell):
    car, D = cell['car'], cell['D']
    ring = cn.ring_profile(car, cn.RING)
    face = np.zeros(30, dtype=np.int64)
    for seam, sign in (((0, 1), 1), ((1, 2), -1), ((0, 2), 1)):
        face[car['index'][seam]] = sign
    single = np.zeros(30, dtype=np.int64)
    single[car['index'][(0, 2)]] = 1
    times = np.arange(2048) * 0.05
    kept = {}
    for name, profile in (('ring', ring), ('equal signs', np.abs(ring)), ('face', face), ('one seam', single)):
        Phi, _, residual = cn.modes(car, profile)
        kept[name] = (residual, float(ins.weight_left(D, Phi, Phi[:, 0], times).min()))
    assert kept['ring'][0] < 1e-12 and kept['ring'][1] > 1 - 1e-9
    for name in ('equal signs', 'face', 'one seam'):
        assert kept[name][0] > 0.1 and kept[name][1] < 0.25


def test_rings_on_twenty_cells():
    fed = cn.federation(0)
    assert cn.line_graph_identity(fed) == 0
    exact = cn.exact_rings(fed)
    assert len(exact) == 92
    assert sorted(set(cn.kept_profiles(fed))) == [6, 7] and sum(cn.kept_profiles(fed)) == 128
    assert cn.pair(fed, cn.ring_profile(fed, exact[0]))[3] < 1e-12
    six_ring, defect = cn.leaking_ring(fed)
    assert defect != 0 and cn.pair(fed, cn.ring_profile(fed, six_ring))[3] > 1e-2
    # numerically: the level 10 of C C^T occurs as often as there are kept profiles
    levels = np.linalg.eigvalsh(cn.seam_laplacian(fed).toarray().astype(float))
    assert int(np.sum(np.abs(levels - 10) < 1e-8)) == 128


def test_the_excitation_the_pulses_and_the_swing():
    raises, kept = ins.excitation_identities()
    assert raises < 1e-12 and kept < 1e-12
    assert max(ins.pulses()) < 1e-12
    times = np.linspace(0, 5, 501)
    computed, closed = ins.swing(times)
    assert np.abs(computed - closed).max() < 1e-12
    step = 1e-4
    values, _ = ins.swing([-step, 0.0, step])
    assert (values[2] - values[0]) / (2 * step) == pytest.approx(0, abs=1e-9)                     # it starts from rest
    assert (values[2] - 2 * values[1] + values[0]) / step ** 2 == pytest.approx(-2 * cn.SIGMA ** 2, rel=1e-5)
    G = np.outer(ins.rotation()[1], ins.rotation()[1])
    assert ins.sectors(G) == pytest.approx(dict(empty=0, reference=0, excitation=1, two=0), abs=1e-12)


def test_three_routes_to_the_law_of_the_two_readings(cell):
    occupation, _ = ins.reference(16)                       # a shorter preparation: any state of independent orbitals serves
    Phi, chi, D, h = cell['Phi'], cell['chi'], cell['D'], cell['h']
    three = np.column_stack([Phi, chi])
    S4 = block_diag(three.T @ occupation @ three, 0.0)
    H4 = block_diag(h, -cn.REST, 0.0)
    waits = np.array([0.0, 0.31, 1.7, 9.125])
    for _, potential, time in ins.PREPARATIONS:
        four = ins.joint_law(H4, np.eye(4)[0], np.eye(4)[2], 3, S4, waits, potential, time)
        assert np.abs(four.sum(axis=(-1, -2)) - 1).max() < 1e-12 and four.min() > -1e-12
        sixteen = ins.fock_joint_law(S4, waits, potential, time)
        every = ins.joint_law(block_diag(D, 0.0), np.r_[Phi[:, 0], 0.0], np.r_[chi, 0.0], 150, block_diag(occupation, 0.0), waits, potential, time)
        assert np.abs(four - sixteen).max() < 1e-12
        assert np.abs(four - every).max() < 1e-10


def test_the_readers_return_the_law_from_its_own_probabilities(cell):
    lower = ins.rotation()[0]
    G = np.outer(lower, lower)
    H4, S4 = block_diag(cell['h'], -cn.REST, 0.0), block_diag(G, 1.0, 0.0)
    shots = 10 ** 6
    for waits in ins.schedules():
        law = ins.joint_law(H4, np.eye(4)[0], np.eye(4)[2], 3, S4, waits, -2 * cn.ENERGY, cn.T_COHERENT)
        phase, record = ins.read_clock(shots * law, shots)
        line = ins.read_line(phase, record)
        assert line['frequency'] == pytest.approx(2 * cn.ENERGY, abs=2e-6) and line['share_of_the_variance'] > 1 - 1e-9
        # the phase of the reference is the wait up to the finite pulses: (g^2 + 1/4) / g^2 = 1 + 6e-5 at g = 64
        assert np.abs(phase - waits).max() < 5e-3
    own = ins.reference_law(H4, np.eye(4)[2], 3, S4, np.arange(ins.REFERENCE_QUERIES) * cn.STEP)
    phases = ins.read_reference(shots * np.stack([1 - own, own], axis=-1), shots)
    dial = ins.read_dial(shots * ins.pointers(ins.dial_law(G)), shots, [phases[q] for q in ins.DIAL_WAITS])
    assert dial['coupling'] == pytest.approx(cn.SIGMA, abs=1e-4)
    for name, potential, time in ins.PREPARATIONS:
        expected = dict(coherent=(1 + cn.REST / cn.ENERGY) / 2, sharp=1.0, none=0.0)[name]
        coordinates = np.array([[replay._pair_before_the_analyzer(cell['h'], G, potential, time, q) for q in ins.WAITS]])
        probabilities = np.einsum('ij,pqj->pqi', ins.effects(cn.SIGMA), coordinates)
        read, condition = ins.read_sectors(shots * ins.pointers(probabilities), shots, cn.SIGMA)
        assert [s['excitation'] for s in read[0]] == pytest.approx([expected] * len(ins.WAITS), abs=1e-9)
        assert 4.6 < condition < 4.7


def test_a_reader_refuses_a_reference_that_may_have_missed_a_turn(cell):
    lower = ins.rotation()[0]
    H4, S4 = block_diag(cell['h'], -cn.REST, 0.0), block_diag(np.outer(lower, lower), 1.0, 0.0)
    slow = ins.joint_law(H4, np.eye(4)[0], np.eye(4)[2], 3, S4, np.arange(8) * 0.01, -2 * cn.ENERGY, cn.T_COHERENT)
    fast = ins.joint_law(H4, np.eye(4)[0], np.eye(4)[2], 3, S4, np.arange(8) * 0.5, -2 * cn.ENERGY, cn.T_COHERENT)
    for law in (slow, fast):
        with pytest.raises(ValueError, match='outside the declared bounds'):
            ins.read_clock(10 ** 6 * law, 10 ** 6)


def dense_pauli(word, qubits):
    phase, x, z = word
    out = np.zeros((1 << qubits, 1 << qubits), complex)
    for column in range(1 << qubits):
        out[column ^ x, column] = 1j ** phase * (-1) ** (column & z).bit_count()
    return out


def test_the_words_of_the_two_operators_against_dense_fermions():
    """Seven modes joined as the tree of the ring: seams 0 to 3, the first fuel 4, the diagonal seam 5, the second fuel 6."""
    graph = compiler.encode_graph(7, ((0, 4), (1, 6), (2, 6), (3, 4), (4, 5), (5, 6)))
    qubits, modes = graph['M'], 7
    identity = np.eye(1 << qubits)
    numbers = [dense_pauli(word, qubits) for word in graph['B']]
    projector = identity.copy()
    for number in numbers[:-1]:
        projector = projector @ ((identity + number) / 2)
    vacuum = projector[:, np.argmax(np.linalg.norm(projector, axis=0))]
    vacuum = vacuum / np.linalg.norm(vacuum)
    flips = [dense_pauli(word, qubits) for word in compiler.logical_pauli_generators(graph)]
    a = ins.fock_operators(modes)
    gammas = [x + x.conj().T for x in a]
    signs = [np.eye(1 << modes) - 2 * x.conj().T @ x for x in a]
    fermion_flips = []
    for i in range(modes - 1):
        operator = -1j * gammas[i] @ gammas[modes - 1]
        for j in range(i + 1, modes - 1):
            operator = operator @ signs[j]
        fermion_flips.append(operator)
    code, fock = [], []
    for word in range(1 << (modes - 1)):
        left, right = vacuum.copy(), np.eye(1 << modes, dtype=complex)[:, 0]
        for bit in range(modes - 1):
            if word & (1 << bit):
                left, right = flips[bit] @ left, fermion_flips[bit] @ right
        code.append(left)
        fock.append(right)
    V, W = np.column_stack(code), np.column_stack(fock)
    assert np.abs(V.conj().T @ V - np.eye(64)).max() < 1e-12 and np.abs(W.conj().T @ W - np.eye(64)).max() < 1e-12
    for i in range(modes):                                   # the map carries numbers to numbers
        assert np.abs(V.conj().T @ numbers[i] @ V - W.conj().T @ signs[i] @ W).max() < 1e-12
    seams, pattern, fuels, diagonal = [0, 1, 2, 3], (1, -1, 1, -1), [4, 6, 6, 4], 5
    dense = lambda poly: sum(complex(*c) * dense_pauli((0, x, z), qubits) for (x, z), c in poly.items())
    profile = sum(s * a[k] for s, k in zip(pattern, seams)) / 2                   # removes an object from e
    Q = registers.profile_number(graph, seams, pattern, fuels, diagonal)
    assert np.abs(V.conj().T @ dense(Q) @ V / 8 - W.conj().T @ profile.conj().T @ profile @ W).max() < 1e-12
    K = registers.converter(graph, seams, pattern, fuels, diagonal, 4)
    exchange = profile.conj().T @ a[4] + a[4].conj().T @ profile
    assert np.abs(V.conj().T @ dense(K) @ V / 4 - W.conj().T @ exchange @ W).max() < 1e-12
    # another pattern of signs is another operator: the words carry the signs
    other = registers.profile_number(graph, seams, (1, 1, 1, 1), fuels, diagonal)
    assert np.abs(dense(other) - dense(Q)).max() > 0.5


def test_the_two_operators_on_the_registers_of_one_cell():
    out = registers.run()
    assert (out['registers'], out['modes'], out['loops']) == (432, 174, 260)
    number, converter = out['profile_number'], out['converter']
    assert (number['words'], number['largest_weight'], number['hermitian']) == (17, 15, True)
    assert number['terms_left_in_Q_squared_minus_8_Q'] == 0 and number['pairs_with_the_loops'] == 4420
    assert number['pairs_with_the_loops_that_do_not_commute'] == 0
    assert number['words_that_do_not_commute_with_a_fuel_or_reference_number'] == 0
    assert number['terms_left_in_the_commutator_with_the_number_of_objects'] == 0
    assert (converter['words'], converter['largest_weight'], converter['hermitian']) == (8, 24, True)
    assert converter['terms_left_in_K_squared_minus_2Q_minus_8R_plus_2QR'] == 0
    assert converter['terms_left_in_the_commutator_with_Q_plus_4R'] == 0 and converter['pairs_with_the_loops'] == 2080
    assert converter['pairs_with_the_loops_that_do_not_commute'] == 0
    for operator, most in ((number, 4), (converter, 3)):
        assert operator['stars_of_all_words'] == [0, 2, 4, 6] and operator['most_stars_in_one_word'] == most


def test_the_stored_receipts(run):
    saved = json.loads(replay.RESULTS.read_text(encoding='utf-8'))
    assert saved['counts_sha256'] == replay.digest(replay.COUNTS)
    assert all(run['gates'].values())
    assert set(run['sim']['residuals']) == set(replay.BOUNDS)
    assert all(run['sim']['residuals'][name] <= bound for name, bound in replay.BOUNDS.items())
    fresh = replay._plain(run['read'])
    assert replay.compare({k: saved[k] for k in fresh}, fresh) == []
    assert replay.compare(saved['construction'], replay._plain(run['sim']['construction'])) == []
    assert replay.compare(saved['object'], replay._plain(run['sim']['object'])) == []
    # the numbers of the manual
    equal, unequal = (saved['clock'][k]['line'] for k in ('steps_of_equal_length', 'steps_of_unequal_length'))
    for line in (equal, unequal):
        assert abs(line['frequency'] - 2 * math.sqrt(6)) < 4 * line['standard_error']
    assert abs(saved['dial']['level_difference'] - 2 * math.sqrt(6)) < 4 * saved['dial']['level_difference_standard_error']
    assert saved['object']['excitation_in_the_state'] == pytest.approx(dict(coherent=0.705931828, sharp=0.999954092, none=0.000004150), abs=5e-10)
    assert [f['rings_on_ports_with_five_seams'] for f in saved['construction']['federations']] == [92, 312, 1272]
    assert [f['profiles_that_keep_the_level'] for f in saved['construction']['federations']] == [128, 504, 2024]


def test_a_changed_result_is_found():
    saved = json.loads(replay.RESULTS.read_text(encoding='utf-8'))
    changed = copy.deepcopy(saved)
    changed['clock']['steps_of_equal_length']['line']['frequency'] += 1e-4
    assert len(replay.compare(saved, changed)) == 1
    changed = copy.deepcopy(saved)
    changed['construction']['federations'][2]['rings_on_ports_with_five_seams'] += 1
    assert len(replay.compare(saved, changed)) == 1
    changed = copy.deepcopy(saved)
    del changed['gates']['the_line_and_the_dial_agree']
    assert replay.compare(saved, changed)


def test_counts_that_count_against_the_claims(run):
    sim, counts = run['sim'], run['counts']
    # the reader is given the counts of the run without the pulse
    swapped = dict(counts, clock_prepared=counts['clock_no_preparation'])
    _, gates = replay.observe(swapped, sim)
    assert not gates['one_line_carries_the_seam_record'] and not gates['the_stored_counts_are_draws_from_the_record_laws']
    # the counts of the law and of the changed law change places
    swapped = dict(counts, persistence=counts['persistence_changed_law'], persistence_changed_law=counts['persistence'])
    _, gates = replay.observe(swapped, sim)
    assert not gates['the_changed_law_is_told_apart_by_its_counts'] and not gates['the_sharp_excitation_is_read_at_every_wait']
    # the changed law's clock counts in place of the law's: no single line
    swapped = dict(counts, clock_prepared=np.stack([counts['clock_changed_law']] * 2))
    _, gates = replay.observe(swapped, sim)
    assert not gates['one_line_carries_the_seam_record']
    # counts that are too regular to be draws: the expected counts themselves
    exact = dict(counts, dial=np.rint(ins.SHOTS['dial'] * sim['laws']['dial']).astype(np.int64))
    _, gates = replay.observe(exact, sim)
    assert not gates['the_stored_counts_are_draws_from_the_record_laws']
    # a record with the converter off
    planted = dict(counts, persistence_no_converter=counts['persistence'])
    _, gates = replay.observe(planted, sim)
    assert not gates['no_record_with_the_converter_off']


def test_a_dial_whose_unit_is_off_moves_the_level_difference(run):
    counts = run['counts']
    own = ins.read_reference(counts['reference'], ins.SHOTS['reference'])
    phases = [float(own[q]) for q in ins.DIAL_WAITS]
    good = ins.read_dial(counts['dial'], ins.SHOTS['dial'], phases)
    off = ins.read_dial(counts['dial'], ins.SHOTS['dial'], phases, dials=1.02 * ins.DIALS)
    assert abs(good['level_difference'] - 2 * cn.ENERGY) < 5e-4 < abs(off['level_difference'] - 2 * cn.ENERGY)
