"""Independent dense/reference checks for the routed measurement instrument."""
import itertools
import numpy as np
import pytest

from .compiler import Tableau, build, fundamental_cycles
from .routed import clifford, controlled_factor, measure_loop, run


def dense_pauli(P, n):
    p, x, z = P
    result = np.zeros((2**n, 2**n), complex)
    for j in range(2**n):
        result[j ^ x, j] = 1j**p * (-1)**((j & z).bit_count())
    return result


def dense_gate(name, a, b, n):
    dim = 2**n
    if name == 'SWAP':
        out = np.zeros((dim, dim), complex)
        for j in range(dim):
            moved = j ^ ((1 << a) | (1 << b)) if ((j >> a) ^ (j >> b)) & 1 else j
            out[moved, j] = 1
        return out
    if name == 'CX':
        out = np.zeros((dim, dim), complex)
        for j in range(dim):
            out[j ^ ((1 << b) if j & (1 << a) else 0), j] = 1
        return out
    if name == 'H':
        return (dense_pauli((0, 1 << a, 0), n) +
                dense_pauli((0, 0, 1 << a), n)) / np.sqrt(2)
    return np.diag([(1j if name == 'S' else -1j) if j & (1 << a) else 1
                    for j in range(dim)])


def tableau_density(tab):
    """Independent joint density from the complete commuting stabilizers."""
    identity = np.eye(2**tab.q)
    rho = identity.copy()
    for row in tab.rows:
        rho = rho @ (identity + dense_pauli(row, tab.q)) / 2
    return rho


@pytest.mark.parametrize('name,a,b', [('H', 1, None), ('S', 0, None),
                                     ('SDG', 2, None), ('CX', 2, 0),
                                     ('SWAP', 0, 2)])
def test_clifford_all_paulis_independent_dense(name, a, b):
    U = dense_gate(name, a, b, 3)
    for x, z in itertools.product(range(8), repeat=2):
        P = ((x & z).bit_count() % 4, x, z)
        tab = Tableau(3)
        tab.rows = [P]
        clifford(tab, name, a, b)
        np.testing.assert_allclose(dense_pauli(tab.rows[0], 3),
                                   U @ dense_pauli(P, 3) @ U.conj().T, atol=1e-14)


def test_local_pointer_kraus_for_all_two_bit_paulis_with_reference():
    # Build the literal gates independently on two data bits and an ancilla.
    # Each resulting Kraus matrix equals (I +/- P)/2; equality as operators
    # also proves equality after tensoring with any inaccessible reference.
    for x, z in itertools.product(range(4), repeat=2):
        for negative in (False, True):
            p = ((x & z).bit_count() + 2 * negative) % 4
            U = dense_gate('H', 2, None, 3)
            if negative:
                U = dense_pauli((0, 0, 4), 3) @ U
            for k in range(2):
                if not (x | z) & (1 << k):
                    continue
                # Extract the implementation's Clifford circuit operator by
                # testing its conjugation against independently formed CY/CZ.
                factor = dense_pauli((((x >> k) & (z >> k) & 1),
                                      ((x >> k) & 1) << k,
                                      ((z >> k) & 1) << k), 2)
                ctrl = np.block([[np.eye(4), np.zeros((4, 4))],
                                 [np.zeros((4, 4)), factor]])
                tab = Tableau(3)
                for xx, zz in itertools.product(range(8), repeat=2):
                    P = ((xx & zz).bit_count() % 4, xx, zz)
                    tab.rows = [P]
                    controlled_factor(tab, 2, k, (x >> k) & 1, (z >> k) & 1)
                    np.testing.assert_allclose(dense_pauli(tab.rows[0], 3),
                                               ctrl @ dense_pauli(P, 3) @ ctrl.conj().T,
                                               atol=1e-14)
                U = ctrl @ U
            U = dense_gate('H', 2, None, 3) @ U
            P = dense_pauli((p, x, z), 2)
            for outcome, sign in enumerate((1, -1)):
                np.testing.assert_allclose(U[outcome*4:(outcome+1)*4, :4],
                                           (np.eye(4) + sign*P)/2, atol=1e-14)


@pytest.mark.parametrize('pauli', [(0, 3, 0), (3, 3, 1)])
def test_routed_pointer_preserves_entangled_reference(pauli):
    # Bit 2 is an untouched reference, entangled with data bit 0. The measured
    # Pauli acts only on bits 0 and 1, across two different stars.
    graph = dict(M=3, N=2, owner=[0, 1, 0],
                 data={'edges': [(0, 1)]})
    initial = dense_gate('CX', 0, 2, 3) @ dense_gate('H', 0, None, 3)[:, 0]
    for sign in (-1, 1):
        tab = Tableau(5)
        clifford(tab, 'H', 0)
        clifford(tab, 'CX', 0, 2)
        rec = measure_loop(tab, graph, 0, pauli, sign, omit_correction=True)
        expected = (np.eye(8) + sign*dense_pauli(pauli, 3)) @ initial / 2
        probability = np.vdot(expected, expected).real
        assert rec['probability'] == pytest.approx(probability)
        assert probability == pytest.approx(.5)
        joint = np.zeros(32, complex)
        joint[:8] = expected / np.sqrt(probability)
        np.testing.assert_allclose(tableau_density(tab), np.outer(joint, joint.conj()),
                                   atol=1e-14)
        assert rec['route'] == [0, 1, 0]


@pytest.mark.parametrize('fault', ['missing_ack', 'wrong_ack'])
def test_failed_ack_retains_partial_entangled_state(fault):
    graph = dict(M=2, N=2, owner=[0, 1], data={'edges': [(0, 1)]})
    tab = Tableau(4)
    clifford(tab, 'H', 0)
    rec = measure_loop(tab, graph, 0, (0, 2, 1), 1, fault=fault)
    # Failure occurs after controlled Z at star 0, before the first transfer.
    initial = dense_gate('H', 2, None, 4) @ dense_gate('H', 0, None, 4)[:, 0]
    controlled_z = np.diag([-1 if (j & 1 and j & 4) else 1 for j in range(16)])
    expected = controlled_z @ initial
    np.testing.assert_allclose(tableau_density(tab), np.outer(expected, expected.conj()),
                               atol=1e-14)
    assert rec['status'] == 'failed' and rec['reason'] == fault
    assert rec['probability'] == 1 and rec['outcome'] is None
    assert not rec['reset'] and not rec['corrected']
    assert rec['route'] == [0]
    assert tab.expectation((0, 0, 4)) == 0  # scratch has not been reset


def test_dirty_entangled_scratch_retains_every_diagnostic_branch():
    graph = dict(M=2, N=2, owner=[0, 1], data={'edges': [(0, 1)]})
    branches = [([-1, 1], 1, .5), ([1, 1], 1, .25), ([1, 1], -1, .25)]
    total = 0
    for diagnostic_signs, sign, probability in branches:
        tab = Tableau(4)
        clifford(tab, 'H', 0)
        clifford(tab, 'CX', 0, 2)
        rec = measure_loop(tab, graph, 0, (0, 1, 0), sign,
                           diagnostic_signs=diagnostic_signs, omit_correction=True)
        expected = np.zeros(16, complex)
        if diagnostic_signs[0] < 0:
            expected[5] = 1  # actual collapsed |data0=1, scratch0=1>
            assert rec['reason'] == 'dirty_ancilla' and not rec['reset']
            assert rec['outcome'] is None and len(rec['diagnostics']) == 1
        else:
            expected[0], expected[1] = 1/np.sqrt(2), sign/np.sqrt(2)
            assert rec['status'] == 'admitted' and len(rec['diagnostics']) == 2
        assert rec['probability'] == probability
        np.testing.assert_allclose(tableau_density(tab), np.outer(expected, expected.conj()),
                                   atol=1e-14)
        total += rec['probability']
    assert total == 1


def test_diagnostic_tree_includes_later_dirty_scratch():
    graph = dict(M=2, N=2, owner=[0, 1], data={'edges': [(0, 1)]})
    branches = [([-1, 1], 1, .5), ([1, -1], 1, .25),
                ([1, 1], 1, .125), ([1, 1], -1, .125)]
    total = 0
    for signs, pointer, expected in branches:
        tab = Tableau(4)
        clifford(tab, 'H', 2)
        clifford(tab, 'H', 3)
        rec = measure_loop(tab, graph, 0, (0, 1, 0), pointer,
                           diagnostic_signs=signs, omit_correction=True)
        assert rec['probability'] == expected
        if -1 in signs:
            first = signs.index(-1)
            assert rec['diagnostics'][-1]['star'] == first
            assert tab.expectation((0, 0, 1 << (2 + first))) == -1
        total += rec['probability']
    assert total == 1


def test_seeded_diagnostics_are_reproducible():
    graph = dict(M=2, N=2, owner=[0, 1], data={'edges': [(0, 1)]})
    results = []
    for _ in range(2):
        tab = Tableau(4)
        clifford(tab, 'H', 2)
        rec = measure_loop(tab, graph, 0, (0, 1, 0), 1, diagnostic_seed=42)
        results.append((rec, tab.rows))
    assert results[0] == results[1]


def test_complete_cell_and_failure_control():
    result = run()
    assert result['retained_sign_bits'] == 259
    assert result['transport_hops'] == 941
    assert result['diagnostic_measurements'] == 12 * 259
    assert result['pointer_measurements'] == 259
    expected_signs = np.random.default_rng(result['seed']).choice((-1, 1), size=259)
    assert [r['outcome'] for r in result['retained_sign_records']] == expected_signs.tolist()
    assert result['no_feedback_control']['detected_sign'] == -1
    assert {r['reason'] for r in result['failure_controls']} == {'missing_ack', 'wrong_ack'}
    assert all(not r['reset'] and r['outcome'] is None for r in result['failure_controls'])
    assert all(len(r['halted_state']['stabilizers']) == 444
               for r in result['failure_controls'])
    assert result['dirty_ancilla_control']['probability'] == 1
    assert result['dirty_ancilla_control']['reason'] == 'dirty_ancilla'
    assert len(result['final_state']['stabilizers']) == 444


def test_routes_do_not_measure_data_directly(monkeypatch):
    graph = build()
    tab = Tableau(graph['M'] + graph['N'])
    chord, _, P = fundamental_cycles(graph)[0][0]
    original = Tableau.measure
    def only_ancilla(self, operator, sign):
        assert operator[1] == 0
        assert operator[2].bit_count() == 1
        assert operator[2].bit_length() - 1 >= graph['M']
        return original(self, operator, sign)
    monkeypatch.setattr(Tableau, 'measure', only_ancilla)
    measure_loop(tab, graph, chord, P, -1)


def test_nonhermitian_measurement_is_rejected():
    graph = build()
    chord, _, P = fundamental_cycles(graph)[0][0]
    with pytest.raises(ValueError, match='Hermitian'):
        measure_loop(Tableau(444), graph, chord, ((P[0]+1) % 4, P[1], P[2]), 1)
