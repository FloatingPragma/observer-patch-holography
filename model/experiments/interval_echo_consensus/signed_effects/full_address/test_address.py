"""Exact local algebra, routed factorization and complete address-instrument tests."""

import json
from itertools import product

import numpy as np
import sympy as s

from model.executable_model import federation
from . import measure


def test_exact_ternary_diagonal_and_common_copy_phase():
    words = tuple(product(range(3),repeat=2))
    swap = s.eye(9)[:,[words.index((b,a)) for a,b in words]]
    v = s.diag(*[8*(a-b)**2+s.Rational(int(a!=b),2) for a,b in words])
    h = (s.eye(9)-swap)/2+s.diag(*[8*(a-b)**2 for a,b in words])
    z = s.diag(*[(-1)**(a*a) for a,b in words])
    t = s.diag(*[s.exp(-s.I*s.pi*a*a/4).expand(complex=True) for a,b in words])
    first = h+z*h*z
    second = t.conjugate().T*first*t
    assert h*(z*h*z) == (z*h*z)*h
    assert s.simplify(first*second-second*first) == s.zeros(9)
    assert s.simplify(first+second-4*v) == s.zeros(9)
    hadamard = s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    x = s.Matrix([[0,1],[1,0]])
    common = s.exp(-s.I*133*s.pi/260)
    for load in range(3):
        cost = lambda p: 8*(load-p)**2+s.Rational(int(load!=p),2)
        action = hadamard*s.diag(*[s.exp(-2*s.I*s.pi*cost(p)/65) for p in (0,2)])*hadamard
        correction = s.exp(-s.I*s.pi*load**2/4-s.I*133*s.pi*cost(1)/2210)
        target = x if load!=1 else s.eye(2)
        assert s.simplify(action*correction/common-target) == s.zeros(2)


def test_actual_routes_reference_degrees_and_tensor_slot_conjugation():
    board = federation(0)
    calls = measure.schedule()
    assert len(calls) == 20
    assert sum(len(call["route"])-1 for call in calls) == 98
    for call in calls:
        path = call["route"]
        assert path[0] == call["queried_address"] and path[-1] == measure.SITES[call["site"]]["source"]
        assert set(path).isdisjoint(measure.RESERVED)
        assert all(b in board.neighbours[a] for a,b in zip(path,path[1:]))
        # Unique tensor-slot labels prove the permutation, rather than sampling
        # one set of physical ternary loads and inferring the operator action.
        slots = np.arange(board.n)
        for a,b in zip(path,path[1:]):
            slots[a],slots[b] = slots[b],slots[a]
        assert slots[path[-1]] == path[0]
        for reserved in measure.RESERVED:
            assert slots[reserved] == reserved
        for a,b in reversed(list(zip(path,path[1:]))):
            slots[a],slots[b] = slots[b],slots[a]
        np.testing.assert_array_equal(slots,np.arange(board.n))
    word = measure.pulse_word(calls)
    assert len(word) == 736
    for kind,a,b,_ in word:
        if kind == "reference":
            assert board.degree[a] == 5
        else:
            assert b in board.neighbours[a]
    local = measure.local_check()
    assert local["local_copy_operator_error"] < 2e-12
    assert local["native_pi_swap_operator_error"] < 2e-12


def test_all_physical_words_restore_and_encode_the_address():
    calls = measure.schedule()
    for sign,address,prior in product((-1,1),range(12),range(16)):
        output,phase = measure.execute_basis(calls,sign,address,prior)
        np.testing.assert_array_equal(output,measure.initial_word(sign,address,prior^address))
        assert measure.read_word(output) == prior^address
        np.testing.assert_allclose(phase,measure.COMMON_PHASE**20,atol=1e-14,rtol=0)
        assert sum(output) == 240


def test_rank_two_branches_references_repeatability_and_coherent_undo():
    encoder,_ = measure.encoder(measure.schedule())
    branches = [encoder[:,label,:] for label in range(17)]
    phase = measure.COMMON_PHASE**20
    np.testing.assert_allclose(sum(k.conj().T@k for k in branches),np.eye(24),atol=3e-13)
    for label in range(12):
        p = np.zeros((24,24),complex)
        p[label,label] = p[label+12,label+12] = 1
        np.testing.assert_allclose(branches[label],phase*p,atol=3e-13,rtol=0)
        assert np.linalg.matrix_rank(branches[label]) == 2
        for other in range(12):
            np.testing.assert_allclose(branches[other]@branches[label],
                                       phase*branches[label] if other==label else 0,
                                       atol=3e-13,rtol=0)
        balanced = np.zeros(24,complex)
        balanced[label],balanced[label+12] = 1/np.sqrt(2),1j/np.sqrt(2)
        np.testing.assert_allclose(branches[label]@balanced,phase*balanced,atol=3e-13,rtol=0)
    assert sum(np.linalg.norm(k) for k in branches[12:]) == 0
    rng = np.random.default_rng(925061)
    state = rng.normal(size=(24,3))+1j*rng.normal(size=(24,3))
    state /= np.linalg.norm(state)
    encoded = np.einsum("dbl,lr->dbr",encoder,state)
    ref = np.einsum("dbr,dbs->rs",encoded,encoded.conj())
    np.testing.assert_allclose(ref,state.T@state.conj(),atol=3e-13,rtol=0)
    # The same isometry's adjoint exactly uncomputes an unread record.
    matrix = encoder.reshape(24*17,24)
    np.testing.assert_allclose(matrix.conj().T@matrix@state,state,atol=3e-13,rtol=0)


def test_complete_banks_controls_and_invalid_reader():
    saved = json.loads(measure.RECEIPT.read_text())
    measure.sequential.compare(saved,measure.report())
    assert len(saved["banks"]) == 24
    assert sum(row["shots"] for row in saved["banks"]) == 98304
    for row in saved["banks"]:
        bins = row["all_17_record_counts"]
        assert len(bins) == 17 and sum(bins) == row["shots"]
        assert bins[row["prepared_address"]] == row["shots"]
        assert sum(bins[12:]) == 0
    for control in ("omit_bit3","wrong_route"):
        assert saved["controls"][control]["vacancy_records"][0] > 1-3e-13
        assert saved["controls"][control]["excess_records"][0] > 1-3e-13
    assert abs(saved["coherent_native_analyzer_return"]-1) < 3e-13
    assert abs(saved["pointer_dephased_native_analyzer_return"]-0.5) < 3e-13
    bad = measure.initial_word(-1,8)
    bad[measure.SITES[0]["host"]] = 0
    assert measure.read_word(bad) == 16
