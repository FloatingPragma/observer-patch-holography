"""Mixed ternary loads, coherent sign inputs and fixed material pointer checks."""

from itertools import product
import json

import numpy as np
from scipy.linalg import expm

from model.executable_model import cell
from . import fixed_pointer as pointer, sequential


def test_native_ternary_primitive_from_full_local_operator():
    law = pointer.candidate.native_law
    a, b = 0, 10
    # One physical source column for every ternary endpoint pair, including 02/20.
    selected = [next(i for i, state in enumerate(pointer.STATES)
                     if state[pointer.PORTS.index(a)] == na
                     and state[pointer.PORTS.index(b)] == nb)
                for na, nb in law.SRC]
    start = np.eye(pointer.DIM)[:, selected]
    h0, h1 = law.local_term_parts()
    full = h0 + pointer.candidate.KAPPA * h1
    exact_local = expm(-1j * np.pi / 37 * full)[:, :9]
    expected = np.zeros_like(start, dtype=complex)
    for col, original in enumerate(selected):
        for out_index, (na, nb) in enumerate(law.SRC):
            amplitude = exact_local[out_index, col]
            if abs(amplitude) < 1e-14:
                continue
            word = list(pointer.STATES[original])
            word[pointer.PORTS.index(a)], word[pointer.PORTS.index(b)] = na, nb
            expected[pointer.INDEX[tuple(word)], col] += amplitude
    actual = pointer.apply([pointer.pulse("contact", a, b, pointer.Fraction(1,37))], start)
    np.testing.assert_allclose(actual, expected, atol=3e-13, rtol=0)
    assert np.linalg.norm(exact_local[9:]) < 3e-13


def test_complete_multiset_and_exact_diagonal_refocusing():
    assert pointer.DIM == 560
    assert all(tuple(np.bincount(state, minlength=3)) == (3,2,3)
               and sum(state) == 8 for state in pointer.STATES)
    assert sum(sum(row) == 8 for row in product(range(3), repeat=8)) == 1107
    word = pointer.diagonal_word(0,10,pointer.Fraction(2,65))
    actual = pointer.apply(word, np.eye(pointer.DIM))
    a, b = pointer.LOADS[:,0], pointer.LOADS[:,pointer.PORTS.index(10)]
    v = 8*(a-b)**2 + (a != b)/2
    expected = np.diag(np.exp(-2j*np.pi*v/65))
    np.testing.assert_allclose(actual, expected, atol=6e-13, rtol=0)
    required = {(0,1),(6,10),(0,10),(4,5),(0,5),(0,7),(0,11)}
    assert all(b in cell().neighbours[a] for a,b in required)


def test_fixed_pointer_copy_all_signs_addresses_and_pointer_words():
    phase = np.exp(-1j * 133 * np.pi / 260)
    for p in (0,1):
        for bits in product((0,1), repeat=2):
            start = pointer.embedding(bits)
            actual = pointer.apply(pointer.copy_word(p), start)
            target = np.zeros_like(start)
            for col in range(4):
                changed = list(bits)
                changed[p] ^= int(col % 2 == 0)
                target[:,col] = phase * pointer.embedding(tuple(changed))[:,col]
            np.testing.assert_allclose(actual,target,atol=7e-13,rtol=0)
    # Pointer preparation is one sharp product state in all four logical inputs.
    for col in range(4):
        state = pointer.STATES[np.argmax(abs(pointer.embedding()[:,col]))]
        assert [state[pointer.PORTS.index(p)] for pair in pointer.POINTERS for p in pair] == [2,0,2,0]
        assert state[pointer.PORTS.index(pointer.HOST)] == 1


def test_balanced_charge_reference_and_cross_sign_instrument_maps():
    # Four-dimensional Choi input includes every cross-sign matrix unit. The
    # physical balancing reference is part of the encoded input, not a pointer.
    bell = np.eye(4).reshape(-1)/2
    phase = np.exp(-1j*133*np.pi/130)
    for swapped in (False,True):
        actual = pointer.logical_maps(swapped)
        expected = pointer.target_maps(swapped)
        for label, k in actual.items():
            a = np.kron(k,np.eye(4)) @ bell
            b = np.kron(phase*expected[label],np.eye(4)) @ bell
            np.testing.assert_allclose(np.outer(a,a.conj()),np.outer(b,b.conj()),
                                       atol=7e-13,rtol=0)
        maps = pointer.branch_maps(swapped)
        np.testing.assert_allclose(sum(k.conj().T@k for k in maps),np.eye(4),
                                   atol=7e-13,rtol=0)
    corrected = pointer.coherence_control(True)
    bare = pointer.coherence_control(False)
    np.testing.assert_allclose(corrected,[1,0,0],atol=7e-13,rtol=0)
    np.testing.assert_allclose(bare,[0,0,1],atol=7e-13,rtol=0)


def test_complete_banks_and_mixed_sign_escape_control():
    saved = json.loads(pointer.RECEIPT.read_text())
    sequential.compare(saved,pointer.report())
    assert len(saved["banks"]) == 20
    assert sum(row["shots"] for row in saved["banks"]) == 81920
    for row in saved["banks"]:
        bins = row["all_243_pointer_host_counts"]
        assert len(bins) == 243 and sum(bins) == row["shots"]
        if row["arm"] == "repeatable":
            assert bins[141]+bins[101] == row["shots"]
    leak = saved["controls"]["incomplete_refocus_escape_by_logical_input"]
    assert min(leak[0],leak[2]) > 1e-5
    assert max(leak[1],leak[3]) < 1e-24
