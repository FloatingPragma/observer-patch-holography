from collections import Counter
from fractions import Fraction as F
from itertools import product

import numpy as np

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from model.experiments.interval_echo_clock.publication import SCHEDULE

from . import archive as a


def test_actual_word_lotteries_preserve_both_singles():
    words = a.posted_words()
    assert tuple(Counter(words)[p] for p in a.PAIRS) == (12, 10, 3, 8)
    hist = {}
    for name, weights in a.address_lotteries(words).items():
        assert len(words) == len(weights) == 33 and sum(weights) == 1 and min(weights) >= 0
        hist[name] = tuple(sum(w for w, x in zip(weights, words, strict=True) if x == p)
                           for p in a.PAIRS)
    assert hist["ordinary"] == tuple(F(n, 33) for n in (4, 18, 11, 0))
    assert hist["held"] == tuple(F(n, 33) for n in (12, 10, 3, 8))
    for law in (hist["ordinary"], hist["held"]):
        assert law[2]+law[3] == F(11, 33)
        assert law[1]+law[3] == F(18, 33)
    assert a.source_densities() == {"blank": (0., 0.), "ordinary": (1/3, 0.),
                                    "held": (1/3, 8/33)}


def test_live_archive_after_both_native_readers(monkeypatch):
    # Tap the *return value* of the old coherent clock's already executed bit
    # read. The old native registers, random stream and gate schedule are not
    # passed to the recorder. Exercise a feedback choice of terminal analyzer.
    car = federation(1)
    old_read = clock.SeamOscillator.on_tick
    for seed, expected in [(31, (12, 10, 3, 8)), (2026112607, (9, 11, 7, 6))]:
        initial = clock.material_replay(car, SCHEDULE, age_bits=3, seed=seed, analyzer=0)
        setting = sum(x["tick"] for x in initial["local_events"][:5]) % 2
        baseline = clock.material_replay(car, SCHEDULE, age_bits=3, seed=seed, analyzer=setting)
        posted = []
        def copy_after(self, *args):
            bit = old_read(self, *args)
            posted.append(bit)
            return bit
        with monkeypatch.context() as patch:
            patch.setattr(clock.SeamOscillator, "on_tick", copy_after)
            observed = clock.material_replay(car, SCHEDULE, age_bits=3, seed=seed, analyzer=setting)
        assert observed == baseline
        words = tuple(zip(posted[::2], posted[1::2], strict=True))
        assert words == tuple((x["tick"], x["home_tick"]) for x in baseline["local_events"])
        assert tuple(Counter(words)[p] for p in a.PAIRS) == expected
        if seed == 31:
            assert words == a.posted_words()


def test_all_transcript_branches_with_an_entangled_untouched_reference():
    # The two-old-qubit blocks retain entanglement and arbitrary phases.
    # After publication, an old observer coherently rotates its own qubit
    # CONDITIONAL on the first transcript, then records a second outcome.
    bell = np.array([1., 0, 0, 1j])/np.sqrt(2)
    phi = np.array([0., 1, 1, 0])/np.sqrt(2)
    original = {0: .4*np.outer(bell, bell.conj()), 1: .6*np.outer(phi, phi)}
    copies = a.refine_published_blocks(original, {0: {"a": F(1, 3), "b": F(2, 3)},
                                                  1: {"a": F(3, 5), "b": F(2, 5)}})
    H = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    for first, last in product(range(2), repeat=2):
        P = np.diag([int(j == last) for j in range(2)])
        A = np.kron(P@(np.eye(2) if first == 0 else H), np.eye(2))
        expected = A@original[first]@A.conj().T
        observed = sum(A@copies[first, label]@A.conj().T for label in ("a", "b"))
        np.testing.assert_allclose(observed, expected, rtol=0, atol=1e-15)


def test_copy_before_coherent_reversal_changes_a_native_bit():
    # Isometric pre-read of an uncommitted binary path in an external register.
    # H followed by a Z read reverses the native path superposition exactly.
    h = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    plus = h[:, 0]
    before_copy = np.diag(plus)  # CNOT from path to fresh |0> yields joint 00+11.
    path = before_copy@before_copy.conj().T
    assert np.isclose(abs((h@plus)[0])**2, 1.)
    assert np.isclose((h@path@h.T)[0, 0], .5)
