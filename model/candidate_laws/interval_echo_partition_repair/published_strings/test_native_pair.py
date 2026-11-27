import itertools
import math

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo_partition_repair.published_strings import native_pair as n
from model.candidate_laws.interval_echo_partition_repair.published_strings.memories import LOG_ODDS, Pair


def codes():
    return (n.NativePair(q, side) for q, side in itertools.product((0, 1), ("left", "right")))


def test_exact_generator_reduction_and_independent_entangled_input_columns():
    full = n.full_generator()
    theta = 0.283
    complete = expm(-2j*theta*full)
    for code in codes():
        expected = np.zeros((4, 4))
        expected[0, 0] = expected[3, 3] = 8.5
        expected[0, 3] = expected[3, 0] = -0.5
        np.testing.assert_allclose(code.generator, expected, atol=3e-14, rtol=0)
        encoding = np.eye(n.LAW.DIM)[:, code.columns]
        np.testing.assert_allclose(complete @ encoding, encoding @ code.gate(theta), rtol=0, atol=3e-14)
        forward = code.gate(theta)[:, 0]
        np.testing.assert_allclose(forward[[0, 3]]*np.exp(17j*theta),
                                   [np.cos(theta), 1j*np.sin(theta)], rtol=0, atol=1e-14)
        num = np.diag([sum(code.physical(*s)) for s in n.LOGICAL])
        rho = code.state(theta)
        np.testing.assert_allclose(num @ rho, rho @ num, atol=1e-14)


def test_local_native_pointers_and_full_spectrum_at_every_held_angle():
    for q in (0, 1):
        for i in (0, 1):
            for effect in n.local_effect(q, i):
                np.testing.assert_allclose(effect, np.diag([int(i == 0), int(i == 1)]),
                                           atol=2e-15, rtol=0)
    for code in codes():
        ref = code.state(n.THETA)
        for d in (-0.8, -0.3, 0, 1e-6, 0.4, 1.7):
            rho = code.state(n.THETA+d)
            target = Pair(d)
            np.testing.assert_allclose(np.linalg.eigvalsh(rho), [1/8]*3+[5/8], atol=1e-14)
            np.testing.assert_allclose(n.relative_entropy(rho, ref), target.whole_modular_change,
                                       atol=2e-14)
            for side in ("left", "right"):
                cut, vacuum = n.endpoint(rho, side), n.endpoint(ref, side)
                np.testing.assert_allclose(cut[1, 1], 3/8+target.probability_change, atol=3e-14)
                np.testing.assert_allclose(n.entropy(cut)-n.entropy(vacuum),
                                           target.cut_entropy_change, atol=2e-14)
                np.testing.assert_allclose(n.relative_entropy(cut, vacuum), target.cut_relative_entropy,
                                           atol=2e-14)
            np.testing.assert_allclose(np.trace(code.generator @ rho), 51/8, atol=4e-13)


def test_held_weighted_work_and_read_timing():
    work = 1e-6
    d = Pair.from_weighted_work(work).angle_increment
    for code in codes():
        ref = code.state(n.THETA)
        cut = n.endpoint(ref, "left")
        transferred = n.endpoint(code.state(n.THETA+d), "left")
        assert abs((n.entropy(transferred)-n.entropy(cut))/work-2*math.pi) < 0.0003
        assert abs(Pair(0).cut_entropy_change) == 0
        assert math.isclose(math.sqrt(3)/4*LOG_ODDS, 0.2211939835426899)
        np.testing.assert_allclose(np.diag(code.undo(ref, n.THETA)).real, [5/8, 1/8, 1/8, 1/8],
                                   atol=1e-14)
        marked = code.undo(np.diag(np.diag(ref)), n.THETA)
        np.testing.assert_allclose(np.diag(marked).real, [7/16, 1/8, 1/8, 5/16], atol=1e-14)
