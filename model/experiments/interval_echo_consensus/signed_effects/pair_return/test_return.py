import json

import numpy as np
import sympy as s

from . import measure as m


def test_native_branch_maps_and_terminal_reader_for_every_reachable_prior():
    calls = m.address.schedule()
    coarse_calls = [call for call in calls if call["bit"] != 0]
    coarse = m.read_branches(coarse_calls, 0)
    fine = m.read_branches(calls, 0)
    eye = np.eye(m.DIM)
    phase = m.address.COMMON_PHASE
    for raw in range(m.LABELS):
        coordinate = np.diag(np.arange(m.DIM) == raw)
        pair = np.diag((np.arange(m.DIM) & ~1) == raw)
        np.testing.assert_allclose(fine[raw], phase**len(calls) * coordinate, atol=2e-13, rtol=0)
        np.testing.assert_allclose(coarse[raw], phase**len(coarse_calls) * pair, atol=2e-13, rtol=0)
    for prior in range(m.DIM):
        second = m.read_branches(calls, prior)
        for port, vector in enumerate(eye):
            k = second[prior ^ port]
            np.testing.assert_allclose(k.conj().T @ k, np.outer(vector, vector), atol=2e-13, rtol=0)


def test_exact_uniform_preparation_analyzer_and_disturbance():
    law = m.pulses.candidate.native_law
    h0, h1 = law.local_term_parts()
    code = np.eye(law.DIM)[:, [law.SRC.index((2, 1)), law.SRC.index((1, 2))]]
    actual_contact = h0 + m.pulses.candidate.KAPPA * h1
    reduced_contact = np.array([[8.5, -0.5], [-0.5, 8.5]])
    np.testing.assert_allclose(actual_contact @ code, code @ reduced_contact, atol=2e-13, rtol=0)
    assert m.pulses.federation(0).degree[0] == 5
    assert [m.pulses.candidate.KAPPA * load**2 for load in (2, 1)] == [32, 8]
    # The excess and the unit host give exchange eigenvalues 8 and 9.
    x = s.Matrix([[0, 1], [1, 0]])
    symmetric = (s.eye(2) + x) / 2
    antisymmetric = s.eye(2) - symmetric
    half = symmetric - s.I * antisymmetric
    # The addressed reference has costs 32 and 8 at port zero.
    first = s.diag(1, -s.I)
    last = s.diag(1, s.I)
    zero = s.diag(1, 0)
    one = s.diag(0, 1)
    probe = first * half
    analyzer = last * half * first
    uniform = s.ones(2) / 2
    assert s.simplify(probe * zero * probe.H) == uniform
    assert s.simplify(analyzer.H * zero * analyzer) == uniform
    for projector in (zero, one):
        assert zero * projector * zero + one * projector * one == projector
    dephased = zero * uniform * zero + one * uniform * one
    assert (dephased * uniform).trace() == s.Rational(1, 2)
    assert (uniform * uniform).trace() == 1


def test_every_held_record_and_original_launch_replays():
    receipt = json.loads(m.RECEIPT.read_text())
    m.compare(receipt, m.report())
    assert len(receipt["banks"]) == 36
    for bank in receipt["banks"]:
        counts = bank["nonzero_joint_counts"]
        assert sum(n for _, _, n in counts) == m.SHOTS
        if bank["return_probability"] == 1:
            assert bank["returned"] == m.SHOTS
        else:
            assert bank["probe"] == "plus" and bank["instrument"] == "address_then_forget"
