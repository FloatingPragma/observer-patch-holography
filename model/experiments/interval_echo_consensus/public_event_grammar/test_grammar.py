"""Native source maps, complete records and independent equivalence checks."""

import copy
from fractions import Fraction
import itertools
import json

import numpy as np
import pytest

from . import grammar as g


def test_full_dirty_collar_native_tree_matches_forward_state_updates():
    p = g.native_program()
    branches = g.compile_program(p, 99)
    rng = np.random.default_rng(8304)
    psi = rng.normal(size=99) + 1j * rng.normal(size=99)
    psi /= np.linalg.norm(psi)
    rho = np.outer(psi, psi.conj())
    law = g.runtime_law(branches, rho)

    def forward(program, vector, prefix=()):
        if isinstance(program, g.Halt):
            return {(*prefix, program.label): np.vdot(vector, vector).real}
        result = {}
        for label, k, next_program in program.branches:
            result.update(forward(next_program, k @ vector, (*prefix, label)))
        return result

    expected = forward(p, psi)
    assert law.keys() == expected.keys()
    np.testing.assert_allclose(
        list(law.values()), list(expected.values()), atol=2e-14, rtol=0
    )
    gram = g.event_effect(branches, lambda _: True)
    np.testing.assert_allclose(gram, np.eye(99), atol=5e-14, rtol=0)
    assert any("refused:empty" in h for h in law)
    assert any("early:unequal" in h for h in law)


def test_labelled_refinements_and_different_postmaps_have_same_event_signature():
    plain = g.compile_program(g.native_program(), 99)
    fine = g.compile_program(g.native_program(refined=True), 99)
    changed = g.compile_program(g.native_program(post=True), 99)
    assert len(plain) == 9 and len(fine) == 17
    map_difference = 0.0
    for r in range(9):
        label = f"material:{g.law.SRC[r][0]}{g.law.SRC[r][1]}"

        def event(history, label=label):
            return label in history

        e = g.event_effect(plain, event)
        for other in (fine, changed):
            np.testing.assert_allclose(
                g.event_effect(other, event), e, atol=5e-14, rtol=0
            )
        for history, k in plain.items():
            if event(history) and "refused:empty" not in history:
                other = next(v for h, v in changed.items() if event(h))
                map_difference = max(map_difference, np.linalg.norm(k - other))
    assert map_difference > 1


def test_independent_complete_signature_not_target_state_equality():
    rng = np.random.default_rng(83)
    a = rng.normal(size=(13, 13)) + 1j * rng.normal(size=(13, 13))
    a = (a + a.conj().T) / 2
    np.testing.assert_allclose(g.reconstruct(g.probes(a), 13), a, atol=1e-14, rtol=0)
    # Identical probabilities on one mixed input fail to certify equivalence.
    e = np.diag([1.0, 0])
    f = np.diag([0.0, 1])
    rho = np.eye(2) / 2
    assert np.trace(rho @ e) == np.trace(rho @ f)
    assert not np.array_equal(g.probes(e), g.probes(f))
    # The signature entries are literal normalized positive pure probe weights.
    eye = np.eye(13, dtype=complex)
    expected = list(a.diagonal().real)
    for i in range(13):
        for j in range(i + 1, 13):
            for vector in (
                (eye[i] + eye[j]) / np.sqrt(2),
                (eye[i] + 1j * eye[j]) / np.sqrt(2),
            ):
                expected.append(np.vdot(vector, a @ vector).real)
    np.testing.assert_allclose(g.probes(a), expected, atol=1e-14, rtol=0)


def test_actual_normal_entrypoint_native_and_classical_contracts():
    audit = g.source_audit()
    assert all(audit["checks"].values())
    _, rates, error = g.clock._primitives()
    assert error < 1e-11
    for sign in (0, 2):
        assert abs(rates[sign] - 0.5) < 1e-13
    # The clock's exact scripted fair stream is separate from a Born inference.
    for q in (1, 2, 3, 8, 16):
        for bits in itertools.product((0, 1), repeat=4):
            tape = itertools.cycle(bits)
            result = g.clock.closing_tick(q, 1, lambda: next(tape))
            assert type(result) is bool
    law = g.bounded_address_reader(3)
    assert sum(map(Fraction, law.values())) == 1
    assert Fraction(law["timeout"]) == Fraction(1, 4) ** 12


def test_context_normalization_alone_fails_same_effect_refinement():
    c = g.exact_counterexample()
    assert c["normalized_each_context"]
    assert Fraction(c["fine_union"]) == Fraction(1, 3)
    assert Fraction(c["coarse_same_effect"]) == Fraction(1, 2)
    assert Fraction(c["discrepancy"]) == Fraction(1, 6)
    # Refusal deletion is another invalid normalization shortcut.
    with pytest.raises(ValueError, match="incomplete primitive"):
        g.compile_program(
            g.Step("missing quiet", (("yes", np.sqrt(0.5) * np.eye(2), g.Halt()),)), 2
        )
    with pytest.raises(ValueError, match="normalized initial"):
        g.runtime_law({("complete",): np.eye(2)}, np.eye(2))


def test_receipt_sources_discrete_records_and_floating_replay():
    saved = json.loads(g.RECEIPT.read_text())
    actual = g.report()
    g.compare(saved, actual)
    damaged = copy.deepcopy(saved)
    damaged["base_histories"] -= 1
    with pytest.raises(ValueError):
        g.compare(damaged, actual)


def test_literal_initial_law_and_actual_overflow_have_distinct_domains():
    from types import SimpleNamespace

    law = g.source_initial_law()["exact_complete_weights"]
    weights = list(map(Fraction, law.values()))
    assert sum(weights) == 1 and min(weights) > 0
    # A host exception must not be silently reported as an observer outcome.
    exhausted = SimpleNamespace(get=lambda digits, owner: 511)
    with pytest.raises(OverflowError, match="elapsed word exhausted"):
        g.clock.Registers.post(exhausted, 0, [None] * 9, None)
