import json

import numpy as np
from scipy.linalg import expm

from . import historical, measure


def test_exact_channels_and_reversible_comparator():
    law = measure.LAW
    assert measure.exact_channel() == dict(
        cell_words=78,
        matrix_units=6084,
        native_matches_E=6084,
        changed_matches_T=6084,
        distinguishing_units=110,
        native_label_identity_all_qutrit_keys=True,
    )
    swap = np.zeros((law.DIM, law.DIM))
    memory = [None, *law.LABELS]
    for m, label in enumerate(memory):
        reversed_memory = (
            0 if label is None else memory.index((1 - label[0], -label[1]))
        )
        for i, key in enumerate(law.SRC):
            swap[
                reversed_memory * law.N_SRC + law.SRC.index(tuple(reversed(key))),
                m * law.N_SRC + i,
            ] = 1
    for u in (law.echo(), measure.total_sensitive()):
        np.testing.assert_allclose(u.T, u, atol=2e-15)
        np.testing.assert_allclose(u @ u, np.eye(law.DIM), atol=2e-15)
        np.testing.assert_allclose(swap @ u @ swap.T, u, atol=2e-15)
    k0, k1 = law.local_term_parts()
    h = k0 + measure.candidate.KAPPA * k1
    # The reciprocal apparatus uses positive-duration complementary native holds.
    for t in (np.pi / 2, np.pi / (4 * (measure.candidate.KAPPA + 0.5))):
        np.testing.assert_allclose(
            expm(-1j * (2 * np.pi - t) * h), expm(1j * t * h), atol=2e-13
        )
    # The registered blank compression commutes with the calibration phases
    # on its spectators and stays inside the sharp two-load word sector.
    np.testing.assert_allclose(h[law.N_SRC :, : law.N_SRC], 0, atol=1e-14)


def test_literal_port_source_and_reciprocal_detector():
    p = historical.preparation()
    source = np.zeros(len(historical.WORDS), dtype=complex)
    source[[historical.X, historical.W]] = 1 / np.sqrt(2)
    np.testing.assert_allclose(abs(source.conj() @ p[:, historical.X]), 1, atol=2e-13)
    np.testing.assert_allclose(p.conj().T @ p, np.eye(len(source)), atol=2e-13)
    # A completely independent partial trace on the 99-state seam, with the
    # spectator configurations retained as two distinct orthogonal kets.
    for name, unitary in [
        ("native", measure.LAW.echo()),
        ("total", measure.total_sensitive()),
    ]:
        records = []
        for key in [(0, 0), (1, 1)]:
            i = measure.LAW.SRC.index(key)
            write = unitary[:, i].reshape(-1, measure.LAW.N_SRC)
            np.testing.assert_allclose(np.delete(write, i, axis=1), 0, atol=2e-15)
            records.append(write[:, i])
        overlap = np.vdot(*records)
        np.testing.assert_allclose(overlap, 1 if name == "native" else 0, atol=2e-15)
        for phi in measure.PHASES:
            odds = np.linalg.norm((records[0] + np.exp(1j * phi) * records[1]) / 2) ** 2
            analytic = (1 + (np.cos(phi) if name == "native" else 0)) / 2
            assert abs(odds - analytic) < 2e-13


def test_original_launch_banks_and_cross_port_predictions():
    a = measure.run()
    b = historical.run(a)
    assert a == json.loads((measure.HERE / "results.json").read_text())
    assert b == json.loads((measure.HERE / "historical.json").read_text())
    for archive, bright in ((a, measure.INITIAL), (b, historical.X)):
        for mode in measure.MODES:
            for phi, row in zip(measure.PHASES, archive["banks"][mode], strict=True):
                counts = row["counts"]
                assert sum(counts) == measure.SHOTS
                assert min(counts) >= 0
                if mode.startswith("only_"):
                    expected = 0.25
                    assert (
                        abs(counts[-1] / measure.SHOTS - 0.5)
                        < b["combined_rate_radius"]
                    )
                else:
                    expected = (
                        (1 + np.cos(phi)) / 2 if mode in ("native", "bypass") else 0.5
                    )
                    assert counts[-1] == 0
                assert (
                    abs(counts[bright] / measure.SHOTS - expected)
                    <= b["combined_rate_radius"]
                )
    assert all(
        r["residual"] <= r["joint_error"]
        for group in b["banks"].values()
        for r in group
    )
    assert all(
        abs(r["forecast"] - r["observed"]) <= r["joint_error"]
        for group in a["held"].values()
        for r in group
    )
    assert a["two_writes_single_load"] == {
        "fresh_memory_return": 0.5,
        "same_memory_return": 1,
    }
