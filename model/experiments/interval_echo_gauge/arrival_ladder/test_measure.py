"""Independent checks of source admission, coherent aperture and remote neutral reads."""

import json
from itertools import product

import numpy as np

from . import measure as m


def test_full_local_phase_and_compensated_charge_step():
    for phi in (0, np.pi / 3, np.pi / 2, np.pi):
        phase = m.selective_phase(phi, m.K)
        expected = np.ones((3, 3), complex)
        expected[0, 1] = np.exp(1j * phi)
        np.testing.assert_allclose(phase / phase[2, 2], expected, atol=8e-14)
        gate = m.interaction(phi, m.K)
        ids = [m.WORDS.index(m.SOURCES[k]) for k in ("negative", "neutral", "positive")]
        prob = abs(gate[np.ix_(ids, ids)]) ** 2
        p = np.sin(phi / 2) ** 2
        np.testing.assert_allclose(
            prob, [[1 - p, p, 0], [p, 1 - p, 0], [0, 0, 1]], atol=4e-14
        )
        charge = np.diag([-1, 0, 1])
        raise_only = np.zeros((3, 3), complex)
        raise_only[1, 0] = gate[ids[1], ids[0]]
        np.testing.assert_allclose(
            charge @ raise_only - raise_only @ charge, raise_only
        )
        for name in ("filled", "vacancy_at_other_port", "donor_depleted"):
            ix = m.WORDS.index(m.SOURCES[name])
            np.testing.assert_allclose(
                abs(gate[:, ix]) ** 2, np.eye(27)[:, ix], atol=4e-14
            )
    full = m.interaction(np.pi, m.K)
    np.testing.assert_allclose(abs(full @ full), np.eye(27), atol=5e-14)


def test_complete_cell_aperture_on_reference_entangled_modes():
    ports = list(range(72, 84))
    maps = m.aperture_map(ports)
    for key in ("charged_to_neutral", "neutral_to_charged", "singlet_survival"):
        transfer = maps[key]
        np.testing.assert_allclose(transfer, transfer[0, 0] * np.eye(12), atol=2e-13)
        np.testing.assert_allclose(transfer.conj().T @ transfer, np.eye(12), atol=2e-13)
    np.testing.assert_allclose(maps["charged_to_neutral"], 1j * np.eye(12), atol=2e-13)
    np.testing.assert_allclose(maps["charged_to_charged"], 0, atol=2e-13)
    # Six independent amplitudes and phases stay correlated with the external register.
    reference = np.fft.fft(np.eye(12), axis=0)[:, :6] / np.sqrt(72)
    np.testing.assert_allclose(
        maps["charged_to_neutral"] @ reference, 1j * reference, atol=2e-13
    )
    assert all(
        path[0] == p and path[-1] == 3
        for p, path in zip(ports, maps["routes"], strict=True)
    )
    np.testing.assert_allclose(
        abs(m.aperture_map(ports, np.pi / 2)["charged_to_neutral"]) ** 2,
        np.eye(12) / 2,
        atol=2e-13,
    )


def test_third_owner_via_independent_full_ternary_path():
    # Include all 243 ternary words, with no sector or permutation reduction.
    site = m.layout()
    path = site["neutral_path"]
    board = m.federation(0)
    words = list(product(range(3), repeat=5))
    index = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    for i, w in enumerate(words):
        h[i, i] = m.K * sum(
            (6 - board.degree[p]) * (v * v - 1) for v, p in zip(w, path, strict=True)
        )
        for j, (a, b) in enumerate(zip(w, w[1:])):
            h[i, i] += m.K * (a - b) ** 2
            if a != b:
                h[i, i] += 0.5
                target = w[:j] + (b, a) + w[j + 2 :]
                h[index[target], i] -= 0.5
    # The complete receiver is the two signed loads, including an undisturbed (1,1).
    freqs, vec = np.linalg.eigh(h)
    for source in ("negative", "neutral", "positive"):
        origin = m.SOURCES[source]
        for t in (0.25, 4, 16):
            launch = index[(origin[1], 1, origin[0], 1, 1)]
            state = vec @ (np.exp(-1j * t * freqs) * vec[launch, :])
            counts = np.zeros(9)
            for w, z in zip(words, state, strict=True):
                counts[3 * w[-2] + w[-1]] += abs(z) ** 2
            np.testing.assert_allclose(
                m.release(site, origin, m.K, t), counts, atol=2e-12
            )
    assert site["source"] // 12 != path[-1] // 12 != site["x"] // 12


def test_archived_bank_completeness_and_falsifiers():
    result = m.measure()
    archive = json.loads((m.HERE / "results.json").read_text())
    for fresh, saved in zip(result["banks"], archive["banks"], strict=True):
        assert sum(saved["missed_and_all_27_counts"]) == m.SHOTS
        np.testing.assert_allclose(
            fresh["per_herald_probabilities"],
            saved["per_herald_probabilities"],
            atol=2e-13,
        )
        assert fresh["missed_and_all_27_counts"] == saved["missed_and_all_27_counts"]
    # Keep rare-count tests statistical across platforms; the tolerances are
    # for the simulated acquisitions. The unitary checks have their own error bound.
    for fresh, saved in zip(
        result["independent_receiver"], archive["independent_receiver"], strict=True
    ):
        assert sum(saved["missed_and_all_9_counts"]) == m.RELEASE_SHOTS
        np.testing.assert_allclose(
            fresh["per_herald_probabilities"],
            saved["per_herald_probabilities"],
            atol=2e-12,
        )
        if (
            fresh["source"] == "negative"
            and fresh["angle_over_pi"] == 1
            and fresh["time"] == 16
        ):
            assert 35 <= fresh["neutral_count"] <= 145
    assert result["read_loads"] == [0, 1, 2]
    active = next(
        r
        for r in archive["banks"]
        if r["source"] == "negative" and r["kappa"] == 4 and r["angle_over_pi"] == 1
    )
    assert (
        abs(
            active["per_herald_probabilities"][m.WORDS.index(m.SOURCES["neutral"])]
            - 0.5
        )
        < 2e-13
    )


def test_static_susceptibility_and_complete_rejections():
    from . import static

    archive = json.loads((m.HERE / "static.json").read_text())
    site = static.geometry()
    assert site["donor"] == 11 and site["comparator"] == 10
    assert archive["site"] == site
    keys = [m.WORDS.index(w) for w in [(0, 2, 1), (0, 1, 2), (2, 2, 1)]]
    for kappa, parity in [(m.K, 0), (m.K, 1), (m.K / 2, 0)]:
        for time in (0, 0.318, np.pi / 2, np.pi, 2 * np.pi):
            for i, source in enumerate(("negative", "neutral", "positive")):
                computed = static.instrument(m.SOURCES[source], time, kappa, parity)
                detune = 4 * kappa + parity if source == "positive" else parity
                predicted = np.sin(np.hypot(1, detune) * time / 2) ** 2 / (
                    1 + detune**2
                )
                np.testing.assert_allclose(computed[keys[i]], predicted, atol=2e-13)
    for row in archive["restored_then_rejected_banks"]:
        prob = static.instrument(
            m.SOURCES[row["source"]],
            np.pi * row["time_over_pi"],
            row["kappa"],
            row["comparator_parity"],
        )
        np.testing.assert_allclose(prob, row["probabilities"], atol=2e-13)
        counts = np.array(row["counts"])
        assert sum(counts) == archive["shots_per_prepared_input"]
        expected = archive["shots_per_prepared_input"] * prob
        assert np.all(
            abs(counts - expected) < 7 * np.sqrt(expected * (1 - prob).clip(0)) + 10
        )
