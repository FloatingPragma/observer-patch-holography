"""Independent local gates, native ports, routing and complete absorption records."""
import numpy as np
import pytest

from . import clock_loader as loader
from . import clock_transport as move
from .. import law


@pytest.mark.parametrize("slot", range(8))
def test_common_input_is_optimal_and_preserves_rail_phase(slot):
    delta, alpha, p = move.forecast()
    matrices = [law.contact_coin(np.eye(8), phi) for phi in law.LAWS["torus"].rail_contact_phases]
    w = np.array([matrix[slot].conj() for matrix in matrices])
    value, vectors = np.linalg.eigh((np.outer(w[0], w[0].conj()) + np.outer(w[1], w[1].conj())) / 2)
    assert value[-1] == pytest.approx(p, abs=2e-15)
    assert p < 1 - 1e-8 and abs(np.vdot(vectors[:, -1], move.common_ray(slot))) > 1 - 2e-15
    coefficients = [matrix[slot] @ move.common_ray(slot) for matrix in matrices]
    np.testing.assert_allclose(np.abs(coefficients) ** 2, p, atol=3e-15, rtol=0)
    assert np.angle(coefficients[0] * coefficients[1].conj()) == pytest.approx(alpha, abs=2e-15)
    assert -delta < alpha < 0
    # The eight slots share the very same two rail phases on all directed links.
    assert np.angle(coefficients[0]) == pytest.approx(0, abs=2e-15)


def test_arrival_at_four_actual_neighbours_and_local_failure_decomposition():
    carrier, root = (5, 5, 3, 3), (4, 2, 1, 1)
    loop = loader.circulation(carrier, root)
    balanced = np.array([loop, 1j * loop]) / np.sqrt(2)
    delta, alpha, p = move.forecast()
    for slot in (0, 1, 2, 3):
        screened, where, failed, selected = move.translate(balanced.copy(), root, slot, True, 128)
        unscreened, _, blank, same = move.translate(balanced.copy(), root, slot, False, 128)
        target = np.roll(loop, 1 if slot % 2 == 0 else -1, axis=slot // 2)
        assert selected == pytest.approx(p) and same == pytest.approx(p)
        np.testing.assert_allclose(screened[0], target * np.sqrt(p / 2), atol=2e-15, rtol=0)
        np.testing.assert_allclose(screened[1], target * 1j * np.exp(-1j * alpha) * np.sqrt(p / 2),
                                   atol=2e-15, rtol=0)
        assert failed.sum() == pytest.approx(1-p, abs=3e-15)
        assert not blank.any()
        assert np.linalg.norm(unscreened - screened) ** 2 == pytest.approx(1-p, abs=3e-15)
        # Independent explicit scattering of all 32 incident slots by each rail coin.
        bins = np.zeros((*carrier, 8))
        ray = move.common_ray(slot)
        for r, phi in enumerate(law.LAWS["torus"].rail_contact_phases):
            scattered = law.contact_coin(ray, phi)
            for site in loader.corners(carrier, root):
                for s in range(8):
                    if s == slot:
                        continue
                    dest = list(site)
                    dest[s // 2] = (dest[s // 2] + (1 if s % 2 == 0 else -1)) % carrier[s // 2]
                    bins[(*dest, s ^ 1)] += abs(scattered[s]) ** 2 / 8
        np.testing.assert_allclose(failed, bins, atol=2e-15, rtol=0)


def test_local_plate_carries_full_complex_ray_without_corner_phase():
    rng = np.random.default_rng(4701)
    for _ in range(10):
        a, b = rng.normal(size=(2, 8)) + 1j * rng.normal(size=(2, 8))
        a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
        gate = move.plate(a, b)
        np.testing.assert_allclose(gate @ a, b, atol=1e-14)
        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(8), atol=1e-14)



def test_registered_transport_and_original_start_counts():
    assert all(move.verify_records().values())


def test_driven_clock_rejects_one_impossible_detection_among_booked_launches():
    # Before the first driven hop no launch can have tripped its absorber.
    with np.load(move.PROTOCOL.with_name("clock_transport_counts.npz"), allow_pickle=False) as saved:
        counts = dict(saved)
    word = counts["forward_1_sieved_0_booked"]
    assert word[2] == 0 and word[0] > 0
    word[0] -= 1
    word[2] += 1
    with pytest.raises(ValueError, match="zero-probability outcome"):
        move.register(counts)


def test_check_command_returns_failure_for_a_missed_gate(monkeypatch):
    monkeypatch.setattr(move, "published_replay",
                        lambda seed_replay: ({"gates": {"complete_count_gate": False}}, {}))
    assert move.main(["--check"]) == 1
