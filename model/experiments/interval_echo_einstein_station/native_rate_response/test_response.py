import numpy as np

from model.experiments.interval_echo_einstein_station import backreaction as seam
from .energy_probe import measured_probe
from .measure import integer_native_parts, rate_experiment


def test_rate_and_native_operator():
    for integer, implementation in zip(
        integer_native_parts(), seam.law.local_term_parts(), strict=True
    ):
        np.testing.assert_allclose(integer / 8, implementation, atol=2e-13, rtol=0)
    rows = rate_experiment()
    native, slowed = next(row[2:] for row in rows if row[:2] == (0.5, np.pi / 2))
    np.testing.assert_allclose(native, [0.25, 0.75], atol=3e-12)
    assert slowed[0] < 0.15 and slowed[1] > 0.85


def test_initial_derivatives_with_entangled_environment():
    # Direct full action on (memory AB, memory BC, load A, load B, load C).
    exchange, contrast = seam.law.local_term_parts()
    h = exchange + seam.KAPPA * contrast
    rng = np.random.default_rng(481)
    env = rng.normal(size=(11, 11, 3)) + 1j * rng.normal(size=(11, 11, 3))
    env /= np.linalg.norm(env)
    blank_weight = np.sum(abs(env[0]) ** 2)

    def act(state):
        result = np.zeros_like(state)
        for axes in ((0, 2, 3), (1, 3, 4)):
            moved = np.moveaxis(state, axes, (0, 1, 2))
            out = (h @ moved.reshape(99, -1)).reshape(moved.shape)
            result += np.moveaxis(out, (0, 1, 2), axes)
        return result

    sharp = np.zeros((11, 11, 3, 3, 3), dtype=complex)
    sharp[:, :, 2, 0] = env
    terminal = act(sharp)[:, :, 0, 2]
    curvature = 2 * np.vdot(terminal, terminal).real
    np.testing.assert_allclose(curvature, blank_weight / 2, atol=2e-14)
    phase = sharp.copy()
    phase[:, :, 0, 2] = 1j * env
    phase /= np.sqrt(2)
    rate = 2 * np.vdot(phase[:, :, 0, 2], -1j * act(phase)[:, :, 0, 2]).real
    np.testing.assert_allclose(rate, blank_weight / 2, atol=2e-14)


def test_dynamic_pointer_and_changed_coupling():
    arms = measured_probe()
    shifts = {
        name: np.array([r["high_yes"] - r["low_yes"] for r in rows])
        for name, rows in arms.items()
    }
    assert shifts["whole"][0] > 0 and shifts["whole"][-1] < 0
    for control in ("contrast", "idle"):
        np.testing.assert_allclose(shifts[control], 0, atol=2e-13)
    for whole, contrast in zip(arms["whole"], arms["contrast"], strict=True):
        np.testing.assert_allclose(
            whole["ordinary_yes"][:2], contrast["ordinary_yes"][:2], atol=2e-13
        )
