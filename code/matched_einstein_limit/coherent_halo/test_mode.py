import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.experiments.seam_program import cubic_resonance as cubic
from . import measure


def test_bound_state_matches_independent_matrix_and_sham():
    critical = 1 / (12 * (cubic.green(0.0) - 1 / 6))
    record = measure.check_finite(8, critical, check_spectrum=True)
    assert record["unitary_eigen_residual"] < 1e-12
    kinetic, shell, _ = cubic.punctured_cube(8)
    np.testing.assert_array_equal(kinetic @ np.ones(len(shell)), 0.0)
    m = record["screening_mass"]
    state, _ = measure.normalized_mode(8, m)
    h = kinetic - sparse.diags(2 * critical * shell)
    after = expm_multiply(-0.4j * h, state)
    np.testing.assert_allclose(after, np.exp(0.4j * m * m) * state, atol=2e-13)
    # The same prepared state is released under the disabled anchor.
    control = expm_multiply(-0.4j * kinetic, state)
    assert np.linalg.norm(abs(control) ** 2 - state**2) > 0.001


def test_number_and_energy_continuity_off_equilibrium():
    kinetic, shell, _ = cubic.punctured_cube(6)
    rng = np.random.default_rng(6971)
    state = rng.normal(size=len(shell)) + 1j * rng.normal(size=len(shell))
    state /= np.linalg.norm(state)
    panel = measure.local_energy(kinetic, shell, state, 0.97, 0.5, 48.0)
    assert panel["max_bond_current"] > 0.001
    assert panel["max_energy_current"] > 0.01
    assert panel["number_error"] < 2e-14
    assert panel["energy_error"] < 2e-13
    # Both time derivatives independently follow from a unitary step of the finite graph.
    h = 0.5 * kinetic + sparse.diags(48.0 - 0.97 * shell)
    from scipy.linalg import expm

    dt = 1e-5
    u = expm(-1j * dt * h.toarray())
    before, after = u.conj().T @ state, u @ state
    z = h @ state
    np.testing.assert_allclose(
        (abs(after) ** 2 - abs(before) ** 2) / (2 * dt),
        2 * (state.conj() * z).imag,
        atol=2e-11,
    )
    secant = (after.conj() * (h @ after) - before.conj() * (h @ before)).real / (2 * dt)
    np.testing.assert_allclose(
        secant, (1j * z.conj() * z - 1j * state.conj() * (h @ z)).real, atol=2e-9
    )


def test_held_radial_ratios_and_native_contact():
    critical = 1 / (12 * (cubic.green(0.0) - 1 / 6))
    for n in (64, 128):
        nominal = measure.profile(n, critical)
        assert 3.85 < nominal["density_ratio"] < 4.15
    assert measure.profile(64, 8.0)["density_ratio"] > 1e8
    hw = measure.hardware(critical)
    assert hw["degree_restoration_error"] < 5e-14
    for key in ("occupied", "sham"):
        assert hw[key]["maximum_column_error"] < 6e-12
