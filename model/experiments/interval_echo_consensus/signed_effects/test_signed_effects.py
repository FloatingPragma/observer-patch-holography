import json

import numpy as np
from scipy.linalg import block_diag

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.measure import (
    generator,
)
from . import measure, refinement


def lie_rank(q1, q2):
    """Independent closure on two three-position sectors using the actual bond matrices."""
    dim = 3
    eye = np.eye(dim)
    raw = []
    for a in range(dim):
        e = np.outer(eye[a], eye[a])
        raw.append(1j * block_diag(8 * eye + q1 * e, 8 * eye + q2 * e))
    for a, b in [(0, 1), (1, 2)]:
        local = np.zeros((dim, dim))
        local[np.ix_([a, b], [a, b])] = [[8.5, -0.5], [-0.5, 8.5]]
        raw.append(1j * block_diag(local, local))
    basis = []

    def append(z):
        for u in basis:
            z = z - np.vdot(u, z).real * u
        norm = np.linalg.norm(z)
        if norm > 1e-7:
            basis.append(z / norm)

    for r in raw:
        append(r)
    for z in basis:
        for r in raw:
            append(z @ r - r @ z)
    return len(basis)


def test_joint_rank_mutation_and_actual_exposures():
    assert lie_rank(-8, 24) == 18
    assert lie_rank(-8, -8) == 9
    receipt = measure.report()
    assert receipt["algebra"]["commutator_determinant"] == 18874368
    assert (
        receipt["joint_counts"]
        == json.loads((measure.HERE / "results.json").read_text())["joint_counts"]
    )
    np.testing.assert_allclose(
        receipt["conditional_remote_probabilities"], [0.75, 0], rtol=0, atol=3e-14
    )
    assert receipt["normalized_fourth_power_spectator_reads"] == ["1/2", "2/3"]
    assert receipt["unverified_prior_remote_probability"] == "3/40"


def test_every_registered_lift_and_free_leakage():
    result = refinement.report()
    assert result == json.loads((measure.HERE / "refinement.json").read_text())
    assert [d["protected_child_encodings"] for d in result["transitions"]] == [20, 80]
    assert [
        d["free_by_sign"]["1"]["nonzero_columns"] for d in result["transitions"]
    ] == [120, 444]


def test_arbitrary_child_charts_cannot_eliminate_exterior_loss():
    # Direct dense computation on two random unitary charts per child.
    rng = np.random.default_rng(920131)
    fine = federation(1)
    charts = []
    for _ in range(fine.cells):
        a = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
        charts.append(np.linalg.qr(a)[0])
    v = [np.vstack(charts[i : i + 4]) / 2 for i in range(0, fine.cells, 4)]
    cols = np.zeros((fine.n, 12), complex)
    cols[:48] = v[0]
    for sign in (-1, 1):
        result = generator(fine, sign, 8) @ cols
        exterior = result[48:].reshape(-1, 48, 12)
        np.testing.assert_allclose(
            np.linalg.norm(exterior) ** 2, 3 / 8, atol=2e-13, rtol=0
        )
        loss = sum(
            np.linalg.norm(a - lift @ (lift.conj().T @ a)) ** 2
            for lift, a in zip(v[1:], exterior)
        )
        assert loss >= 3 / 16 - 1e-13


def test_child_coherence_and_address_coarsening():
    rng = np.random.default_rng(92)
    chart = np.vstack([np.eye(12)] * 4) / 2
    ref = rng.normal(size=(12, 2)) + 1j * rng.normal(size=(12, 2))
    ref /= np.linalg.norm(ref)
    lifted = (chart @ ref).reshape(4, 12, 2)
    a = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
    u = np.linalg.qr(a)[0]
    exact = u @ ref
    after = np.stack([u @ slot for slot in lifted])
    np.testing.assert_allclose(chart.conj().T @ after.reshape(48, 2), exact, atol=6e-15)
    recovered = 0
    for j in range(4):
        branch = np.zeros_like(after)
        branch[j] = after[j]
        recovered += np.linalg.norm(chart.conj().T @ branch.reshape(48, 2)) ** 2
    np.testing.assert_allclose(recovered, 1 / 4, atol=6e-15)


def test_quantum_projection_for_correlated_competitor_and_invisible_entropy_slope():
    rng = np.random.default_rng(791)
    a = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
    sigma = a @ a.conj().T
    sigma /= np.trace(sigma)
    rho = np.trace(sigma.reshape(4, 3, 4, 3), axis1=0, axis2=2)
    tau = np.diag([0.2, 0.3, 0.5])
    reference = np.kron(np.eye(4) / 4, tau)
    maxent = np.kron(np.eye(4) / 4, rho)

    def logm_positive(x):
        v, u = np.linalg.eigh(x)
        return (u * np.log(v)) @ u.conj().T

    def rel(x, y):
        return np.trace(x @ (logm_positive(x) - logm_positive(y))).real

    assert rel(sigma, maxent) > 0
    np.testing.assert_allclose(
        rel(sigma, reference), rel(rho, tau) + rel(sigma, maxent), atol=4e-13, rtol=0
    )

    def child_entropy(coherence):
        p = (1 + 3 * coherence) / 4
        q = (1 - coherence) / 4
        return -p * np.log(p) - 3 * q * np.log(q)

    eps = 1e-5
    derivative = (child_entropy(0.5 + eps) - child_entropy(0.5 - eps)) / (2 * eps)
    np.testing.assert_allclose(derivative, -3 * np.log(5) / 4, rtol=0, atol=1e-9)
