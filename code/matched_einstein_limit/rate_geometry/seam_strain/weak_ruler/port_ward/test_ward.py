import json

import numpy as np
from scipy.linalg import expm
from mpmath.ctx_iv import MPIntervalContext

from . import measure as m


def test_pair_divergence_for_an_entangled_four_port_density():
    # Complex Hermitian hopping with an explicit two-step 0→2 contribution.
    rng = np.random.default_rng(209)
    h = np.zeros((8, 8), complex)
    for p in range(3):
        z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        h[2 * p : 2 * p + 2, 2 * p + 2 : 2 * p + 4] = z
        h[2 * p + 2 : 2 * p + 4, 2 * p : 2 * p + 2] = z.conj().T
    x = rng.normal(size=8) + 1j * rng.normal(size=8)
    rho = np.outer(x, x.conj()) / np.vdot(x, x)
    proj = [np.diag(np.repeat(np.eye(4)[i], 2)) for i in range(4)]
    weights = np.array([0.15, 0.5, 0.7, 0.9])
    velocities = np.array([0.1, -0.3, 0.2, 0.02])
    w = sum(a * p for a, p in zip(weights, proj, strict=True))
    rate = sum(a * p for a, p in zip(velocities, proj, strict=True))
    for action, norm, inventory in (
        (h, 1, lambda p: p),
        (h @ h, 2, lambda p: (p @ h + h @ p) / 2),
    ):
        flow = np.zeros_like(h)
        for i, j in m.combinations(range(4), 2):
            b = proj[i] @ action @ proj[j]
            flow += (weights[i] - weights[j]) * 1j * (b - b.conj().T) / norm
        expected = np.trace((flow - inventory(rate)) @ rho).real

        def endpoint(t):
            u = expm(-1j * t * h)
            return np.trace(inventory(w + t * rate) @ u @ rho @ u.conj().T).real

        dt = 1e-5
        assert abs(expected - (endpoint(-dt) - endpoint(dt)) / (2 * dt)) < 2e-8
    assert np.abs(proj[0] @ (h @ h) @ proj[2]).max() > 0.01


def test_all_native_port_identities_and_finite_nested_accounts():
    run = m.run()
    archive = json.loads((m.HERE / "results.json").read_text())
    native = run["native_operator_checks"]
    assert native["verified_port_operators"] == 960
    assert native["nonzero_regional_commutator_entries"] == 276704
    assert native["graph_bonds_with_number_current"] == 2520
    assert native["energy_pairs_on_graph_bonds"] == 2520
    assert native["energy_pairs_two_hops_apart"] == 3600
    assert native == archive["native_operator_checks"]
    ctx = MPIntervalContext()
    ctx.dps = 35
    for material in ("number", "energy"):
        for arm in m.face.SOURCES:

            def measurements(region):
                return run["regions"][region][material][arm]

            for field in (
                "endpoint_loss",
                "integrated_cut_current",
                "integrated_membership_change",
            ):
                vals = [
                    ctx.mpf(measurements(r)[field])
                    for r in ("inner", "annulus", "exterior", "outer")
                ]
                assert 0 in vals[0] + vals[1] + vals[2]
                assert 0 in vals[0] + vals[1] - vals[3]
            for region in run["regions"]:
                readings = measurements(region)
                assert 0 in ctx.mpf(readings["corrected_residual"])
                for key, measured in readings.items():
                    archival = archive["regions"][region][material][arm][key]
                    assert max(measured[0], archival[0]) <= min(
                        measured[1], archival[1]
                    )
            assert 0 not in ctx.mpf(measurements("outer")["residual_omitting_motion"])
