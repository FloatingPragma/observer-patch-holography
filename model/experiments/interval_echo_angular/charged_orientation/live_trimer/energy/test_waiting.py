import json

import numpy as np

from . import waiting


def test_uniform_wait_matches_independent_native_evolution_and_counts():
    data = json.loads((waiting.energy.HERE / "waiting.json").read_text())
    # Reassemble the full 480-state native sector and propagate explicitly at
    # positive quadrature nodes. This test does not use the line arrays to
    # form the reference channel for the quarter-unit interval.
    law = waiting.energy.trimer
    a, ix = law.sector(1, True, True)
    spec = law.spectrum(a, ix)
    rho = np.zeros((4, 4), complex)
    nodes, weights = np.polynomial.legendre.leggauss(24)
    for z, weight in zip(nodes, weights, strict=True):
        flight = waiting.energy.FLIGHT + (z + 1) * 0.125
        m = law.branch(spec, ix, flight)[:, 1:]
        bell = m.transpose(0, 2, 1).reshape(8, 4) / np.sqrt(2)
        rho += weight / 2 * (bell.T @ bell.conj())
    with np.load(waiting.energy.HERE / "lines.npz") as f:
        expected = waiting.averaged(
            f["energies_1"], f["amplitudes_1"], waiting.energy.FLIGHT, 0.25
        )
    np.testing.assert_allclose(rho, expected, rtol=0, atol=4e-9)
    direct = waiting.run()
    for r, existing in zip(direct["rows"], data["rows"], strict=True):
        assert np.all(np.asarray(r["pauli_counts"]).sum(axis=2) == waiting.SHOTS)
        lo, hi = r["simultaneous_pt_interval"]
        assert lo <= r["theory_pt_min"] <= hi
        assert abs(r["measured_pt_min"] - existing["measured_pt_min"]) < 3e-5
    by_width = {r["width"]: r for r in direct["rows"]}
    assert by_width[2.0]["simultaneous_pt_interval"][1] < 0
    assert by_width[64.0]["accepted_negativity_upper_99"] < 0.00111
    assert by_width[64.0]["slow_pair_visibility"] > 0.9999999
    assert (
        direct["synchronized_center_of_width_64"]["accepted_negativity_lower_99"]
        > 0.027
    )
    assert direct["rows"][-1]["slow_pair_visibility"] < 0.00004
    assert direct["infinite_wait"]["theory_pt_min"] >= 0
    assert direct["no_flight"]["observed_mixed_arrival"] == 0
