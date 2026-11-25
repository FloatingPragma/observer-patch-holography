"""Check state transfer on arbitrary complex inputs and the material-only evidence."""

import json
from pathlib import Path

import numpy as np

from model.experiments.interval_echo_einstein_station.collar_pairing.phase_geometry import (
    measure,
)


def test_output_rails_transfer_the_complete_input_operator_space():
    for q in (1, -1):
        w = measure.Z if q == 1 else measure.Z @ measure.X
        joint = np.kron(w, w)
        for i in range(4):
            for j in range(4):
                dyad = np.zeros((4, 4), complex)
                dyad[i, j] = 1
                np.testing.assert_allclose(
                    measure.transfer(q, 8, dyad),
                    joint @ dyad @ joint.conj().T,
                    atol=1e-13,
                )
                erased = np.zeros_like(dyad)
                erased[0, 0] = int(i == j)
                np.testing.assert_allclose(
                    measure.transfer(q, 0, dyad), erased, atol=1e-13
                )


def test_material_statistics_certify_the_formed_source():
    result = measure.run()
    remote = result["remote"]
    assert [r["rank"] for r in remote["finite_rank_certificates"]] == [16, 16]
    assert [r["highest_order"] for r in remote["finite_rank_certificates"]] == [20, 10]
    assert remote["changed_readout_rank"] == 1
    proof = remote["exact_single_port_minor"]
    assert int(proof["determinant"]) % 1009 == 429
    assert int(proof["determinant"]) % 1013 == 128
    held = remote["held_return"]["settings"]
    assert all(
        row["confidence_interval"][0] < 0 < row["confidence_interval"][1]
        for row in held[:3]
    )
    assert held[-1]["confidence_interval"][1] < 0
    for row in remote["finite_population_inversions"]:
        assert max(row["reconstruction_distance"].values()) < 1e-4
    gated = result["gated"]
    assert gated["family_confidence"] == 0.99
    assert abs(gated["formed_minus_independent_ideal_nats"] + 0.29384843) < 1e-7
    for row in gated["held"]:
        lo, hi = row["confidence_interval"]
        assert lo < gated["formed_minus_independent_ideal_nats"] < hi < 0
        assert row["changed_law_max_population_contrast"] < 1e-13
    archived = json.loads(Path(__file__).with_name("results.json").read_text())
    assert proof == archived["remote"]["exact_single_port_minor"]
    for row, stored in zip(gated["held"], archived["gated"]["held"]):
        np.testing.assert_allclose(
            row["confidence_interval"], stored["confidence_interval"], atol=1e-5
        )
