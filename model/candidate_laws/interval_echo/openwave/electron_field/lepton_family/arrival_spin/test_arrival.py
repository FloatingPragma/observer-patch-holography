"""Guard the independent flux trigger, relative phase and detector eligibility."""

import json
from pathlib import Path
import numpy as np

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    tomography,
)

HERE = Path(__file__).resolve().parent


def test_herald_matches_archived_family_and_spin_can_change_fringe():
    report = json.loads((HERE / "results.json").read_text())
    archive = json.loads((HERE.parent / "dirac_read/results.json").read_text())[
        "physical_degree_six_path"
    ]
    for sector in report["sectors"]:
        for j, source in enumerate(archive["modes"]):
            comparison = {
                x["time"]: x["distal_charge_same_character"]
                for x in source["two_positive_energy_mode_source_reads"]
            }
            for row in sector["records"][f"family_{j}_positive"]:
                p = comparison[row["time"]]
                data = row["fringe"]
                np.testing.assert_allclose(
                    data["charge"]["max_over_axes"], [p, 0, 1 - p], atol=2e-12
                )
                np.testing.assert_allclose(
                    data["supplied_spin"]["min_over_axes"], [0, p, 1 - p], atol=2e-12
                )


def test_count_inversion_keeps_complex_coherences():
    # The Y-positive state has rho[0,1] = -i/2. Integer counts here have no sampling error.
    half = tomography.SHOTS // 2
    record = {
        "sharp_and_miss": [half, half, 0],
        "pair_counts": [[[half, half, 0, 0], [2 * half, 0, 0, 0]]],
    }
    expected = np.outer(tomography.INPUTS[-2], tomography.INPUTS[-2].conj())
    np.testing.assert_array_equal(
        tomography.invert(record, [(0, 1)]), expected.round(14)
    )


def test_detection_accepts_transmitted_quantum_information_and_rejects_misrouted_turns():
    records = json.loads((HERE / "circulation.json").read_text())["sectors"]
    counts = json.loads((HERE / "counts.json").read_text())["panels"]
    for model, panel in zip(records, counts, strict=True):
        assert model["sign"] == panel["sign"] and model["source_spur"] == panel["spur"]
        assert model["rank_preserving_half_turns"] == 0
        assert sum(n > 0 for n in model["closed_address_subspace_sizes"]) == 10
        assert not any(model["order_disagreements_on_closed_address_subspace"])
        read = model["panels"][0]["readings"][2]["source_qubit_to_receiver"]
        assert (
            abs(
                read["minimum_reference_partial_transpose_eigenvalue"]
                - panel["minimum_partial_transpose_eigenvalue"]
            )
            < panel["simultaneous_eigenvalue_error"]
        )
        assert panel["certifies_entanglement"] == bool(model["source_spur"])
        for src in panel["counts"]:
            assert sum(src["sharp_and_miss"]) == tomography.SHOTS
            assert all(
                sum(s) == tomography.SHOTS
                for settings in src["pair_counts"]
                for s in settings
            )
