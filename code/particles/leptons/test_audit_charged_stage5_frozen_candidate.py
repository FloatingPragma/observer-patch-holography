"""Tests for the downstream frozen-candidate comparison audit."""

from __future__ import annotations

import json
from decimal import Decimal

import audit_charged_stage5_frozen_candidate as lane


def test_frozen_candidate_is_accurate_but_not_promoted():
    artifact = lane.build_audit(
        json.loads(lane.CANDIDATE.read_text(encoding="utf-8")),
        json.loads(lane.REFERENCE_JSON.read_text(encoding="utf-8")),
        json.loads(lane.TRUNCATED_D10_PROBE.read_text(encoding="utf-8")),
    )
    assert Decimal(artifact["max_absolute_relative_error_percent"]) < Decimal("0.03")
    assert artifact["candidate_mutation_allowed"] is False
    assert artifact["public_prediction_promotion_allowed"] is False
    assert artifact["status"] == "COMPARE_ONLY_RETRODICTIVE_ACCURACY_AUDIT"
    assert Decimal(
        artifact["pixel_branch_sensitivity"]["truncated_d10_comparison_probe"][
            "max_absolute_relative_error_percent"
        ]
    ) < Decimal("0.1")
    assert artifact["mass_scheme_audit"]["direct_pole_mass_comparison_theorem_validated"] is False
    assert set(artifact["ratios"]) == {
        "muon_over_electron",
        "tau_over_electron",
        "tau_over_muon",
    }
    assert max(
        abs(Decimal(row["relative_error_ppm"])) for row in artifact["ratios"].values()
    ) < Decimal("32")
    assert abs(
        Decimal(artifact["koide_invariant"]["reference_relative_to_two_thirds_ppm"])
    ) < Decimal("4")
