"""Tests for the CP-2 / P4 premise reduction machine check.

Asserts the certified verdicts of cp2_p4_premise_reduction.build(): the RC
identity for the inversion gauge, the named-axiom classification of every
deformation family, the P4 controls, determinism, and the no-discharge status
keys.
"""

from __future__ import annotations

import json
from pathlib import Path

import cp2_p4_premise_reduction as mod


def _families(cert):
    return {f["family"]: f for f in cert["part_A_cp2"]["A2_A3_gauge_class"]["families"]}


def test_build_is_deterministic():
    a = mod.build()
    b = mod.build()
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_rc_identity_for_inversion_gauge():
    cert = mod.build()
    a1 = cert["part_A_cp2"]["A1_read_consistency"]
    assert a1["pass"] is True
    fam = _families(cert)["inversion (phi = id)"]
    assert fam["rc_pass"] is True
    assert fam["p4_cp4_coherent"] is True


def test_every_deformation_family_violates_a_named_axiom_or_rc():
    cert = mod.build()
    fams = _families(cert)
    for name, fam in fams.items():
        if name == "inversion (phi = id)":
            continue
        assert fam["rc_pass"] is False, name
    # power and affine members satisfy every remaining axiom (live gauge
    # freedom) and sit off balance under the CP-4 carrier
    for name in [
        "power gauge, s = 1/2",
        "power gauge, s = exp(-P/24) (poisson reserve)",
        "power gauge, s = 1 - P/24 (presence reserve)",
        "power gauge, s = 0.99",
        "affine capacity gauge, a = 0.3",
        "affine capacity gauge, a = 0.7",
    ]:
        fam = fams[name]
        assert fam["seed_pass"] is True, name
        assert fam["coupled_banach_pass"] is True, name
        assert fam["p4_cp4_coherent"] is False, name
        assert float(fam["separation_from_balance_lo"]) > 0.5, name
    # the shift gauge violates the seed normalization
    assert fams["shift gauge"]["seed_pass"] is False


def test_pinned_gauge_witnesses_harmless_freedom():
    cert = mod.build()
    fam = _families(cert)["balance-pinned sine gauge"]
    assert fam["seed_pass"] is True
    assert fam["coupled_banach_pass"] is True
    assert fam["rc_pass"] is False
    assert float(fam["rc_residual_abs_lo"]) > 8.0
    assert fam["p4_cp4_coherent"] is True


def test_p4_controls():
    cert = mod.build()
    b = cert["part_B_p4"]
    assert b["B1_aggregate_support_collapse"]["pass"] is True
    assert b["B2_positive_control_cp4"]["pass"] is True
    assert b["B3_negative_control_location"]["pass"] is True
    assert b["B4_negative_control_sign_carrier"]["pass"] is True
    assert b["B5_recorded_lattice_control"]["pass"] is True
    assert float(b["B5_recorded_lattice_control"]["min_log10_separation"]) > 121


def test_no_discharge_claimed():
    cert = mod.build()
    assert cert["status"] == "reductions_certified_no_premise_discharged"
    assert cert["moves_cl7"] is False
    assert cert["cl7_status"] == "open"
    assert cert["blindness"]["cone_contains_cl3_bridge_expression"] is True


def test_artifact_matches_build():
    artifact = mod._CANON / "runtime" / "cp2_p4_premise_reduction_check.json"
    if not artifact.exists():
        return
    on_disk = json.loads(artifact.read_text(encoding="utf-8"))
    assert on_disk == mod.build()
