from __future__ import annotations

from fractions import Fraction

import pytest

import carrier_scale_bound_diagnostic as cb


def test_receipt_builds_with_expected_status() -> None:
    receipt = cb.finalize(cb.build_receipt())
    assert receipt["status"] == (
        "EXPOSED_RETROSPECTIVE_CARRIER_SCALE_BOUND__DIAGNOSTIC_ONLY"
    )
    boundary = receipt["comparison_boundary"]
    assert boundary["public_measurement_read"] is True
    assert boundary["comparison_permitted"] is False
    assert boundary["scored"] is False
    assert boundary["fz11_untouched"] is True


def test_bound_value_matches_exact_identification() -> None:
    receipt = cb.finalize(cb.build_receipt())
    assert receipt["bound"]["carrier_scale_upper_bound_m"] == "6.788e-27"
    assert receipt["bound"]["planck_length_headroom"] == "4.199e+08"
    # independent recomputation with floats
    import math

    a_gev = math.sqrt(20) / 1.3e11
    a_m = a_gev * 1.97327e-16
    assert abs(a_m - 6.788e-27) / a_m < 1e-3
    assert abs(a_m / 1.616255e-35 - 4.199e8) / 4.199e8 < 1e-3


def test_no_linear_term_statement_is_true_of_the_stencil() -> None:
    # the stencil expansion carries even powers of (a k . u) only
    import a5_multipole_fixed_point_certificate as base

    verts = base.cartesian_vertices()
    for k in (1, 3, 5, 7):
        assert base.p_is_zero(base.moment_sum(verts, k))


def test_committed_receipt_is_byte_exact() -> None:
    committed = cb.RECEIPT_PATH.read_bytes()
    rebuilt = cb.finalize(cb.build_receipt())
    assert committed == cb.base.canonical_json_bytes(rebuilt)


def test_parent_pin_drift_fails_closed(monkeypatch: pytest.MonkeyPatch, tmp_path) -> None:
    import json

    tampered = json.loads(cb.base.RECEIPT_PATH.read_text())
    tampered["kinetic_stencil_conditional"]["expansion"] = "k^2 - (a^2/10) k^4"
    path = tmp_path / "parent.json"
    path.write_text(json.dumps(tampered))
    monkeypatch.setattr(cb.base, "RECEIPT_PATH", path)
    with pytest.raises(cb.base.FingerprintError):
        cb.build_receipt()
