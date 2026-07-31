from __future__ import annotations

from fractions import Fraction

import pytest

import differentiator_inventory_certificate as di


def test_receipt_builds_with_nine_typed_rows() -> None:
    receipt = di.build_receipt()
    assert receipt["status"] == (
        "DIFFERENTIATOR_INVENTORY_CERTIFIED__NINE_ROWS_TYPED"
    )
    assert len(receipt["rows"]) == 9
    assert receipt["comparison_boundary"]["comparison_permitted"] is False
    names = [row["differentiator"] for row in receipt["rows"]]
    assert "no linear vacuum dispersion" in names
    assert "five-fold azimuthal uniqueness" in names
    assert "global-form fork and charge commensuration" in names


def test_recomputed_facts_are_present_and_exact() -> None:
    receipt = di.build_receipt()
    rows = {row["row"]: row for row in receipt["rows"]}
    assert rows[1]["recomputed_facts"]["m_6"] == 1
    assert rows[2]["recomputed_facts"]["vanishing_odd_moments"] == [1, 3, 5, 7]
    assert rows[4]["recomputed_facts"]["weights"]["6"] == "11/25"
    assert rows[6]["recomputed_facts"]["order"] == 6
    assert rows[7]["recomputed_facts"]["port_branch_su2_su3_ratio"] == "6"
    assert rows[8]["recomputed_facts"]["recomputed_center_in_window"] is True


def test_z6_congruence_detects_a_wrong_row(monkeypatch: pytest.MonkeyPatch) -> None:
    original = di.fact_z6_kernel

    def tampered():
        rows = {
            "Q": (1, 1, Fraction(1, 6)),
            "u_c": (2, 0, Fraction(-1, 3)),
        }
        for name, (t, d, y) in rows.items():
            exponent = (2 * t + 3 * d + 6 * y) % 6
            di.require(exponent == 0, f"Z6 kernel drift on {name}")
        return {}

    monkeypatch.setattr(di, "fact_z6_kernel", tampered)
    with pytest.raises(di.base.FingerprintError):
        di.build_receipt()


def test_committed_receipt_is_byte_exact() -> None:
    committed = di.RECEIPT_PATH.read_bytes()
    assert committed == di.base.canonical_json_bytes(di.build_receipt())
