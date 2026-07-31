from __future__ import annotations

import json
from fractions import Fraction

import pytest

import angular_fingerprint_extension as fx


def test_receipt_builds_with_expected_status() -> None:
    receipt = fx.build_receipt()
    assert receipt["status"] == (
        "EXACT_INVARIANT_SUPPORT_FINGERPRINT__CONDITIONAL_KILL_RULES_TYPED"
    )
    assert receipt["comparison_boundary"]["comparison_permitted"] is False
    support = receipt["certificates"]["support_alignment"]
    assert support["even_zero_levels"] == [2, 4, 8, 14]
    assert support["aligned_on_window"] is True


def test_known_weights_are_exact() -> None:
    sequence = fx.equal_port_comb(fx.WINDOW)
    assert sequence[6] == Fraction(11, 25)
    assert sequence[10] == Fraction(247, 1875)
    assert sequence[12] == Fraction(1071, 3125)
    assert sequence[16] == Fraction(9424, 46875)
    assert sequence[30] == Fraction(3066051913, 18310546875)
    assert all(sequence[level] == 0 for level in range(1, fx.WINDOW + 1, 2))


def test_independent_sympy_recomputation() -> None:
    sympy = pytest.importorskip("sympy")
    x = sympy.sqrt(5) / 5
    sequence = fx.equal_port_comb(fx.WINDOW)
    for level in (6, 14, 16, 22, 30, 40):
        exact = sympy.nsimplify(
            (1 + (-1) ** level) / sympy.Integer(12)
            + sympy.Rational(5, 12)
            * (sympy.legendre(level, x) + sympy.legendre(level, -x)),
            rational=False,
        )
        assert sympy.simplify(exact - sympy.Rational(sequence[level])) == 0


def test_invariant_dimensions_match_generating_function() -> None:
    sympy = pytest.importorskip("sympy")
    t = sympy.symbols("t")
    series = sympy.series(
        (1 + t**15) / ((1 - t**6) * (1 - t**10)), t, 0, fx.WINDOW + 1
    ).removeO()
    poly = sympy.Poly(series, t)
    for level in range(fx.WINDOW + 1):
        expected = int(poly.coeff_monomial(t**level))
        computed = (
            fx.even_invariant_dimension(level)
            if level % 2 == 0
            else fx.odd_invariant_dimension(level)
        )
        assert computed == expected, level


def test_parent_drift_fails_closed(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    tampered = json.loads(fx.PARENT_RECEIPT_PATH.read_text(encoding="utf-8"))
    tampered["equal_port_certificate"]["sequence"][6] = "12/25"
    path = tmp_path / "angular_template_receipt.json"
    path.write_text(json.dumps(tampered), encoding="utf-8")
    monkeypatch.setattr(fx, "PARENT_RECEIPT_PATH", path)
    with pytest.raises(fx.FingerprintError, match="disagrees with the pinned parent"):
        fx.build_receipt()


def test_committed_receipt_is_byte_exact() -> None:
    committed = fx.RECEIPT_PATH.read_bytes()
    assert committed == fx.canonical_json_bytes(fx.build_receipt())
