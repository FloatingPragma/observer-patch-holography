from __future__ import annotations

from fractions import Fraction

import pytest

import spin_six_universality_certificate as su


def test_receipt_builds_with_expected_status() -> None:
    receipt = su.build_receipt()
    assert receipt["status"] == (
        "SPIN_SIX_RESIDUE_UNIVERSAL_ON_INVARIANT_CARRIER_CLASS__AMPLITUDE_OPEN"
    )
    table = receipt["invariant_table"]["table"]
    assert [int(table[str(level)]) for level in range(7)] == [1, 0, 0, 0, 0, 0, 1]
    assert receipt["invariant_table"]["first_odd_invariant_level"] == 15
    boundary = receipt["universality"]["comparison_boundary"]
    assert boundary["comparison_permitted"] is False


def test_orbit_rows_all_land_on_the_same_i6() -> None:
    receipt = su.build_receipt()
    rows = {row["orbit"]: row for row in receipt["orbit_verification"]["orbits"]}
    assert set(rows) == {"vertex_12", "face_20", "edge_30"}
    assert all(row["i6_multiple_nonzero"] for row in rows.values())
    assert rows["vertex_12"]["m6_i6_multiple"] == "64/175+0*sqrt5"


def test_independent_sympy_invariant_dimensions() -> None:
    sympy = pytest.importorskip("sympy")
    theta = sympy.symbols("theta")

    def chi(level, angle):
        return sympy.sin((2 * level + 1) * angle / 2) / sympy.sin(angle / 2)

    for level in range(17):
        total = (
            (2 * level + 1)
            + 15 * chi(level, sympy.pi)
            + 20 * chi(level, 2 * sympy.pi / 3)
            + 12 * chi(level, 2 * sympy.pi / 5)
            + 12 * chi(level, 4 * sympy.pi / 5)
        )
        value = sympy.nsimplify(sympy.simplify(total / 60))
        assert value == su.invariant_dimension(level), level


def test_independent_sympy_dodecahedral_orbit_moments() -> None:
    sympy = pytest.importorskip("sympy")
    sqrt5 = sympy.sqrt(5)
    phi = (1 + sqrt5) / 2
    dirs = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                dirs.append(sympy.Matrix([sx, sy, sz]))
    for s1 in (1, -1):
        for s2 in (1, -1):
            dirs.append(sympy.Matrix([0, s1 / phi, s2 * phi]))
            dirs.append(sympy.Matrix([s1 / phi, s2 * phi, 0]))
            dirs.append(sympy.Matrix([s2 * phi, 0, s1 / phi]))
    assert len(dirs) == 20
    x, y, z = sympy.symbols("x y z", real=True)
    n = sympy.Matrix([x, y, z])
    r2 = x**2 + y**2 + z**2
    m2 = sympy.expand(sum(((d.T * n)[0]) ** 2 for d in dirs) / 3)
    m4 = sympy.expand(sum(((d.T * n)[0]) ** 4 for d in dirs) / 9)
    assert sympy.simplify(m2 - sympy.Rational(20, 3) * r2) == 0
    assert sympy.simplify(m4 - 4 * r2**2) == 0


def test_chebyshev_engine_matches_closed_form() -> None:
    # U_n(cos t) sin t = sin((n+1) t) at the golden angles
    assert su.chebyshev_u(0, su.q5(Fraction(1, 2))) == su.ONE
    assert su.chebyshev_u(2, su.ZERO) == su.q5(-1)
    # U_4(1/2) = -1 + 2*... exact: cos t = 1/2 -> t = pi/3, sin5t/sint = ...
    val = su.chebyshev_u(4, su.q5(Fraction(1, 2)))
    assert val == su.q5(-1)


def test_committed_receipt_is_byte_exact() -> None:
    committed = su.RECEIPT_PATH.read_bytes()
    assert committed == su.base.canonical_json_bytes(su.build_receipt())


def test_parent_pin_matches_fingerprint_receipt() -> None:
    receipt = su.build_receipt()
    pin = receipt["parent_pins"][0]
    payload = su.base.RECEIPT_PATH.read_bytes()
    assert pin["sha256"] == su.base.tagged_sha256(payload)
