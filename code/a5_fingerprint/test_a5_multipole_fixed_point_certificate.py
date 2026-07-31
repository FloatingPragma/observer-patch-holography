from __future__ import annotations

from fractions import Fraction

import pytest

import a5_multipole_fixed_point_certificate as fp


def test_receipt_builds_with_expected_status() -> None:
    receipt = fp.build_receipt()
    assert receipt["status"] == (
        "EXACT_A5_FINGERPRINT_CERTIFICATE__PHYSICAL_MAP_OPEN"
    )
    census = receipt["critical_points"]["census"]
    assert census == {"maxima": 12, "minima": 20, "saddles": 30, "total": 62}
    boundary = receipt["decision_rules_and_ledger"]["comparison_boundary"]
    assert boundary["comparison_permitted"] is False
    assert receipt["fail_closed_controls"]["all_detectors_fired"] is True


def test_orbit_values_and_hessians_exact() -> None:
    receipt = fp.build_receipt()
    orbits = {row["orbit"]: row for row in receipt["critical_points"]["orbits"]}
    assert orbits["vertex_pole"]["value"] == "1+0*sqrt5"
    assert orbits["vertex_ring"]["value"] == "1+0*sqrt5"
    assert orbits["face_high"]["value"] == "-5/9+0*sqrt5"
    assert orbits["face_low"]["value"] == "-5/9+0*sqrt5"
    assert orbits["edge_high"]["value"] == "-5/16+0*sqrt5"
    assert orbits["edge_equator"]["value"] == "-5/16+0*sqrt5"
    assert orbits["face_high"]["hessian_eigenvalues"] == [
        "35/3+0*sqrt5",
        "35/3+0*sqrt5",
    ]
    assert orbits["vertex_pole"]["hessian_eigenvalues"] == [
        "-21+0*sqrt5",
        "-21+0*sqrt5",
    ]
    assert set(orbits["edge_equator"]["hessian_eigenvalues"]) == {
        "105/16+105/16*sqrt5",
        "105/16+-105/16*sqrt5",
    }


def test_independent_sympy_i6_and_orbit_values() -> None:
    sympy = pytest.importorskip("sympy")
    sqrt5 = sympy.sqrt(5)
    phi = (1 + sqrt5) / 2
    raw = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            raw.append(sympy.Matrix([0, s1, s2 * phi]))
            raw.append(sympy.Matrix([s1, s2 * phi, 0]))
            raw.append(sympy.Matrix([s2 * phi, 0, s1]))
    norm = sympy.sqrt(2 + phi)
    verts = [v / norm for v in raw]
    x, y, z = sympy.symbols("x y z", real=True)
    n = sympy.Matrix([x, y, z])
    p6 = sympy.legendre(6, sympy.Symbol("t"))
    total = sum(p6.subs(sympy.Symbol("t"), (v.T * n)[0]) for v in verts)
    i6 = sympy.Rational(25, 132) * total

    vertex_val = sympy.simplify(
        i6.subs({x: verts[0][0], y: verts[0][1], z: verts[0][2]})
    )
    assert vertex_val == 1

    face = (raw[0] + raw[1] + raw[4]) if False else None
    # face: three mutually adjacent vertices (dot 1/sqrt5)
    adj = []
    for i in range(12):
        for j in range(i + 1, 12):
            d = sympy.simplify((verts[i].T * verts[j])[0])
            if sympy.simplify(d - 1 / sqrt5) == 0:
                adj.append((i, j))
    triple = None
    for i, j in adj:
        for k in range(12):
            if k in (i, j):
                continue
            if (min(i, k), max(i, k)) in adj and (min(j, k), max(j, k)) in adj:
                triple = (i, j, k)
                break
        if triple:
            break
    assert triple is not None
    centre = verts[triple[0]] + verts[triple[1]] + verts[triple[2]]
    centre = centre / sympy.sqrt((centre.T * centre)[0])
    face_val = sympy.simplify(
        i6.subs({x: centre[0], y: centre[1], z: centre[2]})
    )
    assert sympy.simplify(face_val + sympy.Rational(5, 9)) == 0

    edge = verts[adj[0][0]] + verts[adj[0][1]]
    edge = edge / sympy.sqrt((edge.T * edge)[0])
    edge_val = sympy.simplify(i6.subs({x: edge[0], y: edge[1], z: edge[2]}))
    assert sympy.simplify(edge_val + sympy.Rational(5, 16)) == 0


def test_independent_sympy_moment_nulls() -> None:
    sympy = pytest.importorskip("sympy")
    sqrt5 = sympy.sqrt(5)
    phi = (1 + sqrt5) / 2
    raw = []
    for s1 in (1, -1):
        for s2 in (1, -1):
            raw.append(sympy.Matrix([0, s1, s2 * phi]))
            raw.append(sympy.Matrix([s1, s2 * phi, 0]))
            raw.append(sympy.Matrix([s2 * phi, 0, s1]))
    norm = sympy.sqrt(2 + phi)
    verts = [v / norm for v in raw]
    x, y, z = sympy.symbols("x y z", real=True)
    n = sympy.Matrix([x, y, z])
    r2 = x**2 + y**2 + z**2
    m2 = sympy.expand(sum(((v.T * n)[0]) ** 2 for v in verts))
    m4 = sympy.expand(sum(((v.T * n)[0]) ** 4 for v in verts))
    assert sympy.simplify(m2 - 4 * r2) == 0
    assert sympy.simplify(m4 - sympy.Rational(12, 5) * r2**2) == 0
    for k in (1, 3, 5):
        mk = sympy.expand(sum(((v.T * n)[0]) ** k for v in verts))
        assert sympy.simplify(mk) == 0


def test_independent_meridian_factorization_and_latitudes() -> None:
    sympy = pytest.importorskip("sympy")
    c = sympy.Symbol("c")
    p6 = sympy.legendre(6, c)
    lhs = sympy.expand(
        sympy.diff(p6, c) ** 2
        - sympy.Rational(441, 64) * (1 - c**2) ** 3 * (1 - 6 * c**2) ** 2
    )
    rhs = sympy.expand(
        sympy.Rational(441, 64)
        * (5 * c**2 - 1)
        * (5 * c**4 - 5 * c**2 + 1)
        * (45 * c**4 - 30 * c**2 + 1)
    )
    assert sympy.simplify(lhs - rhs) == 0
    sqrt5 = sympy.sqrt(5)
    for tval, factor in [
        (sympy.Rational(1, 5), 5 * c**2 - 1),
        ((5 + sqrt5) / 10, 5 * c**4 - 5 * c**2 + 1),
        ((5 - sqrt5) / 10, 5 * c**4 - 5 * c**2 + 1),
        ((5 + 2 * sqrt5) / 15, 45 * c**4 - 30 * c**2 + 1),
        ((5 - 2 * sqrt5) / 15, 45 * c**4 - 30 * c**2 + 1),
    ]:
        assert sympy.simplify(factor.subs(c**2, tval).subs(c**4, tval**2)) == 0


def test_independent_stencil_coefficients() -> None:
    sympy = pytest.importorskip("sympy")
    receipt = fp.build_receipt()
    row = receipt["kinetic_stencil_conditional"]["sixth_moment_relation"]
    # on the sphere: sum (u.n)^6 = alpha + beta I6 with the certified pair
    text = row["on_sphere"]
    assert "12/7" in text and "64/175" in text
    # cross-check: at a vertex the sum is 2 + 10/125, and the stencil
    # anisotropic coefficient is (64/175)/1440 = 2/7875
    assert Fraction(12, 7) + Fraction(64, 175) == Fraction(2) + Fraction(10, 125)
    assert Fraction(64, 175) / 1440 == Fraction(2, 7875)


def test_committed_receipt_is_byte_exact() -> None:
    committed = fp.RECEIPT_PATH.read_bytes()
    assert committed == fp.canonical_json_bytes(fp.build_receipt())


def test_vertex_tamper_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    original = fp.cartesian_vertices

    def tampered():
        verts = original()
        v = verts[0]
        verts[0] = (v[0], v[1], fp.q5_scale(v[2], Fraction(101, 100)))
        return verts

    monkeypatch.setattr(fp, "cartesian_vertices", tampered)
    with pytest.raises(fp.FingerprintError):
        fp.build_receipt()
