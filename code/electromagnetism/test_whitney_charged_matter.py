"""Exact regressions for the dressed scalar action and charged cone data.

Only the registered vertex coordinates are read from Lean. Geometry, element
gradients, integrals and gauge variations are reconstructed here; no producer,
receipt or claimed theorem output supplies the expected equalities. These
tests do not assert a simulated matter trajectory or a spatial error bound.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import re

import pytest
import sympy as s

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = (1 + s.sqrt(5)) / 2


def zero(value):
    return s.simplify(s.expand(value)) == 0


def antisymmetric(values, size=4):
    result = s.zeros(size)
    for (i, j), value in zip(combinations(range(size), 2), values, strict=True):
        result[i, j], result[j, i] = value, -value
    return result


@lru_cache(maxsize=1)
def cone_geometry():
    source = (ROOT / "Lean/Screen/SeamCurrentEdge30Moment.lean").read_text(encoding="utf-8")
    definition = source.split("def portVector", 1)[1].split(
        "theorem portVector_positivePort", 1)[0]
    literals = {"0": s.Integer(0), "1": s.Integer(1), "-1": s.Integer(-1),
                "φ": GOLDEN, "-φ": -GOLDEN}
    boundary = [s.Matrix([literals[token.strip()] for token in row.split(",")])
                for row in re.findall(r"!\[([^\[\]]+)\]", definition)]
    assert len(boundary) == 12 and len({tuple(v) for v in boundary}) == 12
    vertices = [s.zeros(3, 1)] + boundary
    boundary_edges = [pair for pair in combinations(range(1, 13), 2)
                      if zero((vertices[pair[1]]-vertices[pair[0]]).norm()**2 - 4)]
    edge_set = set(boundary_edges)
    triangles = [face for face in combinations(range(1, 13), 3)
                 if all(pair in edge_set for pair in combinations(face, 2))]
    assert len(boundary_edges) == 30 and len(triangles) == 20
    assert all(sum(vertex in face for face in triangles) == 5 for vertex in range(1, 13))
    edges = [(0, vertex) for vertex in range(1, 13)] + boundary_edges
    elements = []
    for face in triangles:
        tet = (0,) + face
        rays = s.Matrix.hstack(*(vertices[v]-vertices[0] for v in face))
        determinant = s.simplify(rays.det())
        volume = s.simplify(abs(determinant)/6)
        inverse = rays.inv().applyfunc(s.simplify)
        gradients = [-(inverse[0, :]+inverse[1, :]+inverse[2, :])]
        gradients.extend(inverse[j, :] for j in range(3))
        assert volume > 0
        elements.append((tet, volume, gradients))
    return vertices, edges, elements


@lru_cache(maxsize=1)
def charged_witness():
    _, edges, elements = cone_geometry()
    charge = s.Symbol("charge", real=True, nonzero=True)
    beta = charge*(2+3*GOLDEN)/10
    potential = [beta] + [s.Integer(0)]*12
    electric = {(i, j): potential[j]-potential[i] for i, j in edges}
    rates = [3] + [-1]*12
    load, gauss = s.zeros(13, 1), s.zeros(13, 1)
    x, y, z = s.symbols("x y z", real=True)
    lam = [1-x-y-z, x, y, z]
    for tet, volume, grad in elements:
        field = s.zeros(1, 3)
        for i, j in combinations(range(4), 2):
            edge = (tet[i], tet[j])
            coefficient = electric[edge] if edge in electric else -electric[edge[::-1]]
            field += coefficient*(lam[i]*grad[j]-lam[j]*grad[i])
        field = field.applyfunc(s.simplify)
        # Reconstruct the field from oriented Whitney edge coefficients,
        # before comparing with the gradient representation.
        assert all(zero(v) for v in field-beta*grad[0])
        for i, vertex in enumerate(tet):
            # Exact degree-two simplex moments integrate 2e lambda_i f_h.
            load[vertex] += 2*charge*sum(
                rates[other]*volume*s.Rational(1+int(i == j), 20)
                for j, other in enumerate(tet))
            # Independently integrate grad(lambda_i) dot E, i.e. D^T M E.
            gauss[vertex] += volume*(grad[i]*field.T)[0]
    return charge, beta, electric, load.applyfunc(s.simplify), gauss.applyfunc(s.simplify)


def test_exact_charged_initial_data_on_all_thirteen_vertices():
    _, _, elements = cone_geometry()
    charge, beta, electric, load, gauss = charged_witness()
    volume = GOLDEN**2/3
    assert all(zero(v-volume) for _, v, _ in elements)
    assert zero(load[0]-6*charge*volume)
    assert all(zero(load[v]+charge*volume/2) for v in range(1, 13))
    assert zero(sum(load)) and any(not zero(v) for v in load)
    assert all(zero(value+beta) for (i, _), value in electric.items() if i == 0)
    assert all(value == 0 for (i, _), value in electric.items() if i != 0)
    assert all(zero(value) for value in gauss-load)


def test_gauss_detects_electric_sign_and_missing_matter_factor():
    _, _, _, load, gauss = charged_witness()
    assert any(not zero(value) for value in -gauss-load)
    assert any(not zero(value) for value in gauss-load/2)


@pytest.mark.parametrize("vertex", range(4))
def test_all_generic_straight_whitney_path_integrals(vertex):
    a = antisymmetric(s.symbols("a01 a02 a03 a12 a13 a23", real=True))
    lam = s.symbols("lambda0:4", real=True)
    t = s.Symbol("t", real=True)
    path = [int(i == vertex)+t*(lam[i]-int(i == vertex)) for i in range(4)]
    derivative = [s.diff(value, t) for value in path]
    pullback = sum(a[i, j]*(path[i]*derivative[j]-path[j]*derivative[i])
                   for i, j in combinations(range(4), 2))
    integral = s.integrate(s.expand(pullback), (t, 0, 1))
    assert zero(integral-sum(a[vertex, j]*lam[j] for j in range(4)))


def test_constant_nodal_field_dressing_time_cancellation_is_special():
    adot = antisymmetric(s.symbols("v01 v02 v03 v12 v13 v23", real=True))
    lam = s.symbols("lambda0:4", real=True)
    phase_dot = adot*s.Matrix(lam)
    assert zero(sum(lam[i]*phase_dot[i] for i in range(4)))
    # The same chain-rule term cannot be dropped for general nodal matter.
    chosen = antisymmetric([1, 0, 0, 0, 0, 0])
    point = s.ones(4, 1)/4
    psi = [1, 2, 3, 4]
    omitted = s.I*sum(point[i]*psi[i]*(chosen*point)[i] for i in range(4))
    assert omitted == -s.I/16


def test_pointwise_covariant_time_derivative_and_chain_rule_control():
    lam = s.Matrix([s.Rational(i, 10) for i in (1, 2, 3, 4)])
    psi = s.Matrix([1+s.I, 2-s.I, 3, -1+2*s.I])
    velocity = s.Matrix([2, -s.I, 1+s.I, -3])
    eta = s.Matrix([1, -2, 3, 4])
    phi = s.Matrix([s.Rational(i, 7) for i in (1, 2, -1, 3)])
    adot = antisymmetric([1, -2, 3, 4, 0, -1])
    transformed_adot = s.Matrix(4, 4, lambda i, j: adot[i, j]+eta[j]-eta[i])
    scalar = (lam.T*psi)[0]

    def covariant(vel, edge_vel, scalar_potential, dressing=True):
        result = (lam.T*vel)[0] + s.I*(lam.T*scalar_potential)[0]*scalar
        if dressing:
            result += s.I*sum(lam[i]*psi[i]*(edge_vel*lam)[i] for i in range(4))
        return s.expand(result)

    # chi(t)=t eta: at t=0 fields coincide, while their velocities change.
    transformed_velocity = velocity+s.I*s.matrix_multiply_elementwise(eta, psi)
    assert zero(covariant(transformed_velocity, transformed_adot, phi-eta)
                - covariant(velocity, adot, phi))
    assert not zero(covariant(transformed_velocity, transformed_adot, phi-eta, False)
                    - covariant(velocity, adot, phi, False))
    assert not zero(covariant(transformed_velocity, transformed_adot, phi+eta)
                    - covariant(velocity, adot, phi))


def test_spatial_gauge_covariance_and_undressed_counterexample():
    x, y, z = s.symbols("x y z", real=True)
    coords, lam = (x, y, z), [1-x-y-z, x, y, z]
    chi, psi = [0, s.pi, 0, 0], [1, 2, 3, 4]
    chi_h = sum(l*c for l, c in zip(lam, chi, strict=True))
    transformed_edges = s.Matrix(4, 4, lambda i, j: chi[j]-chi[i])
    transformed = sum(lam[i]*s.exp(s.I*(transformed_edges*s.Matrix(lam))[i])
                      *s.exp(s.I*chi[i])*psi[i] for i in range(4))
    original = sum(lam[i]*psi[i] for i in range(4))
    assert zero(transformed-s.exp(s.I*chi_h)*original)
    for coordinate in coords:
        derivative = s.diff(transformed, coordinate)-s.I*s.diff(chi_h, coordinate)*transformed
        assert zero(derivative-s.exp(s.I*chi_h)*s.diff(original, coordinate))
    undressed = sum(lam[i]*s.exp(s.I*chi[i])*psi[i] for i in range(4))
    midpoint = {x: s.Rational(1, 2), y: 0, z: 0}
    assert not zero((undressed-transformed).subs(midpoint))


def test_actual_common_face_value_and_tangential_covariant_derivative():
    vertices, _, elements = cone_geometry()
    left, right = next((first, second) for first, second in combinations(elements, 2)
                       if len(set(first[0]) & set(second[0])) == 3)
    face = sorted(set(left[0]) & set(right[0]))
    tangent = vertices[face[1]]-vertices[face[0]]

    def edge(i, j):
        return s.Rational((i+2*j) % 5-2, 7) if i < j else -edge(j, i) if i > j else 0

    def trace(element):
        tet, _, gradients = element
        lam = [s.Rational(1, 3) if vertex in face else 0 for vertex in tet]
        a = s.Matrix(4, 4, lambda i, j: edge(tet[i], tet[j]))
        psi = [vertex+1+s.I*(vertex % 3) for vertex in tet]
        derivatives = [s.simplify((g*tangent)[0]) for g in gradients]
        phases = a*s.Matrix(lam)
        scalar = sum(lam[i]*s.exp(s.I*phases[i])*psi[i] for i in range(4))
        ordinary = sum(s.exp(s.I*phases[i])*psi[i]*(derivatives[i]+s.I*lam[i]
                       *sum(a[i, j]*derivatives[j] for j in range(4))) for i in range(4))
        potential = sum(a[i, j]*(lam[i]*derivatives[j]-lam[j]*derivatives[i])
                        for i, j in combinations(range(4), 2))
        return scalar, ordinary-s.I*potential*scalar

    left_value, left_tangent = trace(left)
    right_value, right_tangent = trace(right)
    assert zero(left_value-right_value)
    assert zero(left_tangent-right_tangent)


def test_coordinate_reconstruction_uses_utf8_with_windows_default(monkeypatch):
    original = Path.read_text
    def windows_default(path, *args, **kwargs):
        if not args and "encoding" not in kwargs:
            kwargs["encoding"] = "cp1252"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", windows_default)
    cone_geometry.cache_clear()
    try:
        vertices, edges, elements = cone_geometry()
        assert (len(vertices), len(edges), len(elements)) == (13, 42, 20)
    finally:
        cone_geometry.cache_clear()
