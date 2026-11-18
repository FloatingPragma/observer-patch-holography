"""Independent controls for uniform refinement and real nonlinear matter."""
from copy import deepcopy
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import whitney_real_continuum as producer
import verify_whitney_real_continuum as verifier
import whitney_interacting_quantum as full_action


@pytest.fixture(scope="module")
def packet():
    return verifier.load()


def test_receipt_replays_independently(packet):
    result = verifier.verify(packet)
    assert result["accepted"] is True
    assert result["conditional_real_sector_trajectory_bound"] is True
    assert result["tetrahedra"] == [20, 160, 1280, 10240]
    for key in ("full_charged_complex_trajectory_bound", "numerical_trajectory_error_certified",
                "physical_source_or_clock_selected", "formalized_in_lean"):
        assert result[key] is False


def test_finite_meshes_retain_scaled_shapes_and_ritz_order(packet):
    rows = packet["mesh_checks"]
    assert all(row["diameter_inradius_max"] < 10 for row in rows)
    assert max(row["h_max_times_n"] for row in rows[1:])-min(row["h_max_times_n"] for row in rows[1:]) < 1e-12
    for coarse, fine in zip(rows, rows[1:]):
        assert fine["ritz_h1_error"] < .7*coarse["ritz_h1_error"]
        assert fine["ritz_l2_error"] < .6*coarse["ritz_l2_error"]


@pytest.mark.parametrize("n", [1, 2])
def test_dyadic_children_lie_in_unique_parent(n):
    vertices, cells, _ = producer.mesh(n)
    fine_vertices, fine_cells, _ = producer.mesh(2*n)
    x = vertices[cells]
    affine = np.concatenate((np.ones(x.shape[:-1]+(1,)), x), axis=-1)
    inverse = np.linalg.inv(affine)
    centers = fine_vertices[fine_cells].mean(axis=1)
    lam = np.einsum("fc,tci->fti", np.column_stack((np.ones(len(centers)), centers)), inverse)
    inside = np.all(lam > 1e-10, axis=2)
    assert np.all(inside.sum(axis=1) == 1)
    parent = inside.argmax(axis=1)
    all_vertices = np.concatenate((np.ones(fine_vertices[fine_cells].shape[:-1]+(1,)), fine_vertices[fine_cells]), axis=-1)
    bary = np.einsum("fic,fcj->fij", all_vertices, inverse[parent])
    assert bary.min() > -1e-12 and bary.max() < 1+1e-12


def test_cubic_force_matches_independent_simplex_tensor_and_variation():
    vertices, cells, _ = producer.mesh(2)
    _, _, volume = producer.geometry(vertices, cells)
    rng = np.random.default_rng(713)
    u, direction = rng.normal(size=len(vertices)), rng.normal(size=len(vertices))
    fourth = np.array([float(verifier.moment(indices)) for indices in product(range(4), repeat=4)]).reshape((4,)*4)
    local = np.einsum("ijkl,tj,tk,tl,t->ti", fourth, u[cells], u[cells], u[cells], volume, optimize=True)
    independent = np.bincount(cells.ravel(), weights=local.ravel(), minlength=len(vertices))
    force = producer.cubic_load(u, cells, volume)
    np.testing.assert_allclose(force, independent, atol=3e-15, rtol=3e-15)
    def potential(t):
        v = u+t*direction
        return v@producer.cubic_load(v, cells, volume)/4
    eps = 1e-4
    derivative = (potential(-2*eps)-8*potential(-eps)+8*potential(eps)-potential(2*eps))/(12*eps)
    assert derivative == pytest.approx(force@direction, abs=3e-11)
    # Vertex lumping loses the nonlocal mixed polynomial moments.
    mass, _ = producer.matrices(vertices, cells)
    assert np.linalg.norm(mass@(u**3)-force) > .1


def test_real_sector_uses_every_full_action_variation():
    mesh = full_action.geometry(4)
    rng = np.random.default_rng(763)
    scalar = rng.normal(size=13)
    velocity = rng.normal(size=13)
    direction_a = rng.normal(size=42)
    direction_v = rng.normal(size=42)
    direction_imaginary = rng.normal(size=13)
    tangent = np.r_[np.zeros(42), velocity, np.zeros(13)]
    def lagrangian(eps):
        metric, potential = full_action.coefficients(eps*direction_a,
            scalar+1j*eps*direction_imaginary, .7, mesh=mesh)
        speed = tangent+eps*np.r_[direction_v, np.zeros(26)]
        return speed@metric@speed/2-potential
    eps = 1e-4
    assert abs((lagrangian(eps)-lagrangian(-eps))/(2*eps)) < 1e-9
    metric, _ = full_action.coefficients(np.zeros(42), scalar, .7, mesh=mesh)
    gauss = full_action.vertical(scalar, .7, mesh).T@metric@tangent
    np.testing.assert_allclose(gauss, 0, atol=2e-14)
    # A complex scalar velocity exits the neutral real sector.
    changed = tangent.copy(); changed[55:] = scalar
    assert np.linalg.norm(full_action.vertical(scalar, .7, mesh).T@metric@changed) > 1


def test_free_boundary_ritz_load_cannot_drop_affine_flux():
    # This is a projection test, not a Neumann continuum-wave reference.
    # An affine function has zero volume Laplacian but nonzero boundary flux.
    vertices, cells, _ = producer.mesh(2)
    mass, stiffness = producer.matrices(vertices, cells)
    affine = .5+vertices[:, 0]
    full_weak_load = (stiffness+mass)@affine
    right = spsolve(stiffness+mass, full_weak_load)
    wrong = spsolve(stiffness+mass, mass@affine)
    np.testing.assert_allclose(right, affine, atol=2e-14, rtol=2e-14)
    assert np.linalg.norm(stiffness@affine) > 1
    error = wrong-affine
    assert error@(stiffness+mass)@error > 1


@pytest.mark.parametrize("keys,value", [
    (("analytic_scope", "full_charged_complex_trajectory_bound"), True),
    (("analytic_scope", "numerical_trajectory_error_certified"), True),
    (("analytic_scope", "physical_source_or_clock_selected"), True),
    (("analytic_scope", "formalized_in_lean"), True),
    (("analytic_scope", "conditional_real_sector_trajectory_bound"), 1),
    (("analytic_scope", "continuum_reference_existence_assumed"), False),
    (("uniform_bounds", "shape_ratio_upper"), "71"),
    (("uniform_bounds", "inradius_times_n_lower"), "2/24"),
    (("mesh_checks", 0, "n"), 1.0),
    (("mesh_checks", 1, "vertices"), 55.0),
    (("mesh_checks", 2, "tetrahedra"), 1280.00000000001),
    (("mesh_checks", 3, "boundary_triangles"), True),
    (("mesh_checks", 0, "mesh_sha256"), "0"*64),
    (("mesh_checks", 0, "ritz_h1_error"), 0.),
])
def test_false_promotions_and_changed_diagnostics_fail(packet, keys, value):
    changed = deepcopy(packet); node = changed
    for key in keys[:-1]:
        node = node[key]
    node[keys[-1]] = value
    with pytest.raises(ValueError):
        verifier.verify(changed)


def test_changed_source_bytes_are_not_cached(packet, monkeypatch):
    original = Path.read_bytes
    target = verifier.ROOT/"code/electromagnetism/whitney_real_continuum.py"
    def changed(path):
        value = original(path)
        return value+b"\n" if path == target else value
    monkeypatch.setattr(Path, "read_bytes", changed)
    with pytest.raises(ValueError, match="source pin"):
        verifier.verify(packet)


def test_changed_wave_endpoint_fails_independent_integration(packet):
    changed = deepcopy(packet)
    changed["trajectories"][0]["final_state"][0] += .001
    with pytest.raises(ValueError, match="independent nonlinear wave endpoint"):
        verifier.verify(changed)


@pytest.mark.parametrize("text", ['{"x":0,"x":1}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e9999}'])
def test_noncanonical_json_rejected(tmp_path, text):
    path = tmp_path/"invalid.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        verifier.load(path)


def test_source_reads_are_explicit_utf8(packet, monkeypatch):
    original = Path.read_text
    def read(path, *args, **kwargs):
        assert kwargs.get("encoding") == "utf-8"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", read)
    assert verifier.load()["scope"] == verifier.SCOPE
    verifier.source_geometry()
