"""Independent checks of released rings, escape, flux and simultaneous charges."""

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.quark_winding.written_junction import (
    backgrounds,
    colour_census,
    coupled_field,
    escape_census,
    holonomy,
    joint_field,
    mass_profile,
    multiplicity,
    reference_scan,
    reservoir,
    ring,
)

HERE = Path(__file__).parent


def saved(name):
    return json.loads((HERE / (name + ".json")).read_text())


def equal(a, b):
    if isinstance(b, dict):
        assert a.keys() == b.keys()
        for key in b:
            equal(a[key], b[key])
    elif isinstance(b, list):
        for x, y in zip(a, b, strict=True):
            equal(x, y)
    elif isinstance(b, float):
        assert a == pytest.approx(b, abs=3e-8, rel=1e-10)
    else:
        assert a == b


@pytest.mark.parametrize(
    "name,module",
    [
        ("backgrounds", backgrounds),
        ("escape_census", escape_census),
        ("joint_field", joint_field),
        ("reference_scan", reference_scan),
        ("holonomy", holonomy),
        ("multiplicity", multiplicity),
        ("colour_census", colour_census),
        ("mass_profile", mass_profile),
        ("reservoir", reservoir),
        ("coupled_field", coupled_field),
    ],
)
def test_archives(name, module):
    args = () if name in ("backgrounds", "escape_census") else (saved("backgrounds"),)
    equal(json.loads(json.dumps(module.run(*args))), saved(name))


def test_all_small_ring_boundary_ranks():
    # Construct the outside-by-inside matrix directly in every binary sector.
    # Counting each sharp word and computing a dense SVD is independent of rotation
    # orbits, the bit test for mobile Gauss sectors and the singleton elimination.
    for n in (3, 4):
        length = 9
        result = escape_census.panel(length, n)
        count, kernel, no_exit = 0, 0, 0
        for sites in combinations(range(length), n):
            g = backgrounds.bits(backgrounds.mask(sites), length)
            data = backgrounds.block(g, np.zeros(length))
            sharp = np.isfinite(data[3])
            direct = ring.generator(data[0], data[2], np.zeros(length), False)
            mat = direct[~sharp][:, sharp].toarray()
            count += int(sharp.sum())
            no_exit += int((abs(mat).sum(axis=0) == 0).sum())
            if mat.size:
                kernel += int(sharp.sum()) - np.linalg.matrix_rank(mat)
            else:
                kernel += int(sharp.sum())
        assert (count, kernel, no_exit) == (
            result["pin_free_words"],
            result["boundary_nullity"],
            result["pin_free_without_exit"],
        )


def test_flat_and_curved_controls():
    d = saved("holonomy")
    for row in d["controls"]:
        assert row["ordinary_spectrum_shift"] > 1
        assert row["ordinary_port_probability_shift"] > 0.15
        assert row["active_face_fourth_moment_difference"] == pytest.approx(
            row["single_triangle_loop_prediction"]
        )
    for rows in d["written"].values():
        assert all(
            x["conjugacy_residual"] < 1e-10 and x["port_probability_residual"] < 1e-10
            for x in rows
        )


def test_held_reference_and_stationary_read():
    d = saved("coupled_field")["extension"][1]
    assert d["hop_components"]["components"] == 1
    assert d["indistinguishability"][-1]["classes"] == 1478
    for tag, src in (("neutral", "reservoir_full"), ("positive", "target_spoke")):
        line = d[tag]
        prob = line["full_read"][src][0]
        assert prob == pytest.approx(0.80895213 if tag == "neutral" else 0.09451105)
        assert sum(sum(x) for x in line["full_read"].values()) == pytest.approx(1)
        for core, q in (("positive", 2 / 3), ("negative", -1 / 3)):
            m = line["field"][tag + "/" + core]
            assert m["phase_units"] == pytest.approx([q, q])
            y = np.asarray(m["signed_joint"])
            np.testing.assert_allclose(y.sum(axis=-1), 1, atol=1e-10)
            np.testing.assert_allclose(y[:, :, 4], 1 - prob, atol=1e-10)
            counts = np.array(m["counts"])
            np.testing.assert_equal(counts.sum(axis=-1), 2**19)
            derived = (counts[:, :, 0] - counts[:, :, 1]) / 2**19
            np.testing.assert_allclose(derived, y[:, :, 0] - y[:, :, 1], atol=0.006)
            np.testing.assert_allclose(m["head_independent_imaginary"], 0, atol=1e-10)


def test_reservoir_initial_evolution():
    s = saved("backgrounds")["paired"]["extension"]
    row = s["selected"][2]
    length = row["length"]
    site = s["additional_site"]
    _, spokes, edges, faces, geo = reservoir.patch(length, site)
    g = backgrounds.bits(row["gauss"], length)
    w, n, entry, tight = reservoir.configurations(g, spokes)
    assert len(faces) == 1 and len(spokes) == 2
    ref = np.array(geo["reference"])
    h = reservoir.generator(w, n, edges, faces, ref, True)
    _, state = backgrounds.terminal(backgrounds.block(g, ref[:-1]), True)
    source = np.pad(state, (0, len(w) - len(state)))
    dist = abs(expm_multiply(-1j * h / 16, source)) ** 2
    rec = reservoir.outcomes(dist, entry, tight, site)
    equal(rec, saved("reservoir")["laws"]["extension"][1]["terminal"][1]["read"])
    # Independent local divergence from oriented physical edges. Sum all tip
    # indicators at each port and compare the difference with its excess load.
    heads = (w[:, None] >> np.arange(len(edges))) & 1
    incoming = np.zeros((len(w), length + 1), int)
    for e, (a, b) in enumerate(edges):
        incoming[:, a] += heads[:, e]
        incoming[:, b] += 1 - heads[:, e]
    gfull = ((n[:, None] >> np.arange(length + 1)) & 1) - incoming
    assert (gfull == gfull[0]).all()


def test_reference_energy_and_radial_charge():
    r = saved("mass_profile")["extension"][1]
    assert r["neutral_minus_positive"] == pytest.approx(76.57137253)
    assert r["difference_slope"] == pytest.approx(0.968926673)
    for name, q in (("neutral", 0), ("positive", 1)):
        line = r[name]
        assert line["measured_slope"] == pytest.approx(
            line["bulk_reservoir_probability"], abs=1e-6
        )
        assert min(x["overlap"] for x in line["tracked"]) > 0.999
        for profile in line["tight_all_centres"]:
            assert sum(profile) == pytest.approx(q, abs=1e-9)


def test_exact_eigenvalue_bound():
    from scipy import sparse
    import sympy as sp

    h = sparse.csr_matrix([[0, 1, 0], [1, 0, 3], [0, 3, 2]], dtype=float)
    for copies in (1, 2, 3):
        m = sparse.block_diag([h] * copies, format="csr")
        assert multiplicity.certificate(m)["distinct_energies_at_least"] == 3
        assert colour_census.capacity(m) == copies
    mat = np.array([[1, 2, 3, 4], [5, 8, 9, 10], [10, 16, 18, 20]])
    assert colour_census.modular_rank(mat) == sp.Matrix(mat).rank() == 2
