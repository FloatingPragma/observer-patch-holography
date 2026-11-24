"""Native operator checks and independent reconstruction of signed wall measurements."""

import json
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.interval_echo.openwave.quark_winding.written_junction import (
    conversion,
    exchange,
    field_probe,
    junction,
    loop_colour,
    orientations,
    relaxation,
    polarization,
    ring,
)
from model.executable_model import carrier

HERE = Path(__file__).parent


def saved(name):
    return json.loads((HERE / (name + ".json")).read_text())


def compare(a, b):
    if isinstance(b, dict):
        assert a.keys() == b.keys()
        for k, v in b.items():
            compare(a[k], v)
    elif isinstance(b, list):
        for u, v in zip(a, b, strict=True):
            compare(u, v)
    elif isinstance(b, float):
        assert a == pytest.approx(b, abs=2e-8, rel=0)
    else:
        assert a == b


@pytest.mark.parametrize(
    "name,module",
    [
        ("polarization", polarization),
        ("field", field_probe),
        ("exchange", exchange),
        ("junction", junction),
        ("conversion", conversion),
    ],
)
def test_record(name, module):
    compare(json.loads(json.dumps(module.run())), saved(name))


def test_direct_probe_channel():
    # Use the exact contact and trace the outgoing sample bits, retaining both detector rails.
    x, c = law.local_term_parts()
    h = x + 8 * c
    basis = [law.SRC.index((i, j)) for i, j in product((0, 1), repeat=2)]
    n = [1, 0, 1, 1]
    b = [0, 0, 1, 0]
    f = np.array([1, 0.2, 0.8, 0.45])
    theta = 2.7
    rho = np.zeros((4, 4), complex)
    rho[1:3, 1:3] = 0.5
    for a, r, w in zip(n, b, f):
        u = expm(-1j * w * theta / 17 * h)[np.ix_(basis, basis)].reshape(2, 2, 2, 2)
        new = np.zeros_like(rho)
        for i, j in product(range(2), repeat=2):
            k = np.kron(u[i, :, a, :], u[j, :, r, :])
            new += k @ rho @ k.conj().T
        rho = new
    predicted = polarization.probe(np.array(n), np.array(b), f, theta, 8, 1)["outcomes"]
    for key, contrast in [("x", rho[2, 1].real), ("y", rho[2, 1].imag)]:
        avg = (rho[1, 1] + rho[2, 2]).real / 2
        np.testing.assert_allclose(
            [avg + contrast, avg - contrast, rho[0, 0].real, rho[3, 3].real],
            predicted[key],
            atol=1e-12,
        )


def test_native_apparatus_and_held_phase():
    j = saved("polarization")
    f = saved("field")
    board = carrier.federation(2)
    g = nx.Graph(board.seams)
    for path in (j["sample_path"], j["bulk_path"]):
        assert all(g.has_edge(a, b) for a, b in zip(path, path[1:]))
        g.remove_edges_from(zip(path, path[1:]))
    assert nx.is_connected(g)
    assert (
        j["cases"]["1"]["source"]["aligned_heads_hex"]
        == j["cases"]["-2"]["source"]["aligned_heads_hex"]
    )
    for shift, q in [(1, -1 / 3), (-2, 2 / 3)]:
        clocks = j["cases"][str(shift)]["source"]["density_clocks"]
        assert (np.argmax(clocks["6"]) - np.argmax(clocks["-6"])) % 3 == 1
        assert (np.argmax(clocks["126"]) - np.argmax(clocks["102"])) % 3 == 2
        for sign in (1, -1):
            fit = f["held"]["samples"][f"{shift}_{sign}"]
            lo, hi = fit["simultaneous_interval"]
            assert lo < q < hi
            assert all(
                sum(counts) == f["held"]["shots"] for counts in fit["counts"].values()
            )
            left = next(
                r
                for r in f["cases"][str(shift)]["reads"]
                if r["width"] == 32
                and r["region"] == "left"
                and r["offset"] == 0
                and r["probe"] == sign
            )
            assert abs(sign * left["quarter_turn"]["phase"] / (np.pi / 2) - q) < 0.0003
            for rd in left["full_turn"]["quadratures"].values():
                assert min(rd) >= 0 and sum(rd) == pytest.approx(1)


def test_ring_against_explicit_local_columns():
    length = 6
    geom = (np.arange(length) % 3 == 0).astype(int)
    g = geom - geom + np.roll(geom, 1)
    w, heads, n = ring.basis(g)
    # Independent enumeration of every possible joint head configuration.
    allowed = [
        sum(b << i for i, b in enumerate(bits))
        for bits in product((0, 1), repeat=length)
        if set(g + bits - np.roll(bits, 1)) <= set((0, 1))
    ]
    assert sorted(allowed) == w.tolist()
    h0, h1 = law.local_term_parts()
    complete = h0 + 8 * h1
    basis = {}
    for coin, a, b in product(range(2), repeat=3):
        v = np.zeros(99)
        for d in (-1, 1):
            v[9 + 9 * law.LABELS.index((coin, d)) + law.SRC.index((1 + a, 1 + b))] = (
                2**-0.5
            )
        basis[coin, a, b] = v
    expected = np.zeros((len(w), len(w)))
    ix = {int(word): i for i, word in enumerate(w)}
    for i, word in enumerate(w):
        for s in range(length):
            initial = basis[heads[i, s], n[i, s], n[i, (s + 1) % length]]
            expected[i, i] += initial @ complete @ initial
            other = ix.get(int(word) ^ 1 << s)
            if other is not None:
                final = basis[heads[other, s], n[other, s], n[other, (s + 1) % length]]
                expected[other, i] += final @ complete @ initial
    actual = ring.generator(
        w, n, np.zeros(length, int), False
    ).toarray() + 8 * length * np.eye(len(w))
    np.testing.assert_allclose(actual, expected, atol=3e-13)


def test_static_rings_and_conversion():
    result = saved("ring")
    assert len(result) == 18
    for r in result:
        if r["length"] == 18 and r["law"] == "native":
            assert r["stationary"]["extrema"]["tight"]["maximum"] < 0.16
        if r["extra_unit_at"] is None:
            for v in r["separation"].values():
                np.testing.assert_allclose(
                    np.diff([p["contact_change"] for p in v]),
                    4 if r["law"] == "native" else -4,
                )
    j = saved("conversion")
    for a, b in zip(j["panels"][:24], j["panels"][24:48], strict=True):
        assert a["input"] == b["input"] and a["reflected"] != b["reflected"]
        assert a["records"] == b["records"]
    j = saved("junction")["cases"]["8"]
    assert [d["commutant_dimension"] for d in j["memory"].values()] == [4, 2, 4]
    assert j["memory"]["1"]["right_material_plus_coin_commutes"]
    costs = {
        n: {r[n]["panels"]["geodesic"]["blank_energy"] for r in j["distances"].values()}
        for n in ("10", "11")
    }
    assert costs == {"10": {488.0}, "11": {521.5}}


def test_joint_loop_words_and_release():
    loops = loop_colour.run()
    compare(loops, saved("loops"))
    for p in loops["panels"].values():
        assert p["without_face_flip"]["components"] == 1
        steps = p["history_refinements"]
        assert len(steps) == 5 and steps[-1]["classes"] == p["joint_states"]
        assert steps[-1]["largest"] == 1 and steps[-2]["largest"] > 1
        histories = p["witness"]["material_masks"]
        a, b = histories["native"]
        assert a[:-1] == b[:-1] and a[-1] != b[-1]
        assert histories["head_independent"][0] == histories["head_independent"][1]
    board, cycles = ring.physical_cycles()
    for r in saved("relaxation"):
        assert r["quadratic"] == pytest.approx([16 / 3, -16 / 3])
        if r["length"] == 12:
            extension = r["law"] == "positive_laplacian"
            calc = relaxation.panel(12, 6 - board.degree[cycles[12]], extension, True)
            compare(calc, r)
    # Verify stationary degeneracy optimization against an independent dense observable in one sector.
    length = 12
    g = np.roll((np.arange(length) % 3 == 0).astype(int), 1)
    w, _, loads = ring.basis(g)
    h = ring.generator(w, loads, 6 - board.degree[cycles[length]], False)
    row = next(
        r
        for r in saved("ring")
        if r["length"] == 12 and r["extra_unit_at"] is None and r["law"] == "native"
    )
    arcs = np.array(
        [ring.gap_arc(tuple(np.flatnonzero(n)), length, 1, 2) for n in loads], float
    )
    e, v = np.linalg.eigh(h.toarray())
    for label, sel in [("any", np.isfinite(arcs)), ("tight", arcs == np.nanmin(arcs))]:
        maxima = []
        for group in np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1e-8) + 1):
            p = v[:, group].T @ (sel[:, None] * v[:, group])
            maxima.append(np.linalg.eigvalsh(p)[-1])
        assert max(maxima) == pytest.approx(
            row["stationary"]["extrema"][label]["maximum"], abs=1e-9
        )


def test_rotated_written_wall_components():
    result = orientations.run()
    compare(result, saved("orientations"))
    rotation = np.array(result["rotation"])
    ports = np.arange(12)
    assert np.array_equal(rotation[rotation[rotation]], ports)
    assert np.all(ports != rotation)
    for panel in result["panels"].values():
        backgrounds = np.array(panel["backgrounds"])
        assert np.array_equal(backgrounds[:, rotation], np.roll(backgrounds, 1, axis=0))
        assert len({tuple(b) for b in backgrounds}) == 3
        assert panel["degree_symmetric_mismatches"] == [0, 0, 0]
        assert min(panel["changed_owner_mismatches"][1:]) > 0
    assert result["largest_diagonal_difference"] == 0
    assert result["determinant_rotated_bright"] == pytest.approx(1)
    assert result["mixture_rotated_bright"] == pytest.approx(0.5)
