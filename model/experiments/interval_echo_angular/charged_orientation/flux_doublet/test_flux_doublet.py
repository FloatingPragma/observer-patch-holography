"""Independent character, pulse, receiver and complete local-law checks."""

import json
import re
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh, expm

from . import bands, measure as m, passive


def test_covariant_rotations_and_public_source():
    basis, projector, (r, s) = m.charged_source()
    np.testing.assert_allclose(r @ s, -s @ r, atol=1e-12)
    np.testing.assert_allclose(projector, basis @ basis.conj().T, atol=1e-13)
    for g in (r, s):
        np.testing.assert_allclose(g @ projector, projector @ g, atol=1e-12)
    for owner in (0, 6):
        frame, _, turns = m.charged_source(owner=owner)
        np.testing.assert_allclose(frame, basis, atol=1e-12)
        for other, base in zip(turns, (r, s), strict=True):
            np.testing.assert_allclose(other, base, atol=2e-12)
    p, q = m.half_pairs()[0]
    at = Path(__file__).resolve().parents[5] / "Lean/SourceLaw/MagneticVacancy.lean"
    text = at.read_text()
    for name, expected in [
        ("turnR", p),
        ("turnS", q),
        ("phaseR", m.compensated_turn(p, m.links(1))),
        ("phaseS", m.compensated_turn(q, m.links(1))),
    ]:
        row = re.search(rf"def {name} .*? := !\[(.*?)\]", text)
        assert list(map(int, row.group(1).split(","))) == list(expected)


def test_inertia_merits_and_different_dual_root():
    assert m.character_census(m.links(1))["commutant_dimension"] == 3
    for q in (-1, 1):
        result = bands.band_counts(m.links(q, 7), m.CELL.seams)
        assert [v["counts_below"] for v in result["inertia"]] == [2, 6]
        h = m.induced_hamiltonian(m.links(q, 7))
        e = eigvalsh(h)
        np.testing.assert_allclose(e, np.repeat(e[[0, 2, 6]], [2, 4, 6]), atol=2e-12)
    levels = m.induced_hamiltonian(m.links(1))
    for k in (8, 32, 128):
        approx = m.split_steps(m.links(1), k)
        np.testing.assert_allclose(approx.conj().T @ approx, np.eye(12), atol=2e-10)
        error = np.linalg.norm(approx - expm(-1j * levels), 2)
        ceiling = 1.1 if k == 8 else 100 / k**2
        assert error < ceiling


def test_passive_finite_pulses_and_complete_acquisitions():
    # A finite crossing on a written seam maps one coin class into the other,
    # even when several intermediate generator entries contribute to the pulse.
    generator = sum(
        v * t for v, t in zip(m.LAW.local_term_parts(), (1, m.K), strict=True)
    )
    u = expm(-0.719j * generator)
    for i, j in [(1, 3), (5, 7), (3, 1), (7, 5)]:
        written = u[j + 9 :: 9, i + 9 :: 9]
        assert np.linalg.norm(written) > 0.1
        np.testing.assert_allclose(written @ written, 0, atol=1e-12)
    x = m.run()
    received = json.loads((m.HERE / "results.json").read_text())
    for panel, arch in zip(x["panels"], received["panels"], strict=True):
        assert panel["commutator_units"] == arch["commutator_units"]
        assert panel["ground_dimension"] == arch["ground_dimension"]
        for read, old in zip(panel["reads"], arch["reads"], strict=True):
            for key in ("bare", "compensated"):
                np.testing.assert_allclose(
                    read[key]["probability"], old[key]["probability"], atol=2e-12
                )
                c = np.array(old[key]["counts"])
                assert c.sum() == m.SHOTS
                np.testing.assert_allclose(
                    c / m.SHOTS, read[key]["probability"], atol=0.005
                )
    for bank in x["coherent_preparations"]["filters"]:
        assert sum(bank["counts"]) == m.SHOTS
        assert bank["band_fidelity"] > 0.99995
        assert max(bank["counts"][12:24]) == 0
    archive = x["coherent_preparations"]
    assert min(archive["fraction_in_ground_band"]) > 0.12
    np.testing.assert_allclose(
        sorted(archive["principal_overlap_squared"])[0], 0, atol=1e-12
    )
    assert np.isclose(sum(archive["principal_overlap_squared"]), 1 / 2)
    exact = passive.exact_columns()
    assert sum(row["checked_nonzero"] for row in exact) == 192
    for bank in x["passive_collar"]["triangle_reads"]:
        if bank["written"]:
            assert bank["spectrum_shift"] < 1e-11
            assert bank["max_port_probability_shift"] < 1e-11
        else:
            assert bank["max_port_probability_shift"] > 0.03
    assert x["aligned_doublet_exchange"]["reads"][0]["counts"] == [m.SHOTS, 0]
    assert x["aligned_doublet_exchange"]["reads"][-1]["counts"] == [0, m.SHOTS]


def test_all_binary_cargo_invariant_spaces():
    from math import comb
    from .cargo import census

    table = census()
    assert table["conjugacy_classes"] == 56
    assert table["constituent_dimensions"] == [1, 10, 44, 110, 165, 132]
    for k, dimensions in enumerate(table["by_double_count"]):
        assert sum(dimensions) == comb(11, k)
        assert dimensions.count(1) == 1 and 2 not in dimensions
