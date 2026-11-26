import json

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import carrier

from .. import material
from . import colour, measure, readback


def test_standing_sources_are_a_compression_of_actual_rail_and_degree_costs():
    n = 6
    board = carrier.federation(0)
    path = np.array(material.route(0)[: 2 * (n + 2)]).reshape(-1, 2)
    words = material.bits(n)
    u = np.pad(words, ((0, 0), (1, 1)), constant_values=0)
    loads = np.where(
        (board.degree[path] == 5)[None, :, :], u[:, :, None], 1 - u[:, :, None]
    )
    marked = 32 * (1 - loads[:, :-1, 1] * loads[:, 1:, 0]).sum(axis=1)
    referenced = 8 * loads[:, board.degree[path] == 5].sum(axis=1)
    diagonal = marked + referenced
    h = np.diag((diagonal - diagonal[0]).astype(float))
    for j in range(n):
        h[np.arange(1 << n), np.arange(1 << n) ^ (1 << j)] = -0.5
    np.testing.assert_allclose(h, material.core(n + 2, 0, 8), atol=1e-13)


def test_registered_temporal_bins_against_single_finite_pointer():
    e = np.array([-0.037, 39.985])
    v = np.array([[0.6, 0.8], [-0.8, 0.6]])
    output = measure.raw_clock(e, v, {"a": 0})
    for b in [0, 1, 638, 639, 640, 8191]:
        ts = 2 * np.pi * np.arange(measure.CLOCK) / measure.SPAN
        amps = np.exp(-1j * np.outer(e - measure.SPAN * b / measure.CLOCK, ts)).mean(
            axis=1
        )
        np.testing.assert_allclose(
            output[0, b], abs(v[0]) ** 2 @ abs(amps) ** 2, atol=3e-13
        )
    ts = 2 * np.pi * np.arange(measure.FILTER) / measure.SPAN
    np.testing.assert_allclose(
        measure.aperture(e, 40),
        np.exp(-1j * np.outer(e - 40, ts)).mean(axis=1) ** 4,
        atol=3e-13,
    )


def test_spectral_comparison_and_complete_banks():
    data = json.loads((measure.HERE / "records.json").read_text())
    with np.load(measure.HERE / "records.npz") as banks:
        rows = readback.decode(data, banks)
    for calc, row in zip(data["panels"], rows, strict=True):
        a, b = calc["band"].values()
        if row["mode"] == "native":
            assert (
                b["vacuum_relative_range"][1]
                < calc["far_pair_conditioned_min"]
                < calc["two_separated_lowest_products"]
            )
            assert b["vacuum_relative_range"][1] / a["vacuum_relative_range"][0] < 1.205
            assert (
                a["single_pair_length_weight_lower"] > 0.98
                and b["single_pair_length_weight_lower"] > 0.969
            )
            assert row["conversion_both"]["1_near"] > 10**7
            assert row["conversion_both"]["1_parked"] == 0
            assert row["flights"]["1"][-2][1] > 14_000_000
        elif row["mode"] == "degree_off":
            assert a["states"] == b["states"] == 0
        else:
            assert all(
                pair[1] == 0 for data in row["flights"].values() for pair in data
            )


def test_conjugate_detector_identity_and_exact_seam_block():
    effect = np.zeros((9, 9), complex)
    for basis in colour.colour_bases():
        for v in basis.T:
            z = np.kron(v, v.conj())
            effect += np.outer(z, z.conj())
    bell = np.eye(3).ravel() / np.sqrt(3)
    np.testing.assert_allclose(effect, np.eye(9) + 3 * np.outer(bell, bell), atol=2e-14)
    w = law._write_pattern()
    rem = 2 * np.eye(90, dtype=int) - w @ w.T
    r, d2 = law._label_ops()
    for g in (np.eye(10, dtype=int) - r, d2):
        assert not np.any(w.T @ np.kron(g, np.eye(9, dtype=int)) @ rem)
    data = json.loads((measure.HERE / "colour.json").read_text())
    obs = readback.read_colour(data)
    assert abs(obs["routed_02"]["0"] - 3 / 64) < data["joint_99_radius"]
    assert (
        max(obs["local_inverse_02"].values()) == max(obs["written_11"].values()) == 0.0
    )
    assert data["colour_checks"]["pairs_checked"] == 45_360
    for rec in data["independent_singlet_tag"]["records"].values():
        assert all(sum(c) == data["shots_per_setting"] for c in rec["counts"])
        assert rec["inferred_singlet"] == 1.0
