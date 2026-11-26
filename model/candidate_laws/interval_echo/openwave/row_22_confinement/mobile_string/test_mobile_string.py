"""Rebuild the material from native ports, then check the held receiver bank."""

import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import commission, material as mat, measure


def physical_sector(ports, board, g, h):
    n = len(ports) // 2
    sites = np.asarray(ports).reshape(n, 2)
    z = mat.bits(n)
    # Multiply the marked connector costs on their actual port occupations.
    load = np.where(
        board.degree[sites][None, :, :] == 5, z[:, :, None], 1 - z[:, :, None]
    )
    cost = 4 * mat.J * (1 - load[:, :-1, 1] * load[:, 1:, 0]).sum(axis=1)
    cost += g * load[:, board.degree[sites] == 5].sum(axis=1)
    idx = np.arange(1 << n)
    generator = sparse.diags(cost, format="csr")
    for i in range(n):
        generator += sparse.csr_matrix(
            (np.full(len(idx), -h), (idx, idx ^ (1 << i))), shape=generator.shape
        )
    return generator


def test_material_and_calibrations():
    recorded = json.loads((measure.HERE / "commission.json").read_text())
    assert json.loads(json.dumps(commission.report())) == recorded
    np.testing.assert_allclose(mat.calibration_energies(), [32, 1, 8], atol=1e-12)
    for lev in range(3):
        board = mat.federation(lev)
        ports = mat.route(lev)[:16]
        graph = physical_sector(ports, board, mat.G, mat.H)
        reduced = mat.release_generator(8, mat.G, mat.H)
        difference = (graph - reduced).toarray()
        np.testing.assert_allclose(difference, difference[0, 0] * np.eye(256))
        words = mat.bits(8)
        for length in range(2, 7):
            readings = []
            # Both external pins lie on actual material, with the same
            # unrestricted outside spins included in the two minima.
            for pin in (0, 1):
                eligible = (words[:, [0, length + 1]] == 0).all(axis=1) & (
                    words[:, [1, length]] == pin
                ).all(axis=1)
                compress = graph[eligible, :][:, eligible].toarray()
                readings.append(np.linalg.eigvalsh(compress)[0])
            np.testing.assert_allclose(
                readings[1] - readings[0], mat.static_laws(length, mat.G)[0], atol=1e-10
            )


def test_unpinned_words_against_physical_short_flight():
    board = mat.federation(0)
    ports = mat.route(0)[:28]
    initial = np.zeros((1 << 14, 2))
    mask = sum(1 << i for i in range(3, 11))
    initial[[mask, ((1 << 14) - 1) ^ mask], [0, 1]] = 1
    with np.load(measure.HERE / "receivers.npz") as bank:
        for name, (g, h) in measure.FREE.items():
            graph = physical_sector(ports, board, g, h)
            psi = expm_multiply(-1j * mat.RELEASE[1] * graph, initial)
            np.testing.assert_allclose(
                (abs(psi) ** 2).T, bank["word_law_" + name][1], rtol=8e-8, atol=2e-12
            )


def test_receiver_enrollment():
    result = measure.decode()
    assert json.loads(json.dumps(result)) == json.loads(
        (measure.HERE / "results.json").read_text()
    )
    for panel in result["levels"].values():
        native, off = panel["native"], panel["degree_off"]
        assert abs(native["slope_2_to_3"] - 8) < 0.01
        for row in (native, off):
            assert (
                max(abs(x) for x in row["held_intact_line_minus_extrapolated"].values())
                < 0.001
            )
            for s in row["static"].values():
                assert (
                    abs(
                        s["calculated_ground_work"]
                        - s["calibrated_zero_swap_prediction"]
                    )
                    < 0.3
                )
        for length, s in native["static"].items():
            reader = s["sharp_intact" if length < 6 else "sharp_screened"]
            assert reader["ground_product_weight"] > 0.989
            assert abs(reader["line"] - s["calculated_ground_work"]) < 0.001
            assert reader["min_visibility"] > 0.97
            assert abs(off["static"][length]["calculated_ground_work"] - 32) < 1e-10
        assert native["static"][6]["sharp_intact"]["ground_product_weight"] < 0.0001
        for name in measure.FREE:
            release = panel[name]["release"]
            joint = [s[1]["count"] / measure.LAUNCHES for s in release[1]]
            if name == "native":
                assert joint[0] - joint[1] > 0.001
            elif name == "degree_off":
                assert abs(joint[0] - joint[1]) < 0.0001
            else:
                assert all(ln[1]["count"] == 0 for t in release for ln in t)


def test_ambient_ruler_all_embeddings():
    from . import ruler

    audit = ruler.audit()
    assert json.loads(json.dumps(audit)) == json.loads(
        (measure.HERE / "ruler.json").read_text()
    )
    assert (
        sum(len(data["positions"]) for level in audit for data in level.values()) == 234
    )
    assert all(
        6 in {s["ambient"] for s in audit[0][length]["positions"]} for length in (4, 5)
    )
    assert all(
        s["ambient"] <= data["material"]
        for level in audit
        for data in level.values()
        for s in data["positions"]
    )
