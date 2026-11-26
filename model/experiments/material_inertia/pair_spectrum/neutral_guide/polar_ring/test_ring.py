"""Independent finite-word and released-record checks."""

from pathlib import Path
import json

import numpy as np
from scipy.linalg import eigh
from . import fringes, measure, source, stock
from model.executable_model import federation

HERE = Path(__file__).parent


def test_static_geometry_and_physical_costs():
    charts = stock.charts()
    assert [(r["cells"], r["requested_owners"]) for r in charts if not r["stock"]] == [
        (20, 6)
    ]
    for row in charts:
        if not row["stock"]:
            continue
        s = row["stock"]
        board = federation({20: 0, 80: 1, 320: 2}[row["cells"]])
        ring = s["ring"]
        assert board.degree[s["exterior_unit"]] == 5 and s["exterior_unit"] not in sum(
            map(list, ring), []
        )
        small = set(s["private"])
        n = len(ring)
        material = []
        for word in range(1 << n):
            loads = {}
            for j, pair in enumerate(ring):
                upper = word >> j & 1
                for p in pair:
                    loads[p] = upper if p in small else 1 - upper
            cost = 8 * sum(loads[p] for p in small) + sum(
                32 * (1 - loads[a] * loads[b]) for a, b in s["connectors"]
            )
            local = (
                8 * word.bit_count()
                + 16 * (word ^ ((word << 1) % (1 << n) | (word >> (n - 1)))).bit_count()
            )
            material.append(cost - local)
        assert len(set(material)) == 1
    m = measure.meter()
    assert {(v["private"], v["gateway"]) for v in m.values()} == set(np.ndindex(3, 3))
    assert [
        (m[k]["signed_change"], m[k]["squared_change"])
        for k in ("01", "10", "00", "02", "11", "12")
    ] == [(0, 0), (0, 0), (-1, 1), (1, 1), (1, -1), (2, 0)]


def test_microscopic_momentum_blocks():
    # Small independent full Hilbert space, directly flipping every dimer.
    n = 6
    h = np.zeros((1 << n, 1 << n))
    for i in range(1 << n):
        b = [i >> j & 1 for j in range(n)]
        h[i, i] = 8 * sum(b) + 16 * sum(x != b[j - 1] for j, x in enumerate(b))
        for j in range(n):
            h[i, i ^ (1 << j)] = -0.5
    direct, vectors = eigh(h)
    initial = np.eye(1 << n)[:, 1]
    raw = vectors @ (np.exp(-13j * direct) * (vectors.T @ initial))
    _, reduced = measure.spectrum(n, 8.0)
    recovered = sum(
        trans @ (v @ (np.exp(-13j * e) * local)) for trans, e, v, local, _ in reduced
    )
    np.testing.assert_allclose(raw, recovered, atol=3e-12)
    np.testing.assert_allclose(
        sorted(np.concatenate([b[1] for b in reduced])), direct, atol=1e-12
    )
    rng = np.random.default_rng(49)
    complex_profile = rng.normal(size=n) + 1j * rng.normal(size=n)
    complex_profile /= np.linalg.norm(complex_profile)
    formed = source.prepared(complex_profile)
    np.testing.assert_allclose(formed[1 << np.arange(n)], complex_profile, atol=2e-14)
    assert np.count_nonzero(formed) == n


def test_recoil_and_original_launch_denominators():
    result = json.loads((HERE / "records.json").read_text())
    clocks = json.loads((HERE / "fringes.json").read_text())
    seeds = sorted([p["prefix"] for p in result["panels"]])
    assert len(set(seeds)) == len(seeds) == 10
    with (
        np.load(HERE / "word_counts.npz") as original,
        np.load(HERE / "fringe_counts.npz") as moving,
    ):
        all_native = []
        for p, clock in zip(result["panels"], clocks["panels"], strict=True):
            n = clock["n"]
            lines = [q["host_gap"] for q in clock["clocks"]]
            lower = clock["sharp_lower_return"]
            for line in clock["clocks"]:
                product = [
                    fringes.amplitude(x) * fringes.amplitude(y).conjugate()
                    for x, y in zip(lower, line["plus_counts"], strict=True)
                ]
                assert (
                    abs(fringes.unwrap(measure.HOLDS, product) - line["host_gap"])
                    < 1e-10
                )
            for j, panel in enumerate(clock["fringes"]):
                count = moving[panel["counts_key"]]
                assert np.all(count.sum(axis=1) == measure.SHOTS)
                read = fringes.recoil_from_words(count, n)
                if j == 0:
                    assert read["max_nonconstant"] < 0.005
                else:
                    assert read["mode"] == 1 and read["min_amplitude"] > 0.2
                    assert (
                        abs(
                            (-1) ** (j - 1) * read["signed_recoil"]
                            - (lines[j] - lines[0])
                        )
                        < 0.0004
                    )
            for flight in p["flight_summary"]:
                counts = original[f"{p['prefix']}_t{flight['time']}"]
                assert np.all(counts.sum(axis=1) == measure.SHOTS)
                remote = counts[:, 1 << np.arange(2, n)].sum(axis=1)
                np.testing.assert_array_equal(
                    remote, flight["unique_upper_other_owner"]
                )
            native = p["kappa"] == 8
            if native:
                assert p["minimum_one_island_weight"] > 0.98
                assert (
                    0.81
                    < p["flight_summary"][3]["unique_upper_other_owner"][0]
                    / measure.SHOTS
                    < 0.82
                )
                assert abs(lines[0] - 39.9132) < 0.0001
                assert 0.998 < clock["ground_line_weight"] < 1
                all_native.append((n, lines[1] - lines[0]))
            else:
                assert abs(lines[0] - 19.8301) < 0.0001
        training = next(x for n, x in all_native if n == 10)
        forecast = (
            training * (1 - np.cos(2 * np.pi / 12)) / (1 - np.cos(2 * np.pi / 10))
        )
        assert all(abs(forecast / r - 1) < 0.15 for n, r in all_native if n == 12)
