import gzip
from io import BytesIO
import json
import math

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from .aperture import flights
from .fragments import regional_meter, assignments
from .material import Assembly, bound, competing_orbit
from .measure import HERE
from .probe import decode, sampler_radius


def test_fragment_floors_and_matching_host():
    assert bound(2, 1) == (112, 119.5)
    assert bound(2, 2) == (136, 144)
    assert bound(4, 3) == (200, 209)
    assert bound(4, 0)[1] + bound(0, 3)[1] == 185.5
    competitor = competing_orbit()
    np.testing.assert_array_equal(competitor["gram"], np.eye(3))
    np.testing.assert_array_equal(competitor["doubled_total_energy"], 371 * np.eye(3))
    a = Assembly(2, 1)
    # A direct swap-difference generator is independent of sparse construction.
    psi = np.arange(len(a.w)) % 7
    direct = a.v * psi + sum((psi - psi[s]) / 2 for s in a.swaps)
    np.testing.assert_array_equal(a.h @ psi, direct)


def test_scalar_clock_against_held_coupled_flights():
    h = np.array([[209, 1.5], [1.5, 328]])
    seed = np.array([[1], [0]], complex)
    e, u = eigh(h)
    times = (0, 1)
    result, _ = flights(sparse.csr_matrix(h), (200, 350), seed, 209, 8, times)
    clock = np.exp(-1j * (e[:, None] - 209) * np.arange(32) / 256).mean(axis=1) ** 8
    for j, t in enumerate(times):
        for selected in (0, 1):
            ideal = u @ (
                (np.exp(-1j * t * e) * clock**selected)[:, None] * (u.T @ seed)
            )
            np.testing.assert_allclose(result[2 * j + selected], ideal, atol=3e-13)


def test_assignment_precedes_sign_and_phase_reads_resolve_all_charges():
    w = np.array([[1, 1, -1, -1, -1, 0, 0, 1, 0, 0, 1, 0]])
    # A sign reversal on identical occupied geometry reverses all regional charges.
    weights = assignments(w)
    np.testing.assert_allclose(weights[:, 0], assignments(-w)[::-1, 0])
    assert np.isclose(weights[:, 0].sum(), 1)
    assert regional_meter() == {str(q): (130 * q) % 31 for q in range(-4, 5)}


def test_sampled_signs_and_contamination_on_the_selected_sources():
    r = json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    g = r["observed_unit_gain"]
    pol = r["cap"]["banks"]
    assert decode(pol["circular_plus"]["centroid_0"], g) > sampler_radius(
        pol["circular_plus"]["centroid_0"], g
    )
    assert decode(pol["axial_0"]["centroid_0"], g) < -sampler_radius(
        pol["axial_0"]["centroid_0"], g
    )
    f = json.loads((HERE / "fragments.json").read_text())
    assert f["regional_bins"] == regional_meter()
    words = np.load(
        BytesIO(gzip.decompress((HERE / "joined_words.npy.gz").read_bytes())),
        allow_pickle=False,
    )
    assert words.shape == (8, 13, 27720)
    for j, runs in enumerate(f["runs"]):
        np.testing.assert_array_equal(
            words[2 * j : 2 * j + 2].sum(axis=-1), np.array(runs["joint_probe"])[:, :-1]
        )
    n = f["original_shots"]
    observed = f["calibration"]["counts"][2]
    g = (observed[0] - observed[1]) / sum(observed)

    def eps(count):
        return math.sqrt(2 * math.log(2 * 4096 / 0.01) / count)

    amplitude = max(
        math.sin(39 / 16) / (32 * math.sin(39 / 512)), 1 / (32 * math.sin(math.pi / 32))
    )
    high = 4.5 / (4.5 + 39**2) * amplitude**16
    assert high < 3.362e-11
    for row in f["runs"]:
        assert sum(row["charge_and_refusals"]) == n
        scale = row["weight_scale"]
        for direction, bank in zip((1, -1), row["joint_probe"]):
            assert sum(bank) == n
            p, m, zero = bank[6:9]
            kept = p + m + zero
            value = 12 * scale * (p - m) / kept / g
            radius = (12 * scale * eps(kept) + abs(value) * eps(sum(observed))) / (
                g - eps(sum(observed))
            )
            if row["contrast"] == 8:
                source = 1 if row["source"] == "sharp" else 1 / 3
                yield_floor = (kept / n - eps(n) / 2 - 1e-10) / source
                assert yield_floor > 0.00053
                drift = 2 * 7 * scale * math.sqrt(high / yield_floor)
                assert drift < 0.0089
                assert direction * value > radius + drift + 1e-6
            else:
                assert direction * value > radius + 0.15
