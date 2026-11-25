import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_30_collective import (
    experiment,
    packets,
)
from model.candidate_laws.interval_echo.openwave.row_30_collective.packet_speed import (
    measure,
)


def test_density_in_physical_dimer_basis():
    n = 5
    words = np.arange(2**n)
    z = 1 - 2 * ((words[:, None] >> np.arange(n)) & 1)
    ends = np.column_stack((np.ones(len(words)), z, -np.ones(len(words))))
    walls = (1 + ends[:, :-1] * ends[:, 1:]) / 2
    single = np.flatnonzero(walls.sum(axis=1) == 1)
    lookup = single[np.argsort(walls[single].argmax(axis=1))]
    x = np.arange(n + 1)
    for j in (8, 4):
        h = np.diag(
            j * (np.sum(z[:, 1:] * z[:, :-1], axis=1) + z[:, 0] - z[:, -1])
        ).astype(float)
        for i in range(n):
            h[words, words ^ (1 << i)] -= 0.5
        u, e, v = experiment.dual(n, j)
        for k in (-1.1, 0.2, 1.5):
            packet = np.exp(-0.3 * (x - 2.2) ** 2 + 1j * k * x)
            packet /= np.linalg.norm(packet)
            state = np.zeros(len(words), complex)
            state[lookup] = packet
            position = walls @ x
            commutator = 1j * (h * position[None, :] - position[:, None] * h)
            initial = np.vdot(state, commutator @ state).real
            np.testing.assert_allclose(
                initial, np.sum((packet[:-1].conj() * packet[1:]).imag), atol=1e-13
            )
            t = 2.17
            probability = abs(expm(-1j * t * h) @ state) ** 2
            np.testing.assert_allclose(
                probability @ walls, packets.occupations(packet, u, e, v, t), atol=2e-13
            )


def test_speed_bound_for_both_coupling_orderings():
    k = np.linspace(-np.pi, np.pi, 999)
    for j, h in [(8.0, 0.5), (4.0, 0.5), (0.5, 3.0)]:
        large, small = max(j, h), min(j, h)
        x = np.append(k, np.arccos(small / large))
        w2 = 4 * (j - h) ** 2 + 8 * j * h * (1 - np.cos(x))
        numer = 4 * j * h * np.sin(x)
        np.testing.assert_allclose(
            (2 * small) ** 2 * w2 - numer**2,
            16 * small**2 * (large * np.cos(x) - small) ** 2,
            atol=2e-13,
        )
        np.testing.assert_allclose(max(abs(numer) / np.sqrt(w2)), 2 * small)


def test_detector_counts():
    report = json.loads((measure.HERE / "results.json").read_text())
    with np.load(measure.HERE / "counts.npz") as raw:
        for r in report:
            base = f"n{r['dimers']}_p{r['kick_index']}"
            key = base + f"_j{int(r['contrast'])}"
            assert np.all(raw[base + "_wave"].sum(axis=-1) == measure.SHOTS)
            times, intervals = measure.clock(raw[key + "_clock"])
            velocity, band = measure.slope(times, intervals, raw[key + "_wall"])
            np.testing.assert_allclose(
                [velocity, *band], [r["velocity"], *r["velocity_band"]]
            )
    native = [r for r in report if r["contrast"] == 8]
    assert all(
        r["velocity_band"][0] < r["envelope_band_velocity"] < r["velocity_band"][1]
        for r in native
    )
    assert max(r["velocity"] for r in native) < 1.0
