"""Independent native full-cell, detector and held-frequency checks."""

import gzip
import json

import numpy as np
import pytest
from scipy.linalg import eigh, expm
from model.executable_model.carrier import cell
from . import core, measure


def raw_hamiltonian(board, seam, k):
    """Assemble each actual seam and port reference before subtracting the free load."""
    h0, h1 = core.LAW.local_term_parts()
    local = h0 + k * h1
    n = board.n
    h = np.diag(np.tile(k * (6 - board.degree) - 6 * k, 11)).astype(float)
    lookup = {None: 0, **{label: 1 + i for i, label in enumerate(core.LAW.LABELS)}}
    for a, b in board.seams:
        if (a, b) != seam:
            for i in range(11):
                ports = [i * n + a, i * n + b]
                h[np.ix_(ports, ports)] += local[np.ix_([3, 1], [3, 1])]
            continue
        for before in lookup:
            for after in lookup:
                m0, m1 = lookup[before], lookup[after]
                outside = local[
                    core.local_index(after, (0, 0)), core.local_index(before, (0, 0))
                ]
                for j in range(n):
                    if j not in seam:
                        h[m1 * n + j, m0 * n + j] += outside
                rows = [core.local_index(after, w) for w in ((1, 0), (0, 1))]
                cols = [core.local_index(before, w) for w in ((1, 0), (0, 1))]
                h[np.ix_([m1 * n + a, m1 * n + b], [m0 * n + a, m0 * n + b])] += local[
                    np.ix_(rows, cols)
                ]
    return h


@pytest.mark.parametrize("k", [0, 8])
def test_all_132_levels_and_native_phase_flights(k):
    board = cell()
    a, b = seam = (0, 1)
    n = board.n
    full = raw_hamiltonian(board, seam, k)
    f = eigh(full, eigvals_only=True)
    blank = core.standing_modes(board, seam, k, "blank", False)[0]
    cut = core.standing_modes(board, seam, k, "vacuum_write", False)[0]
    classified = [blank, cut, cut]
    for odd in (False, True):
        e, sources, read = core.standing_modes(board, seam, k, "occupied", odd)
        mu = k + int(odd)
        classified.extend([cut + 4 * k + int(odd)] * 2 + [e + mu])
        local = core.occupied_source(odd).reshape(11, 9)
        source = np.zeros((11, n), dtype=complex)
        source[:, a] = local[:, 3]
        source[:, b] = local[:, 1]
        for i, phase in enumerate((1, -1)):
            launch = source.copy()
            launch[:, a] *= phase
            for t in (0, 1.327):
                actual = (expm(-1j * t * full) @ launch.ravel()).reshape(11, n)
                expected = read(sources[i] * np.exp(-1j * t * e))
                np.testing.assert_allclose(
                    (abs(actual) ** 2).sum(axis=0), expected, atol=1e-12
                )
    np.testing.assert_allclose(f, np.sort(np.concatenate(classified)), atol=1e-12)
    if k == 8:
        bound = f[f < -1e-9]
        assert len(bound) == 6 and bound.min() >= -8 and bound.max() <= -5.5
    else:
        assert min(f) > -1e-12


def receipt():
    return json.loads(gzip.decompress(measure.OUTPUT.read_bytes()))


def test_reference_readings_and_public_routes():
    archive = receipt()
    values = [0, 1, 8, 9, 48, 0]
    overlaps = core.separated_standards(np.pi / 128)
    for (key, reading), target in zip(
        archive["separated_energy_standards"].items(), values, strict=True
    ):
        assert abs(reading["reading"] - target) < reading["radius"] < 0.123
        np.testing.assert_allclose(
            overlaps[key], np.exp(1j * target * np.pi / 128), atol=1e-12
        )
        assert all(
            sum(pair) == measure.SHOTS for pair in reading["real_and_imag_counts"]
        )
    for level in range(3):
        board, seam, chart = core.chart(level)
        assert board.degree[chart["phase_route"][-1]] == 5
        assert chart["reference_parks"][0][-1] != chart["reference_parks"][1][-1]
        for path in [chart["phase_route"], *chart["reference_parks"]]:
            assert all(
                q in board.neighbours[p] and {p, q} != set(seam)
                for p, q in zip(path, path[1:])
            )
        expected = archive["boards"][str(board.cells)]["terminal"]
        np.testing.assert_equal(
            json.loads(json.dumps(core.terminal_reader(board, seam))), expected
        )
        assert expected["owners"] == board.cells


def test_clock_by_direct_fourier_and_all_registered_outcomes():
    e = np.array([-13.829, -5.89, 0, 6.01])
    v = np.sqrt([[0.3, 0.2, 0.1, 0.4], [0.5, 0.18, 0.32, 0]])
    direct = (
        abs(
            np.fft.ifft(np.exp(-1j * np.outer(e, measure.STEP * measure.INDEX)), axis=1)
        )
        ** 2
    )
    np.testing.assert_allclose(measure.timing(e, v), v**2 @ direct, atol=2e-12)
    result = receipt()
    for size in ("20", "80"):
        panels = result["boards"][size]["panels"]
        for panel in panels:
            for bank in panel["banks"]:
                assert sum(bank["clock"]["counts"]) == measure.SHOTS
                for flight in bank["flight"]:
                    assert (
                        len(flight["counts"]) == 12 * int(size)
                        and sum(flight["counts"]) == measure.SHOTS
                    )
                    assert (
                        abs(
                            flight["endpoint_counts"] / measure.SHOTS
                            - flight["endpoint_probability"]
                        )
                        < 0.002
                    )
            for window in panel["held_clock_windows"]:
                counts = np.array(panel["banks"][1]["clock"]["counts"])
                assert (
                    counts[np.abs(measure.ENERGY - window["target"]) <= 0.1].sum()
                    == window["counts"]
                )
        for native, control in zip(panels[:2], panels[2:4], strict=True):
            levels = native["below_minus_1e_9"]
            assert len(levels) == 4 and levels[0][1] > 0.994 and levels[0][2] < 0.001
            assert levels[2][2] > 0.9 and levels[1][1] + levels[3][1] < 1e-20
            assert all(w["counts"] < 60 for w in native["held_clock_windows"])
            left, right = [
                r["flight"][-1]["endpoint_probability"] for r in native["banks"]
            ]
            assert left - right > 0.027
            assert control["banks"][0]["peak_energy"] > 0
        assert panels[2]["below_minus_1e_9"] == []
        left, right = [
            r["flight"][-1]["endpoint_probability"] for r in panels[2]["banks"]
        ]
        assert abs(left - right) < 1e-12
        assert -0.07 < panels[3]["below_minus_1e_9"][0][0] < -0.067
    with pytest.raises(ValueError):
        measure.verify({"counts": [2, 3]}, {"counts": [3, 2]})
