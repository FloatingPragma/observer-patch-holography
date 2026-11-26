import json

import numpy as np
from scipy.linalg import expm

from . import measure


def test_finite_time_kick_reader_against_direct_schrodinger_paths():
    # An energy-coherent reference checks the phase of the aperture as well as its width.
    e = np.array([-0.26, 1.3, 3.8])
    start = np.array([1, 1j, -1]) / np.sqrt(3)
    flip = np.eye(3)[[2, 1, 0]]
    target, step, length = 1.2, 0.09, 47
    h = np.diag(e)
    direct = (
        sum(
            np.exp(-1j * target * step * s)
            * expm(-1j * h * step * (length - 1 - s))
            @ flip
            @ expm(-1j * h * step * s)
            @ start
            for s in range(length)
        )
        / length
    )
    offsets = e[:, None] - e[None, :] - target
    actual = (flip * measure.aperture(step * offsets, length)) @ start
    np.testing.assert_allclose(
        direct, np.exp(-1j * e * step * (length - 1)) * actual, atol=1e-13
    )


def test_enclosures_and_held_forecasts():
    old = measure.read(measure.HERE.parent / "records.json")
    root = next(
        p for p in measure.HERE.parents if (p / "model/candidate_laws").is_dir()
    )
    path = (
        root
        / "model/candidate_laws/interval_echo/openwave/row_03_clock/polar_island/enclosures.json"
    )
    checks = {(p["n"], p["kappa"]): p for p in json.loads(path.read_text())}
    data = measure.read(measure.HERE / "records.json")
    assert len(data["panels"]) == 4
    for p in data["panels"]:
        n, k = p["n"], p["kappa"]
        previous = next(
            panel
            for panel in old["panels"]
            if panel["cells"] == 80
            and panel["dimer_count"] == n
            and panel["kappa"] == k
        )
        cert = checks[n, k]
        assert cert["ordered_eigenvalue_error"] < 8e-9
        assert previous["ground_energy"] + cert["ground_to_rest"] - 3 * 8e-9 > 17
        assert 3 * k * n + n / 2 < 512 * np.pi
        assert p["refused"] + np.sum(p["pole_and_complement"]) == data["starts"]
        assert all(
            sum(counts) == data["starts"]
            for r in p["local"]
            for arm in ("live", "exchange_off")
            for counts in r[arm]
        )
        assert all(r["owner"] != p["origin_owner"] for r in p["local"])
        report = measure.decode(p)
        assert report["min_scaled_weight_lower"] > 0.8
        assert report["max_pole_error"] < 0.03
        assert report["max_local_response"] > 0.08
        assert report["max_flat_error"] > 0.08
        assert report["max_exchange_off_response"] < 0.001
        if k == 8:
            assert report["long_time_local_power_lower"] > {10: 0.098, 12: 0.0817}[n]


def test_native_four_panels_replay():
    assert measure.rounded(measure.replay()) == measure.read(
        measure.HERE / "records.json"
    )
