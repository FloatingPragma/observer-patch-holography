import gzip
import json

import numpy as np
from scipy.linalg import eigvalsh

from model.executable_model import cell
from . import measure


def test_clock_selects_the_native_energy_even_on_superpositions():
    car = cell()
    b, c = measure.incidence(car)
    p = (b + c / np.sqrt(5)) / 4
    h = measure.mounted_h(car, 1)
    first = eigvalsh(h)[1]
    np.testing.assert_allclose(h @ p, first * p, atol=1e-12)
    bright = np.eye(12)
    for f in measure.bright_filters(car, 1):
        bright = f @ bright
    np.testing.assert_allclose(bright, p, atol=2e-13)
    assert abs(np.trace(p) - 3) < 1e-12


def test_face_projection_and_changed_exchange():
    p = measure.loop_probs(1)
    np.testing.assert_allclose(
        p, [[0, 1, 0], [1, 0, 0], [0.5, 0.5, 0], [0, 0, 1]], atol=2e-13
    )
    control = measure.loop_probs(0.5)
    np.testing.assert_allclose(np.asarray(control)[:3, 0], 5 / 8, atol=2e-13)


def test_archive_replay_and_simultaneous_gateway_bounds():
    data = measure.replay()
    assert data == json.loads(
        gzip.decompress((measure.HERE / "counts.json.gz").read_bytes())
    )
    report = measure.decode(data)
    allowance = np.sqrt(np.log(2 * 7560 / 0.01) / (2 * measure.STARTS))
    assert report["native"]["max_pair_error"] < allowance < 0.003
    assert report["native"]["max_diagonal_error"] < allowance
    assert report["half_exchange"]["max_pair_error"] > 0.4
    assert sum(report["native"]["pair_classes"].values()) == 1260
    for bank in data["boards"]:
        assert len(bank["rooted_faces"]) == 60 * bank["owners"]
        for arm in bank["arms"].values():
            for key in ("loop_counts", "gateway_counts"):
                assert np.all(np.asarray(arm[key]).sum(axis=-1) == measure.STARTS)
