"""Regression checks against independently computed native observables."""

import json
import numpy as np
from scipy.stats import nbinom

from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.row_30.experiment import (
    Guide,
    HERE,
    select_route,
    spatial_inverse,
    spatial_reader,
    temporal_inverse,
    TIMES,
)
from model.candidate_laws.interval_echo.openwave.row_30.patch_time import decoding_curve


def test_recovered_winding_and_local_beat():
    car = federation(0)
    route, _ = select_route(car)
    guide = Guide(car, route, 0.2)
    k = -2 * np.pi / guide.size
    state = guide.uniform(k)
    assert abs(spatial_inverse(spatial_reader(guide, state)) - k) < 1e-11
    clicks = []
    for _ in TIMES:
        clicks.append(abs(state[route[1]]) ** 2)
        guide.step(state)
    assert (
        abs(temporal_inverse(np.array(clicks)) - np.arccos(np.cos(0.2) * np.cos(k)))
        < 1e-7
    )


def test_local_tick_frequency_transfer():
    p, omega = 0.5, 0.32
    for ticks, predicted in enumerate(decoding_curve(omega, p)[1:], start=1):
        waits = np.arange(ticks, 256)
        direct = nbinom.pmf(waits - ticks, ticks, p) @ (np.sin(omega * waits) ** 2)
        assert abs(direct - predicted) < 1e-12


def test_held_accuracy_and_microscopic_challenge():
    rows = json.loads((HERE / "results.json").read_text())["trials"]
    for board in rows:
        for band in board["bands"]:
            if band["detuning"] < 0.8:
                assert band["worst_kg_exact"] < 0.02
                assert band["worst_kg_sampled"] < 0.06
                assert band["worst_parabolic_exact"] > 0.4
            else:
                assert band["worst_kg_exact"] > 0.08
            for packet in band.get("packets", []):
                assert packet["absolute_error"] < 0.08
