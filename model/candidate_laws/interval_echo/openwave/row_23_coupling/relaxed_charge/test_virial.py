"""Check the thermal transform against microscopic heat traces and opposite cuts."""

import gzip
import json
import numpy as np
from scipy.linalg import expm
from model.executable_model import federation
from . import confirmation, readback, virial
from .test_static import archive, direct_words


def load(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def test_independent_microscopic_heat_traces_and_cut():
    r = load(confirmation.OUTPUT)["rows"]["20,3"]
    original = load(virial.OUTPUT)["panels"]["20,8,3"]
    b = federation(0)
    home = archive()["boards"]["20,8"]["source"]
    y = r["probe"]
    for beta in virial.BETAS:
        score = {}
        for k, setting in r["signs"].items():
            s, t = map(int, k.split(","))
            free, h = direct_words(b, 8, {home: s + 1, y: t + 1})
            trace = np.trace(expm(-beta * h))
            observed = setting["readings"][str(beta)]
            assert abs(observed["trace"] - trace) < observed["error"]
            score[k] = trace
        weight = -readback.signed(score) / beta
        held = r["readings"][str(beta)]
        assert weight < 0 and abs(weight - held["value"]) < held["error"]
        assert (
            abs(original["readings"][str(beta)]["cut"])
            <= original["readings"][str(beta)]["cut_error"]
        )


def test_all_calibrated_returns_and_foil_decisions():
    regular = load(virial.OUTPUT)
    confirmation_rows = load(confirmation.OUTPUT)["rows"]
    for row in regular["panels"].values():
        for v in row["readings"].values():
            assert (
                abs(v["intact"] - v["forecast"])
                < v["intact_error"] + v["forecast_error"]
            )
    agreed = []
    for key, row in confirmation_rows.items():
        for beta, v in row["readings"].items():
            assert v["agrees"]["1.0"]
            agreed.append((key, beta, v["agrees"]["0.9"], v["agrees"]["1.1"]))
    assert [(key, beta) for key, beta, lower, upper in agreed if lower or upper] == [
        ("80,3", "0.0625")
    ]


def test_transfer_sign_and_exact_thermal_factorization_in_a_small_host():
    adj = np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]]) * 0.5
    for hop in (True, False):
        h = -adj.copy()
        if not hop:
            h[1, 2] = h[2, 1] = 0
        corners = {}
        for s in (-1, 1):
            for t in (-1, 1):
                hd = h + np.diag([2.0 * s, 0.0, 0.0, 2.0 * t])
                corners[f"{s},{t}"] = np.trace(expm(-0.8 * hd))
        assert (
            readback.signed(corners) > 1e-6
            if hop
            else abs(readback.signed(corners)) < 1e-14
        )
