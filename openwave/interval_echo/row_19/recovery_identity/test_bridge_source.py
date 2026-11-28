from copy import deepcopy
import json

from mpmath import iv
import numpy as np
import pytest

from openwave.interval_echo.row_19.recovery_identity import bridge_dynamics as native
from openwave.interval_echo.row_19.recovery_identity import bridge_radar as work
from openwave.interval_echo.row_19.recovery_identity import bridge_range as ranges
from openwave.interval_echo.row_19.recovery_identity import stored_source as store


def test_contact_projection_and_affine_population_control():
    # A swap-symmetric scalar contact has a distinct outside event and two
    # one-dimensional inside irreducible sectors, even with spectator heads.
    q, _ = np.linalg.qr(np.random.default_rng(72).normal(size=(9, 4)))
    grouped = q.reshape(3, 3, 4)
    effects = work.bridge_effects(grouped, (1, 2))
    for i, projector in enumerate(
        (
            np.outer([0, 1, 1], [0, 1, 1]) / 2,
            np.outer([0, 1, -1], [0, 1, -1]) / 2,
            np.diag([1, 0, 0]),
        )
    ):
        np.testing.assert_allclose(
            effects[i], q.T @ np.kron(projector, np.eye(3)) @ q, atol=5e-16
        )
    a = [iv.mpf(x) for x in (2, -3, 5, 7)]
    assert work.calibrated(a, [4 * x for x in a])["compatible"] == [True] * 3
    s = [iv.mpf(x) for x in (1, 0, 0)]
    b = [iv.mpf(x) for x in (0, 1, 0)]
    assert (
        ranges.balance(s, b, [3 * x + 4 * y for x, y in zip(s, b, strict=True)]) is None
    )
    c = ranges.balance(s, b, [iv.mpf(x) for x in (3, 4, 5)])
    assert c["delta_range"] == [5, 5]


def test_population_work_and_balanced_ranges_against_counts():
    d = json.loads(work.RECEIPT.read_text())
    n = json.loads(native.RECEIPT.read_text())
    r = json.loads(ranges.RECEIPT.read_text())
    assert work.decode(d) == d["results"]
    assert ranges.run() == r
    for arm, v in d["results"].items():
        _, q, _ = arm.split(",")
        assert n["contact"][q]["work_ratio"] in iv.mpf(
            v["fixed_material"]["source_paired_ratio"]
        )
        assert 0 in iv.mpf(v["blank_area"]) and 0 in iv.mpf(v["blank_work"])
        material = n["balance"][arm]
        areas = material["areas"]
        for val, bounds in zip(
            work.slopes(areas, d["steps"]), v["geometry_slopes"], strict=True
        ):
            assert val in iv.mpf(bounds)
        for i, (src, dst) in enumerate(work.radar.ROUTES):
            k = f"{src}>{dst} paired,formed_vacancy,formed_excess"
            cert = r["mixtures"][arm][k]
            if cert is not None:
                for weight, conf in zip(
                    material["weights"], cert["signed_weights"], strict=True
                ):
                    assert weight in iv.mpf(conf)
                assert material["range_differences"][i] in iv.mpf(cert["delta_range"])
    first = r["mixtures"]["0,1,-1"]["0>1 paired,formed_vacancy,formed_excess"]
    second = r["mixtures"]["0,1,1"]["0>1 paired,formed_vacancy,formed_excess"]
    assert first["delta_range"][0] > 0 > second["delta_range"][1]
    assert d["results"]["0,1,1"]["fixed_material"]["compatible"] == [False] * 3
    bad = deepcopy(d)
    bad["contact"]["1"][0][0] += 1
    with pytest.raises(ValueError):
        work.decode(bad)


def test_routed_source_and_local_trace_identity():
    d = json.loads(store.RECEIPT.read_text())
    n = json.loads(native.RECEIPT.read_text())
    assert store.decode(d) == d["measurements"]
    for key, rows in n["archive"].items():
        f = d["failure"][key]
        for arm, v in rows.items():
            report = d["measurements"][key]
            for i, label in enumerate(("ordinary", "held")):
                assert v["source"][i] in iv.mpf(report[label])
                assert v["geometry"][i] in iv.mpf(report["radar"][arm][label])
            expected = f * n["archive"]["all_fail"][arm]["source"][1]
            assert abs(v["source"][1] - expected) < 3e-12
    assert (
        sum(
            v["radar"]["1,1"]["compatible"] is False
            for k, v in d["measurements"].items()
            if "/" in k
        )
        == 5
    )
    # Implement an arbitrary trace-preserving reset on an entangled, discarded
    # three-level subsystem. Check every matrix element on the kept subsystem.
    rng = np.random.default_rng(2419)
    v = rng.normal(size=(6, 3)) + 1j * rng.normal(size=(6, 3))
    rho = v @ v.conj().T
    rho /= rho.trace()
    kets = np.eye(3)
    ops = [np.kron(np.eye(2), np.outer(kets[1], e)) for e in kets]
    changed = sum(k @ rho @ k.conj().T for k in ops)

    def keep(x):
        return np.einsum("iaja->ij", x.reshape(2, 3, 2, 3))

    np.testing.assert_allclose(keep(rho), keep(changed), atol=1e-15)
    bad = deepcopy(d)
    bad["panels"]["all_accept"]["source"][2][1][4][0] -= 1
    with pytest.raises(ValueError):
        store.decode(bad)


def test_continuous_calibration_including_tiny_failure_fractions():
    from openwave.interval_echo.row_19.recovery_identity import continuous_storage as c

    out = json.loads(c.RECEIPT.read_text())
    assert c.run() == out
    r = out["results"]["1,1"]
    assert r["all_fractions_excluded"] and r["lower"] > 2.5e-8
    native_data = json.loads(native.RECEIPT.read_text())["archive"]
    for key, f in store.fractions().items():
        if f == 0:
            assert native_data[key]["1,1"]["geometry"][1] == 0
            continue
        obs = native_data[key]["1,1"]
        so, sh = obs["source"]
        go, gh = obs["geometry"]
        residual = (gh * so - go * sh) / f
        box = min(int(c.BINS * f), c.BINS - 1)
        assert residual in iv.mpf(r["boxes"][box])
    # Direct analytic Heron difference, including the limit at f=0. The
    # baseline is a 3-4-5 triangle. Change its first side by f/2.
    from math import sqrt

    for f in (0, 1e-10, 0.125, 1):
        x = [iv.mpf(v) for v in (3, 4, 5)]
        y = list(x)
        y[0] += iv.mpf(f) / 2
        delta = c.factors([iv.mpf(0.5), iv.mpf(0), iv.mpf(0)])
        measured = c.difference_area(c.factors(y), c.factors(x), delta)
        # A^2=(4*a^2*b^2-(a^2+b^2-c^2)^2)/16. This expanded
        # difference has no cancellation at f=0.
        change_a2_over_f = 3 + f / 4
        predicted = change_a2_over_f * (64 - f * change_a2_over_f) / 16
        value = predicted / (6 + sqrt(36 + f * predicted))
        assert abs(float(measured.mid) - value) < 2e-14
