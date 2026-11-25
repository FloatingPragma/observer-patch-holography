import json

import numpy as np

from . import measure as m


def test_exact_coefficients_and_clock_coordinates():
    report = json.loads((m.HERE / "results.json").read_text())
    weight = m.weights()
    assert weight == report["log_range_weights"]
    assert min(lo for lo, _ in weight["directed"]) > 14.71
    assert len(weight["edges"]) == 6 and weight["total"][1] < 237.55
    ref = m.preparations()["reference"]
    radar = m.station.radar_menus()[6]
    times = m.station.pair(radar, ref)
    record = json.loads((m.ALGEBRA / "results.json").read_text())
    area, volume = m.station.range_gradient(times)
    prices = times * (record["ordinary_eta"] * area - record["price"] * volume)
    upper = dict(zip(zip(*m.station.ORDER), prices, strict=True))
    edge_prices = np.array([upper[i, j] + upper[j, i] for i, j in m.EDGES])
    np.testing.assert_allclose(prices, m.log_length_map().T @ edge_prices, atol=3e-12)
    for state in m.preparations().values():
        delta = state - ref
        z = np.log(m.station.pair(radar, state) / times)
        work = -m.station.pair(m.station.log_blocks(ref), delta)
        divergence = m.station.pair(
            m.station.log_blocks(state) - m.station.log_blocks(ref), state
        )
        np.testing.assert_allclose(work, -prices @ np.expm1(z), atol=2e-10)
        rhs = (
            divergence
            + prices @ (np.expm1(z) - z)
            + edge_prices @ (m.log_length_map() @ z)
        )
        np.testing.assert_allclose(
            m.station.entropy(ref) - m.station.entropy(state), rhs, atol=2e-10
        )


def test_independent_clock_and_tomography_banks():
    report = json.loads((m.HERE / "results.json").read_text())
    readings = report["acquisitions"]
    ref = readings["reference"]
    for name, row in readings.items():
        assert all(sum(v) == m.TICK_BUDGET for v in row["clocks"]["counts"])
        assert all(
            np.sum(v) == m.TICK_BUDGET for v in row["material"]["counts"].values()
        )
        if name != "reference":
            assert row["material"]["entropy"][0] > ref["material"]["entropy"][1]
            assert row["clocks"]["area"][1] < ref["clocks"]["area"][0]
            assert row["clocks"]["volume"][1] < ref["clocks"]["volume"][0]
    data = m.preparations()
    np.testing.assert_allclose(
        np.trace(data["reference"][1]), np.trace(data["independent_right"][1])
    )
    assert readings["ordinary_write"]["material"]["odd_energy"][0] > 0.49
