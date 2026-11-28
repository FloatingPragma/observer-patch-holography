"""Reconstruct the all-phase mass comparison directly from archived interval witnesses."""

import json
from fractions import Fraction as F

import pytest

from openwave.interval_echo.row_19.selected_collar import (
    pointer_mass_equivalence as phase,
)


def load_record(name):
    return json.loads((phase.HERE / (name + ".json")).read_bytes())


def test_charge_selected_flight_measurements():
    data = load_record("charge_selected_flight")
    for window in data["windows"].values():
        assert window["selection_charge"][0] > 0 and window["held_charge"][0] > 0
        assert max(window["held_charge"][0], window["selection_charge"][0]) <= min(
            window["held_charge"][1], window["selection_charge"][1]
        )
        for key in ("area", "volume"):
            estimate = window["geometry_from_source"][key]
            radar = window["geometry_from_radar"][key]
            assert max(estimate[0], radar[0]) <= min(estimate[1], radar[1])
        for label, meter in window["material"].items():
            assert all(
                sum(x) == data["material_starts"] for x in meter["counts"].values()
            )
            forecast = window["forecast"][label]
            measured = meter["direct"]
            assert max(forecast[0], measured[0]) <= min(forecast[1], measured[1])


@pytest.mark.parametrize("name", ("pointer_mass_equivalence", "pointer_shape"))
def test_signed_quadrature_certificates(name):
    data = load_record(name)
    for case, groups in data["certificates"].items():
        times = list(data["ordinary"]) if case == "mixture" else [case]
        for group, tests in groups.items():
            assert len(tests) == len(times) * 4 ** len(times)
            for identifier, certificate in tests.items():
                pattern, base = identifier.split("/")
                keys = [f"{t}/{axis}" for t in times for axis in "RI"] + [
                    t + "/ordinary" for t in times if t != base
                ]

                def interval(key, name):
                    t, axis = key.split("/")
                    if axis == "ordinary":
                        a, b = data["ordinary"][t][name]
                        return -F(b), -F(a)
                    j = "RI".index(axis)
                    a, b = data["quadratures"][t][j][name]
                    return (
                        (F(a), F(b))
                        if pattern[2 * times.index(t) + j] == "+"
                        else (-F(b), -F(a))
                    )

                if certificate["status"] == "excluded":
                    values = dict.fromkeys(keys, F(0))
                    bound = F(0)
                    for name, side, n, d in certificate["coefficients"]:
                        assert name in data["groups"][group] and side in (
                            "lower",
                            "upper",
                        )
                        assert n > 0 and d > 0
                        endpoint = int(side == "upper")
                        weight = F(n, d) if endpoint == 0 else -F(n, d)
                        bound += weight * F(data["ordinary"][base][name][1 - endpoint])
                        for k in keys:
                            values[k] += weight * interval(k, name)[endpoint]
                    assert min(values.values()) >= 0 > bound
                elif certificate["status"] == "compatible_relaxation":
                    weights = {k: F(n, d) for k, n, d in certificate["allocation"]}
                    assert all(w >= 0 for w in weights.values()) and set(
                        weights
                    ) <= set(keys)
                    for name in data["groups"][group]:
                        lo, hi = (
                            sum(w * interval(k, name)[q] for k, w in weights.items())
                            for q in (0, 1)
                        )
                        assert lo <= F(data["ordinary"][base][name][1]) and hi >= F(
                            data["ordinary"][base][name][0]
                        )
                else:
                    assert certificate["status"] == "unresolved"
