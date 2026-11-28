"""Verify the two probes, ordinary distributions, and cancellation on recorded endpoints."""

from fractions import Fraction as F
import json

from openwave.interval_echo.row_19.selected_collar import density_material as experiment


def read(name):
    return json.loads((experiment.HERE / (name + ".json")).read_bytes())


def test_direct_tagged_and_untagged_forecasts():
    for name in (
        "density_material",
        "density_material_active",
        "density_material_impulse",
    ):
        for frame in read(name)["windows"].values():
            assert frame["ordinary_energy"][0] > 0
            assert frame["sparse_energy"][0] <= 0 <= frame["sparse_energy"][1]
            for power, meters in frame["exposures"].items():
                for meter in meters.values():
                    for probe in meter["probes"].values():
                        for a, b, miss in probe["counts"].values():
                            assert min(a, b, miss) >= 0 and a + b + miss == 2 ** int(
                                power
                            )
                        for arm, result in probe["contrasts"].items():
                            forecast = probe["source_forecasts"][arm]
                            assert max(forecast[0], result[0]) <= min(
                                forecast[1], result[1]
                            )
                            if arm == "erased":
                                assert result[0] <= 0 <= result[1]
    earlier = read("density_material_active")["windows"]["2"]["exposures"]["72"][
        "shift"
    ]["probes"]
    assert (
        earlier["-1"]["ordinary_energy_equivalent"][0]
        <= 0
        <= earlier["-1"]["ordinary_energy_equivalent"][1]
    )
    held = read("density_material_impulse")["windows"]["2"]["exposures"]["76"]["shift"][
        "probes"
    ]
    assert (
        held["1"]["ordinary_energy_equivalent"][1]
        < 0
        < held["-1"]["ordinary_energy_equivalent"][0]
    )


def test_every_rational_service_witness():
    for name in ("ordinary_service_distribution", "charge_even_material"):
        for window in read(name)["windows"].values():
            assert not window.get("ordinary_controls_disagree", [])
            columns = window["ordinary_columns"]
            assert len(columns) == 12
            for source, groups in window["comparisons"].items():
                target = window["target"][source]
                for group, certificate in groups.items():
                    names = window["groups"][group]
                    if certificate["status"] == "compatible_relaxation":
                        allocation = {
                            key: F(n, d) for key, n, d in certificate["allocation"]
                        }
                        assert min(allocation.values(), default=0) >= 0
                        for label in names:
                            low, high = (
                                sum(
                                    w * F(columns[k][label][endpoint])
                                    for k, w in allocation.items()
                                )
                                for endpoint in (0, 1)
                            )
                            assert low <= F(target[label][1]) and high >= F(
                                target[label][0]
                            )
                    else:
                        assert certificate["status"] == "excluded"
                        response = dict.fromkeys(columns, F(0))
                        rhs = F(0)
                        for label, side, n, d in certificate["coefficients"]:
                            assert (
                                label in names
                                and side in ("lower", "upper")
                                and n > 0
                                and d > 0
                            )
                            w = F(n, d) if side == "lower" else -F(n, d)
                            index = int(side == "upper")
                            rhs += w * F(target[label][1 - index])
                            for key, z in columns.items():
                                response[key] += w * F(z[label][index])
                        assert min(response.values()) >= 0 > rhs


def test_charge_even_proof_from_original_readouts():
    paired = read("ordinary_service_distribution")["windows"]["2"]
    even = read("charge_even_material")["windows"]["2"]
    for source in ("direct", "source"):
        for meter in experiment.METERS:
            endpoints = list(
                zip(
                    paired["target"][source][meter + "/-1"],
                    paired["target"][source][meter + "/1"],
                    strict=True,
                )
            )
            low, high = [sum(map(F, z)) / 2 for z in endpoints]
            assert (
                F(even["target"][source][meter][0])
                <= low
                <= high
                <= F(even["target"][source][meter][1])
            )
    assert min(z["shift"][0] for z in even["ordinary_columns"].values()) > 5.6493e-10
    assert even["target"]["direct"]["shift"][1] < -6.8106e-12
    cost = even["uniform_energy_equivalents"]["direct"]
    assert cost["marked_force"][0] > 0 > cost["free_rate"][1]
    assert cost["shift"][1] < 0
