"""Cross-check the source-only, stationary and spatial inferences against raw records."""

import json
from itertools import product

import numpy as np

from openwave.interval_echo.row_19.selected_collar import persistent_density as dynamics


def record(name):
    return json.loads((dynamics.HERE / (name + ".json")).read_bytes())


def overlaps(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def test_source_reader_and_independent_instruments():
    data = record("position_source")
    for power, setting in data["exposures"].items():
        for counts in setting["records"].values():
            assert len(counts) == 132
            assert all(
                len(row) == 7 and min(row) >= 0 and sum(row) == 2 ** int(power)
                for row in counts
            )
    high = data["exposures"]["70"]["decoded"]
    held = high["held"]
    for field, independent in (
        ("connected_charge_change", "independent_held_head_charge"),
        ("entropy_change", "independent_held_head_entropy"),
        ("inferred_area_change", "independent_held_radar_area"),
    ):
        assert overlaps(held[field], data[independent])
        assert overlaps(high["sham"][field], [0, 0])
    assert (
        high["reverse"]["connected_charge_change"][1]
        < 0
        < held["connected_charge_change"][0]
    )
    assert (
        0
        < held["minimum_contraction_from_entropy"]
        <= -data["independent_held_radar_area"][1]
    )
    assert overlaps(
        data["exposures"]["66"]["decoded"]["held"]["connected_charge_change"], [0, 0]
    )


def test_link_population_witness_and_volume_bracket():
    data = record("link_source_alternative")
    low, high = data["triangles"]
    assert low["volume_minus_held"][1] < 0 < high["volume_minus_held"][0]
    for triangle in (low, high):
        weights = triangle["barycentric_weights"]
        assert min(p[0] for p in weights) > 0
        assert sum(p[0] for p in weights) <= 1 <= sum(p[1] for p in weights)
        assert max(data["modular_change"][k][1] for k in triangle["sources"]) < 0
    assert data["modular_change"]["held"][0] > 0
    assert (
        data["held_area_minus_reference"][1]
        < 0
        < data["all_alternatives_area_minus_reference"][0]
    )


def test_channel_action_and_measured_phase_dependence():
    # Trace an explicitly attached two-dimensional pointer. The second
    # computation uses the Pauli-basis superoperator of the process reader.
    rng = np.random.default_rng(119)

    def noise():
        return rng.normal(size=(2, 4, 4)) + 1j * rng.normal(size=(2, 4, 4))

    click, miss, prep = noise(), noise(), noise()[0]
    source = prep @ prep.conj().T
    source /= source.trace()
    axes = dynamics.pos.AXES
    coords = np.einsum("kij,ji->k", axes, source).real

    def read(g):
        parts = [miss + g * click, np.sqrt(1 - abs(g) ** 2) * click]
        ops = np.concatenate(parts)
        direct = (ops @ source @ ops.conj().swapaxes(-1, -2)).sum(axis=0)
        decoded = np.tensordot(dynamics.superoperator(ops) @ coords, axes, axes=1)
        np.testing.assert_allclose(direct, decoded, atol=2e-13)
        return direct

    u, v = (read(1) - read(-1)) / 2, (read(1j) - read(-1j)) / 2
    for gamma in (0.3 + 0.4j, -0.4 + 0.3j):
        np.testing.assert_allclose(
            read(gamma) - read(0), gamma.real * u + gamma.imag * v, atol=2e-13
        )
    windows = record("source_interference")["windows"]
    erased = record("observation_density")["windows"]
    for duration, window in windows.items():
        predictions = window["predictions"]
        for arm, panel in predictions.items():
            source_area, radar = panel["source_area"], panel["radar_area"]
            assert overlaps(source_area, radar)
            if arm == "orthogonal":
                assert overlaps(source_area, [0, 0]) and overlaps(radar, [0, 0])
            else:
                assert source_area[0] * source_area[1] > 0 and radar[0] * radar[1] > 0
        assert (predictions["omitted"]["radar_area"][1] < 0) == (duration == "2")
        assert all(
            overlaps(z, [0, 0]) for z in erased[duration]["erased_minus_dense"].values()
        )


def test_material_price_all_coordinates_against_integer_flight():
    import sympy as sp
    from matched_einstein_limit.source_implication.directional_stress.translation import (
        certificate,
    )
    from openwave.interval_echo.row_19.selected_collar import (
        density_material as material,
    )

    bridge = certificate.bridge
    native = certificate.measure
    atlas = native.Atlas(1)
    route = bridge.native.iso.pairing.formation.region(atlas.board)
    launches = bridge.beam.indices(sum(bridge.beam.screens(atlas.board, route), []))
    weights = bridge.native.prior.weights()[7]
    numerical, _, _ = material.beam.flight(material.METERS, include_edges=False)
    for sign, odd in product((-1, 1), (False, True)):
        h, ops = native.op_table(atlas.board, route, native.chosen(atlas), sign, odd)
        indices = np.unique(
            np.concatenate([ops[k][0].nonzero()[0] for k in material.METERS])
        )
        panels = certificate.four_tick_columns(h, indices, launches, weights)
        for name in material.METERS:
            op, unit = ops[name]
            assert abs(op).sum(axis=1).max() / abs(unit) < 200
            fixed = certificate.compress(panels, op[indices, :][:, indices], unit)
            price = numerical[sign, name][int(odd)]
            # The rational L1 sum bounds the operator norm, even for complex
            # entries. The fixed-point propagation allowance is per 4x4 block.
            distance = sum(
                abs(sp.re(z) - sp.Rational(float(v.real)))
                + abs(sp.im(z) - sp.Rational(float(v.imag)))
                for z, v in zip(fixed, price.ravel(), strict=True)
            )
            assert distance + sp.Rational("2e-15") < sp.Rational("1e-13")


def test_local_price_exclusions_from_determinant_signs():
    data = record("regional_density")
    rows = [r for window in data["windows"].values() for r in window["balls"]]
    assert len(rows) == 108
    for window in data["windows"].values():
        assert all(len(z) == 7 for z in window["source_by_receiver"].values())
    assert (
        sum(not overlaps(r["determinants"]["connected"], [0, 0]) for r in rows) == 107
    )
    for model in ("marginal", "energy"):
        for filling in {tuple(r["triangle"]) for r in rows}:
            lower, upper = 0, np.inf
            for row in (r for r in rows if tuple(r["triangle"]) == filling):
                a, b = (row["determinants"][k] for k in (model, "connected"))
                # Necessary inequalities on alpha >= 0, reconstructed without
                # using the decoder's saved feasible interval.
                for intercept, slope in ((a[0], b[0]), (-a[1], -b[1])):
                    if slope > 0:
                        upper = min(upper, -intercept / slope)
                    elif slope < 0:
                        lower = max(lower, -intercept / slope)
                    elif intercept > 0:
                        upper = -np.inf
            assert lower > upper


def test_relaxation_and_finite_gain_at_each_stationary_source():
    schedules = record("persistent_density")["schedules"]
    test = record("stationary_source_test")["works_at"]
    assert set(test) == set(schedules) == {"uniform", "4", "8", "12"}
    for name, panel in schedules.items():
        for arm, sector in product(range(2), repeat=2):
            power = panel["process_calibration"][arm][sector]
            periods, remainder = divmod(panel["rounds"], power["power"])
            assert remainder == 0 and 2 * power["power_norm_bound"] ** periods <= 1e-11
        q = panel["stationary_connected"]
        assert q["ordinary_reference"][1] < 0 < q["dense_reference"][0]
        for geometry in ("pulse_geometry", "terminal_geometry"):
            assert panel[geometry]["area"][0] > 0 > panel[geometry]["adjusted_nats"][1]
        for arm in test[name].values():
            slope = arm.get("followup", {}).get(
                "entropy_slope", arm["entropy_slope_in_frozen_direction"]
            )
            assert (
                min(v * arm["sign"] for v in slope) * arm["step"]
                >= arm["minimum_finite_entropy_gain"]
                > 0
            )
    assert overlaps(
        test["uniform"]["sparse"]["entropy_slope_in_frozen_direction"], [0, 0]
    )
