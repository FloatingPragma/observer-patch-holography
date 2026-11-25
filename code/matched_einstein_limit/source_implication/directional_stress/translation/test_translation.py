import json

import numpy as np

from ...action_response.four_ticks import parse_interval
from . import bridge_limit, counts, measure, smooth, transmitting


def test_two_regional_balances_and_native_flight():
    data = measure.run()
    certified = json.loads((measure.HERE / "certificate.json").read_text())
    assert len(data["cuts"]["domain"]) == 14 and len(data["chart"]["cells"]) == 12
    assert data["finite_survey"][320]["maximum_witnessed_connected_cells"] == 85
    for sign, row in data["signs"].items():
        local = data["local_rates_and_edges"][sign]
        for k in row["free_rate"]:
            assert (
                abs(row["free_rate"][k] - row["independent_free_shift_difference"][k])
                < 1e-8
            )
            for name in ("boundary", "marked_force", "free_rate"):
                assert row[name][k] in parse_interval(certified[sign][name][k])
            for u in data["chart"]["cells"]:
                cuts = sum(
                    v[k] for tag, v in local.items() if tag.startswith(f"edge_{u}_")
                )
                np.testing.assert_allclose(
                    cuts + local[f"body_{u}"][k], local[f"rate_{u}"][k], atol=1e-12
                )


def test_meter_resolves_body_and_boundary():
    data = counts.run()
    for row in data["signs"].values():
        for arm in ("ordinary", "pairing"):

            def read(k):
                return row[k]["slopes"][arm]

            d = (
                read("free_rate")["estimate"]
                - read("boundary")["estimate"]
                - read("marked_force")["estimate"]
            )
            assert abs(d) < sum(read(k)["halfwidth"] for k in row)


def test_trace_obstruction_and_smooth_sector_bounds():
    data = bridge_limit.run()
    for s in data["irreducible_marked_force"].values():
        assert s["integer_traces"] == [0] * 7 + [64]
    for rows in data["charges"].values():
        for v in rows.values():
            assert abs(v["energy_load_covariance"]) < v["energy_covariance_bound"]
            assert abs(v["actual_path_error"]) < v["path_error_bound"]


def test_uniform_cell_transfer_and_native_leakage():
    values = smooth.run()["charge_results"]
    for sign, rows in values.items():
        assert rows["central_position_blocks_proved_one_head_only"] > 0
        assert rows["arbitrary_common_orbital_pairing_response_exact"] == 0
        assert abs(rows["central_flux_at_pi_over_two"]["pairing"]) < 1e-12
        assert abs(rows["central_flux_at_pi_over_two"]["ordinary"]) > 1e-5
        assert (
            0.06
            < rows["original_mixture_uniform_subspace_weight"]["equilibrium"]
            < 0.09
        )
        assert (
            rows["lower_bound_on_squared_leakage_speed"]
            == {"-1": "64/3", "1": "75/4"}[sign]
        )
        shift = rows["readings"]["shift"][2]
        local = rows["readings"]["routed_energy_over_40"][2]
        assert (
            shift["estimate"] - shift["halfwidth"]
            > local["estimate"] + local["halfwidth"]
        )


def test_invariant_transmitting_sector_keeps_the_force_obstruction():
    result = transmitting.run()
    assert result["registered_transmitting_weight_exact"] == "4/5"
    assert result["shift_witness_gateway_ports"] == [1, 229]
    assert all(v["transmitting_dimensions"] == 2912 for v in result["sectors"])
    assert result["native_trace_H7_force_restricted_to_transmitting_sector"] == "1/8"
