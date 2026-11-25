import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply

from . import aperture, current, regional
from .. import live_radar, live_source
from ..source import measure as source

HERE = Path(__file__).parent


def test_region_restriction_and_modular_slopes_from_spent_mobile():
    board = source.formation.federation(1)
    rho, product = source.form(board, source.formation.region(board), -1)
    full = live_source.head_channel()
    for saved in json.loads((HERE / "regional.json").read_text())["regions"]:
        certified = regional.certify(saved["radius"])
        ports = certified["ports"]
        outside = np.setdiff1d(range(board.n), ports)
        w = full[ports].reshape(-1, 8)

        def blocks(head):
            s = np.kron(np.diag(source.RAILS["reference"]), head)
            return w @ s @ w.conj().T, live_source.apply_channel(
                full[outside], head, source.RAILS["reference"]
            )

        r = blocks(rho)
        u = blocks(product)
        # QR on the injected mobile columns identifies their support without
        # dividing by a small reference density eigenvalue. This replay starts
        # with the full local density and uses a different compression from J.
        occupied, _ = np.linalg.qr(w, mode="reduced")
        total = 0
        for a, b, basis in zip(r, u, (occupied, np.eye(8)), strict=True):
            a0, b0 = basis.conj().T @ a @ basis, basis.conj().T @ b @ basis
            np.testing.assert_allclose(
                basis @ a0 @ basis.conj().T, a, atol=1e-13, rtol=0
            )
            np.testing.assert_allclose(
                basis @ b0 @ basis.conj().T, b, atol=1e-13, rtol=0
            )
            val, vec = np.linalg.eigh(a0)
            assert min(val) > 1e-8
            total -= np.log(val) @ (vec.conj().T @ (b0 - a0) @ vec).diagonal().real
        np.testing.assert_allclose(
            [np.trace(x).real for x in r],
            saved["occupied_and_vacuum_trace"],
            atol=1e-13,
        )
        assert abs(sum(np.trace(x).real for x in r) - 1) < 1e-11
        lo, hi = certified["regional_entropy_slopes"]["decorrelated"]
        assert lo < total < hi
        assert certified["fixed_metric_pairing_slope"][0] > 0.19
        assert saved["fixed_metric_pairing_slope"][0] > 0.19


def test_native_energy_current_matches_energy_lost_between_appointments():
    board = source.formation.federation(1)
    route = source.formation.region(board)
    h = source.system(board, route, -1)[0] / 40
    rho, product = source.form(board, route, -1)
    delta = np.kron(np.diag(source.RAILS["reference"]), product - rho)
    v = (
        live_source.head_channel()
        .reshape(-1, 2, 4, 8)
        .transpose(1, 0, 2, 3)
        .reshape(-1, 8)
    )
    dt = 2e-5
    record = json.loads((HERE / "current.json").read_text())
    for row in record["regions"]:
        ports = row["ports"]
        bounds = current.slopes(ports)
        addresses = live_radar.station_rows(board, ports).ravel()
        reads = []
        for t in (-dt, dt):
            w = expm_multiply(-1j * t * h, v, traceA=-1j * t * h.diagonal().sum())
            x = w[addresses]
            e = x.conj().T @ (h @ w)[addresses]
            e = (e + e.conj().T) / 2
            reads.append([np.trace(op @ delta).real for op in (x.conj().T @ x, e)])
        observed = (np.array(reads[0]) - reads[1]) / (2 * dt)
        for i, k in enumerate(("number", "energy")):
            lo, hi = bounds[k]["decorrelated"]
            assert abs(observed[i] - (hi + lo) / 2) < 2e-10
            assert (
                lo < row["derivatives"][k]["decorrelated"][1]
                and row["derivatives"][k]["decorrelated"][0] < hi
            )
    data = json.loads((HERE / "regional.json").read_text())
    assert len(data["fixed_census_transfers"]) == 2
    source_record = json.loads((HERE / "aperture.json").read_text())
    for charge, rows in source_record["panels"].items():
        result = regional.match_all_three(rows[1], data["regions"][0])
        assert result["finite_entropy_difference"][0] > 1.2e-15
        assert result["regional_entropy_slope"][0] > 0.21
        assert data["fixed_census_and_both_currents"][charge] == result
        assert (
            min(
                x["held_residual"][0]
                for x in data["fixed_census_transfers"][charge].values()
            )
            > 0.19
        )
    assert data["regions"][0]["ports"] == aperture.enroll(board, 1)[0]


def test_front_identical_native_lotteries_carry_different_energy_flux_on_both_cuts():
    rows = json.loads((HERE / "current.json").read_text())["regions"]
    assert (
        current.universal_current_gap(rows)["uniform_over_energy_origin_lower"]
        == 3.4e-11
    )
    assert rows[0]["matched_front_lottery_difference"]["energy"][0] > 1.13e-9
    assert rows[1]["matched_front_lottery_difference"]["energy"][1] < -2.28e-11
    assert rows[0]["front_energy_origin_for_zero_response"][1] < -3.32
    assert rows[1]["front_energy_origin_for_zero_response"][0] > 6.34
    bank = json.loads((HERE / "late_read.json").read_text())
    for row in rows:
        raw = current.clock_transfer(row["derivatives"], bank)
        assert raw["coefficient_error_upper"] < 2e-7
        for key, held in raw["held"].items():
            assert abs(held["direct"] - held["from_clock_probabilities"]) < 2e-11
            np.testing.assert_allclose(
                held["direct"],
                row["independent_clock_flux"]["held"][key]["direct"],
                atol=1e-14,
            )


def test_regional_charge_and_two_native_currents_exhaust_the_missing_directions():
    from . import source_rank

    saved = json.loads((HERE / "source_rank.json").read_text())
    for record in (saved, source_rank.run()):
        assert record["inverse_defect_upper"] < 0.012
        assert (
            record["rank_with_regional_entropy_and_both_energy_currents"]
            - record["front_rank"]
            == 3
        )
