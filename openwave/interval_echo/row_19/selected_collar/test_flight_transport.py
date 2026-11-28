"""Check independent signed instruments, exact phase implications and the regional calibration."""

from fractions import Fraction as F
import json
from math import comb, factorial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def test_native_transport_identity_and_measured_current():
    H = np.array([[1, 2 + 0.3j], [2 - 0.3j, -1]])
    P = np.diag([1, 0])
    Q = np.eye(2) - P
    x = np.array([0.6, 0.8j])
    rho = np.outer(x, x.conj())
    gamma = 0.3 + 0.4j
    dense = P @ rho @ P + Q @ rho @ Q
    cross = P @ rho @ Q
    active = dense + gamma * cross + gamma.conjugate() * cross.conj().T
    z = np.trace(H @ cross)
    np.testing.assert_allclose(np.diag(active), np.diag(dense))
    np.testing.assert_allclose(
        np.trace(H @ (active - dense)), 2 * (gamma.real * z.real - gamma.imag * z.imag)
    )
    current = 1j * (H @ P - P @ H)
    np.testing.assert_allclose(
        np.trace(current @ rho), 2 * np.imag(H[0, 1] * rho[1, 0])
    )
    record = read("directed_service")
    for h in record["histories"].values():
        for arm in h.values():
            for table in arm.values():
                assert all(
                    min(tape) >= 0 and sum(tape) == record["starts"] for tape in table
                )
    for v in record["continuity"]["overlap"].values():
        assert v["at_switch"][0] <= 0 <= v["at_switch"][1]
        assert overlap(v["from_directed_bonds"], v["from_position_slope"])
        assert all(
            a * b > 0 for a, b in (v["from_directed_bonds"], v["from_position_slope"])
        )


def test_energy_from_independent_late_fringe_and_early_bonds():
    r = read("flight_energy")
    N = r["starts"]
    # The forward derivative weights reproduce every power through order six.
    weights = [
        sum(-F((-1) ** (i + 1) * comb(6, i), i) for i in range(1, 7)),
        *(F((-1) ** (i + 1) * comb(6, i), i) for i in range(1, 7)),
    ]
    for degree in range(7):
        assert sum(w * i**degree for i, w in enumerate(weights)) == int(degree == 1)
    assert abs(sum(w * i**7 for i, w in enumerate(weights))) / factorial(7) == F(1, 7)
    for history, z in r["histories"].items():
        scores = {"work": F(0), "incoming_current": F(0)}
        for bond in z["bonds"]:
            for label in scores:
                plus, minus, zero = bond[label]
                assert plus + minus + zero == N and min(plus, minus, zero) >= 0
                scores[label] += F(bond["unit"]) * (plus - minus) / N / 256
        a, b = map(F, r["knobs"][history.split("/")[0]])
        results = r["result"]["histories"][history]
        for name, price in (
            ("bare", scores["work"]),
            ("overlap", a * scores["work"] - b * scores["incoming_current"]),
        ):
            lo, hi = results["forecasts"][name]
            assert lo <= price <= hi
            assert overlap((lo, hi), results["late_minus_dense"][name]["total"])
        for panel in z["ramsey"].values():
            for table in panel:
                assert (
                    sum(sum(site) for site in table) == N
                    and min(min(site) for site in table) >= 0
                )
    assert all(
        row["fixed_radius_remainder"][0] > 0 and row["fixed_volume_remainder"][0] > 0
        for row in r["result"]["cap_remainders"]
    )


def test_zero_work_and_source_predictions():
    z = read("zero_flight_work")
    for window in z["windows"].values():
        q = window["quadratures"]
        phase = window["phase"]

        def product(a, b):
            v = [F(x) * F(y) for x in a for y in b]
            return min(v), max(v)

        w = product(q["work"], phase[0])
        j = product(q["current"], phase[1])
        assert w[0] - j[1] <= 0 <= w[1] - j[0]
        assert overlap(window["work_at_center"], window["late_work_at_center"])
        assert (
            overlap(window["direct_charge"], window["modular"]["connected"])
            and window["direct_charge"][0] > 0
        )
        assert overlap(
            window["source_geometry"]["area"], window["direct_geometry"]["area"]
        )
        assert overlap(
            window["area_from_stationarity"], window["direct_geometry"]["area"]
        )
        assert window["area_from_stationarity"][1] < 0
        assert not overlap(*window["ordinary_impulse_equivalents"].values())


def test_cap_price_includes_the_shared_miss_outcome():
    from openwave.interval_echo.row_19.selected_collar.local_cap import local_prices

    block = np.arange(1, 9).reshape(2, 4) / 100
    direction = np.array([[0.001, -0.002, 0.003, 0.001], [0.004, 0.002, -0.001, 0.001]])
    price, _ = local_prices((np.array([np.diag(b) for b in block]), np.zeros(2)))

    def entropy(t):
        p = (block + t * direction).ravel()
        p = np.r_[p, 1 - p.sum()]
        return -p @ np.log(p)

    predicted = sum(
        p.diagonal().real @ v
        for p, v in zip(price["modular"][0], direction, strict=True)
    )
    np.testing.assert_allclose(
        (entropy(0.0001) - entropy(-0.0001)) / 0.0002, predicted, rtol=2e-9
    )


def test_local_energy_reference_certificates():
    evidence = read("cap_energy_zero")
    for arms in evidence["counts"].values():
        for tape in arms.values():
            assert sum(map(sum, tape)) == evidence["starts"]
    for scale, result in evidence["tests"].items():
        target = result["target"]
        for name, orthants in result["groups"].items():
            for signs, entry in orthants.items():
                cols = {
                    k: {n: v if sign == "+" else (-v[1], -v[0]) for n, v in col.items()}
                    for (k, col), sign in zip(
                        evidence["offset_columns"][name].items(), signs, strict=True
                    )
                }
                if scale == "free":
                    cols["flight_energy"] = evidence["energy"]
                witness = entry["strictly_negative_witness"]
                if witness:
                    weights = {k: F(*v) for k, v in witness["weights"].items()}
                    assert min(weights.values()) >= 0
                    margin = min(
                        sum(w * F(cols[k][n][0]) for k, w in weights.items()) - F(v[1])
                        for n, v in target.items()
                    )
                    assert margin == F(*witness["strict_margin"]) > 0
                for case in ("zero", "nonpositive"):
                    cert = entry[case]
                    if case == "zero" or (
                        scale == "fixed" and name in ("common", "parity")
                    ):
                        assert cert["status"] == "excluded"
                    if cert["status"] != "excluded":
                        continue
                    columns = dict(cols)
                    if case == "nonpositive":
                        columns.update(
                            {
                                f"surplus/{n}": {p: (-int(n == p),) * 2 for p in target}
                                for n in target
                            }
                        )
                    lhs = dict.fromkeys(columns, F(0))
                    rhs = F(0)
                    for n, side, num, den in cert["coefficients"]:
                        i = int(side == "upper")
                        assert num > 0 and den > 0
                        w = F(num, den) * (-1 if i else 1)
                        rhs += w * F(target[n][1 - i])
                        for k, col in columns.items():
                            lhs[k] += w * F(col[n][i])
                    assert min(lhs.values()) >= 0 > rhs


def test_exact_volume_work_bracket_and_independent_observers():
    r = read("fixed_volume_flight")
    a, b = r["amplitude_bracket"]
    assert 0 < a < b < 1
    assert r["volume_at_ends"][0][0] > 0 > r["volume_at_ends"][1][1]
    qlo, qhi = map(F, r["full_at_exact_volume"])
    eps = F(r["stationarity_error_upper"])
    eta_min, eta_max = map(F, r["ordinary_gain"])
    lo, hi = r["area_at_exact_volume"]
    assert 0 < eps < qlo and lo <= float((-qhi - eps) / eta_min)
    assert float((-qlo + eps) / eta_max) <= hi < 0
    assert overlap((lo, hi), r["direct_geometry"]["area"])
    assert r["direct_geometry"]["area"][1] < 0
    assert r["direct_geometry"]["volume"][0] <= 0 <= r["direct_geometry"]["volume"][1]
    for name in ("full", "connected"):
        assert overlap(r[name + "_at_exact_volume"], r["direct_modular"][name])
    for tape in r["clocks"].values():
        assert len(tape) == 30
        assert all(sum(counts) == r["starts"] for counts in tape)
    for name in ("equivalent_at_exact_volume", "equivalent_at_programmed_center"):
        assert not overlap(*(z[name] for z in r["material"].values()))


def test_three_radius_obstruction_and_energy_convention():
    result = read("cap_stress_radii")
    for row in result["balls"]:
        qlo, qhi = map(F, row["intercept"])
        assert qlo * qhi > 0
        lo, hi = map(F, row["three_radius_contrast"])
        assert lo <= -qhi / 4 <= -qlo / 4 <= hi
        b = row["fixed_radius"]
        assert b["remainder"][0][0] > 0
        if qlo > 0:
            assert hi < 0 and row["minimum_negative_shell_energy"] > 0
        else:
            assert b["remainder"][1][1] < 0
        # Compare the finite-difference witness to separately decoded radii.
        direct = [(F(1, 2), F(1)), (F(1, 4), F(-5)), (F(1, 8), F(4))]
        x = y = F(0)
        for (r, c), values in zip(direct, b["remainder"], strict=True):
            lower, upper = sorted(r * c * F(v) for v in values)
            x += lower
            y += upper
        assert x <= lo <= hi <= y
        for frame in ("fixed_radius", "fixed_volume"):
            z = row[frame]
            radius = F(z["radius_retaining_half_intercept"])
            slope = max(map(abs, map(F, z["coefficient_of_radius"])))
            assert 2 * radius * slope <= min(abs(qlo), abs(qhi)) * (1 + F(1, 10**13))
    # Extremizing shell is at the middle radius. Continuity and the
    # three linear pieces give the same bound between these exact knots.
    for x in (F(0), F(1, 64), F(1, 16), F(1, 4), F(1)):
        val = sum(c * max(0, r * r - x) for r, c in direct)
        assert 0 <= val <= F(3, 16)


def test_resolved_addressed_shell_weights_and_held_word_readout():
    record = read("resolved_shells")
    for panel in record["levels"]:
        total = panel["cells"]
        h = panel["eccentricity"]
        radius = panel["radius_hops"]
        weights = list(map(F, panel["shell_weights_in_kappa_over_pi"]))
        assert (
            len(weights) == radius
            and min(weights) > 0
            and F(panel["next_weight_in_kappa_over_pi"]) > 0
        )
        count = 0
        for j, n in enumerate(panel["shell_populations"][: radius + 1]):
            count += n
            x = F(j, h)
            expected = x * F(count * (count - 1), total * (total - 1))
            integrated = sum(
                v * max(0, x * x - F(k, h) ** 2) for k, v in enumerate(weights)
            )
            assert expected == integrated
        assert overlap(
            panel["forecast_entropy_tangent"], panel["direct_entropy_tangent"]
        )
        assert panel["subhop_secant_change"][1] < 0
        for arm in panel["word_counts"].values():
            assert min(arm) >= 0 and sum(arm) == record["starts"]
        assert F(panel["eligible_pair_probability"]) == F(
            count * (count - 1), total * (total - 1)
        )


def test_phase_spectroscopy_reads_the_price_and_rejects_a_common_time_change():
    row = read("word_phase_read")
    for schedule, clock in row["clocks"].items():
        for board in clock["boards"]:
            for setting in board["tape"].values():
                assert all(
                    sum(count) == row["starts"] and min(count) >= 0
                    for count in setting.values()
                )
            for q in board["readings"]:
                assert overlap(
                    q["native_time_derivative"], q["measured_time_derivative"]
                )
                assert overlap(q["held"], q["from_measured_modular"])
                assert overlap(q["held"], q["source"])
                if schedule == "2":
                    assert q["held"][0] * q["held"][1] > 0
            if schedule == "2":
                assert not overlap(*(q["time_equivalent"] for q in board["readings"]))
    coarse = row["clocks"]["2"]["boards"][0]["readings"]
    assert not overlap(*(q["held_time_equivalent"] for q in coarse))
