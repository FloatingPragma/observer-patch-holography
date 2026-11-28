"""Offline recalculation from archived records. Run from any directory with Python 3.11+."""

from fractions import Fraction
from itertools import combinations
import json
from math import exp, fsum, isclose, log, pi, sqrt
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parents[3]
OW = ROOT / "model/candidate_laws/interval_echo/openwave"
SPHERE = ROOT / "model/candidate_laws/interval_echo_published_records/spherical_dust"
BAO = ROOT / "code/cosmology/paid_trail_bao"


def read(path):
    return json.loads(path.read_text())


def overlap():
    records = read(
        ROOT / "openwave/interval_echo/row_19/selected_collar/parked_energy_charge.json"
    )
    counts = records["counts"]
    ref, fresh, shared = (counts[k] for k in ("reference", "fresh", "shared"))
    n = sum(ref)
    assert len(ref) == 81 and sum(fresh) == sum(shared) == n
    change = [h - o for h, o in zip(shared, fresh, strict=True)]
    assert all(
        sum(change[9 * i : 9 * i + 9]) == sum(change[i::9]) == 0 for i in range(9)
    )
    q = -fsum(v * log(w / n) for v, w in zip(change, ref, strict=True)) / n
    assert isclose(q, log(169 / 105) / 144, abs_tol=1e-15)
    return q


def finite_balance():
    archive = read(SPHERE / "results.json")
    stocks = archive["ball_surveys"]
    summary = []
    frames = ["rest", "+x", "-x", "+y", "-y", "+z", "-z", "xy", "xz", "yz"]
    assert [row["inverse_radius"] for row in stocks["ordinary"]] == [512, 2048, 8192]
    for ordinary, held, half in zip(
        stocks["ordinary"], stocks["held"], stocks["half_collar_force"], strict=True
    ):
        n = ordinary["inverse_radius"]
        assert n == held["inverse_radius"] == half["inverse_radius"]
        assert [ob["frame"] for ob in ordinary["observers"]] == frames
        assert [ob["frame"] for ob in half["observers"]] == frames
        eta = archive["ordinary_calibration"][str(n)]["eta"]
        for first, second in zip(ordinary["observers"], held["observers"], strict=True):
            assert first["frame"] == second["frame"]
            s0, s1 = first["entropy_over_r4"], second["entropy_over_r4"]
            a0, a1 = first["area_over_r4"], second["area_over_r4"]
            floor = abs(s1 - s0) / 2
            errors = [abs(s0 + eta * a0), abs(s1 + eta * a1)]
            # Use finite residuals for the descriptive comparison. The
            # continuum conclusion uses the separate theorem's o(ell^4) bounds.
            assert s1 - s0 > sum(errors) and a1 < a0 and eta > 0
            if first["frame"] == "rest":
                half_rest = half["observers"][0]
                half_error = abs(
                    half_rest["entropy_over_r4"] + eta * half_rest["area_over_r4"]
                )
                assert abs(half_error - floor) < 2e-6 and max(errors) < 1e-6
                assert isclose(
                    half_error,
                    abs(half_rest["calibrated_balance_over_r4"]),
                    abs_tol=1e-14,
                )
                summary.append(
                    dict(
                        inverse_radius=n,
                        common_area_error_floor=floor,
                        actual_max_residual=max(errors),
                        area_difference=a1 - a0,
                        half_coupling_error=half_error,
                    )
                )
    return summary


def regional_clock_check():
    rows = read(SPHERE / "clock_pairs/regions.json")["held_out"]
    ordinary, held = rows["ordinary"], rows["held"]
    for key in ("frame", "center", "proper_time", "inverse_radius"):
        assert ordinary[key] == held[key]
    floor = abs(held["entropy_over_r4"] - ordinary["entropy_over_r4"]) / 2
    actual = max(abs(row["balance_over_r4"]) for row in (ordinary, held))
    frozen = abs(rows["frozen_pair_clock"]["balance_over_r4"])
    assert actual < floor < frozen
    return {
        "common_area_error_floor": floor,
        "measured_max_residual": actual,
        "frozen_pair_clock_error": frozen,
    }


def releases():
    data = read(OW / "row_26_newton/radial_release.json")
    tick = Fraction(data["circumference_over_2pi_tick"])
    tau = Fraction(data["duration"])
    windows = {}
    assert set(data["experiments"]) == {
        "blank",
        "ordinary",
        "held",
        "changed_outer_profile",
        "doubled_coupling",
    }
    for name, arm in data["experiments"].items():
        assert [r["start_radius_ticks"] * tick for r in arm["radii"]] == [
            2,
            4,
            6,
            8,
            16,
        ]
        assert all(len(r["end_radius_ticks_masses_1_3"]) == 2 for r in arm["radii"])
        bounds = []
        for row in arm["radii"]:
            start = row["start_radius_ticks"] * tick
            for ticket in row["end_radius_ticks_masses_1_3"]:
                end = ticket * tick
                # Both reported radii are conservatively given half a bin.
                loss = start - end
                lo = 2 * max(0, loss - tick) * (end - tick / 2) ** 2 / tau**2
                hi = 2 * (loss + tick) * (start + tick / 2) ** 2 / tau**2
                bounds.append((lo, hi))
        windows[name] = max(v[0] for v in bounds), min(v[1] for v in bounds)
    for name, expected in (
        ("blank", 0),
        ("ordinary", Fraction(1, 12)),
        ("held", Fraction(19, 132)),
        ("doubled_coupling", Fraction(19, 66)),
    ):
        lo, hi = windows[name]
        assert lo <= expected <= hi
    assert windows["changed_outer_profile"][0] > windows["changed_outer_profile"][1]
    return {name: list(map(float, values)) for name, values in windows.items()}


def optics():
    data = read(OW / "row_27_metric/dust_exterior/results.json")
    tick = data["angular_and_time_tick"]
    results = {}
    assert set(data["sources"]) == {"blank", "ordinary", "held"}
    for stock, source in data["sources"].items():
        assert [shot["local_sine"] for shot in source["shots"]] == [0.5, 0.75]
        emitted = source["emitted_inner_ticks"]
        received = source["received_outer_ticks"]
        frequency = (emitted[1] - emitted[0]) / (received[1] - received[0])
        for shot in source["shots"]:
            reading = shot["reading"]
            assert reading["event"] == "returned"
            assert (
                reading["matched_blank_mirror_address_ticks"]
                == reading["mirror_address_ticks"]
            )
            angles = sum(
                reading[k]
                for k in (
                    "mirror_address_ticks",
                    "emission_angle_ticks",
                    "reception_angle_ticks",
                )
            )
            bend = angles * tick - pi
            delay = (
                reading["echo_clock_ticks"] - reading["matched_blank_echo_clock_ticks"]
            ) * tick
            # Absolute tick rounding can put the straight blank infinitesimally below zero.
            for key, value in (
                ("bending", bend),
                ("delay", delay),
                ("clock_ratio_2_to_4", frequency),
            ):
                lo, hi = shot["forecast"][key]
                assert lo - 4 * tick <= value <= hi + 4 * tick
            values = {"bend_rad": bend, "delay": delay, "frequency_ratio": frequency}
            if stock == "held":
                control = shot["angular_light_control"]
                assert control["event"] == "returned"
                assert (
                    control["matched_blank_mirror_address_ticks"]
                    == control["mirror_address_ticks"]
                )
                changed_angles = sum(
                    control[k]
                    for k in (
                        "mirror_address_ticks",
                        "emission_angle_ticks",
                        "reception_angle_ticks",
                    )
                )
                changed_bend = changed_angles * tick - pi
                changed_delay = (
                    control["echo_clock_ticks"]
                    - control["matched_blank_echo_clock_ticks"]
                ) * tick
                assert (
                    changed_bend < shot["forecast"]["bending"][0]
                    and changed_delay < shot["forecast"]["delay"][0]
                )
                values.update(
                    changed_bend_rad=changed_bend, changed_delay=changed_delay
                )
            results[f"{stock}/{shot['local_sine']}"] = values
    return results


def likelihood_score(table, cov, z, d, matter, b):
    def dh(t):
        return 1 / sqrt(matter * (1 + t) ** 3 + b)

    chisq = 0.0
    for zz in sorted(set(z)):
        ids = [i for i in range(len(z)) if z[i] == zz]
        assert all(cov[i][j] == 0 for i in ids for j in range(len(z)) if j not in ids)
        # Independent composite Simpson integration for D_M/r_d.
        steps = 256
        h = zz / steps
        dm = (
            (
                dh(0)
                + dh(zz)
                + fsum((4 if i % 2 else 2) * dh(i * h) for i in range(1, steps))
            )
            * h
            / 3
        )
        pred = {
            "DH_over_rs": dh(zz),
            "DM_over_rs": dm,
            "DV_over_rs": (zz * dm**2 * dh(zz)) ** (1 / 3),
        }
        diff = [d[i] - pred[table[i][2]] for i in ids]
        if len(ids) == 1:
            chisq += diff[0] ** 2 / cov[ids[0]][ids[0]]
        else:
            i, j = ids
            a, c, e = cov[i][i], cov[i][j], cov[j][j]
            u, v = diff
            chisq += (e * u * u - 2 * c * u * v + a * v * v) / (a * e - c * c)
    return chisq


def acoustic():
    prefix = "desi_gaussian_bao_ALL_GCcomb_"
    table = [
        line.split()
        for line in (BAO / (prefix + "mean.txt")).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    cov = [
        list(map(float, line.split()))
        for line in (BAO / (prefix + "cov.txt")).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    z = [float(row[0]) for row in table]
    d = [float(row[1]) for row in table]
    radial = [i for i, row in enumerate(table) if row[2] == "DH_over_rs"]
    assert len(radial) == 6
    k = NormalDist().inv_cdf(1 - 0.04 / (2 * len(radial)))
    scale = (100 * 147.87 / 299792.458) ** 2
    floors = []
    for ids in combinations(radial, 3):
        s = [1 + Fraction(table[i][0]) for i in ids]
        c = [s[(j + 1) % 3] ** 2 - s[(j + 2) % 3] ** 2 for j in range(3)]
        norm = sum(t**3 * a for t, a in zip(s, c))
        w = [a / norm for a in c]
        assert [sum(a * t**p for a, t in zip(w, s)) for p in (0, 2, 3)] == [0, 0, 1]
        total = fsum(
            float(a) / (d[i] + (k if a > 0 else -k) * sqrt(cov[i][i])) ** 2
            for a, i in zip(w, ids)
        )
        if total > 0:
            matter = total / scale - 0.00005 * max(
                0, float(sum(a * t**4 for a, t in zip(w, s)))
            )
            floors.append((float((matter - 0.003) / 0.02383 - 1), [z[i] for i in ids]))
    abundance = max(floors)
    assert 1.23 < abundance[0] < 1.24

    # Score the noisy radial calibration and the full-data matter fit against
    # the same measured covariance. Their predictive reference laws differ.
    budget = 6 * pi / (1.6309720958588974 * 0.041124247441816685)
    x = (budget + sqrt(budget**2 - 24 * pi)) / 2
    mpc = 1e6 * 149597870700 * 648000 / pi
    b = exp(-x) * (147.09 * mpc / 1.616255e-35) ** 2

    anchor = max(radial, key=lambda i: z[i])
    anchor_matter = (d[anchor] ** -2 - b) / (1 + z[anchor]) ** 3
    anchor_chisq = likelihood_score(table, cov, z, d, anchor_matter, b)
    assert 31.33 < anchor_chisq < 31.35
    chisq = likelihood_score(table, cov, z, d, 0.0003542525052, b)
    assert 17.26 < chisq < 17.28
    return dict(
        fixed_capacity_B=b,
        chi_squared_at_profiled_matter=chisq,
        chi_squared_at_single_radial_anchor=anchor_chisq,
        minimum_extra_per_baryon=abundance[0],
        strongest_radial_triple=abundance[1],
    )


def main():
    if not __debug__:
        raise RuntimeError("Run without -O so that the scientific checks execute.")
    print(
        json.dumps(
            dict(
                modular_change_nats=overlap(),
                finite_balance=finite_balance(),
                regional_clock_check=regional_clock_check(),
                common_GM_windows=releases(),
                optical_records=optics(),
                DESI=acoustic(),
            ),
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
