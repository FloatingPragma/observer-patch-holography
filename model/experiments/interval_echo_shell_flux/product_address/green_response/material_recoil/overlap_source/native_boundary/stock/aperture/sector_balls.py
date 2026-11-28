"""Regional source, ordinary marker geometry and fixed-volume stationarity."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np

from . import sector_ball_chart as chart
from . import sector_ball_control as control
from . import sector_ball_meter as meter
from . import sector_ball_local
from . import sector_loaded_ball
from . import sector_area_lab
from . import sector_ball_source as source
from . import sector_ball_survey as survey
from . import sector_sampling
from . import sector_regional

HERE = source.balls.HERE


def bracket_witness(charts, source_at_mid, delta, rates, words, reference):
    p0 = 29 / 388
    bracket = rates["continuous_cadence"]["root_bracket"]
    miss = max(abs(p - p0) / max(p, p0) for p in bracket)
    # Couple reporting coins by one uniform ticket until either observer
    # accepts. The chance of a mismatch then equals |p-p0|/max(p,p0).
    d, err = source_at_mid
    witnesses = []
    full_witness = None
    for (load, desc, r), row in charts.items():
        if load != 2 or desc:
            continue
        fixed = (
            row["fixed_members"]
            & row["fixed_boundary_count"]
            & (row["volumes"] < len(row["volumes"]))
        )
        for j in np.flatnonzero(fixed):
            mask = int(row["records"][0, j])
            if mask == 127 and full_witness is None:
                full_witness = dict(
                    radius=r,
                    center=int(j),
                    volume=int(row["volumes"][j]),
                    area=row["areas"][0, j].tolist(),
                )
            if mask not in (1, 2, 4, 8):
                continue
            k = mask.bit_length() - 1
            levels = words[:, k]
            probs = np.array([reference[levels == v].sum(0) for v in range(3)])
            prices = -np.log(probs[:, ::-1])
            weights = np.bincount(levels, weights=d, minlength=3)
            endpoints = (
                np.maximum(weights, 0) @ prices
                + np.minimum(weights, 0) @ prices[:, ::-1]
            )
            span = prices.max() - prices.min()
            radius = span * (sum(abs(delta)) * miss + err) + 1e-10
            band = endpoints + np.array([-radius, radius])
            margin = max(band[0], -band[1])
            witnesses.append(
                (
                    margin,
                    dict(
                        radius=r,
                        center=int(j),
                        volume=int(row["volumes"][j]),
                        area=row["areas"][0, j].tolist(),
                        home=k,
                        tangent_nats=band.tolist(),
                        conditional_clock_probabilities=probs.tolist(),
                        maximum_coupling_error=miss,
                        range_extension_nats=radius,
                    ),
                )
            )
    assert full_witness is not None
    best = max(witnesses, key=lambda x: x[0])
    return dict(
        source_rate_bracket=bracket,
        witness_rate="29/388",
        complete_ball=full_witness,
        bracket_ball=best[1],
        excludes_fixed_volume_stationarity_at_all_rates=bool(best[0] > 0),
    )


def donors(words, d, error):
    last = words[:, -1]
    singles = (words[:, :4] == 1).sum(1)
    partial = np.isin(last, source.experiment.arrival.source.PARTIAL)
    stats = np.array([last == 0, partial, singles, partial * singles], float)
    a, b, n, c = stats.mean(1)
    A, B, N, C = stats @ (d + 1 / len(words))
    coeff = np.array(
        [(4 * (A - a) + B - b) / b, (B * N - b * n) / b, (C - B * N - c + b * n) / b]
    )
    # Each probability error is at most the entrywise bound on the end
    # distribution. The count of ordinary singles never exceeds four.
    width = np.array([5.0, 8.0, 12.0]) * error / b + 1e-12
    return dict(
        per_pool_port_and_population=coeff.tolist(),
        coefficient_bands=np.column_stack((coeff - width, coeff + width)).tolist(),
    )


def run():
    words, hops, homes, seams, edges, geometry = chart.charts()
    fixed = chart.fixed_volumes(hops, homes, seams, edges, geometry)
    lookup = chart.regional_masks(hops, homes, seams)
    request = {0, 127}
    for masks in lookup.values():
        request.update(map(int, masks))
    rho = source.thermal.bath()
    cal_path = HERE / "sector_calorimetry.json"
    saved = json.loads(cal_path.read_text())
    banks = saved["ordinary_clock_calibration"]
    sharp, coherent, clock = source.meters(words, banks, request)
    rates = json.loads((HERE / "sector_retention.json").read_text())
    prepared, delta, subset, full, mats = source.sources((1, 8, 64, Q(388, 29)))
    bracket = bracket_witness(
        fixed,
        prepared.pop(Q(388, 29)),
        np.array(delta, float),
        rates,
        words,
        sector_regional.calibrate(banks)["retained"],
    )
    witness = bracket["bracket_ball"]
    center, radius, home = witness["center"], witness["radius"], witness["home"]
    witness["markers"] = np.flatnonzero(hops[center] < radius).tolist()
    for load in (1, 2):
        for reverse in (False, True):
            row = fixed[load, reverse, radius]
            assert row["fixed_members"][center] and row["fixed_boundary_count"][center]
            assert np.all(row["records"][:, center] == (1 << home))
            assert np.all(row["areas"][:, center] == witness["area"])
    # Fixing the other three home loads fixes the fourth on the full
    # support. Both one-home quantum marginals are consequently diagonal.
    assert np.ptp(words[:, :4].sum(1)) == 0
    weight = sector_regional.calibrate(banks)["exchanged"]
    witness["exchanged_clock_probabilities"] = [
        weight[words[:, home] == j].sum(0).tolist() for j in range(3)
    ]
    bracket["single_home_all_rates"] = sector_ball_local.certify(
        words, delta, subset, mats, rates, bracket["bracket_ball"]
    )
    local = meter.calibrate()
    _, _, g, blocks, ops = source.thermal.apparatus()
    heat = ops[0].diagonal(axis1=1, axis2=2).sum(0) / 4
    means = np.array([heat[i].mean() for i in blocks])
    thermo = (
        np.array([b["centered_log_band"] for b in banks])
        - source.thermal.BETA * means[:, None]
    )[g]
    pairs = control.contrasts(words, local, geometry, thermo, subset, full, mats)
    rng = np.random.default_rng(853124)
    tables = source.experiment.acquisition.returns.ticket_tables(mats)
    native = np.full(len(full), -1)
    native[subset] = np.arange(len(words))
    observed = {}
    stat = {}
    held = {}
    counts = np.asarray(delta, float)
    for div, copies in ((8, 1024), (64, 2048)):
        init, final, duration = sector_sampling.trajectories(
            tables, subset, native, rng, div, copies
        )
        obs = np.bincount(final, weights=counts[init] / copies, minlength=len(words))
        observed[div] = obs
        stat[div] = math.sqrt(
            (counts @ counts) * math.log(4 * (len(request) + 2) / 0.002) / (2 * copies)
        )
        outcomes = np.bincount(words[:, home], weights=obs, minlength=3)[[0, 2]]
        endpoints = np.column_stack((outcomes - stat[div], outcomes + stat[div]))
        assert np.all(endpoints > 0)

        def gains(probs):
            z = np.array(probs)
            return np.log(z[1] / z[[0, 2], ::-1])

        ret = gains(witness["conditional_clock_probabilities"])
        exc = gains(witness["exchanged_clock_probabilities"])
        modular = ret - exc[:, ::-1]
        assert np.all(modular > 0)
        # Save the signed word counts needed for every regional comparison.
        held[str(div)] = dict(
            trials=len(final),
            copies_per_start=copies,
            maximum_sweeps=duration,
            home_probability_gains=endpoints.tolist(),
            home_retained_entropy_nats=(ret * endpoints).sum(0).tolist(),
            home_relative_modular_nats=(modular * endpoints).sum(0).tolist(),
            terminal_counts=np.bincount(final, minlength=len(words)).tolist(),
            contrast_by_terminal_word=obs.tolist(),
        )
    rate_rows = {}
    keys = {
        "head": words[:, -1, None],
        "three_contact_means": ops[0].diagonal(axis1=1, axis2=2).T,
        "head_and_total_heat": np.column_stack((words[:, -1], heat)),
    }
    for div, (d, err) in prepared.items():
        redistribution = {}
        for name, key in keys.items():
            _, k, n = np.unique(key, axis=0, return_inverse=True, return_counts=True)
            residual = d - source.group_mean(k, n, d)
            redistribution[name] = float(residual @ sharp[127])
        source_heat = float(source.thermal.BETA * (d @ heat))
        rate_rows[str(div)] = dict(
            transported_double_equivalent_from_six_marker_area=-float(
                d @ (words[:, 1] == 0)
            ),
            ordinary_heat_nats=source_heat,
            source_excess=float(d @ sharp[127]),
            retained_tangent=source_heat + float(d @ sharp[127]),
            distribution_error_l1=err,
            residual_after_matched_renewal=redistribution,
            ordinary_relay=donors(words, d, err),
        )
    hot, cold = [
        source.thermal.bath(beta=source.thermal.BETA + step)["retained"].diagonal()
        for step in (-1 / 256, 1 / 256)
    ]
    loaded_stat = {
        div: math.sqrt((counts @ counts) * math.log(2 * 40 / 0.002) / (2 * copies))
        for div, copies in ((8, 1024), (64, 2048))
    }
    assert all(np.ptp((words[group, 1] == 0).astype(int)) == 0 for group in blocks)
    empty = words[:, 1] == 0
    assert sum(empty) == 24 and len(empty) == 448
    board = source.balls.energy.returned(words)[0]
    loaded = sector_loaded_ball.evaluate(
        words,
        homes,
        seams,
        board,
        center,
        prepared,
        rho["retained"].diagonal(),
        sector_regional.calibrate(banks)["retained"],
        hot - cold,
        observed,
        loaded_stat,
        np.array(pairs["ordinary_move_input"]),
    )
    eq0 = np.array(witness["exchanged_clock_probabilities"])[0]
    p_empty = Q(int(sum(empty)), len(empty))
    log_odds = np.log(float((1 - p_empty) / p_empty))
    # The two-outcome local instrument reports empty or occupied from its
    # membership. Derive its modular difference from the independent baths.
    beta_rel = log_odds - np.log((1 - eq0[::-1]) / eq0[::-1])
    beta_rel.sort()
    assert min(beta_rel) > 0
    six = loaded["surveys"]["ascending/6"]
    assert six["binary_identity"]["extra_boundary_links_per_empty_home"] == -4
    assert six["moved_double"]["area_change"] == 4
    assert all(
        loaded["surveys"][key]["distinct_geometries"] == 1
        for key in ("ascending/11", "ascending/13")
    )
    observed_lower = bracket["single_home_all_rates"]["outcomes"][0][
        "probability_gain_per_rate_lower"
    ]
    loaded["six_marker_binary_source"] = dict(
        log_price_difference_from_exchanged_reference=beta_rel.tolist(),
        local_excess_per_unit_reporting_probability_lower=float(
            beta_rel[0] * observed_lower - 1e-12
        ),
        boundary_decrease_per_unit_reporting_probability_lower=4 * observed_lower
        - 1e-12,
        ordinary_doubles_per_unit_empty_probability=-1,
        measured_at_rates={
            k: (np.array(v["home_probability_gains"][0]) * beta_rel).tolist()
            for k, v in held.items()
        },
    )
    return dict(
        independent_area_acquisition=sector_area_lab.read(
            tables, subset, native, words, loaded["six_marker_surface_by_inventory"]
        ),
        occupied_home_radar=loaded,
        retained_empty_home_probability=str(Q(int(sum(empty)), len(empty))),
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                cal_path,
                HERE / "sector_retention.json",
                HERE / "SECTOR_BALLS.md",
            )
        },
        common_chart_markers=len(hops),
        local_thermometer=local,
        controls=pairs,
        first_return_rates=rate_rows,
        held_native=held,
        joint_sharp_measurement_coverage_at_least=0.997,
        fixed_volume_test=bracket,
        fixed_volume_surveys=survey.adaptive(
            words,
            fixed,
            prepared,
            [rho[k].diagonal() for k in ("retained", "exchanged")],
            hot,
            cold,
        ),
        lapse_surveys=survey.lapse(
            hops, lookup, sharp, coherent, clock, observed, prepared, stat, words
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "sector_balls.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_text()), out
        )
    print(out["fixed_volume_test"])
    print(out["controls"]["equal_charge_pair"])
