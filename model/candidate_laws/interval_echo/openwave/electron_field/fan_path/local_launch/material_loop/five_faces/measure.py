"""Replay charge, reference, ring-current, angular and outside-head records."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from random import Random

import numpy as np

from .....fan_recoil.self_start.measure import bound, draw
from .....replay import compare
from ..measure import cohorts, matched_axis
from .apron import response
from .calibration import reported_interval, rotor_check, ruler
from .local_read import read as local_observers
from .spectrum import closed_form, flux_test, lines, readers, sectors, spectral_filter
from .wheel import LAND, RIM, preparation, wheel_basis, wheel_matrices


HERE = Path(__file__).resolve().parent
READS = HERE.parents[1] / "word_reads.json"
TRIALS = 2**37


def inputs():
    case = next(
        c for c in json.loads(READS.read_text())["cases"] if c["fan_count"] == 4
    )
    stock = next(a for a in case["stock_arms"] if a["preparation"] == "clock_bypassed")
    selected = []
    for u in (27.0625, 30.8125):
        f = next(
            f
            for f in stock["flights"]
            if f["appointment"] == u and f["material_bridge_rate"] == 4
        )
        pair = next(
            p
            for p in f["pair_readings"]
            if p["port_face_and_words"] == [1045, 0, 334, 365]
        )
        selected.append(
            (
                u,
                np.array(pair["receiver_reference_pauli_moments"]),
                f["observable_absolute_error_upper_per_pair"],
            )
        )
    axis = matched_axis(selected[-1][1])
    return axis, [
        (u, {s: r for s, (r, _) in cohorts(t, axis).items()}, 4 * pad + 4e-10, t)
        for u, t, pad in selected
    ]


def positional_read(basis, rho):
    prob = np.diag(rho).real
    return np.array(
        [sum(p for p, (x, _) in zip(prob, basis, strict=True) if x == v) for v in RIM]
    )


def angular_movie(basis, source):
    harmonics = np.exp(4j * np.pi * np.arange(5) / 5)
    zero = positional_read(basis, source)
    first = harmonics @ zero / zero.sum()
    assert abs(first) > 1e-3
    h, _ = wheel_matrices(basis, 0.5, 1 / 64)
    spec = sectors(basis, h)
    splitting = spec[1][0][0] - spec[-1][0][0]
    cases = []
    for sign in (-1, 0, 1):
        h, _ = wheel_matrices(basis, 0.5, sign / 64)
        e, v = np.linalg.eigh(h)
        for q in range(5):
            elapsed = q * np.pi / (2 * abs(splitting))
            u = (v * np.exp(-1j * e * elapsed)) @ v.conj().T
            final = positional_read(basis, u @ source @ u.conj().T)
            mode = harmonics @ final / final.sum()
            target = first * np.exp(-1j * sign * splitting * elapsed)
            angle = float(np.angle(mode / target))
            assert abs(angle) < 0.04
            cases.append(
                {
                    "external_flux_sign": sign,
                    "quarter_turns": q,
                    "wait": elapsed,
                    "all_five_charge_bins_per_original": final.tolist(),
                    "azimuth_second_harmonic": [mode.real, mode.imag],
                    "angular_error_from_splitting": angle,
                }
            )
    return cases


def interval_mean(j, p, pad):
    return [
        min(a / b for a in (j - 4 * pad, j + 4 * pad) for b in (p - pad, p + pad)),
        max(a / b for a in (j - 4 * pad, j + 4 * pad) for b in (p - pad, p + pad)),
    ]


def exterior_counts(arms, admitted, pad, rng):
    for arm in arms:
        for dwell in arm["dwells"]:
            bins = []
            for s in ("plus", "minus"):
                z = dwell["records"][s]
                bins.extend(
                    [
                        z["outside_reversed"],
                        z["source_admission"]
                        * z["original_energy_and_outside_word_given_source"],
                        z["other_source_energy_with_outside_unchanged"],
                    ]
                )
            bins.extend([admitted - sum(bins), 1 - admitted])
            counts = draw(bins, TRIALS, rng)
            intervals = []
            for i in (0, 3):
                num = bound(counts[i], TRIALS, pad)
                den = bound(sum(counts[i : i + 3]), TRIALS, pad)
                intervals.append((max(0, num[0] / den[1]), min(1, num[1] / den[0])))
            delta = [
                intervals[0][0] - intervals[1][1],
                intervals[0][1] - intervals[1][0],
            ]
            dwell["joint_original_counts"] = counts
            dwell["conditional_plus_minus_outside_difference_simultaneous_99pct"] = (
                reported_interval(*delta)
            )
            if arm["exterior_face_rate"] == 0.5 and dwell["wait"] in (4, 8):
                assert delta[1] < 0


def run():
    route = preparation()
    basis = wheel_basis()
    ids = [basis.index((LAND, w)) for w in route["wheel_input_words"]]
    embed = np.eye(len(basis))[:, ids]
    axis, arrivals = inputs()
    rng = Random(10471045)
    local_rng = np.random.RandomState(10461048)
    casebook = []
    for g in (0.5, 0.25, 0):
        h, j, m, rows = lines(basis, g)
        counter = readers(basis, j)
        target = next(i for i, (_, _, tags, _) in enumerate(rows) if (1, 0) in tags)
        select = spectral_filter(h, rows, target)
        flux = flux_test(basis, g, m, j)
        for phase in (0, 1 / 32, -1 / 32, 1 / 64, -1 / 64, 2 * np.pi):
            closed_form(basis, wheel_matrices(basis, g, phase)[0], g, phase)
        arrivals_out = []
        for time, bank, pad, joint_moments in arrivals:
            local = {s: embed @ r @ embed.T for s, r in bank.items()}
            admitted = float(bank["unread"].trace().real)
            readings = []
            unread_current = {s: 0.0 for s in local}
            for e, p, tags, bracket in rows:
                after = {s: p @ r @ p for s, r in local.items()}
                measured = {}
                for s, r in after.items():
                    prob = r.trace().real
                    side = np.array(
                        [
                            [np.trace(a @ r).real, np.trace(b @ r).real]
                            for a, b in counter
                        ]
                    )
                    curr = 4 * (side[:, 1] - side[:, 0])
                    np.testing.assert_allclose(curr, curr[0], atol=2e-13, rtol=0)
                    np.testing.assert_allclose(
                        curr, [np.trace(a @ r).real for a in j], atol=2e-13, rtol=0
                    )
                    assert prob > pad
                    if g == 0 or len(tags) == 1:
                        assert abs(curr[0]) < 1e-13
                    unread_current[s] += curr[0]
                    measured[s] = {
                        "line_per_original": prob,
                        "five_side_source_destination_bins_per_original": side.tolist(),
                        "single_side_current_per_original": curr[0],
                        "current_given_pass": curr[0] / prob,
                        "conditional_current_enclosure": interval_mean(
                            curr[0], prob, pad
                        ),
                    }
                readings.append(
                    {
                        "energy": e,
                        "energy_bracket": bracket,
                        "character_and_band": tags,
                        "cohorts": measured,
                    }
                )
            for s in local:
                ps = sum(line["cohorts"][s]["line_per_original"] for line in readings)
                np.testing.assert_allclose(ps, local[s].trace().real, atol=2e-14)
            selected = {s: select @ r @ select for s, r in local.items()}
            out = {
                "source_wait": time,
                "original_face_refusal": 1 - admitted,
                "raw_probability_guard": pad,
                "energy_bank": readings,
                "current_per_original_with_energy_forgotten": unread_current,
                "stationary_charge_and_local_face": local_observers(
                    basis, embed, select, joint_moments, pad, local_rng, g == 0.5
                ),
            }
            if g == 0.5:
                tests = [
                    r["minimum_partial_transpose_interval"]
                    for r in out["stationary_charge_and_local_face"]
                ]
                far = next(
                    r["local_witness"]
                    for r in out["stationary_charge_and_local_face"]
                    if "local_witness" in r
                )
                if time == 30.8125:
                    assert (
                        all(hi < 0 for _, hi in tests)
                        and far["simultaneous_99pct_interval"][1] < 0
                    )
                else:
                    assert (
                        all(lo > 0 for lo, _ in tests)
                        and far["simultaneous_99pct_interval"][0] > 0
                    )
            if time == 30.8125 and g == 0.5:
                true = readings[target]["cohorts"]
                assert true["plus"]["conditional_current_enclosure"][1] < 0
                assert true["minus"]["conditional_current_enclosure"][0] > 0
                normalized = [
                    selected[s] / selected[s].trace() for s in ("plus", "minus")
                ]
                np.testing.assert_allclose(
                    normalized[0].real, normalized[1].real, atol=3e-13, rtol=0
                )
                draws = {}
                for s in ("plus", "minus", "unread"):
                    z = true[s]
                    a, b = z["five_side_source_destination_bins_per_original"][0]
                    line = z["line_per_original"]
                    cohort = np.trace(local[s]).real
                    bins = [
                        a,
                        b,
                        line - a - b,
                        cohort - line,
                        admitted - cohort,
                        1 - admitted,
                    ]
                    counts = draw(bins, TRIALS, rng)
                    aa, bb = (bound(k, TRIALS, pad) for k in counts[:2])
                    ll = bound(sum(counts[:3]), TRIALS, pad)
                    ends = [4 * (b - a) / den for a in aa for b in bb for den in ll]
                    draws[s] = {
                        "seeded_side_counts": counts,
                        "conditional_current_simultaneous_99pct": reported_interval(
                            min(ends),
                            max(ends),
                        ),
                    }
                    lo, hi = draws[s]["conditional_current_simultaneous_99pct"]
                    assert lo * hi > 0
                out["simulated_side_banks"] = draws
                out["flux_precession"] = angular_movie(basis, selected["unread"])
                ext = {
                    f"ideal_character_{k}": np.outer(
                        m[k][1][:, 0], m[k][1][:, 0].conj()
                    )
                    for k in (1, -1)
                }
                out["exterior_head"] = response(basis, route, select, selected, ext)
                exterior_counts(out["exterior_head"], admitted, pad, rng)
                calibration = ruler(
                    basis, ids[0], bank["unread"][0, 0].real, admitted, TRIALS, pad, rng
                )
                out["physical_flux_ruler"] = calibration
                out["first_gap_orbital_comparison"] = rotor_check(m, flux, calibration)
                # A real initial word mixture gives zero directed equilibrium current at every line.
                for _, p, _, _ in rows:
                    mix = (
                        p @ embed @ np.diag(np.diag(bank["unread"]).real) @ embed.T @ p
                    )
                    assert max(abs(np.trace(b @ mix)) for b in j) < 1e-13
            arrivals_out.append(out)
        casebook.append(
            {
                "face_rate": g,
                "selected_line_index": target,
                "flux_calibration": flux,
                "arrival_studies": arrivals_out,
            }
        )
    return {
        "design_sha256": sha256((HERE / "DESIGN.md").read_bytes()).hexdigest(),
        "probe_design_sha256": sha256((HERE / "PROBE.md").read_bytes()).hexdigest(),
        "incoming_receipt_sha256": sha256(READS.read_bytes()).hexdigest(),
        "preparation": route,
        "source_axis": axis.tolist(),
        "original_launches_per_simulated_side_bank": TRIALS,
        "count_order": [
            "source_port",
            "destination_port",
            "other_on_line",
            "other_energy",
            "other_reference",
            "face_refused",
        ],
        "experiments": casebook,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    d = json.loads(json.dumps(run(), allow_nan=False))
    if args.output:
        with args.output.open("x") as f:
            f.write(json.dumps(d, indent=1) + "\n")
    else:
        compare(d, json.loads((HERE / "records.json").read_text()), "five_faces")
        print(
            "Verified native tag transfer, stationary five-side records, field splitting and the exterior probe."
        )


if __name__ == "__main__":
    main()
