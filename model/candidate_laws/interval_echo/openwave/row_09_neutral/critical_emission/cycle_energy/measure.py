"""Read ordinary active-contact energy on the actual periodic daughter."""

import argparse
from hashlib import sha256
import json
from math import log, sin, sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .......executable_model.carrier import federation
from ......contrast_memory_spherical_v1 import law
from ....row_20_electric import experiment as prior
from ....row_20_electric.schedule import DELTA, rotate_pairs
from ...passive_decay.measure import compare_records
from ..apparatus import Guide
from ..emitter import launch, reflected, selected_parent
from ..mode_reader import analyzer
from ..release import LEVELS, SHOTS, calibrate, enroll, local_bins, measure_contact


HERE = Path(__file__).resolve().parent
DRIVES = (("native", 2 / 65), ("half", 1 / 65))
NATIVE_LINK = 0.5


def layers(n, duty, factor):
    """Construct each native Majorana layer and its uniform-scale derivative."""
    j, h = prior.couplings(8)
    a_indices = np.arange(0, 2 * n, 2)
    b_indices = np.arange(1, 2 * n - 1, 2)
    a_angle, b_angle = -2 * h * DELTA, 2 * j * DELTA * duty
    a, b = [np.eye(2 * n) for _ in range(2)]
    rotate_pairs(a, a_indices, a_angle)
    rotate_pairs(b, b_indices[1:], b_angle)
    rotate_pairs(b, b_indices[:1], b_angle * factor)
    ga, gb = np.zeros_like(a), np.zeros_like(b)
    ga[a_indices, a_indices + 1] = a_angle
    gb[b_indices, b_indices + 1] = b_angle
    gb[1, 2] *= factor
    return a, b, ga - ga.T, gb - gb.T


def slope_and_error(g, a, b, ga, gb):
    whole = b @ a
    np.testing.assert_allclose(whole, g.r, atol=2e-14, rtol=0)
    raised = gb @ whole + b @ ga @ a
    v = g.modes[:, 0] / sqrt(2)
    value = np.exp(-1j * g.phi[0])
    # Verify every Schur pair and pairwise phase separation. The real
    # rotation is normal and these guards include the products of the
    # sparse 2-by-2 layers, represented eigenpairs and local contractions.
    if len(g.phi) > 1:
        gap = min(g.phi[1] - g.phi[0], 2 * g.phi[0])
    else:
        gap = 2 * g.phi[0]
    e = np.exp(-1j * g.phi)
    theory = g.modes / sqrt(2)
    residual = np.linalg.norm(whole @ theory - theory * e, ord="fro") + 2**-39
    chord = 2 * sin(gap / 2) - 2**-38
    assert 0.0008 < chord and residual < 2**-38
    # Trace-norm allowance for a rank-two covariance change, plus the
    # frame defect of the precomputed vacuum and two sparse rotations.
    probability_error = 16 * g.n * residual / chord
    assert probability_error < 2e-6
    return -float(np.imag(np.vdot(v, raised @ v) / value)) / g.period, probability_error


def records(ports, board, covariances, admission, rng):
    signal, blank, signal_b, blank_b = covariances
    weights = 6 - board.degree
    edges = set(board.seams)
    for i in range(len(ports) // 2):
        assert weights[ports[2 * i]] == weights[ports[2 * i + 1]]
    rows = []
    for start in range(len(ports) - 1):
        u, v = start, start + 1
        pair = [int(ports[u]), int(ports[v])]
        assert tuple(sorted(pair)) in edges
        p, q = (signal_b, blank_b) if start % 2 else (signal, blank)
        tests = [local_bins(p, u, v), local_bins(q, u, v)]
        d, de, nd = measure_contact(tests[0], admission, rng)
        c, ce, nc = measure_contact(tests[1], (0.0, 0.0, 1.0), rng)
        rows.append(
            {
                "physical_ports": pair,
                "layer": "link" if start % 2 else "intradimer",
                "daughter_and_cold_bins": [d, c],
                "nominal_ordinary_energy_difference": float(
                    LEVELS @ (tests[0] - tests[1])
                ),
                "simulated_mean_energy_difference": de - ce,
                "sum_of_inverse_net_admissions": 1 / nd + 1 / nc,
            }
        )
    return rows


def summary(rows, coefficients, allowance, families=24):
    exact = sum(
        w * r["nominal_ordinary_energy_difference"]
        for w, r in zip(coefficients, rows, strict=True)
    )
    estimate = sum(
        w * r["simulated_mean_energy_difference"]
        for w, r in zip(coefficients, rows, strict=True)
    )
    clocks = [
        (1 if row["layer"] == "intradimer" else 33) * a
        for row, a in zip(rows, coefficients, strict=True)
    ]
    scale = sum(abs(x) for x in clocks)
    q = sum(
        a * a * r["sum_of_inverse_net_admissions"]
        for a, r in zip(clocks, rows, strict=True)
    )
    radius = sqrt(0.5 * log(2 * families / 0.01) * q) + scale * allowance
    return {
        "calculated": exact,
        "sampled": estimate,
        "simultaneous_99pct_interval": [estimate - radius, estimate + radius],
    }


def acquired_clock(panels, g):
    clock = panels["clock"]
    limits = np.array(
        [
            clock["energy"] - clock["simultaneous_error"],
            clock["energy"] + clock["simultaneous_error"],
        ]
    )
    return (0.1 * limits / g.period).tolist()


def terminal_charge_control(signal, duty, admission, rng, tolerance):
    """Compare full versus projected native transfer on the final link."""
    v = len(signal) - 2
    unlike = float(sum(local_bins(signal, v - 1, v)[-2:]))
    indices = [law.SRC.index((0, 2)), law.SRC.index((2, 0))]
    dur = DELTA * duty
    free = expm(-1j * dur * law.seam_term(8))
    changed = abs(free[indices[0], indices[1]]) ** 2
    block = free[np.ix_(indices, indices)]
    np.testing.assert_allclose(block @ block.conj().T, np.eye(2), atol=2e-14, rtol=0)
    np.testing.assert_allclose(changed, sin(dur / 2) ** 2, atol=1e-18, rtol=0)
    # The parity eigenstate gives equal weights to 02 and 20. The other
    # rail of the terminal dimer has no link during this layer.
    assert 0 < unlike < 1
    escaped = changed * unlike
    lower = changed * (unlike - tolerance)
    assert lower > 0
    ready = float(admission[2])
    kept = [*admission[:2], 0.0, ready, 0.0]
    removed = [
        *admission[:2],
        0.5 * escaped * ready,
        (1 - escaped) * ready,
        0.5 * escaped * ready,
    ]
    bins = [rng.multinomial(SHOTS, p).tolist() for p in (kept, removed)]
    assert bins[0][2] == bins[0][4] == 0 and bins[1][2] > 0 and bins[1][4] > 0
    return {
        "count_bin_order": [
            "parent refused",
            "lowest one-mode refused",
            "minus_two",
            "zero",
            "plus_two",
        ],
        "guarded_and_open_guard_counts": bins,
        "motion_per_admitted_exposed_dimer": escaped,
        "motion_lower_over_calibrated_covariance": lower,
        "motion_per_original_launch": ready * escaped,
        "dwell": dur,
        "unlike_terminal_link_probability_at_layer_opening": unlike,
    }


def run():
    rng = np.random.default_rng(90131)
    guard_rng = np.random.default_rng(90132)
    b = federation(2)
    _, physical, _ = prior.guide(b)
    archive = HERE.parent / "emitter.json"
    published = json.loads(archive.read_text())
    results = []
    for n in (23, 45):
        ports = physical[1 : 1 + 2 * n]
        a, rb, *_ = layers(n, DRIVES[0][1], 0)
        disconnected = Guide(n, "native", 0)
        full = rb @ a
        assert not np.any(full[2:, :2]) and not np.any(full[:2, 2:])
        # With the link cut the space of the initial X=- parent is invariant.
        # Other guide excitations can be measured, but this herald cannot pass.
        launch_cut = launch(disconnected, "parent")
        evolved = np.linalg.matrix_power(full, 4 * n)
        pg, _ = selected_parent(evolved @ launch_cut @ evolved.T, 1)
        assert pg == 0
        cut = [SHOTS, 0, *([0] * len(LEVELS))]
        for mode, duty in DRIVES:
            g = Guide(n, mode, NATIVE_LINK)
            a, rb, ga, gb = layers(n, duty, NATIVE_LINK)
            one, admission = enroll(g)
            checker, word, _ = analyzer(g)
            tags = word.copy()
            tags[0] *= -1
            np.testing.assert_allclose(
                (checker @ one @ checker.T)[::2, 1::2],
                np.diag(tags),
                atol=3e-12,
                rtol=0,
            )
            np.testing.assert_allclose(
                reflected(g.c, g.modes[:, 0].real), one, atol=2e-12, rtol=0
            )
            slope, allowance = slope_and_error(g, a, rb, ga, gb)
            rows = records(
                ports, b, (one, g.c, a @ one @ a.T, a @ g.c @ a.T), admission, rng
            )
            coefficients_a = [int(i % 2 == 0) for i in range(len(rows))]
            coefficients_b = [0, *[i % 2 for i in range(1, len(rows))]]
            coefficients_late_b = coefficients_b.copy()
            coefficients_late_b[1] = 0
            weighted = [
                (x + 0.5 * duty * (y + z)) / (1 + duty)
                for x, y, z in zip(
                    coefficients_a, coefficients_b, coefficients_late_b, strict=True
                )
            ]
            work = [y - x for x, y in zip(coefficients_a, coefficients_b, strict=True)]
            measured = {
                name: summary(rows, c, allowance)
                for name, c in (
                    ("intradimer", coefficients_a),
                    ("link_layer_including_first", coefficients_b),
                    ("link_layer_after_first_shutters", coefficients_late_b),
                    ("daughter_relative_initial_switch_work", work),
                    ("cycle_average", weighted),
                )
            }
            np.testing.assert_allclose(
                measured["cycle_average"]["calculated"], slope, atol=2e-11, rtol=0
            )
            threshold = acquired_clock(published, g)
            lo, hi = measured["cycle_average"]["simultaneous_99pct_interval"]
            inherited_emission = None
            if mode == "native":
                assert 0 < lo <= hi < threshold[0]
                assert (
                    g.phi[0] / g.period
                    < threshold[0]
                    < threshold[1]
                    < g.phi[1] / g.period
                )
                original = [x for x in published["panels"] if x["dimers"] == n]
                chosen = next(
                    x for x in original if x["source_link_fraction"] == NATIVE_LINK
                )
                lower = chosen["sources"][0]["joint_lower_99"]
                bright = chosen["bright_only_tag_upper"]
                assert lower > bright
                inherited_emission = {
                    "original_starts": published["shots"],
                    "joint_parent_ground_absorption_vacuum_counts": chosen["sources"][
                        0
                    ]["joint_counts"][0],
                    "lower_99pct_from_original_joint_banks": lower,
                    "bright_only_probability_upper": bright,
                    "lowest_mode_mass_lower_from_coherent_amplitude_triangle": (
                        sqrt(lower) - sqrt(bright)
                    )
                    ** 2,
                    "first_bright_relative_frequency": float(g.phi[1] / g.period),
                    "registered_source_duty_and_joint_lower": [
                        [x["source_link_fraction"], x["sources"][0]["joint_lower_99"]]
                        for x in original
                    ],
                }
            else:
                assert lo > threshold[1]
            results.append(
                {
                    "dimers": n,
                    "schedule": mode,
                    "relative_link_duty": duty,
                    "calibrated_period": g.period,
                    "lowest_relative_quasifrequency": float(g.phi[0] / g.period),
                    "slope_of_lowest_relative_phase_per_original_period": slope,
                    "light_cutoff_from_separate_clock": threshold,
                    "parent_refused_mode_refused_admitted": admission.tolist(),
                    "source_link_substitution_full_refusal": cut,
                    "daughter_minus_cold_contact_records": rows,
                    "observed_energy_and_switch_work": measured,
                    "conservative_covariance_allowance": allowance,
                    "inherited_same_mode_emission": inherited_emission,
                    "guard_only_terminal_charge_control": terminal_charge_control(
                        a @ one @ a.T, duty, admission, guard_rng, allowance
                    )
                    if mode == "native"
                    else None,
                }
            )
    return {
        "protocol_sha256": sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        "source_emission_sha256": sha256(archive.read_bytes()).hexdigest(),
        "seed": 90131,
        "guard_control_seed": 90132,
        "shots_per_contact_and_branch": SHOTS,
        "count_bin_order": [
            "parent refused",
            "lowest one-mode refused",
            *LEVELS.tolist(),
        ],
        "clock_calibration": calibrate(),
        "panels": results,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    data = run()
    if args.acquire:
        with args.acquire.open("x") as out:
            out.write(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        compare_records(json.loads((HERE / "results.json").read_text()), data)
        print(
            "Verified the original daughter, both contact layers and their measured work."
        )


if __name__ == "__main__":
    main()
