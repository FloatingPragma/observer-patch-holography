"""Repeated native shutter firings bound the timing contribution to the held face mean.

Separate equal-degree hands receive the two arms' actual duration commands and
are read at their terminal rail. For a stable command of unknown duration t,
n uninterrupted firings give sin²(nt/2). The endpoint is the only measurement
in a train. The source, random arm and trit draws, incoming sharp records and
exterior holds remain as registered. A pi or 3pi preparation hold uses one or
three repetitions of the full-exchange clock command. Each pair of sensing
holds uses the same firing count and duration. A time-dependent phase offset,
loss of coherence, or device-specific difference between a hand's clock and a
guide's clock requires its own test. Counts are synthetic multinomial marginals.

Given overlapping confidence intervals, all compatible aliases are retained
until a unique interval survives. Observing no failures or no successes gives
an endpoint-adaptive absolute probability bound of 20/N. For other outcomes
the ordinary Hoeffding bound is used. A union of both possible errors, over
all banks, also covers this data-dependent choice. The subtraction of a scalar
from each isolated source exchange is valid since the schedule controls it
classically. Trace distance contracts under every subsequent coherent or
measured operation, including adaptive draining and recording. In a telescoping
replacement of all source pulses, the ideal prior segment has expected lifetime
at most 120 appointments. This gives a total-variation bound on the entire
stopped record of one release, including the final cutoff at 2048.
"""

from fractions import Fraction as F
from math import ceil, exp, floor, pi, log, sqrt
import json
from pathlib import Path

import numpy as np
from mpmath import iv

from . import measure, receiver

HERE = Path(__file__).resolve().parent
REPEATS = (1, 1021, 1021**2)
SHOTS = 2**33
GRID = 2**44
MARGIN = F(1, 32768)
NAMES = (
    "sensing subpulse",
    "full release",
    "half release",
    "half exchange",
    "full exchange",
)
ARMS = ("near", "far")


def standards():
    """Equal-degree material hands outside all routed registers and their neighbors."""
    board, guides, flag = receiver.original.rails.layout()
    active = set(flag).union(*guides)
    for alternate in (False, True):
        left, right, pad, ref, paths = receiver.apparatus(alternate)
        active.update((left, right, pad, ref))
        for path in paths:
            active.update(path)
    used = active.union(*(set(board.neighbours[a]) for a in active))
    choices = sorted(
        (int(a), int(b))
        for a, b in board.seams
        if a // 12 == b // 12 and board.degree[a] == board.degree[b]
    )
    clocks = []
    for arm in ARMS:
        left, right = next(
            (a, b) for a, b in choices if a not in used and b not in used
        )
        used.update([left, right, *board.neighbours[left], *board.neighbours[right]])
        for stiffness in (4, 8):
            diagonal, hop = receiver.original.terms(1, stiffness)
            # Blank neighbors stay held. Equal degrees give equal complete site
            # energies, including the original host term and every incident contact.
            energies = [
                stiffness * (6 - int(board.degree[site]))
                + int(board.degree[site]) * diagonal[1, 0]
                for site in (left, right)
            ]
            assert energies[0] == energies[1] and hop[0, 1] == -0.5 == hop[1, 0]
        clocks.append(
            dict(arm=arm, rails=[left, right], degree=int(board.degree[left]))
        )
    return clocks


def design():
    return dict(
        zip(
            NAMES,
            (iv.pi / 3, 2 * iv.pi / iv.sqrt(5), iv.pi / iv.sqrt(5), iv.pi / 2, iv.pi),
            strict=True,
        )
    )


def bounds(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def canonical_probability(duration, repeats):
    exact = iv.sin(duration * repeats / 2) ** 2
    low, high = bounds(exact * GRID + iv.mpf(1) / 2)
    a, b = floor(low), floor(high)
    if a != b:
        raise ArithmeticError(
            "Native dimer probability straddles the sampling-grid midpoint"
        )
    return a


def interval_probability(hits, shots):
    p = F(hits, shots)
    margin = F(20, shots) if hits in (0, shots) else MARGIN
    width = margin + F(1, 2 * GRID)
    return max(F(), p - width), min(F(1), p + width)


def angle(q):
    p = iv.mpf(q.numerator) / q.denominator
    if q == 0:
        return iv.mpf(0)
    if q == 1:
        return iv.pi
    return 2 * iv.atan2(iv.sqrt(p), iv.sqrt(1 - p))


def refine(prior, repeats, hits, shots):
    lo, hi = prior
    low, high = (angle(q) for q in interval_probability(hits, shots))
    # Both inverse-sine branches and every compatible turn. The prior alone
    # fixes which aliases survive; the design duration is not consulted.
    images = []
    for turn in range(
        floor(lo * repeats / (2 * pi)) - 2, ceil(hi * repeats / (2 * pi)) + 3
    ):
        base = 2 * turn * iv.pi
        for left, right in ((low, high), (-high, -low)):
            a = max(lo, bounds((base + left) / repeats)[0])
            b = min(hi, bounds((base + right) / repeats)[1])
            if a <= b:
                images.append((a, b))
    merged = []
    for a, b in sorted(images):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    if len(merged) != 1:
        raise ArithmeticError(
            f"Clock observations leave {len(merged)} compatible branches"
        )
    return merged[0]


def run():
    iv.dps = 55
    rng = np.random.default_rng(20019667)
    calibration = {}
    recipes = design()
    hands = standards()
    for arm in ARMS:
        clocks = {}
        for name, tau in recipes.items():
            # A single full-exchange service can be on either side of pi.
            # The other commands use the original first monotone [0,pi] branch.
            prior = [0.0, bounds((2 if name == "full exchange" else 1) * iv.pi)[1]]
            stages = []
            sequence = (*REPEATS, 1021**3) if name == "full exchange" else REPEATS
            for n in sequence:
                shots = 2**39 if n == 1021**3 else SHOTS
                ticket = canonical_probability(tau, n)
                hits = int(rng.binomial(shots, ticket / GRID))
                prior = refine(prior, n, hits, shots)
                stages.append(
                    dict(
                        repeats=n,
                        right=hits,
                        left=shots - hits,
                        probability_tickets=ticket,
                        duration_interval=prior,
                    )
                )
            target = bounds(tau)
            error = max(abs(prior[0] - target[1]), abs(prior[1] - target[0]))
            clocks[name] = dict(
                stages=stages, duration_error_upper=float(np.nextafter(error, np.inf))
            )
        calibration[arm] = clocks
    bank_count = len(ARMS) * (len(NAMES) * len(REPEATS) + 1)
    alpha = bank_count * 2 * (exp(-2 * SHOTS * float(MARGIN) ** 2) + exp(-20))
    # Share the original 1% budget with the Markov experiment by allocating
    # 1/(panels+1) of it to all duration readings combined.
    parents = list(measure.settings())
    assert alpha < 0.01 / (1 + len(parents))
    u = measure.visits_before_drain()
    panels = []
    raw = json.loads((HERE / "results.json").read_text())["panels"]
    for setting, record in zip(parents, raw, strict=True):
        # Actual material routing uses one forward and one restored endpoint
        # transposition for each guide port and for the independent reference.
        # One further swap advances the selected guide on a live appointment.
        *_, paths = receiver.apparatus(setting.equal_degree)

        def length(path):
            return max(0, 2 * len(path) - 3)

        per_visit = 1 + 2 * (
            sum(max(map(length, paths[i : i + 12])) for i in (0, 12, 24))
            + 3 * length(paths[-1])
        )
        max_time_error = max(
            calibration[a]["full exchange"]["duration_error_upper"] for a in ARMS
        )
        # On the nine sharp rails the extreme eigenvalues are 0 and 4*kappa+1.
        centered_norm = (4 * setting.stiffness + 1) / 2
        law = receiver.original.rails.candidate.native_law
        h = law.seam_term(setting.stiffness)
        assert np.max(abs(h - centered_norm * np.eye(len(h))).sum(1)) <= centered_norm
        # Compare full individual records by telescoping with the ideal sharp
        # prefix. Dead streams use no gates, and the ideal mean lifetime <=120.
        # The initial corner uses one flag pulse and up to three guide pulses.
        total_variation = (4 + 120 * per_visit) * centered_norm * max_time_error
        horizon = measure.PLAN["maximum_appointments"]
        route_bias = horizon * total_variation
        errors = []
        for hand in calibration.values():

            def err(name):
                return hand[name]["duration_error_upper"]

            # Three independent phase channels, trit counts 0,4,8. A general
            # subsequent right-port effect has absolute off-diagonal <= 1/2.
            energy = receiver.native_terms(
                setting.probe, setting.stiffness, setting.equal_degree
            )[0]
            gaps = energy[:, 0] - energy[:, 1]
            phase_gap = float(max(gaps) - min(gaps))
            sense = 6 * phase_gap * err("sensing subpulse")
            flight = 0.5 * err("half release" if setting.short else "full release")
            # The one-load standing preparation has isolated exchange norm 1/2
            # after phase removal, and blank rail detuning 1/2 in its holds.
            # Six isolated rail interchanges invert/restore three comparisons.
            prep = 0.5 * err("half exchange") + 0.75 * err("full exchange")
            swaps = 3 * err("full exchange")
            errors.append(sense + flight + prep + swaps)
        bias = 120 * max(errors)
        reward = measure.signed_readout(setting)
        sampling_width = measure.bound(reward, u, len(parents) + 1)
        # Rounding merit for the altered per-visit means belongs in the path
        # second moment too. Each perturbed reward magnitude is <= |r|+max(errors).
        variance_width = measure.bound(abs(reward) + max(errors), u, len(parents) + 1)
        # |S|<=H/2, so the route perturbation changes E S by <=H*TV and
        # E S^2 by <=H^2*TV/4. Add its root using sqrt(a+b)<=sqrt(a)+sqrt(b).
        deviation_log = log(2 * (len(parents) + 1) / 0.01)
        extra_width = sqrt(
            2
            * (horizon * horizon * total_variation / 4)
            * deviation_log
            / measure.PLAN["releases"]
        )
        overall = bias + route_bias + variance_width + extra_width
        residual = abs(record["mean_per_release"] - record["predicted"])
        assert residual <= overall
        panels.append(
            dict(
                name=setting.name,
                face=setting.face,
                source=setting.source,
                probe=setting.probe,
                exact_source_timing_bias_per_release_upper=bias,
                max_source_routing_pulses_per_visit=per_visit,
                source_record_total_variation_upper=total_variation,
                source_routing_bias_per_release_upper=route_bias,
                timing_bias_per_release_upper=bias + route_bias,
                extra_source_sampling_width=extra_width,
                ideal_resampling_half_width=sampling_width,
                perturbed_resampling_half_width=variance_width,
                combined99_half_width=overall,
                observed_deviation=residual,
            )
        )
    return dict(
        repetitions=REPEATS,
        ordinary_clock_starts_per_bank=SHOTS,
        last_full_exchange_starts_per_bank=2**39,
        grid=GRID,
        material_hands=hands,
        joint_clock_failure_probability_upper=alpha,
        duration_readings=calibration,
        panels=panels,
        largest_timing_bias_per_release=max(
            x["timing_bias_per_release_upper"] for x in panels
        ),
        largest_combined_half_width=max(x["combined99_half_width"] for x in panels),
    )


if __name__ == "__main__":
    HERE.joinpath("clock_results.json").write_text(
        json.dumps(run(), allow_nan=False, indent=2) + "\n"
    )
