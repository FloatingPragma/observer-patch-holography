"""Signed three-corner correlations, complete drain records, and four-point ranges."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import log, prod
from pathlib import Path

import numpy as np

from .. import measure as apparatus

HERE = Path(__file__).parent
SHOTS = 65536
CORNERS = tuple(product(apparatus.CENTRES, repeat=3))
SIGNS = (1, 2, 3)  # numbers of positive-parity coin outcomes out of four


def angular_score(address):
    return prod(2 * x - 11 for x in address)


def source_tags(positive_tags):
    """Sixteen equiprobable settings physically dispatch to the eight corners."""
    addresses = []
    for tag in range(16):
        sx = 1 if tag & 1 else -1
        sy = 1 if tag & 2 else -1
        parity = 1 if tag // 4 < positive_tags else -1
        corner = tuple((11 + v) // 2 for v in (sx, sy, parity * sx * sy))
        addresses.append(apparatus.release_corner(corner))
    return addresses


def active_range(start, end, parked=False):
    """Actual parallel contact pulses after the source flag has fired."""
    _, axes, _ = apparatus.layout()
    address = list(start)
    for ticks in range(12):
        if tuple(address) == tuple(end):
            return ticks
        before = address.copy()
        for i, (a, b) in enumerate(zip(before, end, strict=True)):
            if a != b:
                dest = a + (1 if b > a else -1)
                edge = (axes[i][a], axes[i][dest])
                address[i] = axes[i].index(apparatus.pulse(edge, edge[0], parked))
        if before == address:
            break
    return None


def ranges():
    panels = []
    for k in range(1, 6):
        points = [(5 + k, 5, 5), (5, 5 + k, 5), (5 - k, 5, 5), (5, 5 - k, 5)]
        distances = {
            f"{a}{b}": active_range(points[a], points[b])
            for a, b in combinations(range(4), 2)
        }
        diagonals = distances["02"] * distances["13"]
        sides = distances["01"] * distances["23"] + distances["03"] * distances["12"]
        if (
            Fraction(diagonals, sides) != 2
            or active_range(points[0], points[1], True) is not None
        ):
            raise ArithmeticError(
                "unexpected native four-point range or parked arrival"
            )
        panels.append(
            dict(
                step=k,
                points=points,
                outward_times=distances,
                ptolemy_ratio="2",
                parked_crossing_arrives=False,
            )
        )
    return panels


def one_step(address, tag):
    if tag >= 6:
        return address
    _, axes, _ = apparatus.layout()
    axis = tag // 2
    old = address[axis]
    dest = old + (1 if tag % 2 else -1)
    path = axes[axis]
    received = apparatus.pulse((path[old], path[dest]), path[old])
    result = list(address)
    result[axis] = path.index(received)
    return tuple(result)


def harmonic_check():
    """Check the stopped-score equation at every native interior address."""
    counts = []
    for drain in range(2, 7):
        interior = [
            p for p in product(range(12), repeat=3) if apparatus.layer(p) < drain
        ]
        for p in interior:
            if sum(
                angular_score(one_step(p, tag)) for tag in range(8)
            ) != 8 * angular_score(p):
                raise ArithmeticError("native transition breaks the angular martingale")
        counts.append(
            dict(
                drain=drain,
                interior_addresses=len(interior),
                checked_directed_offers=8 * len(interior),
            )
        )
    return counts


def nearest_exit(positive_tags, seed):
    rng = np.random.default_rng(seed)
    releases = source_tags(positive_tags)
    tag_ids = rng.integers(16, size=SHOTS)
    ids = np.array([CORNERS.index(p) for p in releases])[tag_ids]
    positions = list(CORNERS)
    transitions = np.empty((len(CORNERS), 8), dtype=int)
    for i, corner in enumerate(CORNERS):
        for tag in range(8):
            p = one_step(corner, tag)
            if p not in positions:
                positions.append(p)
            transitions[i, tag] = positions.index(p)
    active = np.arange(SHOTS)
    arrivals = []
    while len(active):
        ids[active] = transitions[ids[active], rng.integers(8, size=len(active))]
        finished = ids[active] >= len(CORNERS)
        arrivals.append(
            np.bincount(ids[active[finished]], minlength=len(positions))[8:].tolist()
        )
        active = active[~finished]
    s = Fraction(positive_tags, 2) - 1
    score = np.array([angular_score(x) for x in positions[8:]])
    banks = np.array(arrivals, dtype=np.int64)
    if not np.all(abs(score) == 3) or int(banks.sum()) != SHOTS:
        raise ArithmeticError("lost a launch or observed an unexpected final address")
    counts = Counter(releases)
    marginal = {
        str(axes): sorted(Counter(tuple(p[k] for k in axes) for p in releases).values())
        for axes in combinations(range(3), 2)
    }
    return dict(
        seed=seed,
        shots=SHOTS,
        s=str(s),
        source_counts_per_16=[counts[p] for p in CORNERS],
        each_two_register_marginal_per_16=marginal,
        observed_launch_counts=np.bincount(tag_ids, minlength=16).tolist(),
        collectors=positions[8:],
        exit_counts_by_tick_and_collector=arrivals,
        predicted_positive_score=str((1 + s / 3) / 2),
        positive_score_count=int(banks[:, score > 0].sum()),
        observed_mean_score=float(banks.sum(axis=0) @ score / SHOTS),
        predicted_mean_score=str(s),
        predicted_conditional_first_two_means=[str(3 * s), str(-3 * s / 5)],
    )


def reference_modular_meter(banks):
    # Sample the reference on a fresh bank and read the three actual addresses.
    seed = 87371
    reference = source_tags(1)
    indices = np.array([CORNERS.index(p) for p in reference])
    outcomes = indices[np.random.default_rng(seed).integers(len(indices), size=SHOTS)]
    counts = np.bincount(outcomes, minlength=len(CORNERS))
    signs = np.array([angular_score(p) for p in CORNERS])
    nplus = int(counts[signs > 0].sum())
    if not 0 < nplus < SHOTS:
        raise ArithmeticError("missing a modular reference eigenspace")
    coefficient = log((SHOTS - nplus) / nplus) / 2
    first, last = banks[0], banks[-1]
    score_change = last["observed_mean_score"] - first["observed_mean_score"]
    return dict(
        reference_tag_seed=seed,
        reference_s="-1/2",
        reference_shots=SHOTS,
        reference_counts_by_corner=counts.tolist(),
        predicted_K0_parity_coefficient=log(3) / 2,
        observed_K0_parity_coefficient=coefficient,
        held_initial_s=first["s"],
        held_final_s=last["s"],
        observed_terminal_score_change=score_change,
        observed_fixed_reference_modular_change_nats=coefficient * score_change,
        predicted_fixed_reference_modular_change_nats=log(3) / 2,
        predicted_held_entropy_change_nats=0,
    )


def population_clock():
    # These complete native release-to-target times were stored independently
    # of either launch bias. Include the dispatch appointment from READY.
    archive = json.loads((apparatus.HERE / "records.json").read_text())
    clocks = np.asarray(archive["native_range_reads"], dtype=np.int64)
    if clocks.shape != (8, 12**3):
        raise ArithmeticError("incomplete source-to-target survey")
    signs = np.array([angular_score(p) for p in CORNERS])
    xsign = np.array([2 * p[0] - 11 for p in CORNERS])
    corr = signs @ clocks
    bias = xsign @ clocks
    shells = {}
    resolved = []
    diagonal_slopes = {}
    for i, target in enumerate(product(range(12), repeat=3)):
        offsets = [min(abs(x - 5), abs(x - 6)) for x in target]
        axes = [j for j, d in enumerate(offsets) if d == max(offsets)]
        near = [1 if x >= 6 else -1 for x in target]
        expected_s = -prod(near) if len(axes) == 3 else 0
        expected_b = -near[0] * (8 // 2 ** len(axes)) if 0 in axes else 0
        if (corr[i], bias[i]) != (expected_s, expected_b):
            raise ArithmeticError("source-weighted native clock predictions differ")
        shell = 1 + max(offsets)
        total = shells.setdefault(shell, [0, 0, 0])
        total[0] += 1
        total[1] += int(corr[i])
        total[2] += int(bias[i])
        if corr[i]:
            resolved.append([*target, int(corr[i])])
            diagonal_slopes.setdefault(shell, {})[tuple(near)] = (
                int(corr[i]),
                int(bias[i]),
            )
    parity_tests = []
    monomials = [axes for n in (0, 1, 2) for axes in combinations(range(3), n)]
    for shell, slopes in sorted(diagonal_slopes.items()):
        if len(slopes) != 8:
            raise ArithmeticError("incomplete diagonal direction panel")
        tensor_parities = [
            sum(prod(t) * prod(t[i] for i in axes) for t in slopes)
            for axes in monomials
        ]
        if any(tensor_parities):
            raise ArithmeticError("quadratic directional response has cubic parity")
        corr_contrast = Fraction(sum(prod(t) * s for t, (s, _) in slopes.items()), 64)
        if corr_contrast != Fraction(-1, 8) or any(
            b != -t[0] or s != -prod(t) for t, (s, b) in slopes.items()
        ):
            raise ArithmeticError("clock tensor contrasts disagree with native reads")
        parity_tests.append(
            dict(
                shell=shell,
                targets=len(slopes),
                parity_sums_for_tensor_basis=tensor_parities,
                correlation_parity_projection=str(corr_contrast),
                best_uniform_tensor_error=str(abs(corr_contrast)),
                ordinary_dipole_uniform_error="0",
            )
        )
    return dict(
        checked_targets=clocks.shape[1],
        shell_totals_targets_and_eight_times_s_and_b_slopes=[
            shells[k] for k in sorted(shells)
        ],
        targets_and_eight_times_nonzero_correlation_slope=resolved,
        eight_diagonal_tensor_tests=parity_tests,
    )


def energy_census():
    board, axes, flag = apparatus.layout()
    hop, contrast = apparatus.candidate.native_law.local_term_parts()
    h = hop + apparatus.candidate.KAPPA * contrast
    src = apparatus.candidate.native_law.SRC
    for left, right in product((0, 2), repeat=2):
        i, j = src.index((left, right)), src.index((right, left))
        expected = np.zeros(len(h))
        if left != right:
            expected[i], expected[j] = 32.5, -0.5
        np.testing.assert_allclose(h[:, i], expected, rtol=0, atol=2e-13)
    # A native half-exposure measures the occupation of guide index five on
    # fresh releases. The new index-four count gives the seam mean without
    # reading the tag used in preparation.
    start, end = src.index((0, 2)), src.index((2, 0))
    pulse = apparatus.expm(-0.5j * np.pi * h)
    parked = apparatus.expm(-0.5j * np.pi * apparatus.candidate.KAPPA * contrast)
    np.testing.assert_allclose(
        [abs(pulse[end, start]) ** 2, abs(parked[end, start]) ** 2],
        [0.5, 0.0],
        atol=1e-13,
        rtol=0,
    )
    occupied = [
        {axis[x] for axis, x in zip(axes, p, strict=True)} | {flag[1]} for p in CORNERS
    ]
    # The ready flag has fired on every branch. Only the four sharp doubles
    # contribute. Unequal endpoints have equal odd/even spectral weights.
    unlike = np.array(
        [[(u in sites) != (v in sites) for sites in occupied] for u, v in board.seams],
        dtype=np.int64,
    )
    ordinary_edge = tuple(sorted((axes[0][4], axes[0][5])))
    edge_number = next(
        i for i, e in enumerate(board.seams) if tuple(sorted(e)) == ordinary_edge
    )
    arms = []
    for s, b in (
        (Fraction(-1, 2), 0),
        (0, 0),
        (Fraction(1, 2), 0),
        (Fraction(-1, 2), Fraction(-1, 4)),
        (Fraction(-1, 2), Fraction(1, 4)),
    ):
        weights = np.array(
            [int(4 * (1 + s * angular_score(p) + b * (2 * p[0] - 11))) for p in CORNERS]
        )
        if min(weights) <= 0 or sum(weights) != 32:
            raise ArithmeticError("invalid full-support ordinary calibration")
        signed = unlike @ weights
        seed = 173015 + len(arms)
        rng = np.random.default_rng(seed)
        branches = np.repeat(np.arange(8), weights)[rng.integers(32, size=SHOTS)]
        split = rng.integers(2, size=SHOTS)
        at_five = np.array([p[0] == 5 for p in CORNERS])[branches]
        arrived = int(np.sum(at_five & (split == 1)))
        arms.append(
            dict(
                s=str(s),
                b=str(b),
                weights_per_32=weights.tolist(),
                prob_each_energy_32_and_33_per_64=signed.tolist(),
                ordinary_edge_mean=str(Fraction(65 * int(signed[edge_number]), 64)),
                host_sum_mean=str(Fraction(65 * int(signed.sum()), 64)),
                half_pulse_meter=dict(
                    seed=seed,
                    launches=SHOTS,
                    arrivals=arrived,
                    refusals=SHOTS - arrived,
                    parked_arrivals=0,
                    predicted_arrival=str((1 - Fraction(b)) / 4),
                ),
            )
        )
    return dict(
        host_seams=len(board.seams),
        ordinary_edge=ordinary_edge,
        native_unlike_seams_per_corner=unlike.sum(axis=0).tolist(),
        arms=arms,
    )


def run():
    banks = [nearest_exit(a, 92753 + a) for a in SIGNS]
    return dict(
        ranges=ranges(),
        harmonic_rows=harmonic_check(),
        angular_banks=banks,
        reference_modular=reference_modular_meter(banks),
        released_energy=energy_census(),
        source_weighted_clock=population_clock(),
        released_entropy_gain_from_minus_half_to_zero_nats=3 * log(3) / 4 - log(2),
        entropy_slope_at_minus_half_nats=log(3) / 2,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), separators=(",", ":"), allow_nan=False) + "\n"
    path = HERE / "records.json"
    if args.write:
        path.write_text(out)
    elif json.loads(out) != json.loads(path.read_text()):
        raise SystemExit("three-guide native source records differ")
    else:
        print(
            "Replayed all corner sources, clocks, stopped-score checks and terminal records."
        )
