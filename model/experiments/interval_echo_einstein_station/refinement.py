"""Finite-screen transfer of the occupied-source material-clock radar apparatus."""

import argparse
import json
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from scipy.special import jv

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.contact import radar
from . import volume

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / "refinement.json"
PAIRS = ((0, 1), (2, 3))
ARMS = ("vacant_reference", "occupied_reference", "occupied_minus", "occupied_plus")


def station(level, source, tag):
    board = federation(level)
    a, b, c = source["source_ports"]
    hand = tuple(source["hand_seam"])
    fixed = {b, *tag}
    marked = {(a, b), (b, c), hand}
    seams = set(board.seams)
    if not ({tuple(e) for e in cell().seams} == {e for e in seams if max(e) < 12}):
        raise ValueError("Home cell incidence changed")
    if not (marked | {tag} <= seams):
        raise ValueError("Source, tag or hand seam is missing")
    idle = {e for e in board.seams if fixed.intersection(e) and e not in marked}
    gateways = sorted(e for e in board.external if e[0] < 12)
    if not ([a for a, _ in gateways] == [1, 3, 10]):
        raise ValueError("Home gateway arrangement changed")
    (bridge,) = [e for e in gateways if e not in idle]
    remote = bridge[1]
    remote_gateways = sorted(
        (a % 12, b % 12)
        for x, y in board.external
        for a, b in ((x, y), (y, x))
        if a // 12 == remote // 12
    )
    if not (
        bridge[0] == 3
        and remote % 12 == 0
        and remote_gateways == [(0, 3), (2, 1), (8, 10)]
    ):
        raise ValueError("Bridge or neighbouring gateway arrangement changed")
    sites = [7, 8, remote, remote + 5]
    if (remote, remote + 5) not in seams:
        raise ValueError("Remote control seam is missing")
    # The two metering loads can still reach the same hand without disturbing
    # the source, tag or parked probe, in either direction on both dyads.
    for origin, (p, q) in zip(sites, [(4, 8), (4, 9), (4, 8), (4, 8)], strict=True):
        if (min(p, hand[0]), max(p, hand[0])) not in seams:
            raise ValueError("First hand rail is unreachable")
        if (min(q, hand[1]), max(q, hand[1])) not in seams:
            raise ValueError("Second hand rail is unreachable")
        if not (not {p, q}.intersection({a, b, c, *hand, *tag, origin})):
            raise ValueError("Meter spares overlap a protected register")

    law = candidate.native_law
    full = sum(
        weight * term
        for weight, term in zip(
            (1, candidate.KAPPA), law.local_term_parts(), strict=True
        )
    )
    words = [law.SRC.index(w) for w in ((0, 0), (1, 0), (0, 1), (1, 1))]
    blank = full[np.ix_(words[1:3], words[1:3])]
    energies = source["source_word_energies"]
    for e, mark in zip(
        energies,
        (volume.label(0, 1), volume.label(2, 1), volume.label(2, -1)),
        strict=True,
    ):
        cols = np.kron(mark[:, None], np.eye(law.N_SRC)[:, words])
        np.testing.assert_allclose(full @ cols, cols * [e, e, e, 0], atol=3e-13, rtol=0)
    base = candidate.generator(board).tolil()
    for edge in sorted(marked | idle):
        for i in (0, 1):
            for j in (0, 1):
                base[edge[i], edge[j]] -= blank[i, j]
    base = base.tocsr()
    free = sorted(set(range(board.n)) - fixed)
    coupled = base[free][:, sorted(fixed)]
    if not (max(abs(coupled.data), default=0) < 1e-14):
        raise ValueError("A protected site is coupled to a free flight")
    return (
        base,
        sites,
        dict(
            level=level,
            cells=board.cells,
            ports=board.n,
            home_gateways=gateways,
            serviced_bridge=bridge,
            remote_gateway_local_ends=remote_gateways,
            sites=sites,
        ),
    )


def sparse_hit_laws(h, sites, pairs=PAIRS, *, return_intensity=False):
    """Unmeasured basis flights for all native appointments, then renewal."""
    count = radar.TICKS
    initial = np.zeros((h.shape[0], len(sites)))
    initial[sites, np.arange(len(sites))] = 1
    flight = expm_multiply(
        -1j * h,
        initial,
        start=0,
        stop=count * radar.DT,
        num=count + 1,
        traceA=-1j * h.diagonal().sum(),
    )[:, sites, :]
    laws = np.zeros((len(sites), len(sites), count + 1))
    for left, right in pairs:
        for start, end in ((left, right), (right, left)):
            hits = np.zeros(count + 1, dtype=complex)
            for k in range(1, count + 1):
                hits[k] = (
                    flight[k, end, start] - flight[1:k, end, end] @ hits[k - 1 : 0 : -1]
                )
            laws[start, end, :-1] = abs(hits[1:]) ** 2
            laws[start, end, -1] = 1 - laws[start, end, :-1].sum()
    if not (np.all(laws >= 0)):
        raise ValueError("Negative first-hit probability")
    if return_intensity:
        # Fresh launch at each read time, with no earlier detector appointment.
        return laws, np.abs(flight[1:].transpose(2, 1, 0)) ** 2
    return laws


def projector_check(h, sites, pairs=PAIRS):
    """Independent half-unit Chebyshev evolution and actual first-hit removals."""
    diagonal = h.diagonal()
    off = np.asarray(abs(h).sum(axis=1)).ravel() - abs(diagonal)
    lo, hi = min(diagonal - off), max(diagonal + off)
    center, radius = (lo + hi) / 2, (hi - lo) / 2
    scaled = (h - center * sparse.eye(h.shape[0], format="csr")) / radius
    x = radar.DT * radius
    order = int(np.ceil(x))

    # |J_k(x)| <= (x/2)^k/k!. Uniform operator error in exact arithmetic.
    def tail(k):
        return 2 * (x / 2) ** (k + 1) / factorial(k + 1) / (1 - x / (2 * (k + 2)))

    while tail(order) > 1e-18:
        order += 1
    orders = np.arange(order + 1)
    coeff = (2 * (-1j) ** orders) * jv(orders, x)
    coeff[0] /= 2
    starts = [i for pair in pairs for i in pair]
    ends = [i for pair in pairs for i in reversed(pair)]
    state = np.zeros((h.shape[0], len(starts)), dtype=complex)
    columns = np.arange(len(starts))
    state[np.array(sites)[starts], columns] = 1
    targets = np.array(sites)[ends]
    laws = np.zeros((len(sites), len(sites), radar.TICKS + 1))
    for n in range(radar.TICKS):
        prev, curr = state, scaled @ state
        result = coeff[0] * prev + coeff[1] * curr
        for c in coeff[2:]:
            prev, curr = curr, 2 * (scaled @ curr) - prev
            result += c * curr
        state = result
        laws[starts, ends, n] = abs(state[targets, columns]) ** 2
        state[targets, columns] = 0
    laws[starts, ends, -1] = (abs(state) ** 2).sum(axis=0)
    return laws, tail(order)


def summarize(joint, half_hand):
    rows = []
    for left, right in PAIRS:
        one = joint[left, right, :-1, :-1]
        two = joint[right, left, :-1, :-1]
        np.testing.assert_allclose(one, two.T, rtol=0, atol=5e-15)
        returned = (one + two) / 2
        success = returned.sum()
        rows.append(
            dict(
                return_probability=float(success),
                conditional_half_hand=float((returned * half_hand).sum() / success),
            )
        )
    return rows


def clocked_ranges(joint, hand_yes, gap, pairs, rng):
    """Fresh return-count enrollment preserves both timeout records at rare dyads."""
    returns = 2**21  # per direction and quadrature, set before drawing
    bins = np.add.outer(np.arange(radar.TICKS), np.arange(radar.TICKS)).ravel()
    durations = len(hand_yes[0])
    readings = np.zeros((2, durations, 2), dtype=np.int64)
    counts = np.zeros((len(pairs), durations), dtype=np.int64)
    enrolled, timeouts, return_rates = [], [], []
    for k, (i, j) in enumerate(pairs):
        # For each persistent word a completed trip has weight A(m)B(n).
        # The reverse has B(n)A(m). The identity survives mixing the words
        # and ensures equal success rates and duration laws in both directions.
        np.testing.assert_allclose(
            joint[i, j, :-1, :-1], joint[j, i, :-1, :-1].T, atol=0, rtol=3e-15
        )
        return_rates.append(float(joint[i, j, :-1, :-1].sum()))
        timeouts.append([])
        for start, end in ((i, j), (j, i)):
            law = joint[start, end]
            ok = np.bincount(bins, weights=law[:-1, :-1].ravel(), minlength=durations)
            success = ok.sum()
            quiet = law[-1, -1] + law[:-1, -1].sum()
            np.testing.assert_allclose(success + quiet, 1, atol=1e-14, rtol=0)
            # Keeping relative precision here is essential on remote dyads.
            # Rounding this enrollment rate affects quiet counts only. Durations
            # have their own conditional lattice, independent of the stop count.
            rate = float(f"{success:.9g}")
            outward_quiet = float(f"{law[-1, -1] / quiet:.12g}")
            duration_law = radar.acquisition_weights(ok / success) / radar.LATTICE
            for axis, yes in enumerate(hand_yes):
                missed = int(rng.negative_binomial(returns, rate))
                outward = int(rng.binomial(missed, outward_quiet))
                timeouts[k].extend((outward, missed - outward))
                enrolled.append(returns + missed)
                by_duration = rng.multinomial(returns, duration_law)
                yes = np.round(yes, 12)
                seen = rng.binomial(by_duration, yes)
                readings[axis, :, 0] += seen
                readings[axis, :, 1] += by_duration - seen
                counts[k] += by_duration
    exposures = readings.sum(axis=2)
    if np.any(exposures == 0):
        raise ValueError("A duration lacks independent material-hand reads")
    # Union bound for four arms, two clock axes per bin and each range.
    log_factor = np.log(2 * 4 * (2 * durations + len(pairs)) / 0.01)
    xy = 2 * readings[:, :, 0] / exposures - 1
    phase = np.unwrap(np.arctan2(xy[1], xy[0]))
    phase_disk = np.linalg.norm(np.sqrt(2 * log_factor / exposures) + 1e-12, axis=0)
    if np.max(phase_disk) >= 1:
        raise ValueError("Clock phases are unresolved")
    phase_radius = np.arcsin(phase_disk)
    if abs(phase[0]) + phase_radius[0] >= np.pi or np.any(
        abs(np.diff(phase)) + phase_radius[:-1] + phase_radius[1:] >= np.pi
    ):
        raise ValueError("Clock winding is unresolved")
    np.testing.assert_array_equal(counts.sum(axis=1), 4 * returns)
    np.testing.assert_array_equal(exposures.sum(axis=0), counts.sum(axis=0))
    if 4 * len(pairs) * returns + sum(map(sum, timeouts)) != sum(enrolled):
        raise ValueError(
            "Retained and quiet terminal records do not exhaust enrollments"
        )
    proportions = counts / (4 * returns)
    distance = proportions @ phase / (2 * gap)
    width = radar.DT * (radar.TICKS - 1)
    radius = proportions @ phase_radius / (2 * gap) + width * (
        np.sqrt(log_factor / (8 * returns)) + durations / radar.LATTICE
    )

    # Decimal rows keep the receipt readable and huge enrollment counts portable.
    def decimal_rows(rows):
        return [",".join(map(str, row)) for row in rows]

    return (
        distance,
        radius,
        dict(
            returns_per_direction_and_axis=returns,
            enrollments_by_menu_min_max=[str(min(enrolled)), str(max(enrolled))],
            outward_timeouts=str(sum(sum(row[::2]) for row in timeouts)),
            return_timeouts=str(sum(sum(row[1::2]) for row in timeouts)),
            native_round_trip_probability_by_dyad=return_rates,
            timeout_pairs_by_dyad_and_menu=decimal_rows(timeouts),
            returned_duration_counts_by_dyad=decimal_rows(counts),
            hand_exposures_by_axis_duration=decimal_rows(exposures),
            hand_yes_by_axis_duration=decimal_rows(readings[:, :, 0]),
        ),
    )


def space_form_gram(u, x):
    """Five-point tangent Gram, u = K L² and x = measured distance / L."""
    r, ij = x[0, 1:], x[1:, 1:]
    if u == 0:
        return (r[:, None] ** 2 + r[None, :] ** 2 - ij**2) / 2
    c = np.cos if u > 0 else np.cosh
    scale = np.sqrt(abs(u))
    return (c(scale * ij) - np.outer(c(scale * r), c(scale * r))) / u


def angular_comparison_bound(x, errors):
    """Unit tangent Gram interval for |sectional curvature| L² ≤ 1.

    Dyer, Vegter & Wintraecken, Theorem 33, doi:10.1007/s10711-015-0069-5.
    Assumes minimizing triangles in a common convex comparison neighborhood.
    All three lengths are below L, ensuring the positive perimeter domain.
    """
    low, high = x - errors, x + errors
    a, b = low[0, 1:], high[0, 1:]
    if min(a) <= 0 or high.max() > 1:
        raise ValueError("Hinge distances exceed the enrolled comparison domain")
    # Spherical comparison gives the lowest possible true cosine, hyperbolic
    # the highest. Both positive denominator endpoints are used for either sign.
    lower_numerator = np.cos(high[1:, 1:]) - np.outer(np.cos(a), np.cos(a))
    upper_numerator = np.outer(np.cosh(b), np.cosh(b)) - np.cosh(low[1:, 1:])
    lower = np.minimum(
        *[lower_numerator / np.outer(np.sin(r), np.sin(r)) for r in (a, b)]
    )
    upper = np.maximum(
        *[upper_numerator / np.outer(np.sinh(r), np.sinh(r)) for r in (a, b)]
    )
    np.fill_diagonal(lower, 1)
    np.fill_diagonal(upper, 1)
    if np.any(lower > upper) or np.any(lower > 1) or np.any(upper < -1):
        raise ValueError("Geodesic comparison angles have an empty interval")
    lower, upper = np.maximum(lower, -1), np.minimum(upper, 1)
    midpoint, radius = (lower + upper) / 2, (upper - lower) / 2
    # Symmetric perturbation: the spectral norm is bounded by the max row sum.
    return np.linalg.eigvalsh(midpoint)[0] - np.linalg.norm(radius, np.inf)


def candidate_geometry():
    """Measure clocked ranges among six stations on the physically matched two-write station."""
    from scipy.optimize import brentq
    from . import shared_variance as shared

    sites = [4, 9, 11, 48, 49, 53]
    base, old_sites, hand_yes, source = shared.shared_station(extra_sites=(11, 49))
    if old_sites != [4, 9, 48, 53]:
        raise ValueError("Enrolled six-site geometry placement changed")
    pairs = list(combinations(range(len(sites)), 2))
    subsets = list(combinations(range(len(sites)), 5))
    weights, _ = shared.source_words(source)
    laws, independent, tail = [], [], 0.0
    for h in shared.word_generators(base, source["source_word_energies"]):
        check, bound = projector_check(h, sites, pairs)
        laws.append(radar.round_trip_law(check))
        independent.append(radar.round_trip_law(sparse_hit_laws(h, sites, pairs)))
        tail = max(tail, bound)
    joint = np.asarray(laws)
    mixed = np.einsum("sw,wijab->sijab", weights.reshape(3, 81), joint)
    independent = np.asarray(independent)
    alternate = np.einsum("sw,wijab->sijab", weights.reshape(3, 81), independent)
    length = radar.DT * radar.TICKS
    grid = np.linspace(-1, 1, 1001)
    # For x<=1, cosine product-to-sum bounds an entry's curvature derivative by
    # sum_{n>=2}(n-1)(1+4**n)/(2n)! < 2. Its three distance derivatives sum
    # to at most sinh(1)+2*sinh(1)*cosh(1)<5. Row sums bound the 4x4 norm.
    between_grid = 8 * (grid[1] - grid[0]) / 2
    rng = np.random.default_rng(2026112261)
    rows = []
    for name, law, check in zip(
        ("vacant_reference", *shared.ARMS),
        (joint[0], *mixed),
        (independent[0], *alternate),
        strict=True,
    ):
        expected, total = radar.distance_statistics(law, radar.DT)
        check_distance, _ = radar.distance_statistics(check, radar.DT)
        disagreement = np.max(abs(expected - check_distance))
        if disagreement > 1e-9:
            raise ValueError("Independent six-station propagations disagree")
        # Full two-leg record mass error from telescoping projected Chebyshev steps.
        floor = np.min(total[np.triu_indices(len(sites), 1)]) / 2
        mass_error = 4 * radar.TICKS * tail * np.exp(4 * radar.TICKS * tail)
        if floor <= mass_error:
            raise ValueError(
                "Rare-return conditioning exceeds the propagated numerical budget"
            )
        truncation = (radar.TICKS - 1) * radar.DT * mass_error / (floor - mass_error)
        d, errors, counts = clocked_ranges(
            law, hand_yes, source["hand_gap"], pairs, rng
        )
        errors += truncation
        if np.any(abs(d - expected[np.triu_indices(len(sites), 1)]) > errors):
            raise ValueError(
                "Ideal distances fall outside the acquired simultaneous interval"
            )
        if max(d + errors) >= length or min(d - errors) < 0:
            raise ValueError(
                "Distance confidence intervals exceed the derivative domain"
            )
        distance = np.zeros((len(sites), len(sites)))
        ranges_error = np.zeros_like(distance)
        for (i, j), value, error in zip(pairs, d, errors, strict=True):
            ranges_error[i, j] = ranges_error[j, i] = error / length
            distance[i, j] = distance[j, i] = value / length
        floors, angular_floors, roots = [], [], []
        for subset in subsets:
            selected = np.ix_(subset, subset)
            sub = distance[selected]
            angular_floors.append(angular_comparison_bound(sub, ranges_error[selected]))

            def least_eigenvalue(u):
                return np.linalg.eigvalsh(space_form_gram(u, sub))[0]

            floors.append(
                min(least_eigenvalue(u) for u in grid)
                - between_grid
                - 20 * max(errors) / length
            )
            # This sign-changing positive branch is an exploratory fit. Each
            # five-station calibration is separate, with no extra fit to site 6.
            roots.append(
                brentq(least_eigenvalue, 1, np.pi**2 / sub.max() ** 2) / length**2
            )
        rows.append(
            dict(
                name=name,
                material_ranges=d.tolist(),
                simultaneous_99pct_radius=errors.tolist(),
                weakest_eigenvalue_lower_bound=min(floors),
                weakest_variable_curvature_angular_bound=min(angular_floors),
                independent_range_max_disagreement=disagreement,
                conditional_truncation_bound_in_exact_arithmetic=truncation,
                isolated_spherical_roots=roots,
                terminal_counts=counts,
            )
        )
    # Analytic positive/flat/negative controls in 3D built independently from
    # tangent vectors and the standard S^3 / hyperboloid embeddings.
    v = np.array(
        [[0, 0, 0], [0.3, 0, 0], [0, 0.32, 0], [0, 0, 0.27], [0.23, 0.25, 0.1]]
    )
    norms = np.linalg.norm(v, axis=1)
    controls = []
    angular_controls = []
    for u in (-0.5, 0, 0.5):
        if u == 0:
            x = np.linalg.norm(v[:, None] - v[None, :], axis=2)
        else:
            a = np.sqrt(abs(u))
            c, s = (np.cos, np.sin) if u > 0 else (np.cosh, np.sinh)
            p = s(a * norms)[:, None] * v / np.where(norms == 0, 1, norms)[:, None]
            dots = np.outer(c(a * norms), c(a * norms)) + np.sign(u) * (p @ p.T)
            x = (
                np.arccos(np.clip(dots, -1, 1))
                if u > 0
                else np.arccosh(np.maximum(1, dots))
            ) / a
            np.fill_diagonal(x, 0)
        controls.append(float(np.linalg.eigvalsh(space_form_gram(u, x))[0]))
        angular_controls.append(angular_comparison_bound(x, np.zeros_like(x)))
    np.testing.assert_allclose(controls, 0, atol=1e-13, rtol=0)
    if max(angular_controls) > 1e-13:
        raise ValueError("Comparison test rejected a known 3D geometry")
    return dict(
        sites=sites,
        dyads=[[sites[i], sites[j]] for i, j in pairs],
        five_site_subsets=[[sites[i] for i in idx] for idx in subsets],
        length_unit=length,
        dimensionless_curvature_window=[-1, 1],
        controls_at_dimensionless_curvatures=[-0.5, 0, 0.5],
        control_kernel=controls,
        angular_comparison_control_lower_bounds=angular_controls,
        between_grid_eigenvalue_allowance=between_grid,
        chebyshev_half_step_truncation_in_exact_arithmetic=tail,
        maximum_each_direction_record_mass_error_in_exact_arithmetic=mass_error,
        count_layout="Decimal count rows are comma-delimited. Dyads follow dyads above. Each timeout row lists outward then return failures for the four menus in forward-cosine, forward-sine, reverse-cosine, reverse-sine order. Every menu stops at the fixed number of completed returns. Each duration row spans 2 through 24 half-unit appointments. The two hand rows are cosine and sine.",
        seed=2026112261,
        arms=rows,
        scope="Conditional material-clock half-trips of a stationary six-port apparatus. The constant-curvature Gram and the variable-curvature angle enclosure both reject these ranges as geodesic distances on a three-dimensional Riemannian comparison neighborhood with |sectional curvature| L² ≤ 1. The angle test needs minimizing triangles within the convexity radius. Isolated positive roots are exploratory five-site fits. Other radar instruments, coarse graining and spacetime Ricci response remain untested.",
    )


def run():
    coarse = json.loads((HERE / "volume.json").read_text())["occupied_source_radar"]
    _, _, hand_yes, source = volume.radar_source_and_hand(federation(0))
    tag = tuple(coarse["protected_tag_seam"])
    weights = np.array(source["source_units"]).reshape(3, 9) / 36
    energies = source["source_word_energies"]
    durations = np.arange(2, 2 * radar.TICKS + 1)
    clock = (
        np.unwrap(np.arctan2(2 * hand_yes[1] - 1, 2 * hand_yes[0] - 1))
        / source["hand_gap"]
    )
    np.testing.assert_allclose(clock, durations * radar.DT, atol=5e-13, rtol=0)
    half_hand = clock[np.add.outer(np.arange(radar.TICKS), np.arange(radar.TICKS))] / 2
    levels = []
    for level in (0, 1, 2):
        base, sites, row = station(level, source, tag)
        laws, checked = [], []
        check_error = dense_error = tail_bound = 0.0
        for first in energies:
            for second in energies:
                diag = np.zeros(base.shape[0])
                a, _, c = source["source_ports"]
                diag[[a, c]] = [-first, -second]
                h = base + sparse.diags(diag)
                renewed = sparse_hit_laws(h, sites)
                direct, tail = projector_check(h, sites)
                tail_bound = max(tail_bound, tail)
                check_error = max(check_error, abs(renewed - direct).max())
                if level == 0:
                    dense = radar.first_hit(
                        radar._unitary(h.toarray(), radar.DT), sites, radar.TICKS
                    )
                    for i, j in PAIRS:
                        dense_error = max(
                            dense_error,
                            abs(renewed[i, j] - dense[i, j]).max(),
                            abs(renewed[j, i] - dense[j, i]).max(),
                        )
                laws.append(radar.round_trip_law(renewed))
                checked.append(radar.round_trip_law(direct))
        if not (max(check_error, dense_error) < 1e-11):
            raise ValueError("Independent propagation algorithms disagree")

        def arms(word_laws):
            mixture = np.einsum("sw,wijab->sijab", weights, word_laws)
            if not (abs(mixture[1] - (mixture[0] + mixture[2]) / 2).max() < 3e-15):
                raise ValueError("Source mixture fails its midpoint identity")
            return [word_laws[0], mixture[1], mixture[0], mixture[2]]

        values = {
            name: summarize(joint, half_hand)
            for name, joint in zip(ARMS, arms(laws), strict=True)
        }
        controls = [summarize(joint, half_hand) for joint in arms(checked)]
        distance_error = max(
            abs(row["conditional_half_hand"] - ref["conditional_half_hand"])
            for measured, check in zip(values.values(), controls, strict=True)
            for row, ref in zip(measured, check, strict=True)
        )
        # Embed every stopping record in its own orthogonal branch. Each of
        # the at most 2*T clock steps has exact norm one, and its Chebyshev
        # replacement is within e. Telescoping the joint record amplitudes
        # and taking the squared norms bounds total L1 probability error by
        # 4*T*e*exp(4*T*e). Classical source mixing preserves this bound.
        steps = 2 * radar.TICKS
        mass_error = 2 * steps * tail_bound * np.exp(2 * steps * tail_bound)
        floor = min(r["return_probability"] for arm in controls for r in arm)
        truncation_in_mean = (
            (radar.TICKS - 1) * radar.DT * mass_error / (floor - mass_error)
        )
        if level == 0:
            for name, acq in zip(ARMS, coarse["acquisitions"], strict=True):
                if not (name == acq["name"]):
                    raise ValueError("Coarse calibration has the wrong preparation")
                for dyad, result in zip(([7, 8], [48, 53]), values[name], strict=True):
                    k = coarse["dyads"].index(dyad)
                    np.testing.assert_allclose(
                        result["conditional_half_hand"],
                        acq["ideal_distance"][k],
                        rtol=0,
                        atol=1e-11,
                    )
        contrasts = []
        for k, (left, right) in enumerate(PAIRS):

            def diff(key, baseline, probe):
                return values[probe][k][key] - values[baseline][k][key]

            contrasts.append(
                dict(
                    sites=[sites[left], sites[right]],
                    added_unit_half_hand=diff("conditional_half_hand", *ARMS[:2]),
                    signed_half_hand=diff("conditional_half_hand", *ARMS[2:]),
                    added_unit_return_probability=diff("return_probability", *ARMS[:2]),
                    signed_return_probability=diff("return_probability", *ARMS[2:]),
                )
            )
        row.update(
            arms=values,
            contrasts=contrasts,
            maximum_leg_absolute_disagreement=float(check_error),
            maximum_mean_absolute_disagreement=float(distance_error),
            dense_coarse_leg_disagreement=float(dense_error) if level == 0 else None,
            chebyshev_half_step_uniform_truncation=float(tail_bound),
            conditional_mean_truncation_bound_in_exact_arithmetic=float(
                truncation_in_mean
            ),
        )
        if levels:
            previous = levels[-1]["arms"]
            row["largest_absolute_mean_change_from_preceding_screen_by_dyad"] = [
                max(
                    abs(
                        values[name][k]["conditional_half_hand"]
                        - previous[name][k]["conditional_half_hand"]
                    )
                    for name in ARMS
                )
                for k in range(len(PAIRS))
            ]

        levels.append(row)
    return dict(
        source_ports=source["source_ports"],
        hand_seam=source["hand_seam"],
        protected_tag_seam=tag,
        energies=energies,
        source_units=source["source_units"],
        hand_gap=source["hand_gap"],
        period=radar.DT,
        appointments=radar.TICKS,
        method="Full native sparse flights and renewal per persistent word, independently checked by sparse projected Chebyshev evolutions. Ensemble material phase from both local quadratures; all timeout mass retained.",
        scope="Three supplied screens, one fixed native home, unchanged local address/clock and fixed microscopic service duration. The remote dyad follows the actual surviving bridge. Physical scale matching, proper volume and observer-independent curvature require other measurements.",
        levels=levels,
        six_station_radar_geometry=candidate_geometry(),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = json.loads(json.dumps(run(), allow_nan=False))
    if args.write:
        RECEIPT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        recorded = json.loads(RECEIPT.read_text())
        volume.cap.compare(result, recorded)
        for row, expected in zip(result["levels"], recorded["levels"], strict=True):
            for measured, target in zip(
                row["contrasts"], expected["contrasts"], strict=True
            ):
                for key in (
                    "signed_half_hand",
                    "added_unit_half_hand",
                    "signed_return_probability",
                    "added_unit_return_probability",
                ):
                    np.testing.assert_allclose(
                        measured[key], target[key], rtol=1e-8, atol=1e-13
                    )
    for row in result["levels"]:
        print(row["cells"], row["contrasts"])
