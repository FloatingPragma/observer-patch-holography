"""Finite-screen transfer of the occupied-source material-clock radar apparatus."""

import argparse
import json
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
    if not ((remote, remote + 5) in seams):
        raise ValueError("Remote control seam is missing")
    # The two metering loads can still reach the same hand without disturbing
    # the source, tag or parked probe, in either direction on both dyads.
    for origin, (p, q) in zip(sites, [(4, 8), (4, 9), (4, 8), (4, 8)], strict=True):
        if not ((min(p, hand[0]), max(p, hand[0])) in seams):
            raise ValueError("First hand rail is unreachable")
        if not ((min(q, hand[1]), max(q, hand[1])) in seams):
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


def sparse_hit_laws(h, sites):
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
    for left, right in PAIRS:
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
    return laws


def projector_check(h, sites):
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
    starts = [i for pair in PAIRS for i in pair]
    ends = [i for pair in PAIRS for i in reversed(pair)]
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
