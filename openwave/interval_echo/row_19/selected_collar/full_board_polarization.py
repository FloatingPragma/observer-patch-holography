"""Certify steady ordinary polarization with all external contacts open throughout."""

from fractions import Fraction as F
from math import isqrt

import numpy as np
import sympy as sp

from . import stationary_polarization as steady

motion = steady.motion
exp = steady.exp
local = steady.local
SEED = 28122916
DEN = 2**50


def gram(q):
    # Three signed-base digits keep every integer product in signed int64.
    base = 2**20
    parts = []
    for _ in range(2):
        part = q % base
        parts.append(part)
        q = (q - part) // base
    parts.append(q)
    return sum(
        (a.T @ b).astype(object) * base ** (i + j)
        for i, a in enumerate(parts)
        for j, b in enumerate(parts)
    )


def certified_blocks(h):
    n = len(h)
    values, v = np.linalg.eigh(h / 2)
    cp = sp.Matrix(h.tolist()).to_DM().charpoly()
    t = sp.Symbol("t")
    distinct = sp.Poly.from_list([int(z) for z in cp], t).sqf_part().degree()
    cuts = sorted(np.argsort(np.diff(values))[-(distinct - 1) :] + 1)
    groups = np.split(np.arange(n), cuts)
    lab = np.empty(n, dtype=np.int64)
    centers = []
    for g in groups:
        mid = round((values[g[0]] + values[g[-1]]) * DEN / 2)
        centers.append(mid)
        lab[g] = mid
    q = np.rint(v * DEN).astype(np.int64)
    orth = F(
        int(np.max(np.sum(abs(gram(q) - np.eye(n, dtype=object) * DEN**2), axis=1))),
        DEN**2,
    )
    assert orth < F(1, 10000)
    # h is twice H and the absolute row sum bounds prevent int64 overflow.
    assert max(np.sum(abs(h), axis=1)) * DEN < 2**63
    residual = (h @ q).astype(object) * DEN - 2 * q.astype(object) * lab.astype(object)
    sq = int(np.sum(residual**2))
    residual_bound = F(isqrt(sq) + 1, 2 * DEN**2)
    gap = F(min(y - x for x, y in zip(centers, centers[1:])), DEN)
    assert gap > 4 * residual_bound
    # The exact root count and the residual isolate every distinct energy.
    # Project columns onto the associated true eigenspaces. Their removed
    # Frobenius norm is at most residual_bound/(gap-2*residual_bound).
    # Polar orthogonalization there gives a blockwise exact eigenbasis.
    distance = 2 * residual_bound / (gap - 2 * residual_bound) + orth
    err = 5 * distance
    return (
        q / DEN,
        np.equal.outer(lab, lab),
        {
            "levels": int(distinct),
            "gram_error": str(orth),
            "residual_frobenius": str(residual_bound),
            "minimum_physical_gap_lower": str(gap - 4 * residual_bound),
            "pinching_trace_error": str(err),
        },
    )


def pinch(q, mask, rho):
    return q @ (mask * (q.T @ rho @ q)) @ q.T


def empty_table(setup):
    _, _, ports, ab, _, degree, edges, _ = setup
    pair = [
        local.spectrum(ports, edges, degree, (0, *ab), 0, q, motion.FLIGHTS)
        for q in (1, 2)
    ]
    table = np.zeros((len(motion.FLIGHTS), 8, 4))
    for k in range(len(motion.FLIGHTS)):
        for i in range(2):
            table[k, 2 * i : 2 * i + 2, :3] = pair[i]["panels"][k][0]
        for side in (0, 1):
            a, b = [r["panels"][k][1][side] for r in pair]
            retained = (abs(a) ** 2 + abs(b) ** 2) / 2
            vis = a.conjugate() * b
            for j, val in enumerate((vis.real, -vis.imag)):
                table[k, 4 + 2 * side + j, :3] = [
                    (retained + val) / 2,
                    (retained - val) / 2,
                    1 - retained,
                ]
    return table


def native_layouts(board, setup):
    ports, second = setup[2], setup[7]
    match = set(local.survey.routes.codebook()[0][0])
    active = [
        e
        for e in board.seams
        if e[0] // 12 != e[1] // 12 or (e[0] % 12, e[1] % 12) not in match
    ]
    marks = [{(ports[i], ports[j]) for i, j in bank} for bank in second]
    n = board.n
    rho = np.ones((n, n)) / n
    probs = {}
    kinetic = {}
    certificates = {}
    for flags in steady.FLAGS:
        h = np.diag(16 * (6 - board.degree) - 80)
        for a, b in active:
            if any(flag and (a, b) in m for flag, m in zip(flags, marks, strict=True)):
                continue
            h[a, a] += 17
            h[b, b] += 17
            h[a, b] -= 1
            h[b, a] -= 1
        q, mask, cert = certified_blocks(h)
        if flags == (0, 0):
            restored = q, mask
            kinetic_operator = h / 2
        full = pinch(*restored, pinch(q, mask, rho))
        z = full.diagonal()
        kinetic[flags] = float(np.sum(kinetic_operator * full.T) / z.sum())
        assert min(z) > 0
        probs[flags] = z / z.sum()
        certificates[str(flags)] = cert
    return active, probs, kinetic, certificates


def run():
    rng = np.random.default_rng(SEED)
    board = local.survey.bulk.federation(0)
    setup, _ = steady.symmetric_instrument(board)
    _, cells, ports, ab, _, degree, edges, _ = setup
    gates = [ports[z] for z in ab]
    n = board.n
    active, probs, kinetic, certificates = native_layouts(board, setup)
    # Bound approximation of the two energy averages in trace norm.
    # For double products the Frobenius roundoff is at most gamma_n times
    # the two Frobenius norms. Each density has Frobenius norm at most 1
    # before rounding, and the near-isometries have norm sqrt(n*(1+orth)).
    roundoff = F(100 * n * n, 2**52)
    first = max(F(z["pinching_trace_error"]) for z in certificates.values())
    final = F(certificates[str((0, 0))]["pinching_trace_error"])
    error = first + final + first * final + roundoff
    # A finite implementation uses four independent uniform integer waits
    # per stage. Each suppresses off-diagonal blocks by a Dirichlet factor.
    gap = min(F(z["minimum_physical_gap_lower"]) for z in certificates.values())
    bits = 1
    while 2 * (n - 1) * (F(22, 7) / (2**bits * steady.STEP * gap)) ** 4 > F(1, 10**11):
        bits += 1
    wait_err = 2 * (n - 1) * (F(22, 7) / (2**bits * steady.STEP * gap)) ** 4
    # score and the local finite reader conditional on the measured site.
    score = [
        int(j not in gates and tuple(sorted((gates[0], j))) in active)
        - int(j not in gates and tuple(sorted((gates[1], j))) in active)
        for j in range(n)
    ]
    subpanels = motion.second_instruments(ports, ab, degree, edges)
    panels = np.tile(empty_table(setup)[:, None, :, :], (1, n, 1, 1))
    panels[:, ports] = subpanels
    near = np.array([int(j // 12 == cells[0]) for j in range(n)])
    fwd, rev, dispatch = steady.parity(board, cells)
    out = {}
    for label, law in (("forward", fwd), ("reversed", rev)):
        accepted = sum(float(w) * probs[f] for f, w in law.items())
        dp = (accepted - probs[0, 0]) / 144
        target = float(dp @ score)
        arms = {}
        for name, prob in (("fresh", 25 / 36), ("shared", 101 / 144)):
            dist = (1 - prob) * probs[0, 0] + prob * accepted
            dist /= dist.sum()
            arms[name] = motion.observe(
                dist, panels, near, score, gates, rng, starts=2**60
            )
        a, b = arms["shared"], arms["fresh"]
        out[label] = {
            "steady_available": target,
            "steady_near": float(dp @ near),
            "steady_mass_by_cell": dp.reshape(-1, 12).sum(axis=1).tolist(),
            "available_from_counts": local.difference(a["available"], b["available"]),
            "predicted_release": [v * target for v in (8.5, 16.25, -15.5)],
            "arms": arms,
            "held_flight": [
                [
                    local.difference(x, y)
                    for x, y in zip(r["scores"], s["scores"], strict=True)
                ]
                for r, s in zip(a["held_flights"], b["held_flights"], strict=True)
            ],
        }
    # Masked cells keep all accepted and unaccepted routes. At either tag
    # the same conditional ensembles are used, so the contrast error is
    # at most twice the branch error times the terminal excess probability.
    for name, law in (("forward", fwd), ("reversed", rev)):
        weighted = sum(float(w) * kinetic[f] for f, w in law.items())
        energies = {}
        for tag, prob in (("fresh", 25 / 36), ("shared", 101 / 144)):
            chance = ((1 - prob) * kinetic[0, 0] + prob * weighted) / len(active)
            energies[tag] = local.acquisition(
                [1 - chance, chance], [0, len(active)], rng, 2**60
            )
        out[name]["kinetic_meter"] = {
            "active_seams": len(active),
            "arms": energies,
            "contrast_counts": local.difference(
                energies["shared"]["mean"], energies["fresh"]["mean"]
            ),
            "predicted_energy": (weighted - kinetic[0, 0]) / 144,
        }
    # Normalize each approximate branch before mixing. Trace error expands
    # by at most twice its bound divided by one minus that bound.
    contrast_error = float(2 * (2 * error / (1 - error) + wait_err) / 144) + 1e-14
    pos_score = np.array([0, 1, 2, 0])
    clock_score = np.array([1, -1, 0, 0])
    for row in out.values():
        row["held_flight_enclosed"] = []
        for raw, meter, dt in zip(
            row["held_flight"], panels, motion.FLIGHTS, strict=True
        ):
            observables = [
                (meter[:, i] - meter[:, i + 1]) @ pos_score / dt**2 for i in (0, 2)
            ]
            observables.append((meter[:, 5] - meter[:, 7]) @ clock_score / dt)
            row["held_flight_enclosed"].append(
                [
                    [
                        lo - contrast_error * max(abs(obs)) - 1e-11,
                        hi + contrast_error * max(abs(obs)) + 1e-11,
                    ]
                    for (lo, hi), obs in zip(raw, observables, strict=True)
                ]
            )
    return {
        "seed": SEED,
        "cells": list(cells),
        "gates": gates,
        "opposite_dispatch": dispatch,
        "certificates": certificates,
        "per_branch_trace_error_upper": float(error),
        "secondary_starts": 2**60,
        "source_contrast_score_error_upper": contrast_error,
        "unnormalized_contrast_error": float(
            2 * (error + wait_err) / 144 + F(1, 10**14)
        ),
        "rest_energy": 40,
        "kinetic_bandwidth_upper": 5,
        "arrival_excess": "1/144",
        "homogeneous_rest_fraction": 0,
        "maximum_extra_energy_over_ordinary_rest": "1/1152",
        "waits_per_stage": 4,
        "bits_each": bits,
        "two_stage_wait_error_upper": str(wait_err),
        "outcomes": out,
        "inputs": exp.save_inputs(
            "FULL_BOARD_POLARIZATION.md",
            "stationary_polarization.json",
            "routed_contacts.json",
        ),
    }


if __name__ == "__main__":
    r = run()
    exp.emit("full_board_polarization", r)
    print(
        r["source_contrast_score_error_upper"],
        [
            (k, v["steady_available"], v["available_from_counts"])
            for k, v in r["outcomes"].items()
        ],
    )
