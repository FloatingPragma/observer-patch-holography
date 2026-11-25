"""Certified finite-time energy and number balances on nested moving regions."""

import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from math import factorial, isqrt
from pathlib import Path

import numpy as np
from scipy import sparse
import sympy as sp

from .. import current, flight, momentum_screen as face
from ... import live_radar, live_source
from ...source import certificate as cert

HERE = Path(__file__).parent
SCALE = 1 << 48
HALF = F(1, 1 << 16)
GROUPS = ("left_core", "left_rest", "shared_left", "shared_right", "right")


def memberships(board):
    left, shared, right = face.zones(board)
    left_hops, right_hops = (face.distances(board, p) for p in face.POLES)
    core = left & (left_hops <= 3)
    masks = (
        core,
        left & ~core,
        shared & (left_hops <= right_hops),
        shared & (left_hops > right_hops),
        right,
    )
    assert np.all(np.sum(masks, axis=0) == 1) and all(mask.any() for mask in masks)
    # Numerators of the linear membership w(t)=w0+(t-4)v, in units of 1/16.
    inner = np.array([16, 12, 2, 0, 0])
    outer = np.array([16, 16, 8, 8, 0])
    iv = np.array([0, -2, 1, 0, 0])
    ov = np.array([0, 0, 4, -2, 0])
    profiles = {
        "inner": (inner, iv),
        "annulus": (outer - inner, ov - iv),
        "exterior": (16 - outer, -ov),
        "outer": (outer, ov),
    }
    for w, v in profiles.values():
        assert min(w - abs(v) * HALF) >= 0 and max(w + abs(v) * HALF) <= 16
    return masks, profiles


def all_port_balance(board, h):
    """One sparse row per port, one column per Hamiltonian matrix entry.

    Enumerate signed contributions to [A,{P_p,A}] for all p and compare
    each resulting operator to [A²,P_p]. Here A=40H is integral.
    """
    d = h.shape[0]
    port_rows = live_radar.station_rows(board, list(range(board.n)))
    owner = np.empty(d, dtype=np.int32)
    owner[port_rows] = np.arange(board.n)[:, None]
    a = h.tocoo()
    counts = np.diff(h.indptr)[a.col]
    i = np.repeat(a.row, counts)
    k = np.repeat(a.col, counts)
    start = np.repeat(h.indptr[a.col] - np.cumsum(counts) + counts, counts)
    legs = start + np.arange(len(k))
    j = h.indices[legs]
    coeff = np.repeat(a.data, counts) * h.data[legs]
    cell = i.astype(np.int64) * d + j
    shape = (board.n, d * d)

    def packet(ids, slots, values):
        return sparse.coo_matrix(
            (np.concatenate(values), (np.concatenate(ids), np.concatenate(slots))),
            shape=shape,
        ).tocsr()

    direct = packet(
        [owner[j], owner[k], owner[k], owner[i]],
        [cell] * 4,
        [coeff, coeff, -coeff, -coeff],
    )
    a2 = (h @ h).tocoo()
    spot = a2.row.astype(np.int64) * d + a2.col
    reduced = packet([owner[a2.col], owner[a2.row]], [spot] * 2, [a2.data, -a2.data])
    assert (direct - reduced).count_nonzero() == 0
    reduced.eliminate_zeros()
    edges = {tuple(sorted((int(i), int(j)))) for i, j in board.seams}

    def pairs(x):
        return {
            tuple(sorted((int(i), int(j))))
            for i, j in zip(owner[x.row], owner[x.col], strict=True)
            if i != j
        }

    first, second = pairs(a), pairs(a2)
    assert first <= edges and all(
        set(board.neighbours[i]).intersection(board.neighbours[j])
        for i, j in second - edges
    )
    return dict(
        verified_port_operators=board.n,
        nonzero_regional_commutator_entries=int(reduced.nnz),
        graph_bonds_with_number_current=len(first),
        energy_pairs_on_graph_bonds=len(second & edges),
        energy_pairs_two_hops_apart=len(second - edges),
    )


def sparse_checks(h, groups, profiles):
    """Check all three bin assemblies as exact integer operator identities."""
    labels = np.empty(h.shape[0], dtype=int)
    for i, g in enumerate(groups):
        labels[g] = i
    h2 = h @ h
    mat = (h, h2)
    stats = []
    for a in mat:
        coo = a.tocoo()
        foreign = labels[coo.row] != labels[coo.col]
        stats.append(int(np.count_nonzero(foreign)))
    for weights, rates in profiles.values():
        for w in (weights, rates):
            d = sparse.diags(w[labels], format="csr", dtype=np.int64)
            e = d @ h + h @ d
            assert ((h @ e - e @ h) - (h2 @ d - d @ h2)).count_nonzero() == 0
            for a in mat:
                coo = a.tocoo()
                delta = w[labels[coo.row]] - w[labels[coo.col]]
                z = sparse.coo_matrix(
                    (coo.data * delta, (coo.row, coo.col)), shape=h.shape
                ).tocsr()
                assert (z - (d @ a - a @ d)).count_nonzero() == 0
    for i in (0, 1):
        b = [profiles[x][i] for x in ("inner", "annulus", "exterior")]
        assert np.all(sum(b) == (16 if i == 0 else 0))
    return dict(
        hilbert_dimension=h.shape[0],
        ports=len(labels) // 8,
        cross_group_oriented_matrix_entries=dict(number=stats[0], energy=stats[1]),
    )


def center_waves(h):
    hi, lo, e, _ = live_source.pump_center(True)
    shift = (1 << 82) // SCALE
    center = np.asarray(
        ((hi.astype(object) * (1 << 36) + lo + shift // 2) // shift), dtype=np.int64
    )
    # All 8 columns together are bounded by their Frobenius norm.
    rounding = F(1 + isqrt(2 * center.shape[0] * 8), 2 * SCALE)
    base = 3 * F(e) + rounding
    hbound = F(int(abs(h).sum(axis=1).max()), 40)

    def endpoint(sign):
        degree = 4
        wave, term = center.copy(), center.copy()
        error = 0
        total = 0
        for j in range(1, degree + 1):
            denom = 40 * HALF.denominator * j
            assert int(abs(h).sum(axis=1).max()) * int(abs(term).max()) + denom < 2**62
            re, im = np.hsplit(h @ term, 2)
            term = (sign * np.hstack([im, -re]) + denom // 2) // denom
            wave += term
            error = error * hbound * HALF / j + rounding
            total += error
        tail = (HALF * hbound) ** (degree + 1) / factorial(degree + 1)
        return wave, total + tail + (1 + tail) * base

    return endpoint(-1), (center, base), endpoint(1)


def gram(left, right):
    a, b = np.hsplit(left.astype(object), 2)
    c, d = np.hsplit(right.astype(object), 2)
    return (
        sp.Matrix(a.T @ c + b.T @ d) + sp.I * sp.Matrix(a.T @ d - b.T @ c)
    ) / SCALE**2


def enclose_effect(op, norm, wave_error, preparations):
    radius = 2 * wave_error + wave_error**2
    out = {}
    for name, (state, err) in preparations.items():
        value = sp.re(sp.trace(op * state))
        e = norm * (err + (1 + err) * radius)
        # Enclose the elementary result on a binary grid before summing bonds.
        mid = round(value * (1 << 70))
        width = int(
            sp.ceiling((abs(value - sp.Rational(mid, 1 << 70)) + e) * (1 << 70))
        )
        out[name] = (F(mid, 1 << 70), F(width, 1 << 70))
    return out


def inventories(h, wave, error, groups, bound, states):
    assert int(abs(h).sum(axis=1).max()) * int(abs(wave).max()) < 2**62
    hw = h @ wave
    out = {"number": [], "energy": []}
    for g in groups:
        n = gram(wave[g], wave[g])
        ph = gram(wave[g], hw[g])
        out["number"].append(enclose_effect(n, 1, error, states))
        out["energy"].append(enclose_effect((ph + ph.H) / 80, bound, error, states))
    return out


def bond_currents(h, wave, error, groups, bound, states):
    # Convert the 48-bit center to the two-limb 82-bit layout of the native
    # cross-current contraction, exactly, with no extra rounding.
    hi = wave >> 2
    lo = (wave & 3) << 34
    out = {"number": {}, "energy": {}}
    for tag, m, den, lim in (
        ("number", h, 40, 2 * bound),
        ("energy", h @ h, 3200, bound**2),
    ):
        for i, j in combinations(range(len(groups)), 2):
            inside = groups[i]
            keep = np.zeros(h.shape[0], dtype=np.int64)
            keep[groups[j]] = 1
            cross = m[inside].multiply(keep)
            cross.eliminate_zeros()
            op = current.compressed(cross, inside, hi, lo, den)
            out[tag][i, j] = enclose_effect(op, lim, error, states)
    return out


def sum_read(weights, readings):
    v = e = 0
    for w, (z, r) in zip(weights, readings, strict=True):
        v += w * z
        e += abs(w) * r
    return v, e


def band(v):
    return cert.enclosed(*v)


def panel(arms, profile, left, mid, right, bonds, bound):
    w, velocity = (tuple(F(int(x), 16) for x in a) for a in profile)
    duration = 2 * HALF
    d = max(map(abs, velocity))
    out = {}
    for mat, lim, flux_bound in (("number", 1, 2 * bound), ("energy", bound, bound**2)):
        rows = {}
        for arm in arms:

            def stock(table, weights):
                return sum_read(weights, [z[arm] for z in table[mat]])

            initial = stock(left, [a - HALF * b for a, b in zip(w, velocity)])
            final = stock(right, [a + HALF * b for a, b in zip(w, velocity)])
            loss = sum_read([1, -1], [initial, final])
            curr = bonds[mat]
            cut = sum_read(
                [w[i] - w[j] for i, j in curr], [z[arm] for z in curr.values()]
            )
            swept = stock(mid, velocity)
            # For F=Tr(rho J_w), J linear in weights, ||F''|| <=
            # (2||H||)^2 ||J_w|| + 4||H|| ||J_w'||. For sweep w' is constant.
            cm, cw = cut
            sm, sw = swept
            curv = (2 * bound) ** 2 * flux_bound + 4 * bound * d * flux_bound
            motion = (2 * bound) ** 2 * d * lim
            integrated = (duration * cm, duration * cw + duration**3 * curv / 24)
            transport = (duration * sm, duration * sw + duration**3 * motion / 24)
            prediction = sum_read([1, -1], [integrated, transport])
            gap = sum_read([1, -1], [loss, prediction])
            assert abs(gap[0]) <= gap[1]
            frozen = sum_read([1, -1], [loss, integrated])
            rows[arm] = dict(
                endpoint_loss=band(loss),
                integrated_cut_current=band(integrated),
                integrated_membership_change=band(transport),
                corrected_residual=band(gap),
                residual_omitting_motion=band(frozen),
            )
        out[mat] = rows
    return out


@lru_cache(maxsize=1)
def run():
    h = flight.hardware(-1)[0]
    bound = F(int(abs(h).sum(axis=1).max()), 40)
    board = live_source.source.formation.federation(1)
    masks, profiles = memberships(board)
    groups = [
        live_radar.station_rows(board, np.flatnonzero(mask).tolist()).ravel()
        for mask in masks
    ]
    checks = sparse_checks(h, groups, profiles)
    checks.update(all_port_balance(board, h))
    states = {
        s: (sp.kronecker_product(sp.diag(*rail) / 3, head), F(e))
        for s, (head, rail, e) in face.source_states().items()
    }
    waves = center_waves(h)
    tables = [inventories(h, v, e, groups, bound, states) for v, e in waves]
    cut = bond_currents(h, *waves[1], groups, bound, states)
    return dict(
        appointments=[float(4 - HALF), 4, float(4 + HALF)],
        port_groups=dict(zip(GROUPS, [int(mask.sum()) for mask in masks], strict=True)),
        native_operator_checks=checks,
        region_weights_and_rates_over_16={
            k: [a.tolist(), b.tolist()] for k, (a, b) in profiles.items()
        },
        signed_group_currents={
            k: {
                GROUPS[i] + ":" + GROUPS[j]: {s: band(x) for s, x in vals.items()}
                for (i, j), vals in pairs.items()
            }
            for k, pairs in cut.items()
        },
        regions={k: panel(states, v, *tables, cut, bound) for k, v in profiles.items()},
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    payload = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(payload)
    else:
        print(payload, end="")
