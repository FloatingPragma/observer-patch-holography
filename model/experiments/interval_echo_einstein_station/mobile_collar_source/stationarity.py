"""Measure local stationarity of the four generated sharp records in exposure time."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import betaincinv, xlogy

from . import measure as formation

HERE = Path(__file__).resolve().parent
STEP = 1 / 64
DELTAS = np.array([-STEP, -STEP / 2, 0.0, STEP / 2, STEP])
SHOTS = formation.SHOTS
SEED = 763017
TAIL = 0.01 / (2 * 1000)


def bin_intervals(counts):
    """Clopper–Pearson interval endpoints, including zero or unanimous counts."""
    n = np.asarray(counts)
    lo = np.zeros(n.shape)
    hi = np.ones(n.shape)
    pos = n > 0
    missing = n < SHOTS
    lo[pos] = betaincinv(n[pos], SHOTS - n[pos] + 1, TAIL)
    hi[missing] = betaincinv(n[missing] + 1, SHOTS - n[missing], 1 - TAIL)
    return np.stack([lo, hi], axis=-1)


def algebras(joint):
    table = joint.reshape(2, 2)
    return dict(joint=joint.ravel(), left=table.sum(axis=1), right=table.sum(axis=0))


def entropy(p):
    return float(-xlogy(p, p).sum())


def entropy_read(counts):
    iv = bin_intervals(counts)
    a, b = iv.T
    bounds = -xlogy(iv, iv)
    maximum = np.maximum(bounds[:, 0], bounds[:, 1])
    maximum[(a <= 1 / np.e) & (b >= 1 / np.e)] = 1 / np.e
    return dict(
        nats=entropy(counts / SHOTS),
        interval=[float(bounds.min(axis=1).sum()), float(maximum.sum())],
    )


def draw(p, level, q, arm, index):
    rng = np.random.default_rng([SEED, level, q + 2, arm, index])
    counts = rng.multinomial(SHOTS, p)
    return dict(
        counts=counts.tolist(),
        entropy={name: entropy_read(ns) for name, ns in algebras(counts).items()},
    )


def panel(level, q, arm, appointments, states, reference):
    pilot = draw(reference, level, q, arm + 10, 0)
    r = bin_intervals(pilot["counts"])
    rows = []
    for i, (t, p) in enumerate(zip(appointments, states, strict=True)):
        read = draw(p, level, q, arm, i)
        v = bin_intervals(read["counts"])
        distances = np.maximum(0, np.maximum(v[:, 0] - r[:, 1], r[:, 0] - v[:, 1]))
        tv = np.abs(np.array(read["counts"]) - np.array(pilot["counts"])).sum() / (
            2 * SHOTS
        )
        rows.append(
            dict(
                appointment=float(t),
                joint_probabilities=p.tolist(),
                reference_tv=float(tv),
                reference_tv_lower_99=float(distances.sum() / 2),
                **read,
            )
        )
    slopes = []
    for first_index, last_index in ((0, 4), (1, 3)):
        dt = appointments[last_index] - appointments[first_index]
        values = {}
        for name in ("joint", "left", "right"):
            first, last = (
                rows[first_index]["entropy"][name],
                rows[last_index]["entropy"][name],
            )
            interval = [
                (last["interval"][0] - first["interval"][1]) / dt,
                (last["interval"][1] - first["interval"][0]) / dt,
            ]
            values[name] = dict(
                nats_per_exposure=float((last["nats"] - first["nats"]) / dt),
                interval=interval,
            )
        slopes.append(dict(half_step=dt / 2, entropies=values))
    return dict(reference_counts=pilot["counts"], reads=rows, symmetric_slopes=slopes)


def run():
    rows = []
    k = formation.emission.K
    native = formation.emission.four_h(k)
    labels = np.repeat(
        [-1, *(abs(d) for _, d in formation.emission.LAW.LABELS)],
        formation.emission.LAW.N_SRC,
    )
    blank = (labels == -1).astype(int)
    blank_entries = int(np.count_nonzero((blank[:, None] - blank[None, :]) * native))
    assert blank_entries == 0
    for level in (0, 1):
        board = formation.federation(level)
        route = formation.region(board)
        for q in (-1, 1):
            coins = formation.starting_coins(route, q)
            h = formation.hamiltonian(board, route, q, k)
            start = np.zeros(h.shape[0], complex)
            start[4 * route[0] + 2 * coins[0] + coins[1]] = 1
            states = formation.flow(h, start, 4 + DELTAS)
            distributions = [np.diag(formation.reduce(w)).real for w in states]
            middle = states[2]
            rate = 2 * np.real(middle.conj() * (-1j * (h @ middle))).reshape(-1, 4).sum(
                axis=0
            )
            np.testing.assert_allclose(rate.sum(), 0, rtol=0, atol=4e-13)
            ref = distributions[2]
            assert min(ref) > 0
            slopes = {
                name: float(-deriv @ np.log(algebras(ref)[name]))
                for name, deriv in algebras(rate).items()
            }
            contact = panel(level, q, 0, 4 + DELTAS, distributions, ref)
            block = formation.local_parts(q, k)[0, 0]
            uniform = np.kron(block, formation.COIN_IDENTITY) + np.kron(formation.COIN_IDENTITY, block)
            np.testing.assert_array_equal(uniform, 2 * k * np.eye(4))
            rho = formation.reduce(middle)
            delays = STEP + DELTAS
            stored = []
            for t in delays:
                flight = expm(-1j * t * uniform)
                stored.append(np.diag(flight @ rho @ flight.conj().T).real)
            np.testing.assert_allclose(
                stored, np.broadcast_to(ref, (len(delays), 4)), rtol=0, atol=1e-14
            )
            removed = panel(level, q, 1, delays, stored, ref)
            sharp = np.eye(4)[2 * coins[0] + coins[1]]
            vacuum = draw(sharp, level, q, 2, 0)
            rows.append(
                dict(
                    cells=board.cells,
                    charge=q,
                    route=route,
                    exact_contact_probability_derivative=rate.tolist(),
                    exact_contact_entropy_derivatives=slopes,
                    contact=contact,
                    removed=removed,
                    restored_material_head_energy=float(uniform[0, 0].real),
                    uniform_armed_sharp_control=vacuum,
                )
            )
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "STATIONARITY_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        shots=SHOTS,
        seed=SEED,
        per_tail_probability=TAIL,
        blank_commutator_nonzero_entries=blank_entries,
        arms=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as stream:
        json.dump(run(), stream, indent=2, allow_nan=False)
        stream.write("\n")
