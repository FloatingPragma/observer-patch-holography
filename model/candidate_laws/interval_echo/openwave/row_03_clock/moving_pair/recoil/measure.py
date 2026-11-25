"""A kinetic parity forecast of outward flow for the unfiltered moving ordered pair."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import jv

from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.counting import (
    binomial_limits,
)
from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_03_clock.moving_pair.measure import (
    fixtures,
)
from model.executable_model import federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index

HERE = Path(__file__).resolve().parent
SHOTS = 2**36
ANGLE = 2 * np.pi / 3
AGES = (0.0, 1.0, 4.0)
DELAYS = (0.125, 0.5)
# Fewer than 40,000 scalar bin comparisons in this acquisition.
TAIL = 0.01 / (2 * 40000)


def fly(hamiltonian, states, times, depth):
    center, radius = (12 - depth) / 2, (12 + depth) / 2
    scaled = (
        hamiltonian - center * sparse.eye(hamiltonian.shape[0], format="csr")
    ) / radius
    arg = radius * max(times)
    order = int(np.ceil(arg + 20 * (arg + 1) ** (1 / 3) + 40))
    degrees = np.arange(order + 1)
    coeff = 2 * (-1j) ** degrees * jv(degrees, radius * np.array(times)[:, None])
    coeff[:, 0] /= 2
    previous, curr = states, scaled @ states
    out = [a[0] * previous + a[1] * curr for a in coeff]
    for n in range(2, order + 1):
        previous, curr = curr, 2 * (scaled @ curr) - previous
        for i in range(len(times)):
            out[i] += coeff[i, n] * curr
    for t, v in zip(times, out, strict=True):
        v *= np.exp(-1j * center * t)
        np.testing.assert_allclose(
            np.sum(abs(v) ** 2, axis=0),
            np.sum(abs(states) ** 2, axis=0),
            rtol=0,
            atol=3e-12,
        )
    return out


def cut_indices(board, owner):
    cuts = []
    for a, b in board.external:
        if b // 12 == owner:
            a, b = b, a
        if a // 12 != owner:
            continue
        other = np.array([i for i in range(board.n) if i not in (a, b)])
        left = np.r_[ordered_index(a, other, board.n), ordered_index(other, a, board.n)]
        right = np.r_[
            ordered_index(b, other, board.n), ordered_index(other, b, board.n)
        ]
        cuts.append((left, right))
    assert len(cuts) == 3
    return cuts


def boundary_probabilities(state, cuts, instrument):
    """Randomly select one of the three bridges; score zero outside the one-occupied flag."""
    plus, minus = 0.0, 0.0
    for left, right in cuts:
        inside, outside = state[left], state[right]
        if instrument == "current":
            # Full local propagator at pi/2, up to a common phase on each pair.
            plus += np.sum(abs(outside + 1j * inside) ** 2) / 2
            minus += np.sum(abs(inside + 1j * outside) ** 2) / 2
        elif instrument == "parity":
            plus += np.sum(abs(inside + outside) ** 2) / 2
            minus += np.sum(abs(inside - outside) ** 2) / 2
        else:
            raise ValueError(instrument)
    return np.array([plus / 3, minus / 3, 1 - (plus + minus) / 3])


def counts(rng, probs):
    p = np.round(np.asarray(probs), 14)
    if np.min(p) < -1e-13 or abs(p.sum() - 1) > 2e-12:
        raise ValueError(probs)
    p = np.maximum(p, 0)
    p[np.argmax(p)] += 1 - p.sum()
    return rng.multinomial(SHOTS, p).tolist()


def interval(record, weights):
    k = np.array(record)
    assert k.sum() == SHOTS
    lo, hi = binomial_limits(k, SHOTS, TAIL)
    w = np.array(weights)
    return np.array([w @ np.where(w >= 0, lo, hi), w @ np.where(w >= 0, hi, lo)])


def interval_product(a, b):
    products = np.multiply.outer(a, b)
    return np.array([products.min(), products.max()])


def calibrated_sine(meter):
    # The pi/2 analyzer gives stay - transfer = sin(s theta) for each load.
    bands = []
    for load in (1, 2):
        for sign in (1, -1):
            bands.append(interval(meter[f"{load}/{sign}"][1], [-sign, sign]))
    common = np.array([max(x[0] for x in bands), min(x[1] for x in bands)])
    if common[0] > common[1]:
        raise ValueError("Kick calibrations on the two species disagree")
    return common


def minus(a, b):
    return a - b[::-1]


def half_states(state, rev):
    for name, phase in [("sharp", 1), ("quarter", -1j), ("reverse_quarter", 1j)]:
        yield name, ((1 + phase) * state + (1 - phase) * state[rev]) / 2


def deposits(board, owner):
    return [
        [
            p,
            p
            if board.degree[p] == 5
            else min(
                q
                for q in board.neighbours[p]
                if q // 12 == owner and board.degree[q] == 5
            ),
        ]
        for p in range(12 * owner, 12 * (owner + 1))
    ]


def evaluate(row):
    parity = interval(row["parity"], [1.5, -1.5, 0])
    j = [
        interval(row["current"][name], [1.5, -1.5, 0])
        for name in ("plus", "minus", "sham")
    ]
    diff = minus(j[0], j[1]) / 2
    predicted = np.sin(ANGLE) * parity
    result = {
        "kinetic_coherence": parity.tolist(),
        "opposite_kick_current": diff.tolist(),
        "kinetic_forecast": predicted.tolist(),
        "forecast_residual": minus(diff, predicted).tolist(),
    }
    # For this two-sided phase step the mean of the kicked currents is cos(theta) times sham.
    residual = (j[0] + j[1]) / 2 - np.cos(ANGLE) * j[2]
    result["cosine_residual"] = residual.tolist()
    even_work = (1 - np.cos(ANGLE)) * parity
    odd_work = np.sin(ANGLE) * j[2]
    result["postkick_work"] = {
        "plus": (even_work + odd_work).tolist(),
        "minus": minus(even_work, odd_work).tolist(),
    }
    result["departure"] = []
    for i in range(len(DELAYS)):
        p, m = (row["flight"][s][i] for s in ("plus", "minus"))
        result["departure"].append(
            {
                "centre_difference": minus(
                    interval(p, [0, 0.5, 1, 1, 1]), interval(m, [0, 0.5, 1, 1, 1])
                ).tolist(),
                "receiver_intact_difference": minus(
                    interval(p, [0, 0, 0, 0, 1]), interval(m, [0, 0, 0, 0, 1])
                ).tolist(),
            }
        )
    return result


def acquire():
    rng = np.random.default_rng(171291)
    board = federation(0)
    x, y, adjacent, kinetic, _, _ = generators(board)
    rev = ordered_index(y, x, board.n)
    routes = fixtures(board)
    starts = np.zeros((len(x), len(routes)))
    for i, (entry, _) in enumerate(routes):
        starts[ordered_index(*entry, board.n), i] = 1
    arms = {}
    for law, depth in [("native", 4 * candidate.KAPPA), ("exchange_only", 0.0)]:
        h = (kinetic - sparse.diags(depth * adjacent)).tocsr()
        expanding = fly(h, starts, AGES, depth)
        for it, age in enumerate(AGES):
            for i, (entry, receiver) in enumerate(routes):
                owner, dest = entry[0] // 12, receiver[0] // 12
                cuts = cut_indices(board, owner)
                outside = (x // 12 != owner).astype(int) + (y // 12 != owner)
                amount_square = (x // 12 != owner).astype(int) + 4 * (y // 12 != owner)
                both_at = (x // 12 == dest) & (y // 12 == dest)
                bins = np.select(
                    [outside == 0, outside == 1, ~both_at, ~adjacent],
                    [0, 1, 2, 3],
                    default=4,
                )
                state = expanding[it][:, i]
                signs = {"plus": 1, "minus": -1, "sham": 0}
                # The kick depends only on the number of occupied ports in the owner.
                # It commutes with global exchange of the two distinct species.
                base = np.column_stack(
                    [state * np.exp(1j * s * ANGLE * outside) for s in signs.values()]
                )
                flight = fly(h, base, DELAYS, depth)
                for kind, w in half_states(state, rev):
                    row = {
                        "parity": counts(
                            rng, boundary_probabilities(w, cuts, "parity")
                        ),
                        "current": {},
                        "flight": {},
                    }
                    for j, (sign, s) in enumerate(signs.items()):
                        v = w * np.exp(1j * s * ANGLE * outside)
                        row["current"][sign] = counts(
                            rng, boundary_probabilities(v, cuts, "current")
                        )
                        # Final categories are symmetric in the two species, hence
                        # insensitive to their initial relative parity phase.
                        row["flight"][sign] = [
                            counts(
                                rng, np.bincount(bins, abs(f[:, j]) ** 2, minlength=5)
                            )
                            for f in flight
                        ]
                    selective = w * np.exp(0.5j * np.pi * amount_square)
                    row["species_selective_current"] = counts(
                        rng, boundary_probabilities(selective, cuts, "current")
                    )
                    row["readings"] = evaluate(row)
                    if law == "native":
                        highest_work = max(
                            band[1]
                            for band in row["readings"]["postkick_work"].values()
                        )
                        row["adjacent_lower_after_either_kick"] = (
                            27 - highest_work
                        ) / depth
                    arms[f"{law}/{i}/{age:g}/{kind}"] = row
    # Isolated blank dimer split, phase one arm, and finish with a second
    # quarter exchange. Scanning an analyzer phase gives independent sine/cosine.
    meter = {}
    for load in (1, 2):
        for s in (0, 1, -1):
            phase = s * ANGLE * load**2
            # Analyzer phases 0,pi/2 are supplied calibrated references to the two paths.
            meter[f"{load}/{s}"] = [
                counts(
                    rng, [0.5 * (1 + np.cos(phase + b)), 0.5 * (1 - np.cos(phase + b))]
                )
                for b in (0, np.pi / 2)
            ]
    observed_sine = calibrated_sine(meter)
    for row in arms.values():
        readings = row["readings"]
        predicted = interval_product(observed_sine, readings["kinetic_coherence"])
        observed = np.array(readings["opposite_kick_current"])
        readings["ramsey_forecast"] = predicted.tolist()
        readings["ramsey_residual"] = minus(observed, predicted).tolist()
    clock = {
        str(t): counts(rng, [np.sin(t / 16) ** 2, np.cos(t / 16) ** 2])
        for t in AGES + DELAYS
    }
    return {
        "shots": SHOTS,
        "scalar_tail": TAIL,
        "routes": routes,
        "routed_ports": [deposits(board, r[0][0] // 12) for r in routes],
        "ages": AGES,
        "delays": DELAYS,
        "clock_swap_stay": clock,
        "phase_dimer_transfer_stay": meter,
        "observed_kick_sine": observed_sine.tolist(),
        "pointer_labels": ["plus", "minus", "unsuited"],
        "address_labels": [
            "both_inside",
            "one_out",
            "both_out_other",
            "receiver_apart",
            "receiver_adjacent",
        ],
        "settings": arms,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    payload = json.dumps(acquire(), allow_nan=False, indent=2) + "\n"
    if args.output:
        with args.output.open("x") as f:
            f.write(payload)
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(
            json.loads(payload),
            json.loads((HERE / "records.json").read_text()),
            "moving-pair recoil",
        )
    print(
        "Verified the clock, kick and kinetic/position readings on all three boundaries."
    )
