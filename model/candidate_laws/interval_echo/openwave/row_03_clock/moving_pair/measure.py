"""Four source phases and two local order reads across an open cell boundary."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import gammaln, jv

from .counting import binomial_limits

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.material_inertia.neutral_defect import generators, ordered_index
from model.experiments.material_inertia.pair_recoil.measure import PLACEMENTS

HERE = Path(__file__).parent
SHOTS = 2**32
TIMER_SHOTS = 2**26
TIMES = np.arange(49) * (np.pi / 4)
PHASES = np.arange(4) * (np.pi / 2)
TAIL = 2e-15
BIN_TAIL = 1e-15


def fixtures(board):
    routes = []
    for bridge_id in PLACEMENTS:
        ends = []
        for gateway in board.external[bridge_id]:
            regular = min(
                j
                for j in board.neighbours[gateway]
                if j // 12 == gateway // 12 and board.degree[j] == 5
            )
            ends.append([regular, gateway])
        routes.append(ends)
    return routes


def amplitudes(h, source, arrival, lower, upper, times):
    """Return only the two sharp receiving addresses from a full Chebyshev flight."""
    center, radius = (lower + upper) / 2, (upper - lower) / 2
    scaled = (h - center * sparse.eye(h.shape[0], format="csr")) / radius
    argument = radius * np.asarray(times)
    largest = argument.max()
    order = max(2, math.ceil(largest))
    while True:
        n = order + 1
        ratio = largest / (2 * (n + 1))
        log_bound = (
            np.log(2) + n * np.log(largest / 2) - gammaln(n + 1) - np.log1p(-ratio)
        )
        if log_bound < np.log(TAIL):
            break
        order += 1
    columns = np.arange(len(source))
    previous = np.zeros((h.shape[0], len(source)))
    previous[source, columns] = 1
    current = scaled @ previous
    moments = np.empty((order + 1, len(arrival), 2))
    moments[0] = previous[arrival, columns[:, None]]
    moments[1] = current[arrival, columns[:, None]]
    for n in range(2, order + 1):
        previous, current = current, 2 * (scaled @ current) - previous
        moments[n] = current[arrival, columns[:, None]]
    n = np.arange(order + 1)
    coefficients = 2 * (-1j) ** n * jv(n, argument[:, None])
    coefficients[:, 0] /= 2
    waves = np.einsum("tn,nrs->rts", coefficients, moments)
    return waves * np.exp(-1j * center * times)[None, :, None], order


def acquire(even, odd, rng, erased=False, shots=SHOTS):
    moving = odd[:, None] * np.exp(-1j * PHASES)
    p = np.stack([abs((even[:, None] + s * moving) / 2) ** 2 for s in (1, -1)], axis=-1)
    if erased:
        p = (p + np.roll(p, 2, axis=1)) / 2
    if np.any(p < 0) or np.any(p.sum(axis=-1) > 1 + 1e-13):
        raise ValueError("Invalid terminal probability")
    p = np.round(p, 14)
    counts = []
    for row in p:
        cohort = []
        for forward, reverse in row:
            # Sample the miss first; the final category is its exact count complement.
            total = forward + reverse
            if total == 0:
                cohort.append([0, 0, shots])
            else:
                hits = rng.binomial(shots, min(1.0, total))
                ordered = rng.binomial(hits, forward / total)
                cohort.append([ordered, hits - ordered, shots - hits])
        counts.append(cohort)
    return counts


def decode(banks, shots=SHOTS):
    cells = np.asarray(banks)
    if cells.shape[1:] != (4, 3) or not np.all(cells.sum(axis=-1) == shots):
        raise ValueError("Incomplete receipt")
    counts = cells[:, :, :2]
    rates = counts / shots
    # Exact marginal binomial intervals, combined by the union and triangle bounds.
    lo, hi = binomial_limits(counts, shots, BIN_TAIL)
    errors = np.maximum(rates - lo, hi - rates).sum(axis=-1)
    d = rates[:, :, 0] - rates[:, :, 1]
    z = (d[:, 0] - d[:, 2] + 1j * (d[:, 1] - d[:, 3])) / 2
    dz = np.hypot(errors[:, 0] + errors[:, 2], errors[:, 1] + errors[:, 3]) / 2
    return rates.sum(axis=-1).mean(axis=1), z, dz


def phase_witness(z, dz, bench_z, bench_dz, i, j):
    u = z * bench_z.conj()
    du = abs(z) * bench_dz + abs(bench_z) * dz + dz * bench_dz
    witness = abs((u[i] * u[j].conj()).imag)
    allowance = abs(u[i]) * du[j] + abs(u[j]) * du[i] + du[i] * du[j]
    w = np.argmax(witness / allowance)
    return dict(
        appointments=[int(i[w]), int(j[w])],
        imaginary_product=float(witness[w]),
        simultaneous_allowance=float(allowance[w]),
        exceeds_allowance=bool(witness[w] > allowance[w]),
    )


def assess(banks, bench):
    arrival, z, dz = decode(banks)
    _, bench_z, bench_dz = decode(bench)
    train, held = slice(1, 25), slice(25, 49)
    x = arrival * bench_z
    factor = np.vdot(x[train], z[train]) / np.vdot(x[train], x[train]).real
    retained = np.arange(len(z))[held]
    i, j = retained[np.array(np.triu_indices(len(retained), k=1))]
    best = np.argmax(abs(z) / dz)
    return dict(
        peak_arrival=float(max(arrival)),
        peak_hand=float(max(abs(z))),
        best_resolved_hand=dict(
            j=int(best),
            clicks_per_launch=float(abs(z[best])),
            simultaneous_radius=float(dz[best]),
            excludes_zero=bool(abs(z[best]) > dz[best]),
        ),
        first_half_complex_fit=[float(factor.real), float(factor.imag)],
        held_rms_per_launch=float(
            np.sqrt(np.mean(abs(z[held] - factor * x[held]) ** 2))
        ),
        phase_only_witness=phase_witness(z, dz, bench_z, bench_dz, i, j),
    )


def energy_meter(board, routes):
    """Sample actual mean exchange parity after the four local preparations."""
    rng = np.random.default_rng(843197)
    results = []
    for start, _ in routes:
        a, b = start
        edges = [e for e in board.seams if a in e or b in e]
        kinetic = len(edges)
        if kinetic != 10:
            raise ValueError("Source is not the enrolled degree-five/gateway contact")
        arms = []
        for phi in PHASES:
            # The pulse preserves equal exchange-parity weights on the source bond.
            # On every other incident contact exactly one endpoint is empty. Its
            # swapped word is orthogonal to both prepared words.
            amplitudes = np.array([1 + np.exp(-1j * phi), 1 - np.exp(-1j * phi)]) / 2
            source_odd = abs(amplitudes[0] - amplitudes[1]) ** 2 / 2
            odd_probability = (source_odd + (len(edges) - 1) / 2) / len(edges)
            exposures = 2**20
            odd = rng.binomial(exposures, odd_probability)
            static = 6 * candidate.KAPPA * (1**2 + 2**2) - 2 * candidate.KAPPA * 1 * 2
            arms.append(
                dict(
                    odd=odd,
                    even=exposures - odd,
                    static=static,
                    mean_energy=static + kinetic * odd / exposures,
                )
            )
        results.append(arms)
    return results


def measure():
    rng = np.random.default_rng(3217536)
    board = federation(0)
    _, _, adjacent, kinetic, _, _ = generators(board)
    routes = fixtures(board)
    source = [ordered_index(*s, board.n) for s, _ in routes]
    arrival = np.array(
        [
            [ordered_index(*r, board.n), ordered_index(*r[::-1], board.n)]
            for _, r in routes
        ]
    )
    bench = acquire(np.ones(len(TIMES)), np.exp(-1j * TIMES), rng)
    timer = rng.binomial(TIMER_SHOTS, np.sin(TIMES / 64) ** 2)
    out = dict(
        cells=board.cells,
        configurations=kinetic.shape[0],
        shots=SHOTS,
        timer_shots=TIMER_SHOTS,
        timer_swaps=timer.tolist(),
        isolated_bond_orders=bench,
        bounds_per_marginal_tail=BIN_TAIL,
        routes=routes,
        settings={},
    )
    for label, depth in (("registered", 4 * candidate.KAPPA), ("exchange_only", 0)):
        h = (kinetic - sparse.diags(depth * adjacent, dtype=float)).tocsr()
        waves, order = amplitudes(h, source, arrival, -depth, 12, TIMES)
        scans = []
        for address, (forward, reverse) in zip(
            source, waves.transpose(0, 2, 1), strict=True
        ):
            even, odd = forward + reverse, forward - reverse
            bank = acquire(even, odd, rng)
            erased = acquire(even, odd, rng, erased=True)
            scans.append(
                dict(
                    orders=bank,
                    erased_orders=erased,
                    measured=assess(bank, bench),
                    erased=assess(erased, bench),
                    mean_shifted_energy=float(h[address, address]),
                )
            )
        out["settings"][label] = dict(degree=order, analytic_tail=TAIL, scans=scans)
        if label == "registered":
            fresh_rng = np.random.default_rng(809328)
            fresh = []
            shots = 2**38
            for (forward, reverse), js in zip(
                waves.transpose(0, 2, 1), ([34, 44], [34, 44], [42, 44]), strict=True
            ):
                even = (forward + reverse)[js]
                odd = (forward - reverse)[js]
                bank = acquire(even, odd, fresh_rng, shots=shots)
                isolated = acquire(
                    np.ones(2), np.exp(-1j * TIMES[js]), fresh_rng, shots=shots
                )
                _, z, dz = decode(bank, shots)
                _, b, db = decode(isolated, shots)
                witness = phase_witness(z, dz, b, db, np.array([0]), np.array([1]))
                witness["appointments"] = js
                fresh.append(dict(orders=bank, bench=isolated, witness=witness))
            out["independent_followup"] = dict(shots=shots, scans=fresh)
    out["initial_energy_readings"] = energy_meter(board, routes)
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    payload = json.dumps(measure(), separators=(",", ":"), allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(payload)
    elif json.loads(payload) != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("Reproduced data differ from the saved receipt")
    print("Verified four-phase flights on three registered routes.")
