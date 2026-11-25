"""Compare attracting head preparations under one fixed four-station observer survey."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import root
from scipy.sparse.linalg import expm_multiply

from .. import measure as reuse
from ... import native_precision, station
from ....directional_stress import bridge

ODDS = np.array([0.8657, 0.1343])
UNIFORM = np.full(12, 1 / 12)
BASIS = station.iso.basis(False)
TANGENT = null_space(np.ones((1, 12)))
AXES = (
    "IX",
    "IY",
    "IZ",
    "XI",
    "XX",
    "XY",
    "XZ",
    "YI",
    "YX",
    "YY",
    "YZ",
    "ZI",
    "ZX",
    "ZY",
    "ZZ",
)
HERE = Path(__file__).resolve().parent


def edges(ranges):
    ordered = station.square(ranges)
    np.fill_diagonal(ordered, 1)
    logged = np.log(ordered)
    gain = (logged - logged.T).mean(axis=1)
    corrected = logged - gain[:, None]
    return np.exp((corrected + corrected.T)[np.triu_indices(4, 1)] / 2)


def entropy(blocks):
    eigen = np.concatenate([np.linalg.eigvalsh(b) for b in blocks])
    if eigen.min() <= 0:
        raise ValueError("non-faithful endpoint")
    return -float(eigen @ np.log(eigen))


def prepare(channel, weights):
    basis = BASIS.reshape(15, 16).T
    origin = np.eye(4).ravel() / 4
    states, contractions = [], []
    for phi in np.einsum("d,bdij->bij", weights, channel):
        jac = (basis.conj().T @ phi @ basis).real
        shift = (basis.conj().T @ phi @ origin).real
        coord = np.linalg.solve(np.eye(15) - jac, shift)
        states.append((origin + basis @ coord).reshape(4, 4))
        contractions.append(np.linalg.norm(jac, 2))
    return np.array(states) * ODDS[:, None, None], max(contractions)


def record(state, effects):
    ranges = station.pair(effects, state)
    return ranges, edges(ranges), station.geometry(ranges), entropy(state)


def edge_or_geometry(state, effects, mode):
    ranges = station.pair(effects, state)
    return edges(ranges) if mode == "six_edges" else station.geometry(ranges)


def difference_check(channel, effects, mode, spacing):
    # Each column is a direction on the simplex; order has no fitted scale.
    def evaluate(weights):
        state, _ = prepare(channel, weights)
        return np.r_[edge_or_geometry(state, effects, mode), entropy(state)]

    return np.stack(
        [
            (evaluate(UNIFORM + spacing * v) - evaluate(UNIFORM - spacing * v))
            / (2 * spacing)
            for v in TANGENT.T
        ],
        axis=1,
    )


def endpoints(channel, effects, mode):
    state, _ = prepare(channel, UNIFORM)
    basis = BASIS.reshape(15, 16).T
    changes = []
    for b, rho in enumerate(state):
        avg = channel[b].mean(axis=0)
        jac = (basis.conj().T @ avg @ basis).real
        shifted = np.einsum("dij,j->di", channel[b], rho.ravel()).T @ TANGENT
        coords = np.linalg.solve(np.eye(15) - jac, (basis.conj().T @ shifted).real)
        changes.append((basis @ coords).T.reshape(-1, 4, 4))
    variations = np.stack(changes, axis=1)
    ranges = station.pair(effects, state)
    geometry_read = edges if mode == "six_edges" else station.geometry
    geo_jac = np.stack(
        [geometry_read(ranges + 1e-17j * e).imag / 1e-17 for e in np.eye(12)], axis=1
    )
    range_slopes = np.einsum("pbij,dbji->pd", effects, variations).real
    rows = np.vstack(
        (geo_jac @ range_slopes, -station.pair(variations, station.log_blocks(state)))
    )
    confirm = difference_check(channel, effects, mode, 0.00002)
    geometry, score = rows[:-1], rows[-1]
    _, singular, vt = np.linalg.svd(geometry, full_matrices=False)
    normal = vt.T
    hidden = null_space(geometry)
    direction = hidden @ (hidden.T @ score)
    derivative = np.linalg.norm(direction)
    direction /= derivative
    goal = edge_or_geometry(prepare(channel, UNIFORM)[0], effects, mode)
    runs = []
    for radius in (1 / 64, 1 / 32):
        for sign in (-1, 1):

            def weights(v):
                return UNIFORM + TANGENT @ (sign * radius * direction + normal @ v)

            def residual(v):
                return (
                    edge_or_geometry(prepare(channel, weights(v))[0], effects, mode)
                    - goal
                )

            fit = root(residual, np.zeros(len(goal)), tol=1e-10)
            # Deposit an exactly normalized route lottery reproducible on integer tickets.
            tickets = np.rint(weights(fit.x) * 2**45).astype(np.int64)
            tickets[-1] += 2**45 - int(tickets.sum())
            w = tickets / 2**45
            error = np.max(abs(residual(fit.x)))
            if error > 1e-11 or w.min() < 0.005 or abs(w.sum() - 1) > 1e-14:
                raise ArithmeticError(f"failed {mode} match: {error}, {w.min()}")
            runs.append((sign * radius, w))
    return runs, dict(
        projected_entropy_slope=derivative,
        geometry_singular_values=singular.tolist(),
        derivative_step_agreement=float(abs(rows - confirm).max()),
    )


def appointment():
    """Fixed four-unit unmeasured flight, four station loads + heads and central energy current."""
    board = station.iso.pairing.formation.federation(1)
    route = station.iso.pairing.formation.region(board)
    launch = bridge.beam.indices(sum(bridge.beam.screens(board, route), []))
    u0 = np.zeros((4 * board.n, 40))
    u0[launch] = np.eye(40)
    maps = []
    for odd in (False, True):
        local = []
        for charge in (-1, 1):
            h = native_precision.generator(board, route, charge, odd)
            # Short subintervals keep the unitary Gram residual below the audit tolerance.
            op = -0.125j * h
            u = u0
            for _ in range(16):
                u = expm_multiply(op, u, traceA=op.diagonal().sum())
            np.testing.assert_allclose(u.conj().T @ u, np.eye(40), atol=2e-13, rtol=0)
            cells = bridge.operators(board, route, charge, odd)
            bridge.audit(board, cells, route, h)
            flux = 1j * bridge.current(board, cells, route[1:3]) / 8
            uu = u.reshape(board.n, 4, 10, 4)
            inside = uu[list(route)].reshape(16, 10, 4).transpose(1, 0, 2)
            outside = (
                np.delete(uu, route, axis=0).transpose(0, 2, 1, 3).reshape(-1, 4, 4)
            )
            current = (u.conj().T @ (flux @ u)).reshape(10, 4, 10, 4)
            current = np.einsum("aiaj->ij", current) / 10
            local.append((inside, outside, current))
        maps.append(local)
    return maps


def local_readout(maps, state):
    blocks = []
    flux = 0.0
    for local, rho in zip(maps, state, strict=True):
        erased = np.zeros((4, 4), complex)
        for inside, outside, current in local:

            def reduce(amps):
                return (
                    np.einsum("aik,kl,ajl->ij", amps, rho, amps.conj(), optimize=True)
                    / 20
                )

            blocks.append(reduce(inside))
            erased += reduce(outside)
            flux += np.trace(current @ rho).real / 2
        blocks.append(erased)
    if abs(sum(np.trace(b).real for b in blocks) - 1) > 1e-12:
        raise ArithmeticError("incomplete regional memory record")
    return entropy(blocks), float(flux)


def acquire(channel, effects, maps, weight):
    state, contraction = prepare(channel, weight)
    ranges, lengths, geometry, total = record(state, effects)
    marg = 2 * np.einsum("kij,bji->k", BASIS, state).real
    regional, current = local_readout(maps, state)
    # Replay from opposite pure input words independently of the fixed point solve.
    phi = np.einsum("d,bdij->bij", weight, channel)
    residual = 0.0
    for block, f, p in zip(state, phi, ODDS, strict=True):
        for i in (0, 3):
            start = np.diag(np.eye(4)[i]).ravel()
            for _ in range(320):
                start = f @ start
            residual = max(
                residual, float(np.linalg.norm(start.reshape(4, 4) - block / p))
            )
    return dict(
        routes=weight.tolist(),
        conditional_contraction=contraction,
        pure_start_320_cycle_residual=residual,
        ranges=ranges.tolist(),
        edges=lengths.tolist(),
        area=float(geometry[0]),
        volume=float(geometry[1]),
        closed_flight_entropy=total + float(np.log(20)),
        four_station_memory_entropy=regional,
        central_bridge_energy_current=current,
        parity_blind_paulis=marg.tolist(),
    )


def load_native():
    bank, _ = reuse.channels()
    nums, unit, _ = native_precision.exact_effects()
    pair = np.array(nums, float) / float(unit)
    kernels = pair[:, :, 0] + 1j * pair[:, :, 1]
    effects = 3 * np.eye(4) + (kernels + kernels.conj().transpose(0, 1, 3, 2)) / 2
    return bank.mean(axis=0), effects


def study(channel, effects, maps):
    reference = acquire(channel, effects, maps, UNIFORM)
    results = {}
    for mode in ("area_volume", "six_edges"):
        candidates, report = endpoints(channel, effects, mode)
        report["endpoints"] = [
            dict(route_coordinate=radius, **acquire(channel, effects, maps, w))
            for radius, w in candidates
        ]
        low, high = report["endpoints"][-2:]
        # The highest expected contrast selects the axis in advance of drawing the two cohorts.
        axis = int(
            np.argmax(
                abs(np.array(high["parity_blind_paulis"]) - low["parity_blind_paulis"])
            )
        )
        n = 2**35
        probabilities = [(1 + s["parity_blind_paulis"][axis]) / 2 for s in (low, high)]
        rng = np.random.default_rng(2026112501 if mode == "six_edges" else 2026112500)
        counts = rng.binomial(n, probabilities).tolist()
        report["independent_pauli_read"] = dict(
            axis=AXES[axis],
            trials_per_source=n,
            plus_counts=counts,
            plus_probabilities=probabilities,
            simultaneous_99pct_hoeffding_halfwidth=float(
                np.sqrt(np.log(400) / (2 * n))
            ),
        )
        results[mode] = report
    return dict(
        ports=[0, 1, 38, 39],
        odd_weight=float(ODDS[1]),
        ordinary_energy=(ODDS[0] - ODDS[1]) / 2,
        charges=[-1, 1],
        directed_station_indices=[list(p) for p in reuse.ROUTES],
        axes=list(AXES),
        reference=reference,
        comparisons=results,
    )


def run():
    channel, effects = load_native()
    return study(channel, effects, appointment())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(receipt)
    else:
        print(receipt, end="")
