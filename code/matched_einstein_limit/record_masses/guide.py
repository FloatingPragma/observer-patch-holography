"""Rest gap and longitudinal kinetic mass from two local clock and range tapes."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_21_boost.moving_clock import (
    apparatus as lab,
    current,
    preparation,
    readback,
)
from model.candidate_laws.interval_echo.openwave.row_21_boost.moving_clock.refinement import (
    measure as guide,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**44
SETTINGS = (0, 1, 2, -1, -2)


def acquire(n, route, rng):
    holds = (np.pi, np.pi - 2 * (4.5 * np.pi / n))
    layers = lab.gates(1, 8, holds=holds)
    alternate = lab.gates(1, 4, holds=holds)
    gates_error = 0.0
    for h, u, v in zip(holds, layers, alternate, strict=True):
        gates_error = max(gates_error, float(abs(v - np.exp(16j * h) * u).max()))
    common = -np.exp(-1j * 32.5 * sum(holds))
    times = n // 5
    stations = [n // 3, 2 * n // 3]
    assert min(min(j, n - 1 - j) for j in stations) > times
    tapes = []
    for i, site in enumerate(stations):
        turn = 2
        echo = lab.radar(n, site, site + 1, layers, SHOTS, rng, turn=turn)
        echo_hand = lab.reference_hand(
            n, [echo[k] for k in ("emission", "turn", "receipt")]
        )
        hand = lab.reference_hand(n, list(range(times + 1)))
        offset = 71 * (i + 1)
        tapes.append(
            dict(
                hand=[t + offset for t in hand],
                echo_hand=[t + offset for t in echo_hand],
                echo_hits=echo["forward_then_return"],
                early=echo["early_clicks"],
                parked=echo["parked_turn"],
                shots=SHOTS,
                reads={},
            )
        )
    prep_error = 0.0
    for tag in SETTINGS:
        k = tag * 4.5 * np.pi / (8 * n)
        phase = np.exp(1j * k * np.arange(n))[:, None] / np.sqrt(n)
        arms = np.stack(
            [phase * lab.spin(layers, k, sign, 0, common) for sign in (1, -1)]
        )
        np.testing.assert_allclose(np.vdot(*arms), 0, atol=1e-13)
        for v in [*arms, *[(arms[0] + q * arms[1]) / np.sqrt(2) for q in guide.PHASES]]:
            prep_error = max(prep_error, preparation.build(v, 1)[1])
        # Read adjacent-site momentum before propagation. Roll only chooses the
        # local addresses for the unchanged route meter.
        for tape, site in zip(tapes, stations, strict=True):
            p = current.route_meter(np.roll(arms[0], -site, axis=0), 0, 1)
            tape["reads"][str(tag)] = dict(
                spatial=[lab.draw(np.array(row), SHOTS, rng).tolist() for row in p],
                temporal=[],
            )
        for _ in range(times + 1):
            for tape, site in zip(tapes, stations, strict=True):
                a, b = arms[:, site, 0]
                tape["reads"][str(tag)]["temporal"].append(
                    [
                        int(rng.binomial(SHOTS, abs(a + q * b) ** 2 / 2))
                        for q in guide.PHASES
                    ]
                )
            arms = lab.step(arms, layers)
    return tapes, dict(
        max_preparation_error=prep_error,
        max_gate_phase_residual=gates_error,
        route_ports=len(route["ports"]),
        analytic_gap=4.5 * np.pi / n,
    )


def interval(fun, *pairs):
    vals = [fun(*args) for args in product(*pairs)]
    return [float(min(vals)), float(max(vals))]


def decode(tape):
    shot = tape["shots"]
    if (
        sum(tape["echo_hits"]) != shot
        or sum(tape["parked"]) != shot
        or tape["echo_hits"][-1] <= 0
        or tape["early"]
        or tape["parked"][-1]
    ):
        raise ValueError("Unresolved neighboring direct echo")
    ticks = np.asarray(tape["hand"])
    start, turn, end = tape["echo_hand"]
    if (
        not (0 < turn - start == end - turn)
        or not np.all(np.diff(ticks) == np.diff(ticks)[0])
        or ticks[-1] <= ticks[0]
    ):
        raise ValueError("Invalid clock appointments")
    span = ticks[-1] - ticks[0]
    length = (end - start) / 2
    reads = {}
    for tag, record in tape["reads"].items():
        if len(record["temporal"]) != len(ticks):
            raise ValueError("Missing temporal appointment")
        time = readback.phase_bins(record["temporal"], shot)
        if not time["resolved"]:
            raise ValueError("Temporal fringe is unresolved")
        angle, da = current.spatial(record["spatial"], shot)
        E, dE = time["elapsed"] / (2 * span), time["phase_error"] / (2 * span)
        reads[tag] = dict(
            E=E,
            e=[E - dE, E + dE],
            p=angle / length,
            p_interval=[(angle - da) / length, (angle + da) / length],
        )
    zero = reads["0"]
    if not zero["p_interval"][0] <= 0 <= zero["p_interval"][1] or zero["e"][0] <= 0:
        raise ValueError("No positive rest half-gap at zero resolved momentum")
    branches = {}
    for sign in (1, -1):
        rows = [reads[str(i * sign)] for i in (1, 2)]
        k2 = [sorted(v**2 for v in row["p_interval"]) for row in rows]
        if (
            any(lo <= 0 <= hi for lo, hi in (row["p_interval"] for row in rows))
            or k2[0][1] >= k2[1][0]
        ):
            raise ValueError("Spatial fringe does not distinguish the two kicks")
        # The finite secants use E(k)-E(0). The extrapolation is affine in
        # each energy. Simple interval arithmetic also covers the common E(0).
        sec = [
            interval(lambda a, b, x: 2 * (a - b) / x, row["e"], zero["e"], sq)
            for row, sq in zip(rows, k2, strict=True)
        ]
        inv = interval(
            lambda small, big, x, y: small + (small - big) * x / (y - x), *sec, *k2
        )
        if inv[0] <= 0:
            raise ValueError("Inertial response unresolved")
        mass = [1 / inv[1], 1 / inv[0]]
        branches[str(sign)] = dict(
            kinetic_mass=mass,
            kinetic_over_gap=interval(lambda m, e: m / e, mass, zero["e"]),
        )
    return dict(rest_half_gap=zero["e"], readings=reads, senses=branches)


def run():
    rng = np.random.default_rng(193071)
    board = lab.federation(2)
    routes = guide.guides(board)
    out = {}
    for n in guide.SIZES:
        route = routes[n]
        assert all(
            b in board.neighbours[a]
            for a, b in zip(
                route["ports"], route["ports"][1:] + route["ports"][:1], strict=True
            )
        )
        tapes, audit = acquire(n, route, rng)
        d = audit["analytic_gap"]
        k = np.array([d / 8, d / 4])
        e = np.arccos(np.cos(d) * np.cos(k))
        sec = 2 * (e - d) / k**2
        central = sec[0] + (sec[0] - sec[1]) / 3
        audit.update(
            exact_curvature_mass=float(np.tan(d)),
            finite_kinetic_mass=float(1 / central),
            relative_bias=float(1 / (central * np.tan(d)) - 1),
        )
        decoded = [decode(t) for t in tapes]
        out[str(n)] = dict(tapes=tapes, decoded=decoded, audit=audit)
    return dict(shots=SHOTS, seed=193071, guides=out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = run()
    path = HERE / "guide.json"
    if args.write:
        path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    for k, v in data["guides"].items():
        print(k, v["decoded"][0]["senses"], v["audit"])
