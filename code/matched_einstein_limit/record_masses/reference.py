"""Spectroscopy of the guide's two rest branches against one parked native load."""

import argparse
import json
from pathlib import Path

import numpy as np
from model.candidate_laws.interval_echo.openwave.row_21_boost.moving_clock import (
    apparatus as lab,
    preparation,
)
from model.candidate_laws.interval_echo.openwave.row_21_boost.moving_clock.refinement import (
    measure as guide,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.matched_probe.signal_time import (
    signal,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**36


def decode(tape):
    temporal = np.asarray(tape["ports"])
    if not (temporal.sum(axis=-1) == SHOTS).all():
        raise ValueError("Missing output of the recombination")
    # Per-bank bounds budgeted for this complete run with a simple Hoeffding
    # rectangle. The real and imaginary differences each use two yes counts.
    eps = np.sqrt(np.log(2 * 24000 / 0.01) / (2 * SHOTS))
    p = temporal[:, :, 0] / SHOTS
    z = p[:, 0] - p[:, 2] + 1j * (p[:, 3] - p[:, 1])
    disk = 2 * np.sqrt(2) * eps
    if min(abs(z)) <= disk:
        return dict(resolved=False)
    err = np.arcsin(disk / abs(z))
    increments = np.angle(z[1:] / z[:-1])
    if np.any(abs(increments) + err[1:] + err[:-1] >= np.pi):
        return dict(resolved=False)
    total = float(sum(increments))
    dt = tape["hand"][1] - tape["hand"][0]
    if dt <= 0:
        raise ValueError("Invalid local hand")
    phase = [total - err[0] - err[-1], total + err[0] + err[-1]]
    clock = np.array(tape["clock_yes"]) / SHOTS
    lo = 64 * np.arccos(np.clip(2 * (clock + eps) - 1, -1, 1)).sum()
    hi = 64 * np.arccos(np.clip(2 * (clock - eps) - 1, -1, 1)).sum()
    return dict(
        resolved=True,
        phase_per_guide_tick=(np.array(phase) / dt).tolist(),
        principal_phase=float(np.angle(z[-1] / z[0])),
        turns=round((total - np.angle(z[-1] / z[0])) / (2 * np.pi)),
        rate_per_empty_clock=[phase[0] / hi, phase[1] / lo],
    )


def run():
    rng = np.random.default_rng(319208)
    board = lab.federation(2)
    routes = guide.guides(board)
    empty, _ = signal.model(0, 0)
    odd, _ = signal.model(0, 2)
    np.testing.assert_array_equal(odd - empty, 4 * np.eye(len(odd)))
    reader = lab.gates(1.0, 8.0, holds=(np.pi / 2, np.pi / 2))[0]
    out = {}
    for n in guide.SIZES:
        origin = routes[n]["ports"][0]
        spare = routes[n]["vacant"][0]
        assert spare not in routes[n]["ports"] and spare in board.neighbours[origin]
        dwell = np.array([np.pi, np.pi - 9 * np.pi / n])
        rows = {}
        for contrast in (8, 4):
            micro = lab.gates(1, contrast, holds=dwell / 64)
            full = lab.gates(1, contrast, holds=dwell)
            common = -np.exp(-1j * (4 * contrast + 0.5) * sum(dwell))
            for sign in (1, -1):
                spin = lab.spin(full, 0, sign, 0, common)
                target = np.tile(spin, (n, 1)) / np.sqrt(n)
                assert preparation.build(target, 1)[1] < 1e-12
                # Identical quarter-split on the spare and emission port.
                # The known constant launch phase disappears in phase increments.
                split = reader[:, 0]
                target *= split[0]
                parked = split[1]
                banks = []
                for i in range(129):
                    bins = []
                    for q in guide.PHASES:
                        a, b = abs(reader @ [target[0, 0], q * parked]) ** 2
                        bins.append(
                            lab.draw(np.array([a, b, 1 - a - b]), SHOTS, rng).tolist()
                        )
                    banks.append(bins)
                    if i < 64:
                        target = target @ micro[0].T
                    elif i < 128:
                        target = lab.step(target, [np.eye(2), micro[1]])
                clock = rng.binomial(SHOTS, (1 + np.cos(dwell / 64)) / 2)
                tape = dict(
                    ports=banks,
                    clock_yes=clock.tolist(),
                    hand=lab.reference_hand(n, [0, 1]),
                )
                rows[f"{contrast}/{sign}"] = dict(tape=tape, reading=decode(tape))
        # Outside/empty outcomes are retained when the reference is removed.
        bare = 1 / (2 * n)
        control = [
            lab.draw(np.array([bare / 2, bare / 2, 1 - bare]), SHOTS, rng).tolist()
            for q in guide.PHASES
        ]
        out[str(n)] = dict(rows=rows, no_reference_control=control)
    return dict(shots=SHOTS, seed=319208, guides=out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = run()
    if parser.parse_args().write:
        (HERE / "reference.json").write_text(
            json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n"
        )
    for n, x in data["guides"].items():
        print(n, {key: z["reading"] for key, z in x["rows"].items()})
