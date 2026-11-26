"""Evaluate the native signed-core evolution, terminal Born probabilities and sampler."""

import argparse
import gzip
import json
from math import sqrt, log
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from ...row_31_orbitals.occupied_write import core
from ...row_31_orbitals.written_transitions.measure import compare
from . import flight

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "records.json.gz"
SEED, SHOTS = 377831, 2**20
RADIUS = sqrt(log(2 * 256 / 0.01) / (2 * SHOTS))


class Sampler:
    def __init__(self, reference=None):
        self.rng = np.random.Generator(np.random.PCG64(SEED))
        self.reference = reference
        self.inputs = {}

    def read(self, key, probabilities):
        p = np.asarray(probabilities, dtype=float)
        np.testing.assert_allclose(p.sum(), 1, atol=2e-9, rtol=0)
        if p.min() < -2e-12 or key in self.inputs:
            raise ArithmeticError("Invalid probability or duplicate acquisition")
        p = np.clip(p, 0, 1)
        total = 2**34
        weights = np.rint(p * total).astype(np.int64)
        weights[np.argmax(p)] += total - weights.sum()
        if self.reference is not None:
            frozen = np.array(self.reference[key])
            np.testing.assert_allclose(frozen / total, p, atol=2e-9, rtol=0)
            weights = frozen
        self.inputs[key] = weights.tolist()
        counts = self.rng.multinomial(SHOTS, weights / total)
        return dict(calculated=p.tolist(), counts=counts.tolist())


def timer(sampler):
    h0, h1 = core.LAW.local_term_parts()
    u, v = [core.LAW.SRC.index(w) for w in ((1, 0), (0, 1))]
    basis = np.eye(99)[:, u]
    out = []
    for j, t in enumerate([0, np.pi / 2, np.pi]):
        state = expm(-1j * t * (h0 + 8 * h1)) @ basis
        probs = abs(state[[u, v]]) ** 2
        np.testing.assert_allclose(probs, [[1, 0], [0.5, 0.5], [0, 1]][j], atol=2e-13)
        out.append(sampler.read(f"timer_{j}", [*probs, 0]))
    # A separate shuttered copy supplies the local timer. Preserve its sharp
    # endpoint history as the record from which the flight appointments are read.
    cols = expm(-1j * np.pi * (h0 + 8 * h1))[:, [u, v]]
    expected = np.eye(99)[:, [v, u]]
    np.testing.assert_allclose(abs(cols) ** 2, expected, rtol=0, atol=3e-13)
    pointer = np.asarray([u, v])
    events = [0]
    for _ in range(max(flight.TICKS)):
        events.append(int(np.argmax(abs(cols[pointer, events[-1]]) ** 2)))
    passes = np.r_[0, np.cumsum(np.diff(events) != 0)]
    appointments = [int(passes[i]) for i in flight.TICKS]
    if appointments != flight.TICKS.tolist():
        raise ArithmeticError("Local exchange pointer missed a half-turn")
    return out, dict(
        endpoint_after_each_half_turn=events, ticks_at_appointments=appointments
    )


def component_probabilities(rho, kappa):
    """Native inverse echo, written magnitude, then the reversal-parity instrument."""
    law = core.LAW
    echo = law.echo()
    r, _ = law._label_ops()
    record = echo.T @ rho @ echo
    outcomes = {0: float(np.trace(record[:9, :9]))}
    data = record[9:, 9:]
    for mag in (0, 1, 2):
        mask = np.diag([abs(d) == mag for c, d in law.LABELS]).astype(float)
        for odd in (False, True):
            energy = kappa * mag**2 + int(odd)
            proj = np.kron(mask @ (np.eye(10) + (-1 if odd else 1) * r) / 2, np.eye(9))
            outcomes[energy] = outcomes.get(energy, 0) + np.sum(proj.T * data)
    h0, h1 = law.local_term_parts()
    np.testing.assert_allclose(
        sum(v * p for v, p in outcomes.items()),
        np.sum(rho.T * (h0 + kappa * h1)),
        atol=1e-12,
    )
    return outcomes


def formation(board, seam, state, sampler):
    # Trace the other material positions and the written coin before the
    # neighbouring seam read. They are orthogonal occupation alternatives.
    law = core.LAW
    core_p = component_probabilities(np.outer(state, state), 8)
    blank = {}
    for j in (0, 1):
        word = [0, 0]
        word[j] = 1
        rho = np.zeros((99, 99))
        rho[0, 0] = rho[law.SRC.index(tuple(word)), law.SRC.index(tuple(word))] = 0.5
        blank[j] = component_probabilities(rho, 8)
    empty = np.zeros((99, 99))
    empty[0, 0] = 1
    pempty = component_probabilities(empty, 8)
    if any(abs(p) > 1e-12 for e, p in pempty.items() if e):
        raise ArithmeticError("Vacant host has a nonzero energy read")
    local = []
    for u, v in board.seams:
        if (u, v) == seam:
            local.append(([u, v], core_p))
        elif u in seam or v in seam:
            local.append(([u, v], blank[0 if u in seam else 1]))
    for p in seam:
        if board.degree[p] < 6:
            local.append(([p], {0: 0.5, 8: 0.5}))
    rows = []
    max_error = 0.0
    measured = 0.0
    for i, (contact, probs) in enumerate(local):
        resolved = sorted(v for v in probs if probs[v] > 1e-12)
        returned = sampler.read(
            f"formation_{board.cells}_{i}", [probs[v] for v in resolved]
        )
        scores = np.array(resolved)
        measured += float(scores @ returned["counts"] / SHOTS)
        # For a three-outcome 0,8,9 score use a union bound on two tails.
        max_error += float(scores[-1] - scores[0]) * RADIUS
        rows.append(dict(contact=contact, energy=resolved, **returned))
    expected = 5 * 8 + (sum(int(board.degree[p]) for p in seam) - 2) / 4
    if abs(measured - expected) > max_error:
        raise ArithmeticError("Independent formation sample failed")
    return dict(
        components=rows,
        host_outcome=0,
        reading=measured,
        simultaneous_sampler_halfwidth=max_error,
        operator_mean=expected,
    )


def fit(panel, appointments):
    x = np.pi * np.asarray(appointments)
    p = np.array([r["counts"][0] for r in panel]) / SHOTS
    train = np.isin(flight.TICKS, [0, 80, 160, 240, 400])

    def profile(freq):
        design = np.array([np.ones(len(x)), np.cos(freq * x), np.sin(freq * x)]).T
        a = np.linalg.lstsq(design[train], p[train], rcond=None)[0]
        forecast = design @ a
        return np.linalg.norm((p - forecast)[train]), forecast

    grid = np.linspace(0.003, 0.008, 101)
    idx = int(np.argmin([profile(f)[0] for f in grid]))
    freq = minimize_scalar(
        lambda f: profile(f)[0],
        bounds=(grid[max(0, idx - 1)], grid[min(100, idx + 1)]),
        method="bounded",
        options={"xatol": 1e-14},
    ).x
    prediction = profile(freq)[1]
    error = abs((p - prediction)[~train])
    return dict(
        angular_frequency=float(freq),
        held_ticks=flight.TICKS[~train].tolist(),
        held_forecasts=prediction[~train].tolist(),
        held_max_error=float(max(error)),
        held_pass=bool(max(error) < 0.025),
    )


def run(reference=None):
    core.integer_blocks()
    state = flight.preparation()
    sampler = Sampler(reference)
    records = {}
    calibrations, local_clock = timer(sampler)
    appointments = local_clock["ticks_at_appointments"]
    for level in (0, 2):
        board, seam, _ = core.chart(level)
        energy = formation(board, seam, state, sampler)
        rows = {}
        for tag, kappa in [("native", 8), ("no_contrast", 0)]:
            positive, audit = flight.clock_trace(
                board, seam, kappa, np.pi * np.asarray(appointments)
            )
            banks = {}
            for source, probabilities in [
                ("positive", positive),
                ("negative", 1 - positive),
                ("dephased", np.full(len(positive), 0.5)),
            ]:
                banks[source] = [
                    sampler.read(f"{board.cells}_{tag}_{source}_{tick}", [p, 1 - p, 0])
                    for tick, p in zip(flight.TICKS, probabilities, strict=True)
                ]
            rows[tag] = dict(banks=banks, ground_audit=audit)
        native = rows["native"]["banks"]["positive"]
        readings = np.array([r["counts"][0] / SHOTS for r in native])
        low = np.flatnonzero(readings < 0.1)
        swing = readings[0] > 0.9 and bool(any(np.max(readings[i:]) > 0.9 for i in low))
        records[str(board.cells)] = dict(
            seam=list(seam),
            formation=energy,
            flight=rows,
            fitted=fit(native, appointments),
            full_swing=bool(swing),
        )
    planted = []
    for j, t in enumerate([0, np.pi / 2, np.pi]):
        h = np.array([[0.5, -0.5], [-0.5, 0.5]])
        state = expm(-1j * t * h)[:, 0]
        planted.append(sampler.read(f"planted_{j}", [*(abs(state) ** 2), 0]))
    if reference is not None and reference.keys() != sampler.inputs.keys():
        raise ArithmeticError("Stored sampler bank has unmatched settings")
    return dict(
        provenance="Calculated native-law probabilities and PCG64-seeded multinomial samples from their deposited dyadic probabilities.",
        seed=SEED,
        shots_per_setting=SHOTS,
        sampler_halfwidth=RADIUS,
        ticks=flight.TICKS.tolist(),
        timer_calibration=calibrations,
        local_clock=local_clock,
        planted_signed_record=planted,
        boards=records,
        sampling_weights_over_2_to_34=sampler.inputs,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    opts = parser.parse_args()
    saved = None if opts.output else json.loads(gzip.decompress(OUTPUT.read_bytes()))
    result = run(saved["sampling_weights_over_2_to_34"] if saved else None)
    if opts.output:
        with opts.output.open("xb") as f:
            f.write(
                gzip.compress(
                    json.dumps(result, separators=(",", ":")).encode(), mtime=0
                )
            )
    else:
        compare(result, saved, "released signed core")
        print(
            "Verified signed-core probabilities, formation components and all sampled records."
        )


if __name__ == "__main__":
    main()
