"""Transfer the measured modular source to the extensive equilibrium energy."""

import argparse
from hashlib import sha256
from itertools import product
import json

from mpmath import iv
import numpy as np

from . import thermal_pair as exp

HERE = exp.HERE
kernel = exp.kernel
STEP = 0.00001


def box(value, error):
    return iv.mpf([value - error, value + error])


def trace_data(h, sites, beta, heat, energy=False, archived=None):
    if archived is None:
        w, err = kernel.baseline(h, beta, energy=energy)
        updates, cert = kernel.single_updates(h, sites, beta=beta, energy=energy)
        baseline = float(w.sum())
        totals = [{u: float(v.sum()) for u, v in row.items()} for row in updates]
    else:
        baseline = sum(archived["ordinary_base_shells"])
        err = archived["baseline_per_address_error"]
        totals = [
            {int(u): sum(v) for u, v in row.items()}
            for row in archived["single_increments"]
        ]
        cert = archived["increment_error"]
    n = h.shape[0]
    base = box(baseline, err * n + abs(baseline) * 1e-10)
    marks = [
        {
            u: box(v, n * cert["tail"] + abs(v) * 2 * cert["rounding_relative"])
            for u, v in row.items()
        }
        for row in totals
    ]
    tail = n * heat["upper_each_diagonal"]
    if energy:
        # A positive power series bounds its derivative by its value at beta+h
        # divided by h. Remove the uniform diagonal shift used in that series.
        step = beta / 16
        large = kernel.mixed_bounds(h, sites, heat["distance"], beta + step)
        tail = (
            cert["shift"] * tail
            + n * np.exp(cert["shift"] * step) * large["upper_each_diagonal"] / step
        )
    residual = box(0.0, tail) if energy else iv.mpf([0, tail])
    return base, marks, residual


def log_mix(trace):
    z, marks, r = trace
    return sum(
        iv.log1p((z * r - x * y) / ((z + x) * (z + y))) / 16
        for x, y in product(*[row.values() for row in marks])
    )


def energy_mix(trace, read):
    z, marks, r = trace
    e, increments, er = read
    total = 0
    for a, b in product((32, 33), repeat=2):
        x, y = marks[0][a], marks[1][b]
        u, v = increments[0][a], increments[1][b]
        # Separately normalized aperture identity, applied to energy numerators.
        val = (
            x
            * y
            / (z + x + y)
            * ((e / z - u / x) / (z + x) + (e / z - v / y) / (z + y))
        )
        # Add the genuine connected heat coefficient, including its derivative.
        val += (er - r * (e + u + v) / (z + x + y)) / (z + x + y + r)
        total += val / 16
    return total


def meter(trace, energy, lower, upper, delta):
    z, marks, r = trace
    e, increments, er = energy
    outcomes = {}
    error = 0.0
    for a, b in product(exp.STRENGTHS, repeat=2):
        prob = (
            (
                e
                + increments[0].get(a, 0)
                + increments[1].get(b, 0)
                + (er if a * b else 0)
            )
            / (z + marks[0].get(a, 0) + marks[1].get(b, 0) + (r if a * b else 0))
            - lower
        ) / (upper - lower)
        lo, hi = exp.interval_output(prob)
        assert 0 <= lo <= hi <= 1
        error = max(error, (hi - lo) / 2)
        outcomes[f"{-a},{-b}"] = (lo + hi) / 2
    bright = exp.motion.binary_responses(outcomes)
    table = np.array([[48, 176], [112, 112]])
    native = table.ravel() / table.sum()
    unrelated = np.outer(table.sum(1), table.sum(0)).ravel() / float(table.sum()) ** 2
    probs = dict(
        zip(exp.LABELS, bright, strict=True),
        native=native @ bright,
        product=unrelated @ bright,
    )
    rng = np.random.default_rng(256193)
    counts = {
        tag: int(rng.binomial(exp.BATCH, p, size=exp.BATCHES).sum())
        for tag, p in probs.items()
    }
    band = np.sqrt(exp.contact.RISK / (2 * exp.SHOTS)) + error
    mixed = (counts["00"] + counts["11"] - counts["01"] - counts["10"]) / exp.SHOTS
    held = (counts["product"] - counts["native"]) / exp.SHOTS
    return dict(
        energy_range=[lower, upper],
        numerical_probability_error=error,
        bright=counts,
        from_charge=(
            exp.probe.mix_bands(delta, mixed + np.array([-4 * band, 4 * band]))
            * (upper - lower)
        ).tolist(),
        held=((held + np.array([-2 * band, 2 * band])) * (upper - lower)).tolist(),
    )


def run():
    iv.dps = 40
    p = HERE / "thermal_pair.json"
    pair = json.loads(p.read_bytes())
    h, sites, *_ = exp.motion.hamiltonian(exp.motion.cold.federation(0), pair["edges"])
    d = pair["heat"]["distance"]
    beta = kernel.BETA
    delta = np.array(pair["source_delta"])
    trace = trace_data(h, sites, beta, pair["heat"], archived=pair)
    energy = trace_data(h, sites, beta, pair["heat"], energy=True)
    q = iv.mpf(delta.tolist())
    F = -q * log_mix(trace) / beta
    U = q * energy_mix(trace, energy)
    lower, upper = kernel.bounds(h)
    lower -= 33
    direct = meter(trace, energy, lower, upper, delta)
    logs = []
    for b in (beta - STEP, beta + STEP):
        heat = kernel.mixed_bounds(h, sites, d, b)
        logs.append(log_mix(trace_data(h, sites, b, heat)))
    # Third cumulant of each bounded energy is bounded by the cube of its span.
    # Four parity pairs of weight 1/16 and four terms per mixed difference.
    derivative = -(logs[1] - logs[0]) / (2 * STEP) + box(
        0.0, (upper - lower) ** 3 * STEP**2 / 6
    )
    check = q * derivative
    a, b = exp.interval_output(U)
    for lo, hi in (exp.interval_output(check), direct["held"], direct["from_charge"]):
        assert lo <= b and a <= hi
    return dict(
        sha256={
            k.name: sha256(k.read_bytes()).hexdigest()
            for k in (p, HERE / "THERMAL_ENERGY.md", HERE / "pair_contact.json")
        },
        isothermal_free_energy=exp.interval_output(F),
        ordinary_internal_energy=exp.interval_output(U),
        reversible_heat_into_ordinary=exp.interval_output(U - F),
        from_temperature_difference=exp.interval_output(check),
        log_trace_mixed=exp.interval_output(log_mix(trace)),
        energy_mixed=exp.interval_output(energy_mix(trace, energy)),
        shots=exp.SHOTS,
        energy_detector=direct,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "thermal_energy.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        exp.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    print({k: v for k, v in out.items() if k != "sha256"})
