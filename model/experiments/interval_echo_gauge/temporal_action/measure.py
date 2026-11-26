"""Evaluate the native-law source and action, then sample the calculated fringes."""

import argparse
import gzip
import json
from math import log, sqrt
from pathlib import Path

import numpy as np

from . import action
from .material import geometry, reader
from .rails import Meter, WIDTH, fringe, hardware

HERE = Path(__file__).resolve().parent
SHOTS = 2**20
SEED = 662901


def constraints(hist):
    d, c = geometry()[:2]
    e, b, a, phi, j, rho = (hist[k] for k in ("E", "B", "A", "phi", "J", "rho"))
    magnet = a @ c.T
    u = np.sum(e**2, axis=1) + np.sum(magnet[:-1] * magnet[1:], axis=1)
    values = dict(
        electric=e + 2 * np.diff(a, axis=0) + phi @ d.T,
        magnetic=b - magnet[:-1],
        gauss=e @ d - rho,
        continuity=2 * np.diff(rho, axis=0) + j[:-1] @ d,
        faraday=2 * np.diff(magnet, axis=0) + e @ c.T,
        ampere=2 * np.diff(e, axis=0) - magnet[1:-1] @ c + j[:-1],
        work=2 * np.diff(u) + np.sum((e[:-1] + e[1:]) * j[:-1], axis=1),
    )
    return {k: int(np.max(abs(v))) for k, v in values.items()}


def audit(meter, hist, settings):
    packet = action.Action(meter, hist, **settings)
    measured = packet.value()
    if action.independent(hist, **settings) != measured:
        raise ArithmeticError(
            "Measured full action differs from its independent definition"
        )
    directions = packet.directions()
    for tag, grid in directions.items():
        expected = action.independently_vary(hist, tag, **settings)
        if expected != grid:
            raise ArithmeticError(f"Action interferometer failed on {tag}")
    return dict(
        action=measured,
        differences=directions,
        nonzero={k: int(np.count_nonzero(v)) for k, v in directions.items()},
    )


def face_layers():
    c = geometry()[1]
    adjacency = abs(c @ c.T)
    distance = np.full(20, -1)
    distance[0] = 0
    for i in range(20):
        for f in np.flatnonzero(distance == i):
            distance[(distance < 0) & (adjacency[f] > 0)] = i + 1
    return [
        [int(f) for f in np.flatnonzero(distance == d)]
        for d in range(max(distance) + 1)
    ]


def counts(histories, rng, reference=None):
    # The fringe uses returned field digits from independent source histories;
    # compose their inverse pulses to reference the idle on the same channel.
    layers = face_layers()
    reads = []
    probabilities = []
    for label, bank in histories.items():
        idle = bank["idle"]
        for case in ("loop", "reverse", "open_quiet", "open_swap"):
            row = bank[case]
            for altered in (False, True):
                for n in range(5):
                    for layer, faces in enumerate(layers):
                        f = faces[0]
                        # Use the data words separately: the receiver never gets a
                        # precomputed subtraction as a phase setting.
                        a = 2 * fringe(row[n, f], action.L, altered) - 1
                        b = 2 * fringe(idle[n, f], action.L, altered) - 1
                        z = complex(*a) * complex(*b).conjugate()
                        p = np.clip((1 + np.array([z.real, z.imag])) / 2, 0, 1)
                        if reference is not None:
                            recorded = reference[len(reads)]
                            np.testing.assert_allclose(p, recorded, rtol=0, atol=2e-13)
                            p = recorded
                        probabilities.append(list(p))
                        bright = rng.binomial(SHOTS, p)
                        reads.append(
                            [label, case, int(altered), n, layer, f, bright.tolist()]
                        )
    if reference is not None and len(reference) != len(reads):
        raise ValueError("Incomplete quadrature calibration bank")
    return reads, probabilities


def run(reference=None):
    meter = Meter()
    d, c, seams, pin, face, hole = geometry()
    source = {}
    scenarios = {}
    flux = {}
    for label, s, m in (
        ("native", 1, 1),
        ("reversed_current", -1, 1),
        ("magnetic_omitted", 1, 0),
    ):
        flux[label] = {}
        for case in action.CASES:
            hist = action.update(meter, case, s, m)
            source[case] = hist["record"]
            resid = constraints(hist)
            if label == "native" and any(resid.values()):
                raise ArithmeticError(
                    "Native returned histories break a field identity"
                )
            row = {
                tag: hist[tag].tolist() for tag in ("A", "E", "B", "J", "rho", "phi")
            }
            row["residuals"] = resid
            if label == "native":
                row["stationarity"] = audit(meter, hist, {})
                if any(row["stationarity"]["nonzero"].values()):
                    raise ArithmeticError("A native source failed stationarity")
            flux[label][case] = hist["B"]
            scenarios[label + "_" + case] = row
    held = action.update(meter, "open_swap")
    gauged = action.gauge(meter, held)
    gauged_packet = audit(meter, gauged, {})
    gchecks = constraints(gauged)
    original = scenarios["native_open_swap"]["stationarity"]["action"]
    if (
        any(gchecks.values())
        or gauged_packet["action"] != original
        or any(gauged_packet["nonzero"].values())
    ):
        raise ArithmeticError("Endpoint-fixed gauge change affected action or field")
    # Off-shell history and the three independent omissions each have the same declared domains.
    altered = held["A"].copy()
    altered[2, 0] += action.L
    checks = {
        "off_shell": audit(meter, {**held, "A": altered}, {}),
        "no_kinetic": audit(meter, gauged, dict(kinetic=0)),
        "no_magnetic": audit(meter, gauged, dict(magnetic=0)),
        "reversed_current": audit(meter, gauged, dict(source_sign=-1)),
        "no_port_charge": audit(meter, gauged, dict(port_sign=0)),
    }
    rng = np.random.Generator(np.random.PCG64(SEED))
    split = int(rng.binomial(SHOTS, source["open_swap"]["half_probabilities"][1]))
    fringes, probabilities = counts(flux, rng, reference)
    num_sampled_rates = 1 + 2 * len(fringes)
    return dict(
        provenance=dict(
            field_words="Deterministic signed-word evaluations of the calibrated reduced native dimer law.",
            source_records="Conditional ternary contact and census words calculated from the native seam law.",
            half_source_original_starts="Seeded binomial sample from the native half-exchange crossing probability.",
            fringes="Seeded binomial bright counts from the calculated X and Y probabilities, including the declared altered-contact controls.",
            sampler="NumPy Generator.binomial with PCG64",
            seed=SEED,
            hoeffding=dict(
                scope="Monte Carlo sampling resolution for the bright and branch frequencies conditional on the calculated probabilities.",
                simultaneous_coverage=0.99,
                number_of_frequencies=num_sampled_rates,
                half_width=sqrt(log(2 * num_sampled_rates / 0.01) / (2 * SHOTS)),
            ),
        ),
        width=WIDTH,
        units=action.L,
        field_words=meter.words,
        max_word=meter.max_magnitude,
        largest_absolute_term_budget=meter.max_budget,
        chart=pin,
        source_face=face,
        hole=hole,
        accumulator=list(hardware()),
        receiver=reader(),
        source_records=source,
        half_source_original_starts=[SHOTS - split, split],
        settings=scenarios,
        gauge=dict(
            packet=gauged_packet,
            residuals=gchecks,
            phi=gauged["phi"].tolist(),
            A=gauged["A"].tolist(),
        ),
        challenges=checks,
        face_layers=face_layers(),
        fringes=fringes,
        sampling_probabilities=probabilities,
        shots_per_quadrature=SHOTS,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    deposited = (
        None
        if args.output
        else json.loads(gzip.decompress((HERE / "records.json.gz").read_bytes()))
    )
    result = run(deposited["sampling_probabilities"] if deposited else None)

    if args.output:
        with args.output.open("xb") as f:
            f.write(
                gzip.compress(
                    json.dumps(result, separators=(",", ":")).encode(), mtime=0
                )
            )
    else:
        if result != deposited:
            raise ArithmeticError("Maxwell source, action or receiver archive differs")
        print(
            "Verified calculated native-law words and probabilities, and replayed the seeded binomial counts."
        )


if __name__ == "__main__":
    main()
