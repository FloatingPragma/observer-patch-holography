"""Detect signed common-mode and gradient responses of incident and outgoing walls."""

import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix, eye
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
from .. import wall_flight, walls
from ..contact_flight.measure import contact_masks
from . import pickup

HERE = Path(__file__).resolve().parent
SHOTS = 2**26
CLOCK = 2**26
QUANTA = 2**34
TIMES = np.arange(41)
CASES = ((8, False), (8, True), (0, False))
BITS = np.array(list(product((0, 1), repeat=10)))
LOADS = np.stack((2 * BITS, 2 - 2 * BITS), axis=2).reshape(1024, 20)
CONTACTS = ((3, 4), (7, 0), (9, 10), (15, 16), (11, 12), (17, 18))


def pickups():
    q = np.stack([LOADS[:, a] + LOADS[:, b] - 2 for a, b in CONTACTS], axis=1)
    d = np.stack([LOADS[:, b] - LOADS[:, a] for a, b in CONTACTS], axis=1)
    if not np.array_equal(q * q + d * d, np.full(q.shape, 4)):
        raise ArithmeticError("A dipole or common-mode code fails")
    return q, d


def initial(anchors):
    states = np.zeros((1024, 2, 2, 2))
    for r, g, wall in product(range(2), repeat=3):
        ref = walls.reference(r)
        seed = tuple(bit ^ (wall and i >= 2) for i, bit in enumerate(ref))
        sample = walls.loads(walls.reference(g))
        modes = [walls.wing_ref(8, sample[a], False, 3)[2] for a in anchors]
        indices = np.all(BITS[:, :4] == seed, axis=1)
        states[indices, r, g, wall] = np.kron(*modes)
    return states.reshape(1024, 8)


def flight(discriminator=False):
    _, anchors, _ = wall_flight.extended_geometry()
    states = initial(anchors)
    output = np.zeros((3, 2, 2, 2, len(TIMES), 1024), dtype=np.int64)
    for branch, (k, cut) in enumerate(CASES):
        h = walls.operators(k, cut, 3)
        vals, vec = eigh(h, driver="evd")
        coefficients = vec.T @ states
        phases = np.exp(-1j * np.outer(TIMES, vals))
        for col in range(8):
            amplitudes = (phases * coefficients[:, col]) @ vec.T
            prob = abs(amplitudes) ** 2
            if np.max(abs(prob.sum(axis=-1) - 1)) > 2e-11:
                raise ArithmeticError("Original-start probability loss")
            cdf = np.rint(np.cumsum(prob, axis=-1) * QUANTA).astype(np.int64)
            if np.max(abs(cdf[:, -1] - QUANTA)) != 0:
                raise ArithmeticError("Acquisition binary roundoff")
            output[branch].reshape((8, len(TIMES), 1024))[col] = cdf
        if discriminator and branch == 0:
            # Independent sparse propagation on both original starts and the
            # held inverted-ring input at fixed guide occupation.
            centered = csr_matrix(h) - np.mean(h.diagonal()) * eye(1024, format="csr")
            direct = expm_multiply(-32j * centered, states[:, [0, 1, 5]])
            nominal = output[branch].reshape(8, len(TIMES), 1024)[[0, 1, 5], 32]
            np.testing.assert_allclose(
                abs(direct.T) ** 2,
                np.diff(nominal, axis=-1, prepend=0) / QUANTA,
                rtol=0,
                atol=1.5 / QUANTA,
            )
    return output


def records(cdf):
    rng = np.random.default_rng(831208)
    probabilities = np.diff(cdf, axis=-1, prepend=0) / QUANTA
    counts = np.array(
        [rng.multinomial(SHOTS, p) for p in probabilities.reshape(-1, 1024)],
        dtype=np.int32,
    ).reshape(cdf.shape)
    cs = np.stack((np.cos(TIMES), np.sin(TIMES)), axis=-1)
    hand = rng.binomial(CLOCK, (1 + cs) / 2).astype(np.int32)
    return counts, hand


def decode(counts, hand, cdf):
    _, anchors, _ = wall_flight.extended_geometry()
    _, _, _, distal = contact_masks(anchors)
    q, d = pickups()
    # Values in [-1,1]; 20000 enrolled probability or signed-mean comparisons.
    error = np.sqrt(2 * np.log(2 * 20000 / 0.01) / SHOTS) + 1024 / QUANTA
    t = np.unwrap(np.arctan2(*(2 * hand.astype(float) / CLOCK - 1).T[::-1]))
    times = dict(read=t.tolist(), bound=float(np.arcsin(np.sqrt(2) * error)))
    source = {}
    for r, wall in product(range(2), repeat=2):
        bits = tuple(z ^ (wall and i >= 2) for i, z in enumerate(walls.reference(r)))
        idx = np.flatnonzero(np.all(BITS[:, :4] == bits, axis=1))[0]
        source[f"r{r}_wall{wall}"] = dict(
            common_mode=q[idx, :2].tolist(), gradient=d[idx, :2].tolist()
        )
    for row in source.values():
        # BRIGHT (20) on first read is first-stage index five. Conditional
        # second-stage BRIGHT and DARK are indices 14 and 10.
        raw = {}
        for mode in ("common_mode", "gradient"):
            events = np.zeros((2, 17), dtype=np.int64)
            for j, charge in enumerate(row[mode]):
                events[j, {0: 5, 2: 14, -2: 10}[charge]] = SHOTS
            raw[mode] = events.tolist()
        row["source_pointer_counts_per_cohort"] = raw
    result = dict(
        times=times, source=source, single_mean_error=float(error), settings={}, at32={}
    )
    # Supplementary energy accounting uses the terminal-load basis at release.
    collision_edges = [
        (1, 2),
        (3, 4),
        (5, 6),
        (7, 0),
        (anchors[0], 8),
        (anchors[1], 14),
    ]
    contact_energy = 32.5 * sum(
        (LOADS[:, a] != LOADS[:, b]).astype(int) for a, b in collision_edges
    )
    attachment_energy = []
    for g in (0, 1):
        expected = (
            np.diff(cdf[0, :, g, :, 0, :], axis=-1, prepend=0) @ contact_energy / QUANTA
        )
        recorded = counts[0, :, g, :, 0, :] @ contact_energy / SHOTS
        attachment_energy.append(
            dict(
                guide=g,
                relative_expected=(expected - expected[g, 0]).tolist(),
                relative_sample=(recorded - recorded[g, 0]).tolist(),
            )
        )
    result["launch_energy_relative_to_matching_sham"] = attachment_energy
    for i, (k, cut) in enumerate(CASES):
        for g in (0, 1):
            key = f"k{k}_cut{int(cut)}_g{g}"
            row = {}
            for j, name in ((2, "near"), (4, "far")):
                positive = (d[:, j] == 2) & (d[:, j + 1] == -2)
                negative = (d[:, j] == -2) & (d[:, j + 1] == 2)

                def mean(banks):
                    return banks[..., positive].sum(axis=-1) - banks[..., negative].sum(
                        axis=-1
                    )

                bank = counts[i, :, g].astype(np.int64)
                m = mean(bank) / SHOTS
                # ring-0 wall-minus-sham minus ring-1 wall-minus-sham
                change = m[0, 1] - m[0, 0] - m[1, 1] + m[1, 0]
                radius = 4 * error
                peak = int(np.argmax(abs(change)))
                found = np.flatnonzero(abs(change) > radius)
                law = mean(np.diff(cdf[i, :, g], axis=-1, prepend=0)) / QUANTA
                at32 = []
                for r, wall in product(range(2), repeat=2):
                    c = bank[r, wall, 32]
                    at32.append(
                        dict(
                            ring=r,
                            wall=wall,
                            positive=int(c[positive].sum()),
                            negative=int(c[negative].sum()),
                            dark_or_other=int(c[~(positive | negative)].sum()),
                        )
                    )
                row[name] = dict(
                    earliest=int(found[0]) if found.size else None,
                    peak=peak,
                    interval=[
                        float(change[peak] - radius),
                        float(change[peak] + radius),
                    ],
                    measured_peak=float(change[peak]),
                    peak_law=float(
                        law[0, 1, peak]
                        - law[0, 0, peak]
                        - law[1, 1, peak]
                        + law[1, 0, peak]
                    ),
                    signs_at_32=at32,
                )
            result["settings"][key] = row
            for r, wall in product(range(2), repeat=2):
                bank = counts[i, r, g, wall, 32]
                # Every category, including equal-sign distal events and refusals.
                signed = [
                    int(bank[distal & (d[:, 4] == a) & (d[:, 5] == b)].sum())
                    for a, b in product((-2, 2), repeat=2)
                ]
                result["at32"][f"{key}_r{r}_wall{wall}"] = dict(
                    distal=int(bank[distal].sum()),
                    far_gradient_counts=signed,
                    other=int(bank[~distal].sum()),
                )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--independent", action="store_true")
    args = parser.parse_args()
    cdf = flight(args.independent)
    calibration = dict(
        routes=pickup.routing(),
        standards=pickup.calibration(),
        full_dimer=pickup.full_dimer_calibration(),
    )
    if args.output:
        args.output.mkdir(exist_ok=False)
        counts, hand = records(cdf)
        np.savez_compressed(
            args.output / "counts.npz", counts=counts, hand=hand, cdf=cdf
        )
        material = decode(counts, hand, cdf)
        (args.output / "results.json").write_text(json.dumps(material, indent=2) + "\n")
        (args.output / "calibration.json").write_text(
            json.dumps(calibration, separators=(",", ":")) + "\n"
        )
        return
    with np.load(HERE / "counts.npz") as stored:
        fixed = stored["cdf"]
        np.testing.assert_allclose(
            cdf / QUANTA, fixed / QUANTA, atol=2 / QUANTA, rtol=0
        )
        counts, hand = records(fixed)
        np.testing.assert_array_equal(counts, stored["counts"])
        np.testing.assert_array_equal(hand, stored["hand"])
        report = decode(counts, hand, fixed)
    compare(
        report, json.loads((HERE / "results.json").read_text()), "signed wall flight"
    )
    compare(
        calibration,
        json.loads((HERE / "calibration.json").read_text()),
        "native pickup",
    )
    print(
        "Verified ternary charge standards, signed launch and remote collision records."
    )


if __name__ == "__main__":
    main()
