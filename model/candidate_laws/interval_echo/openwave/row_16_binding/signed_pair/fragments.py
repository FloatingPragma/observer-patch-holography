"""Sign-blind regional allocation followed by simultaneous physical charge reads."""

import argparse
import gzip
from io import BytesIO
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from model.executable_model import carrier
from .material import Assembly
from .aperture import flights
from .measure import HERE, normalize
from .probe import Quadrupole, draw, calibrate
from .standards import phase_standards

SHOTS = 2**38


def regional_meter():
    """Relative phase from one pointer flips the sector marker into a Fourier bin."""
    q = np.arange(-4, 5)
    _, phase_code = phase_standards()
    k = np.arange(31)
    delta = (130 * q[:, None] - k[None, :]) % 31
    frequencies = phase_code[2] ** q[:, None] * np.exp(-2j * np.pi * k / 31)[None, :]
    phase = (frequencies[:, :, None] ** k[None, None, :]).mean(axis=2)
    expect = (delta == 0).astype(float)
    if np.max(abs(phase - expect)) > 1e-12:
        raise ArithmeticError("Regional pointer loses a signed charge")
    return {str(s): int(np.argmax(abs(row))) for s, row in zip(q, phase)}


def assignments(w):
    total = int(w[0].sum())
    assert np.all(w.sum(axis=1) == total)
    signed = range(-4, 5)
    calibration = regional_meter()
    bins = np.array([calibration[str(q)] for q in signed])
    charge_for_bin = {code: q for code, q in zip(bins, signed)}
    bits = 1 << np.arange(12)
    mask = abs(w) @ bits
    neighbours = carrier.cell().neighbours

    def connected(indices):
        seen = {indices[0]}
        for _ in indices:
            seen.update(v for v in indices if any(p in seen for p in neighbours[v]))
        return len(seen) == len(indices)

    triads = {sum(bits[list(t)]): t for t in combinations(range(12), 3) if connected(t)}
    tetrads = {
        sum(bits[list(t)]): t for t in combinations(range(12), 4) if connected(t)
    }
    weight = np.zeros((4, len(w)))
    for m in np.unique(mask):
        ids = np.flatnonzero(mask == m)
        partitions = [
            (t, tetrads[m ^ u])
            for u, t in triads.items()
            if (u & m) == u and (m ^ u) in tetrads
        ]
        for trimer, tetramer in partitions:
            # Two independent pointer circuits contact their own addresses.
            left = bins[4 + w[ids][:, trimer].sum(axis=1)]
            right = bins[4 + w[ids][:, tetramer].sum(axis=1)]
            for lc in np.unique(left):
                for rc in np.unique(right[left == lc]):
                    q1, q2 = charge_for_bin[lc], charge_for_bin[rc]
                    if q1 + q2 != total:
                        raise ArithmeticError(
                            "Independent charge readings lost total material"
                        )
                    weight[(q1 + 3) // 2, ids[(left == lc) & (right == rc)]] += 1 / len(
                        partitions
                    )
    return weight


def run():
    # Independent draws of the recorded words condition on the original measured category counts.
    word_rng = np.random.default_rng(161404)
    full_banks = []
    rng = np.random.default_rng(161402)
    gain, cal = calibrate(rng)
    result = {
        "original_shots": SHOTS,
        "regional_bins": regional_meter(),
        "charge_pairs": [[q, 1 - q] for q in (-3, -1, 1, 3)],
        "calibration": cal,
        "runs": [],
    }
    for k in (8, 0):
        a = Assembly(4, 3, k)
        seeds = a.sources()[:, [0, 2]]
        out, _ = flights(a.h, a.bracket, seeds, 209 if k else 0, 8, (1,))
        w = assignments(a.w)
        reader = Quadrupole(a.w + 1, abs(a.w), a.axes[:1], 0, gain)
        for i, name in enumerate(("sharp", "residue_1")):
            p = abs(out[1, :, i]) ** 2
            q = w @ p
            passed = p.sum()
            counts = draw(rng, [*q, passed - q.sum(), 1 - passed], SHOTS)
            # Each row of the joined bank shares one selected regional charge
            # category with the random-address shape probe from that launch.
            banks = []
            scale, accept, active, host = reader.effects["centroid_0"]
            assert np.max(abs(host)) == 0
            for polarity in (1, -1):
                probabilities = []
                word_weights = []
                for f in w:
                    v = p * f
                    a_gate = v @ accept
                    s = polarity * (v @ active)
                    word_weights.extend(
                        (
                            v * (accept + polarity * active) / 2,
                            v * (accept - polarity * active) / 2,
                            v * (1 - accept),
                        )
                    )
                    probabilities.extend(
                        [(a_gate + s) / 2, (a_gate - s) / 2, v.sum() - a_gate]
                    )
                probabilities += [passed - q.sum(), 1 - passed]
                banks.append(draw(rng, probabilities, SHOTS))
                word_weights.append(p * np.maximum(1 - w.sum(axis=0), 0))
                words = []
                for cohort, dist in zip(banks[-1][:-1], word_weights):
                    if cohort:
                        dist = np.maximum(dist, 0)
                        words.append(draw(word_rng, dist / dist.sum(), cohort))
                    else:
                        words.append(np.zeros(len(p), dtype=int))
                full_banks.append(words)
            result["runs"].append(
                {
                    "contrast": k,
                    "source": name,
                    "charge_and_refusals": counts,
                    "weight_scale": float(scale),
                    "joint_probe": banks,
                }
            )
    return normalize(result), np.array(full_banks, dtype="<i8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result, words = run()
    archive = BytesIO()
    np.save(archive, words, allow_pickle=False)
    packed = archive.getvalue()
    if args.output:
        args.output.with_name("joined_words.npy.gz").write_bytes(
            gzip.compress(packed, compresslevel=4, mtime=0)
        )
        args.output.write_text(json.dumps(result, indent=1) + "\n")
    else:
        assert gzip.decompress((HERE / "joined_words.npy.gz").read_bytes()) == packed
        assert result == json.loads((HERE / "fragments.json").read_text())
        print(
            "Replayed both regional charge readings, signed quadrupoles and all refused assignments."
        )
