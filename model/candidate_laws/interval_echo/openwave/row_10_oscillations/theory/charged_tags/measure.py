"""Acquire all sharp words on the charged daughter, donor and neutral input bench."""

import argparse
import json
from functools import cache
from itertools import product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from model.candidate_laws.interval_echo.openwave.row_10_oscillations.joint_receivers.measure import (
    sample,
)
from model.candidate_laws.interval_echo.openwave.row_10_oscillations.joint_receivers.readback import (
    interval,
)
from .apparatus import gate, inverse, ladder, stations, WORDS

HERE = Path(__file__).resolve().parent
LAW = HERE / "acquisition_law.json"
STARTS = 2**18
SEED = 26112694
READS = ("sharp", "local_half_swap", "joint", "donor_read", "order_read")
SOURCES = ("02", "20", "phase02", "phase20")
ANGLES = (0, np.pi / 2, np.pi)


def separate(rho, labels):
    return rho * np.equal.outer(labels, labels)


@cache
def instruments(k, step):
    blank = ladder.free_seam(k)
    half = np.kron(expm(-0.5j * np.pi * blank), np.eye(3))
    angle = ANGLES[step]
    convert = gate(angle, 0, k) @ gate(angle, 2, k)
    restore = inverse(angle, 2, k) @ inverse(angle, 0, k)
    return half, convert, restore


def evolve(rho, u):
    return u @ rho @ u.conj().T


def trace_distance(first, second):
    return float(sum(abs(np.linalg.eigvalsh(first - second))) / 2)


def predictions(k, step):
    half, convert, restore = instruments(k, step)
    donors = [w[2] for w in WORDS]
    first = [w[0] for w in WORDS]
    readings = {}
    daughter = []
    for label in SOURCES:
        word = (0, 2, 1) if label.endswith("02") else (2, 0, 1)
        ket = np.eye(27, dtype=complex)[:, WORDS.index(word)]
        if label.startswith("phase"):
            ket = half @ ket
        prepared = np.outer(ket, ket.conj())
        full = evolve(prepared, convert)
        if label.startswith("phase"):
            daughter.append(np.trace(full.reshape(9, 3, 9, 3), axis1=1, axis2=3))
        read = {
            "sharp": full,
            "local_half_swap": evolve(full, half),
            "joint": evolve(full, half @ restore),
            "donor_read": evolve(separate(full, donors), half @ restore),
            "order_read": evolve(separate(prepared, first), half @ restore @ convert),
        }
        for key, rho in read.items():
            probabilities = np.real(np.diag(rho))
            assert min(probabilities) > -1e-12 and abs(sum(probabilities) - 1) < 1e-11
            readings[label, key] = np.maximum(probabilities, 0)
    return readings, trace_distance(*daughter)


def internal_tag(sites, draw):
    # The two carrier rails now share the same signed stock in both inputs.
    # Restore the charged daughter to the blank entrance before the two pulses.
    half, convert, _ = instruments(8, 2)
    swap = half @ half
    banks = []
    for site, sign, reverse, panel in product(
        sites, (-1, 1), (False, True), ("sharp", "half_swap")
    ):
        entrance = (1 + sign, 1 - sign, 1)
        col = convert[:, WORDS.index(entrance)]
        if reverse:
            col = swap @ col
        prepared = half @ col
        detected = half @ prepared if panel == "half_swap" else prepared
        p = abs(detected) ** 2
        assert (
            sum(
                prob
                for w, prob in zip(WORDS, p, strict=True)
                if w[0] + w[1] - 2 == sign
            )
            > 1 - 1e-12
        )
        key = [site["sensor"], sign, reverse, panel]
        banks.append(dict(key=key, counts=draw(["internal", *key], [*p, 0])))
    return banks


def acquire(law=None, used=None):
    """Draw every bank; a stored acquisition law fixes the probabilities of the draws.

    Each freshly computed outcome probability must agree with the stored binary64
    value within 1e-12. The seeded draws then use the stored values, so a different
    roundoff in the last bit cannot reorder near-equal bins or shift the stream.
    """
    rows = None if law is None else iter(law)

    def draw(key, probabilities):
        p = [float(x) for x in probabilities]
        if rows is not None:
            row = next(rows, None)
            if (
                row is None
                or row["key"] != key
                or len(row["probabilities"]) != len(p)
                or max(abs(a - b) for a, b in zip(p, row["probabilities"])) > 1e-12
            ):
                raise ArithmeticError(f"{key} differs from the acquisition law")
            p = row["probabilities"]
        if used is not None:
            used.append(dict(key=key, probabilities=p))
        return sample(key, p, STARTS, seed=SEED)

    sites = stations()
    banks = []
    distances = []
    for k, step in product((8, 4), range(3)):
        readings, dist = predictions(k, step)
        distances.append(dict(strength=k, step=step, daughter_trace_distance=dist))
        for site, (source, read) in product(sites, readings):
            key = [site["sensor"], k, step, source, read]
            banks.append(dict(key=key, counts=draw(key, [*readings[source, read], 0])))
    control = internal_tag(sites, draw)
    if rows is not None and next(rows, None) is not None:
        raise ArithmeticError("Unused rows in the acquisition law")
    return dict(
        seed=SEED,
        starts=STARTS,
        stations=sites,
        words=WORDS,
        final_bin="refused local bench setup",
        banks=banks,
        computed=distances,
        same_charge_control=control,
    )


def read(record):
    banks = {tuple(r["key"]): r["counts"] for r in record["banks"]}
    expected = set(
        product(("opposite", "second_owner"), (8, 4), range(3), SOURCES, READS)
    )
    if (
        set(banks) != expected
        or len(record["banks"]) != len(expected)
        or any(len(x) != 28 or sum(x) != STARTS or min(x) < 0 for x in banks.values())
    ):
        raise ValueError("incomplete calibration")
    rows = []
    for site, k, step in product(("opposite", "second_owner"), (8, 4), range(3)):

        def cases(label, read, keep):
            c = banks[site, k, step, label, read]
            return sum(n for w, n in zip(WORDS, c[:-1], strict=True) if keep(w))

        def score(read):
            return [
                cases(label, read, lambda w: w == (2, 0, 1)) for label in SOURCES[2:]
            ]

        def contrast(c):
            lo, hi = interval(c[0], STARTS)
            l2, h2 = interval(c[1], STARTS)
            return [c, (c[0] - c[1]) / STARTS, [lo - h2, hi - l2]]

        rows.append(
            dict(
                site=site,
                k=k,
                step=step,
                charge_clicks=[
                    cases(label, "sharp", lambda w: w[0] + w[1] != 2)
                    for label in SOURCES
                ],
                phase={r: contrast(score(r)) for r in READS[1:]},
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    p.add_argument(
        "--fresh-law",
        type=Path,
        help="draw from freshly computed probabilities and write them to this path",
    )
    args = p.parse_args()
    if args.acquire:
        used = []
        law = None if args.fresh_law else json.loads(LAW.read_text())
        args.acquire.write_text(
            json.dumps(acquire(law, used), separators=(",", ":"), allow_nan=False)
            + "\n"
        )
        if args.fresh_law:
            args.fresh_law.write_text(
                json.dumps(used, separators=(",", ":"), allow_nan=False) + "\n"
            )
    else:
        print(
            json.dumps(read(json.loads((HERE / "counts.json").read_text())), indent=2)
        )
