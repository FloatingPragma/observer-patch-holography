"""Compare late parity sorting to early native material preparation and a bare detector."""

import argparse
import json
from math import asin, atan2, ceil, log, sqrt

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from . import apparatus as lab
from .measure import HERE, SHOTS, sample


def acquire():
    rng = np.random.default_rng(271217)
    banks = []
    for level, owner in [(-1, 0), (0, 0), (0, 13), (1, 0), (1, 13)]:
        car = lab.board(level)
        for m in (1, 2):
            front = lab.source(level, owner, 0, m)
            back = np.zeros_like(front)
            local = slice(12 * owner, 12 * (owner + 1))
            back[local] = front[local][list(lab.geometry()[1])]
            launch = lab.prepare_load(car, (front - back) / np.sqrt(2))
            for t in (1, 3):
                h = -1j * t * lab.free(level, 1)
                wave = expm_multiply(h, launch)[local]
                raw = [
                    sample(
                        rng, lab.endpoint_probabilities(level, owner, wave, e, False)
                    )
                    for e in lab.cell().seams
                ]
                banks.append(
                    dict(
                        level=level,
                        owner=owner,
                        winding=m,
                        time=t,
                        source="odd",
                        edges=raw,
                    )
                )
                if level == -1 and t == 1:
                    other = lab.prepare_load(car, back)
                    f = expm_multiply(h, front)
                    b = expm_multiply(h, other)
                    raw = []
                    for edge in lab.cell().seams:
                        fp = lab.endpoint_probabilities(level, owner, f, edge, False)
                        bp = lab.endpoint_probabilities(level, owner, b, edge, False)
                        raw.append(sample(rng, (fp + bp) / 2))
                    banks.append(
                        dict(
                            level=level,
                            owner=owner,
                            winding=m,
                            time=t,
                            source="mixture",
                            edges=raw,
                        )
                    )
    clocks = []
    a, b = lab.cell().seams[0]
    inputs = []
    for name, amps in (("sharp", [1, 0]), ("phase", np.array([1, -1j]) / sqrt(2))):
        goal = np.zeros(12, complex)
        goal[[a, b]] = amps
        inputs.append((name, lab.prepare_load(lab.cell(), goal, root=a)[[a, b]]))
    for rate in (1, 0.5):
        _, h = lab.local_source(rate, lab.candidate.KAPPA)
        for t in (0, 0.5, 1, 2, 3, 4, 5, 6):
            u = expm(-1j * t * h)
            # Left-minus-right for sharp left gives cos(rate*t).
            # Plus quadrature left-minus-right gives sin(rate*t).
            for input_name, v in inputs:
                clocks.append(
                    dict(
                        rate=rate,
                        time=t,
                        preparation=input_name,
                        counts=sample(rng, abs(u @ v) ** 2),
                    )
                )
    return dict(seed=271217, shots=SHOTS, banks=banks, clocks=clocks)


def compare(reference, check):
    from .readback import decode

    decode(reference)
    if (
        check["seed"] != 271217
        or check["shots"] != SHOTS
        or len(check["banks"]) != 22
        or len(check["clocks"]) != 32
    ):
        raise ValueError("Missing source or clock acquisition")
    # In each of the two files, include all individual endpoint percentages
    # for all panels as well as all ring means, clock quadratures and sharp reads.
    num = 47 * (3 * 30 + 12) + 72 * 4 + 60 * 2
    extra = 22 * (3 * 30 + 12) + 32

    def width(terms, n, range_size):
        return (
            ceil(
                1e8 * (range_size * sqrt(log(2 * terms / 0.005) / (2 * n)) + 8 / 2**36)
            )
            / 1e8
        )

    p_tol = width(num, SHOTS, 1) + width(extra, SHOTS, 1)
    j_tol = width(num, 5 * SHOTS, 2) + width(extra, 5 * SHOTS, 2)
    hand_tol = width(extra, SHOTS, 2)
    source = {
        (r["level"], r["owner"], r["winding"], r["time"]): r
        for r in reference["circulation"]
        if r["parity"] and r["source"] == 0
    }
    tests = []
    seen_keys = set()
    for r in check["banks"]:
        key = r["level"], r["owner"], r["winding"], r["time"]
        key_source = (*key, r["source"])
        if (
            key_source in seen_keys
            or (
                r["source"] != "odd"
                and (r["source"] != "mixture" or key[0] != -1 or key[3] != 1)
            )
            or key[2] not in (1, 2)
            or key[3] not in (1, 3)
        ):
            raise ValueError("Unexpected source enrollment")
        seen_keys.add(key_source)
        ref = source[key]
        x = np.array(r["edges"])
        y = np.array(ref["edges"])
        if (
            x.shape != (30, 3)
            or not np.issubdtype(x.dtype, np.integer)
            or min(x.flat) < 0
            or not np.all(x.sum(axis=1) == SHOTS)
        ):
            raise ValueError("Incomplete detector table")
        diff = (x - y) / SHOTS
        directed = {}
        current = {}
        for i, (a, b) in enumerate(reference["bond_labels"]):
            directed[a, b] = diff[i, 1] - diff[i, 0]
            directed[b, a] = -directed[a, b]
            current[a, b] = (x[i, 1] - x[i, 0]) / SHOTS
            current[b, a] = -current[a, b]

        def circ(table):
            return np.array(
                [
                    sum(table[a, b] for a, b in zip(c, c[1:] + c[:1], strict=True)) / 5
                    for c in reference["cycles"]
                ]
            )

        gaps = abs(circ(directed))
        j = circ(current)
        tests.append(
            [
                *key,
                r["source"],
                float(j[0]),
                float(j[1]),
                int((abs(diff) > p_tol).sum()),
                int((gaps > j_tol).sum()),
                float(gaps.max()),
            ]
        )
    hands = {}
    for rate in (1, 0.5):
        ordered = []
        for t in (0, 0.5, 1, 2, 3, 4, 5, 6):
            cs = []
            for q in ("sharp", "phase"):
                (row,) = [
                    r
                    for r in check["clocks"]
                    if (r["rate"], r["time"], r["preparation"]) == (rate, t, q)
                ]
                a, b = row["counts"]
                if (
                    type(a) is not int
                    or type(b) is not int
                    or min(a, b) < 0
                    or a + b != SHOTS
                ):
                    raise ValueError("Invalid clock")
                cs.append((a - b) / SHOTS)
            ordered.append(atan2(cs[1], cs[0]))
        phases = np.unwrap(ordered)
        hands[str(rate)] = phases.tolist()
    return dict(
        joint_endpoint_width=p_tol,
        joint_circulation_width=j_tol,
        individual_hand_width=hand_tol,
        hand_angular_radius=asin(sqrt(2) * hand_tol),
        clocks=hands,
        detector_comparisons=tests,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = HERE / "source_counts.json"
    result = acquire()
    if args.write:
        path.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    elif result != json.loads(path.read_text()):
        raise ValueError("Source check records differ")
    print(
        json.dumps(
            compare(json.loads((HERE / "records.json").read_text()), result), indent=1
        )
    )


if __name__ == "__main__":
    main()
