"""Local stock calibration and matching of wall patches in the archived terminal records."""

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "words.npz"
WIDTHS = (2, 3, 4, 5, 6)


def digits(code, width):
    return tuple((code >> k) & 1 for k in reversed(range(width)))


def window_codes(bits, width):
    n = len(bits)
    return [
        sum(bits[(i + j) % n] << (width - 1 - j) for j in range(width))
        for i in range(n)
    ]


def hist(codes, width):
    return np.bincount(codes, minlength=2**width)


def calibrations(n, width):
    """Integer identity per directed edge; tests all real lookup tables at once."""
    m, size = width - 1, 2**width
    if n < 2 * m + 1:
        raise ValueError("not enough zero padding for the two closed stocks")
    unit = np.eye(size, dtype=int)
    vertices = size // 2
    stocks = []
    targets = []
    prefix = []
    for v in range(vertices):
        held = digits(v, m)
        register = 0
        walk = []
        for b in held:
            register = register * 2 + b
            walk.append(register)
        prefix.append(hist(walk, width))
    for e in range(size):
        bits = digits(e, width) + (0,) * (n - width)
        stocks.append(hist(window_codes(bits, width), width))
        targets.append(2 * sum(bits) - n)
    for e in range(size):
        u, v = e >> 1, e & (vertices - 1)
        np.testing.assert_array_equal(
            unit[e] - (prefix[v] - prefix[u]), stocks[e] - stocks[v] + unit[0]
        )

    # Tables below use twice the ordinary load unit.
    canonical = np.array([sum(digits(e, width)[:2]) - 1 for e in range(size)])
    gradient = np.array([(e & 1) - ((e >> 1) & 1) for e in range(size)])
    nonlinear = (
        np.array([int((e & 3) == 3) - int(((e >> 1) & 3) == 3) for e in range(size)])
        if width >= 3
        else np.zeros(size, dtype=int)
    )
    alternating = [window_codes((0, 1), width)[j] for j in (0, 1)]
    trials = {}
    for label, table in (
        ("normalized_two_end", canonical),
        ("redistributed_local", canonical + 3 * nonlinear),
        ("raw_two_end", 2 * canonical),
        ("alternating_bulk_response", canonical + gradient),
    ):
        totals = np.array(stocks) @ table
        trials[label] = dict(
            failed_total_stocks=int(np.count_nonzero(totals != targets)),
            alternating_outputs_twice=table[alternating].tolist(),
        )
    return dict(size=size, edge_identities=size, diagnostics=trials)


def masks(words, n, sign):
    """Histogram equality for each pair is checked on physical, cyclic port words."""
    left = {n - 5, n - 4, n - 3}
    right = {2, 3, 4}
    selected = np.zeros((len(WIDTHS), len(words)), bool)
    gated = np.zeros(len(words), bool)
    for i, code in enumerate(words):
        bits = tuple(((int(code) >> j) & 1) ^ int(sign < 0) for j in range(n))
        aliases = [j for j in range(n) if bits[j] == bits[(j + 1) % n]]
        if len(aliases) != 2:
            continue
        a, b = aliases
        if not ((a in left and b in right) or (b in left and a in right)):
            continue
        if any(2 * bits[j] - 1 != sign for j in aliases):
            raise ValueError("same-sign gate has the wrong wall")
        gated[i] = True
        for wi, width in enumerate(WIDTHS):
            patches = [{(j - s) % n for s in range(width - 1)} for j in aliases]
            if patches[0] & patches[1]:
                continue
            codes = window_codes(bits, width)
            lists = [[codes[j] for j in p] for p in patches]
            if sorted(lists[0]) != sorted(lists[1]):
                raise AssertionError("isolated walls have different local words")
            alt = set(window_codes((0, 1), width))
            if any(
                codes[j] not in alt for j in set(range(n)) - patches[0] - patches[1]
            ):
                raise AssertionError("nonalternating bulk escaped both patches")
            selected[wi, i] = True
    return gated, selected


def decode(path=DATA):
    out = {}
    with np.load(path, allow_pickle=False) as banks:
        for n in (14, 16):
            for sign in (-1, 1):
                gate, selected = masks(banks[f"words_{n}"], n, sign)
                for level in (0, 1, 2):
                    for amount in (1, 2):
                        for kind in ("native", "contrast_off"):
                            name = f"counts_{level}_{n}_{amount}_{sign}_{kind}"
                            if name not in banks:
                                continue
                            counts = banks[name]
                            totals = counts.sum(axis=1).tolist()
                            if totals != [2**22] * 7 or not np.issubdtype(
                                counts.dtype, np.unsignedinteger
                            ):
                                raise ValueError("invalid enrollment")
                            out[name] = dict(
                                gated=counts[:, gate].sum(axis=1).tolist(),
                                qualified_by_width=[
                                    counts[:, g].sum(axis=1).tolist() for g in selected
                                ],
                            )
        keys = {k for k in banks if k.startswith("counts_")}
        if set(out) != keys or len(keys) != 32:
            raise ValueError("incomplete released-word cohorts")
    return dict(
        widths=WIDTHS,
        original_starts_per_time=2**22,
        calibration_by_size={
            n: {w: calibrations(n, w) for w in WIDTHS} for n in (14, 16)
        },
        released_banks=out,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--output", type=Path, help="write the reduction of archived word counts"
    )
    args = cli.parse_args()
    result = json.dumps(decode(), separators=(",", ":")) + "\n"
    if args.output:
        args.output.write_text(result)
    elif result != (HERE / "counts.json").read_text():
        raise SystemExit("Deposited reduction differs")
    else:
        print("All 224 release cohorts and 248 calibration identities reproduced.")
