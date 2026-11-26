"""All-start phase discrimination across two stations with three finite resources."""

import argparse
from functools import cache
import json
from math import ceil, log, pi, sqrt
from pathlib import Path

import numpy as np

from ...signed_effects.sequential import grid_counts
from . import analyzer
from . import apparatus as a

HERE = Path(__file__).resolve().parent
SHOTS = 131072
SEED = 26032615
MENUS = (
    ("native", 0),
    ("native", 1),
    ("native", 2),
    ("path_read", 2),
    ("exchange_only", 2),
)


def decode_word(bits, size=None):
    if bits is None:
        return None  # An explicitly recorded detector erasure.
    if (
        not isinstance(bits, str)
        or not bits
        or len(bits) % 2
        or (size is not None and len(bits) != 2 * size)
        or any(c not in "01" for c in bits)
    ):
        raise ValueError("Malformed detector word")
    pairs = [bits[j : j + 2] for j in range(0, len(bits), 2)]
    if any(p not in ("00", "01") for p in pairs):
        return None
    return sum(1 << i for i, p in enumerate(pairs) if p == "00")


def received_word(payload, size):
    if payload is None:
        return None
    if not isinstance(payload, list) or len(payload) != size:
        raise ValueError("Malformed message packet")
    symbols = [decode_word(raw, 2) for raw in payload]
    if any(value not in (1, 2) for value in symbols):
        return None  # A valid classical bit requires exactly one occupied pointer.
    return sum((value - 1) << i for i, value in enumerate(symbols))


@cache
def bob_program(m, word, kappa):
    target, choices = analyzer.receiver(m, word)
    return a.programmed(target, kappa), choices


def bank(menu, refs, angle, rng):
    m = refs + 1
    dim = 1 << m
    kappa = 0 if menu == "exchange_only" else a.candidate.KAPPA
    # The source-route permutation has been checked on every physical port.
    # Contract the three disjoint pair splits only after that routing check.
    split = a.contact(2, 0, 1, pi / 2, kappa)[:, 1]
    state = np.zeros((dim, dim), complex)
    for x in range(dim):
        state[x, dim - 1 - x] = split[1] ** x.bit_count() * split[2] ** (
            m - x.bit_count()
        )
    state = a.phase(m, 0, angle, kappa) @ state
    fourier = a.programmed(analyzer.fourier(m), kappa)
    records = []
    probabilities = []
    if menu == "path_read":
        # A failed signal-path read is an enrolled refusal. Keep its raw word.
        branches = {}
        for x in range(dim):
            raw_tag = a.word_read(0, 1, x & 1)
            tag = decode_word(raw_tag, 1)
            branches.setdefault(tag if tag is not None else raw_tag, []).append(x)
    else:
        branches = {None: range(dim)}
    for bit, paths in branches.items():
        projected = np.zeros_like(state)
        projected[paths] = state[paths]
        if menu == "path_read" and bit not in (0, 1):
            records.append([bit, None, None, None, None])
            probabilities.append(float(np.vdot(projected, projected).real))
            continue
        left = fourier @ projected
        for alice in range(dim):
            near = a.word_read(0, m, alice)
            sent = decode_word(near, m)
            payload = None
            if sent is not None:
                payload = [a.announce((sent >> i) & 1) for i in range(m)]
            message = received_word(payload, m)
            if message is None:
                records.append([bit, near, payload, None, None])
                probabilities.append(float(np.vdot(left[alice], left[alice]).real))
                continue
            # Bob acts on his received word, including a misdelivered valid bit.
            read, choices = bob_program(m, message, kappa)
            output = read @ left[alice]
            for port in range(dim):
                far = a.word_read(1, m, port)
                terminal = decode_word(far, m)
                guess = choices[terminal] if terminal is not None else None
                records.append([bit, near, payload, far, guess])
                probabilities.append(float(abs(output[port]) ** 2))
    counts = grid_counts(rng, SHOTS, probabilities)
    return dict(
        menu=menu,
        references=refs,
        source_sign=int(angle == 0),
        records=[[*record, int(n)] for record, n in zip(records, counts, strict=True)],
    )


def acquire():
    rng = np.random.default_rng(SEED)
    return dict(
        seed=SEED,
        starts_per_bank=SHOTS,
        geometry=a.layout(),
        banks=[bank(menu, ref, angle, rng) for menu, ref in MENUS for angle in (0, pi)],
    )


class ForecastMismatch(ValueError):
    """Complete enrolled data disagree with a frozen forecast."""


def count_bank(b):
    menu, r = b["menu"], b["references"]
    sign = b["source_sign"]
    if (
        type(r) is not int
        or (menu, r) not in MENUS
        or type(sign) is not int
        or sign not in (0, 1)
    ):
        raise ValueError("Unknown preparation or detector menu")
    size = r + 1
    total = hits = refusals = 0
    hist = np.zeros(1 << size, dtype=int)
    for row in b["records"]:
        if not isinstance(row, (list, tuple)) or len(row) != 6:
            raise ValueError("Incomplete transcript row")
        tag, near, payload, far, reported, count = row
        if type(count) is not int or count < 0 or total + count > SHOTS:
            raise ValueError("Invalid enrolled count")
        total += count
        if menu == "path_read":
            if type(tag) is not int or tag not in (0, 1):
                # A path refusal retains the offending raw word or explicit erasure.
                if decode_word(tag, 1) is not None or any(
                    item is not None for item in (near, payload, far)
                ):
                    raise ValueError("Invalid path refusal transcript")
        elif tag is not None:
            raise ValueError("Unexpected path tag")
        alice = decode_word(near, size)
        if alice is not None:
            hist[alice] += count
        elif payload is not None:
            raise ValueError("Message following a refused source word")
        message = received_word(payload, size)
        if message is None and far is not None:
            raise ValueError("Terminal read without a received setting")
        terminal = decode_word(far, size)
        guess = None
        if terminal is not None:
            # Reconstruct the decision independently from the enrolled rule.
            # Its second one-particle (or second single-hole) exit means minus.
            bob_load = size - message.bit_count()
            if bob_load in (0, size):
                guess = int(bob_load == size)
            else:
                minus = 2 if bob_load == 1 else ((1 << size) - 1) ^ 2
                guess = int(terminal != minus)
        if reported is not None and type(reported) is not int or reported != guess:
            raise ValueError("Published guess differs from the recorded decision")
        if guess is None:
            refusals += count
        elif guess == sign:
            hits += count
    if total != SHOTS:
        raise ValueError("Incomplete enrolled bank")
    return hits, refusals, hist


def forecast(menu, r, sign):
    if menu == "exchange_only":
        return 5 / 8 if sign else 3 / 8
    if menu == "path_read":
        return 0.5
    return (0.5, 0.75, (1 + 1 / sqrt(2)) / 2)[r]


def score(receipt):
    if receipt["starts_per_bank"] != SHOTS:
        raise ValueError("Missing enrolled starts")
    banks = {}
    for b in receipt["banks"]:
        stats = count_bank(b)
        key = b["menu"], b["references"], b["source_sign"]
        if key in banks:
            raise ValueError("Duplicate bank")
        banks[key] = stats
    if set(banks) != {(menu, r, s) for menu, r in MENUS for s in (0, 1)}:
        raise ValueError("Missing preparation or detector menu")
    # One interval for each per-sign accuracy and each word probability. The
    # resulting coverage also covers their equal-prior means and differences.
    bins = 2 * sum(1 + (1 << (r + 1)) for _, r in MENUS)
    width = ceil(1e9 * (sqrt(log(2 * bins / 0.01) / (2 * SHOTS)) + 512 / 2**40)) / 1e9
    scores = {}
    for key, (hits, _, hist) in banks.items():
        if abs(hits / SHOTS - forecast(*key)) > width:
            raise ForecastMismatch(f"Per-sign accuracy outside forecast for {key}")
        if np.max(abs(hist / SHOTS - 1 / len(hist))) > width:
            raise ForecastMismatch(f"Alice's marginal outside forecast for {key}")
    for menu, r in MENUS:
        keys = [(menu, r, s) for s in (0, 1)]
        scores[menu, r] = sum(banks[k][0] for k in keys) / (2 * SHOTS)
        target = sum(forecast(*k) for k in keys) / 2
        if abs(scores[menu, r] - target) > width:
            raise ForecastMismatch(
                f"Equal-prior accuracy outside forecast for {(menu, r)}"
            )
    return scores, width


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    out = json.loads(json.dumps(acquire()))
    path = HERE / "records.json"
    if args.write:
        with path.open("x") as stream:
            stream.write(json.dumps(out, indent=1) + "\n")
    elif json.loads(path.read_text()) != out:
        raise SystemExit("Phase comparison differs from archive")
    scores, width = score(out)
    for b in out["banks"]:
        key = b["menu"], b["references"], b["source_sign"]
        hits, refusals, _ = count_bank(b)
        print(
            f"{key}: {hits}/{SHOTS} correct, {refusals} refusals, target {forecast(*key):.8f}"
        )
    for key, measured in scores.items():
        print(f"{key} pooled: {measured:.8f} (±{width:.8f})")


if __name__ == "__main__":
    main()
