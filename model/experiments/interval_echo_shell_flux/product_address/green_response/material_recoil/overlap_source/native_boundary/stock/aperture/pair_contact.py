"""Measure the connected head logarithm and the native ordinary-pair interaction."""

import argparse
from fractions import Fraction
from hashlib import sha256
import json

import numpy as np
from scipy.linalg import expm
from scipy.sparse.csgraph import connected_components

from . import contact_timing as timing

HERE = timing.HERE
stock = timing.stock
SHOTS = 2**38
TAU = 1 / 32
RISK = timing.lab.lab.RISK


def source_head_population():
    W, V, _, _ = stock.kernels()
    _, seams, _, e, p, ix, words, co = stock.carrier()
    cc = connected_components(
        stock.graph(W, e, p, ix, co) + stock.graph(V, e, p, ix, co), directed=False
    )[1]
    init = ix[(words == [2, 2, 0, 0, 0, 0, 0]).all(axis=1)].item()
    mem = words[cc == cc[init], 4:6]
    return [list(e) for e in seams[:2]], mem


def source_counts(full=False):
    seams, mem = source_head_population()
    mag = stock.MAG[mem]
    keep = ((mem != 0) & (mag != 1)).all(axis=1)
    if full:
        labels = np.flatnonzero((stock.MAG != 1) & (np.arange(11) != 0))
        codes = np.searchsorted(labels, mem[keep])
        size = len(labels)
    else:
        codes = mag[keep] // 2
        size = 2
    table = np.bincount(codes @ [size, 1], minlength=size**2).reshape(size, size)
    return seams, table, int((~keep).sum()), len(mem)


def blocks():
    _, _, H, H0 = stock.kernels()
    law = stock.LAW
    mats = {}
    for label, operator in (("exchange", H0), ("full", H)):
        result = {}
        kappa = 0 if label == "exchange" else 8
        for word, (a, b) in enumerate(law.SRC):
            expected = np.zeros(99, dtype=int)
            expected[word] += 2 + 4 * kappa * (a - b) ** 2
            expected[law.SRC.index((b, a))] -= 2
            np.testing.assert_array_equal(operator[:, word], expected)
        for d in (0, 2):
            ids = np.array(
                [i + 1 for i, (_, x) in enumerate(law.LABELS) if abs(x) == d]
            )
            vacuum = 9 * ids
            B = operator[np.ix_(vacuum, vacuum)]
            delta = operator[np.ix_(vacuum + 4, vacuum + 4)] - B
            for w in range(9):
                cols = vacuum + w
                expected = np.zeros_like(operator[:, cols])
                expected[cols] = B + (delta if w == 4 else 0)
                np.testing.assert_array_equal(operator[:, cols], expected)
            np.testing.assert_array_equal(B @ delta, delta @ B)
            result[d] = (B, delta)
        mats[label] = result
    return mats


def clock_banks(mats, rng):
    B, D = mats["full"][2]
    phase = expm(1j * TAU * B / 4) @ expm(-1j * TAU * (B + D) / 4)
    out = {}
    for coin in (1, -1):
        for sign in (1, -1):
            vec = np.array([1, sign, coin, sign * coin])
            lam = -(32 + (1 - sign) // 2) if coin == 1 else 0
            np.testing.assert_array_equal(D @ vec, 4 * lam * vec)
            z = vec @ phase @ vec / 4
            out[f"{coin},{sign}"] = [
                int(rng.binomial(SHOTS, np.clip((1 + p) / 2, 0, 1)))
                for p in (z.real, z.imag)
            ]
    return out


def acquire_joint(rng, table, reject, total):
    nat = table.ravel() / table.sum()
    prod = np.outer(table.sum(1), table.sum(0)).ravel() / table.sum() ** 2
    source = {}
    for key in ("reference", "native"):
        source[key] = rng.multinomial(
            SHOTS, np.r_[table.ravel(), reject] / total
        ).tolist()
    # Product preparation draws one retained head from each of two independent accepted copies.
    source["product"] = rng.multinomial(SHOTS, prod).tolist()
    source["product_input_rejects"] = [
        int(rng.negative_binomial(SHOTS, table.sum() / total)) for _ in range(2)
    ]
    return source, nat, prod


def overlap_bands(raw, size=2):
    cells = size**2
    means = {}
    lo, hi = {}, {}
    for tag in ("reference", "native", "product"):
        bank = raw[tag]
        assert (
            sum(bank) == SHOTS
            and min(bank) >= 0
            and len(bank) == (cells if tag == "product" else cells + 1)
        )
        n = sum(bank[:cells])
        p = np.array(bank[:cells]).reshape(size, size) / n
        width = np.sqrt(RISK / (2 * n))
        means[tag] = p
        lo[tag], hi[tag] = p - width, p + width
    p = means["reference"]

    def logproduct(a, b):
        return np.log(a.sum(1)[:, None] * a.sum(0)[None, :] / b)

    assert np.min(lo["reference"]) > 0
    low = logproduct(lo["reference"], hi["reference"])
    high = logproduct(hi["reference"], lo["reference"])
    left = lo["product"] - hi["native"]
    right = hi["product"] - lo["native"]
    box = np.array([x * y for x in (left, right) for y in (low, high)])
    mc = logproduct(p, p)
    if size != 2:
        return dict(
            connected_charge_nats=[float(box.min(0).sum()), float(box.max(0).sum())]
        )
    return dict(
        connected_charge_nats=[float(box.min(0).sum()), float(box.max(0).sum())],
        reference_log_odds=[
            float(np.log(a[0, 0] * a[1, 1] / (b[0, 1] * b[1, 0])))
            for a, b in (
                (lo["reference"], hi["reference"]),
                (hi["reference"], lo["reference"]),
            )
        ],
        estimate=float(np.sum((means["product"] - means["native"]) * mc)),
        product_minus_native_both_two=[float(left[-1, -1]), float(right[-1, -1])],
    )


def decode(out):
    result = overlap_bands(out["source"])
    result["full_support_nats"] = overlap_bands(out["complete_head_banks"], 6)[
        "connected_charge_nats"
    ]
    eps = np.sqrt(2 * RISK / SHOTS)
    zg = []
    ranges = {}
    for name, counts in out["clock"].items():
        assert len(counts) == 2 and min(counts) >= 0 and max(counts) <= SHOTS
        z = complex(*(2 * np.asarray(counts) / SHOTS - 1))
        # The unit circle point lies in a disk of radius sqrt(2)*eps around z.
        dz = np.arcsin(np.sqrt(2) * eps / abs(z))
        angle = np.angle(z)
        assert abs(angle) + dz < np.pi
        ranges[name] = [float((-angle - dz) / TAU), float((-angle + dz) / TAU)]
        zg.append(z)
    # Each sharp magnitude-two head has equal weights in all four joint parity sectors.
    g = sum(zg) / 4
    d = abs(g - 1)
    low = max(0, d - np.sqrt(2) * eps) ** 2 / 2
    high = (d + np.sqrt(2) * eps) ** 2 / 2
    a, b = result["product_minus_native_both_two"]
    result["ordinary_energy_bands"] = ranges
    result["from_source_fringe_shift"] = [
        min(a * low, a * high, b * low, b * high),
        max(a * low, a * high, b * low, b * high),
    ]
    assert all(0 <= n <= SHOTS for n in out["receiver"].values())
    change = (out["receiver"]["product"] - out["receiver"]["native"]) / SHOTS
    result["held_fringe_shift"] = [change - eps, change + eps]
    return result


def run():
    edges, table, reject, total = source_counts()
    mats = blocks()
    rng = np.random.default_rng(537098)
    source, p, q = acquire_joint(rng, table, reject, total)
    clocks = clock_banks(mats, rng)
    B, D = mats["full"][2]
    z = expm(1j * TAU * B / 4) @ expm(-1j * TAU * (B + D) / 4)
    # Independent receiver simulation on both four-dimensional native heads.
    # The common ordinary parking phase is removed by calibrating the two written zeros.
    psi = np.array([1, 0, 0, 0])
    phase = [np.eye(4), z]
    probabilities = []
    for a in range(2):
        for b in range(2):
            joint = np.kron(phase[a], phase[b].conj().T)
            vector = np.kron(psi, psi)
            probabilities.append((1 + (vector @ joint @ vector).real) / 2)
    rec = {
        name: int(rng.binomial(SHOTS, np.clip(dist @ probabilities, 0, 1)))
        for name, dist in (("native", p), ("product", q))
    }
    true_q = float(np.sum((q - p) * np.log(q / p)))
    _, full, full_reject, full_total = source_counts(full=True)
    assert full_reject == reject and full_total == total and np.min(full) > 0
    labels = np.flatnonzero((stock.MAG != 1) & (np.arange(11) != 0))
    bits = stock.MAG[labels] // 2
    totals = np.outer(full.sum(1), full.sum(0))
    coarse_totals = np.outer(table.sum(1), table.sum(0))
    # P/Q is constant within each of the four magnitude blocks, entry by entry.
    np.testing.assert_array_equal(
        full * coarse_totals[np.ix_(bits, bits)], totals * table[np.ix_(bits, bits)]
    )
    complete_banks, _, _ = acquire_joint(
        np.random.default_rng(537099), full, reject, total
    )
    prior = HERE / "mediated_release.json"
    gain = timing.lab.ordinary_gain(json.loads(prior.read_bytes()))
    out = dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "PAIR_CONTACT.md",
                HERE / "FULL_PAIR_REFERENCE.md",
                stock.HERE / "results.json",
                prior,
            )
        },
        physical_seams=edges,
        class_words=total,
        accepted_weight=table.tolist(),
        rejected_weight=reject,
        covariance=str(
            Fraction(int(table[1, 1]), int(table.sum()))
            - Fraction(int(table[1].sum() * table[:, 1].sum()), int(table.sum()) ** 2)
        ),
        full_sharp_labels=labels.tolist(),
        full_sharp_weights=full.tolist(),
        complete_head_banks=complete_banks,
        ordinary_four_kappa=gain.tolist(),
        exact_connected_nats=true_q,
        source=source,
        clock=clocks,
        receiver=rec,
    )
    out["decoded"] = decode(out)
    diff = (q - p) @ probabilities
    for key, target in (
        ("connected_charge_nats", true_q),
        ("full_support_nats", true_q),
        ("from_source_fringe_shift", diff),
        ("held_fringe_shift", diff),
    ):
        a, b = out["decoded"][key]
        assert a <= target <= b
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    mode = cli.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--counts-only", action="store_true")
    args = cli.parse_args()
    path = HERE / "pair_contact.json"
    if args.counts_only:
        out = json.loads(path.read_bytes())
        assert out["decoded"] == decode(out)
    else:
        out = run()
        if args.write:
            path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
        else:
            timing.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["covariance"], out["decoded"])
