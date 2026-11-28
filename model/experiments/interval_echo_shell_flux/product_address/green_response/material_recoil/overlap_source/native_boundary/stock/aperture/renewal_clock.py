"""An ordinary Ramsey flag clock with finite, source-dependent residence times."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import mpmath as mp
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from . import pair_work as lab

HERE = lab.HERE
N = 2**18


def read_clock(raw, neighbors, load, source=20):
    """Full native blank pulse, checks of neighboring loads and dark-port read."""
    h = raw[:9, :9] / 4
    split = expm(-0.5j * np.pi * h)
    undo = expm(-1.5j * np.pi * h)
    gate = [load, 3 * load]
    prep = split[np.ix_(gate, gate)][:, 0]
    combine = undo[np.ix_(gate, gate)]
    np.testing.assert_allclose(combine @ prep, [1, 0], atol=1e-13)
    empty_phase = len(neighbors) * h[3 * load, 3 * load]
    dwell = np.pi / (2 * empty_phase)
    contact = expm(-1j * dwell * h)
    banks = []
    for n in range(3):
        amp = prep.copy()
        losses = []
        for site in neighbors:
            x = n if site == source else 0
            codes = [x, 3 * load + x]
            stay = contact[codes, codes]
            # Any changed outside load records which arm exchanged, for either
            # unequal input. A full blank SWAP restores that sharp branch.
            for t, code in zip((0, load), codes, strict=True):
                assert set(np.flatnonzero(abs(contact[:, code]) > 1e-13)) <= {
                    code,
                    3 * x + t,
                }
                if x != t:
                    swapped = 3 * x + t
                    assert np.argmax(abs(expm(-1j * np.pi * h)[:, swapped])) == code
            losses.append(float(np.sum(abs(amp) ** 2 * (1 - abs(stay) ** 2))))
            amp *= stay
        bright, dark = abs(combine @ amp) ** 2
        probs = np.maximum([dark, bright, *losses], 0)
        np.testing.assert_allclose(sum(probs), 1, atol=2e-14)
        banks.append(probs / sum(probs))
    return dwell, np.array(banks)


def exact_bands(load, degree):
    # Interval elementary functions independently enclose the native amplitude.
    iv = mp.iv
    iv.dps = 35
    a = iv.pi / (2 * degree * (8 * load**2 + iv.mpf(1) / 2))

    def stay(t, x):
        if t == x:
            return iv.mpc(1)
        phase = (8 * (t - x) ** 2 + iv.mpf(1) / 2) * a
        return iv.cos(a / 2) * iv.mpc(iv.cos(phase), -iv.sin(phase))

    flags = []
    for n in range(3):
        amp = stay(load, 0) ** (degree - 1) * stay(load, n) - stay(0, n)
        flags.append((amp.real**2 + amp.imag**2) / 4)

    def half(x):
        return (x[0] + x[2]) / 2 - x[1]

    delay = half([1 / p for p in flags])
    ordinary = 1 / flags[1] - 1 / flags[0]

    def ends(x):
        return [
            float(np.nextafter(float(x.a), -np.inf)),
            float(np.nextafter(float(x.b), np.inf)),
        ]

    return dict(
        flags=[ends(x) for x in flags],
        curvature=ends(delay),
        flag_curvature=ends(half(flags)),
        extra_per_ordinary=ends(delay / (14 * ordinary)),
    )


def source_checks(w, codes, h):
    edges = lab.energy.pair.stock.carrier()[3]
    comparisons = 0
    for j, (a, b) in enumerate(edges):
        for i, word in enumerate(w):
            for k in np.flatnonzero(h[:, codes[j, i]]):
                loads = word[:4].copy()
                loads[a], loads[b] = divmod(k % 9, 3)
                np.testing.assert_array_equal(loads[:2], word[:2])
                comparisons += 1
    return comparisons


def renewal_bank(rng, counts, outcomes):
    # Every failed dark-port read consumes its whole slot. Failures, including
    # changed-neighbor records, are kept. No rejection or tail truncation.
    total = np.zeros((3, outcomes.shape[1]), dtype=np.int64)
    for n, p in enumerate(outcomes):
        active = int(counts[n])
        while active:
            flags = rng.multinomial(active, p)
            total[n] += flags
            active -= int(flags[0])
    assert np.array_equal(total[:, 0], counts)
    # Exact negative-binomial intervals conditional on the number of ticks
    # issued to this n. Four arms x three strata; nonunit counts are symmetric
    # in the designated ensembles. Simultaneously bracket their source rates.
    tail = 0.005 / (2 * 12)
    attempts = total.sum(1)
    p_lo = beta.ppf(tail, counts, attempts - counts + 1)
    p_hi = beta.ppf(1 - tail, counts, attempts - counts)
    return total, 1 / p_hi, 1 / p_lo


def timed_bank(rng, words, weights, banks, trials, pilot=False):
    if pilot:
        counts = np.bincount(
            words[rng.choice(len(words), trials, p=np.array(weights, float)), 1],
            minlength=3,
        )
    else:
        word_counts = rng.multinomial(trials, np.array(weights, float))
        counts = np.array([word_counts[words[:, 1] == k].sum() for k in range(3)])
    total, lo, hi = renewal_bank(rng, counts, banks)
    f_observed = (counts[0] + counts[2]) / trials
    ferr = np.sqrt(np.log(2 * 4 / 0.005) / (2 * trials))
    endpoints = [
        np.array([f / 2, 1 - f, f / 2])
        for f in [max(0, f_observed - ferr), min(1, f_observed + ferr)]
    ]
    return dict(
        completed_by_load=counts.tolist(),
        flags=total.tolist(),
        observed_slots_per_tick=float(total.sum() / trials),
        slots_interval_99=[
            float(min(d @ lo for d in endpoints)),
            float(max(d @ hi for d in endpoints)),
        ],
    )


def run():
    _, w, cells, _, codes, h, *_ = lab.energy.pair.setup()
    board, homes, marks, *_ = lab.energy.returned(w)
    near, split = 13, 12
    neighbors = sorted(map(int, board.neighbours[near]))
    assert set(neighbors) & set(homes) == {20} and split in neighbors
    forbidden = set(neighbors) | set(homes) | set(np.ravel(marks)) | {near}
    far = min(set(board.neighbours[split]) - forbidden)
    assert all(near not in m and split not in m and far not in m for m in marks)
    column_checks = source_checks(w, codes, h)
    sources = {
        "original": np.full(len(w), F(1, len(w)), object),
        "repaired": lab.energy.least_change(cells) / len(w),
    }
    means = {}
    bounds = {}
    rows = {}
    rng = np.random.default_rng(522132)
    validation = np.random.default_rng(44828931)
    refresh_rng = np.random.default_rng(249101)
    pilots = {}
    fresh_size = 2**26
    for t in (1, 2):
        dwell, banks = read_clock(h, neighbors, t)
        p = banks[:, 0]
        cert = exact_bands(t, len(neighbors))
        np.testing.assert_allclose(
            p, np.mean(cert["flags"], axis=1), rtol=2e-13, atol=2e-14
        )
        assert cert["curvature"][0] > 0 and cert["flag_curvature"][0] > 0
        row = {}
        for k, weights in sources.items():
            f = weights @ ((w[:, 1] - 1) ** 2)
            expected = [f / 2, 1 - f, f / 2]
            assert [sum(weights[w[:, 1] == n]) for n in range(3)] == expected
            dist = np.array(expected, float)
            pilots[f"{t}/{k}"] = timed_bank(rng, w, weights, banks, N, pilot=True)
            records = timed_bank(validation, w, weights, banks, 2**22)
            band = records["slots_interval_99"]
            word_counts = refresh_rng.multinomial(fresh_size, np.array(weights, float))
            counts = [int(word_counts[w[:, 1] == n].sum()) for n in range(3)]
            fresh = [
                refresh_rng.multinomial(num, probs).tolist()
                for num, probs in zip(counts, banks, strict=True)
            ]
            fraction = sum(v[0] for v in fresh) / fresh_size
            err = np.sqrt(np.log(8 / 0.01) / (2 * fresh_size))
            assert abs(fraction - dist @ p) <= err
            tick = float(dist @ (1 / p))
            r = 1 / tick
            time_weights = dist / p * r
            cv2 = float(time_weights @ ((p / r - 1) ** 2))
            np.testing.assert_allclose(tick * (dist @ p), 1 + cv2)
            assert band[0] <= tick <= band[1]
            v = weights @ ((w[:, 0] - 1) ** 2)
            assert v < F(1, 2)
            row[k] = dict(
                nonunit=str(f),
                **records,
                predicted_slots_per_tick=tick,
                renew_each_attempt_slots=float(1 / (dist @ p)),
                fresh_attempt_flags=fresh,
                fresh_attempt_slots_interval_99=[
                    1 / (fraction + err),
                    1 / (fraction - err),
                ],
                time_sampled_fractional_variance=cv2,
                signed_phase_visibility_lower=str(1 - 2 * v),
            )
        interval = [
            row["repaired"]["slots_interval_99"][0]
            - row["original"]["slots_interval_99"][1],
            row["repaired"]["slots_interval_99"][1]
            - row["original"]["slots_interval_99"][0],
        ]
        assert interval[0] > 0
        assert (
            row["repaired"]["fresh_attempt_slots_interval_99"][1]
            < row["original"]["fresh_attempt_slots_interval_99"][0]
        )
        means[t] = (
            row["repaired"]["predicted_slots_per_tick"]
            - row["original"]["predicted_slots_per_tick"]
        )
        bounds[t] = cert["extra_per_ordinary"]
        rows[str(t)] = dict(
            dwell=dwell,
            exact=cert,
            outcome_order=["tick", "bright", *[f"changed_{v}" for v in neighbors]],
            arms=row,
            delay_increment_99=interval,
            required_delay=means[t],
            ordinary_donor_delay=1 / p[1] - 1 / p[0],
        )
    assert bounds[1][0] > bounds[2][1]
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "RENEWAL_CLOCK.md", HERE / "pair_field.json")
        },
        near=near,
        split=split,
        reference=int(far),
        neighbors=neighbors,
        compact_column_entries=column_checks,
        ticks_per_arm=2**22,
        pilot_ticks_per_arm=N,
        pilot=pilots,
        fresh_attempts_per_arm=fresh_size,
        readers=rows,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "renewal_clock.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.check(json.loads(dst.read_bytes()), result)
    for t, row in result["readers"].items():
        print(
            t,
            row["required_delay"],
            row["delay_increment_99"],
            row["exact"]["extra_per_ordinary"],
        )
