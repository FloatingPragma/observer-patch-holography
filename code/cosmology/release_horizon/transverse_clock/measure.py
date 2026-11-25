"""Replay the clock records and a two-layer material visibility measurement."""

import argparse
import json
from math import log, sqrt
import numpy as np

if __package__:
    from .oscillator import HERE, Instrument, oscillator, rounded
else:
    from oscillator import HERE, Instrument, oscillator, rounded
from angular_transfer.measure import run as angular_run
from relay import apply, transpose
from scipy.special import eval_legendre

SEED = 210022
STARTS = 2**16
CAL = 2**20


def gate_probability(inst, owner, exposure):
    # Borrow the same nominee port to query the pulse on a physically
    # disjoint owner. Parking the other local rails leaves a clean reference.
    pol, net = inst.pol, inst.net
    pair = pol.givens(0, 1, exposure / 2)[:, 0]
    for selected in (2, 3):
        net.word[12 * owner + selected] = 2
        assert net.read(owner, 3) == int(selected == 3)
        net.word[12 * owner + selected] = 0
    return abs(pair[1]) ** 2


def shutter_messages(net, arrivals):
    # Borrow the actual sharp flag and return it, so a second read of the
    # same dimer can carry its first decision to a different detector.
    targets = sorted({0, *(end for _, end in arrivals)})
    receipts = []
    for source in (260, 261):
        net.compile([(source, t) for t in targets], 8)
        base = 12 * source
        lift = [(base + a, base + b) for a, b in transpose(3, 8)]
        for t in targets:
            send = lift + net.ops[base + 8, 12 * t + 8]
            outcomes = []
            for bit in (0, 1):
                load = base + 2 + bit
                net.word[load] = 2
                before = net.word.copy()
                cut = net.read(t, 8)
                apply(net.word, send)
                arrived = net.read(t, 8)
                apply(net.word, send[::-1])
                assert net.word == before and cut == 0 and arrived == bit
                net.word[load] = 0
                outcomes.append([cut, arrived])
            receipts.append([source, t, len(send), outcomes])
    return receipts


def visibility(rng):
    inst = Instrument()
    # Independent material gas shutters on two owners, separate from the
    # angular payload (which occupies only owners 0 through 124 at readout).
    rates = [
        gate_probability(inst, owner, t)
        for owner, t in ((260, np.pi / 3), (261, 2 * np.pi / 3), (260, np.pi / 2))
    ]
    calibration = [rng.binomial(CAL, p) for p in rates]
    parent, data = angular_run()
    messages = shutter_messages(inst.net, parent["endpoints"]["arrivals"])
    directions = np.asarray(parent["directions"])
    h6 = eval_legendre(6, directions[:12] @ directions.T).mean(axis=0)
    contrast = np.r_[np.full(12, 1 / 12), np.full(20, -1 / 20)]
    q = eval_legendre(2, directions @ directions[0])
    single = data["tomography"][:, :, :2, 0, :2].sum(axis=-1).mean(axis=2) / 32e9
    double = data["second_kernels"][0] / 1e9
    identity = np.eye(32)
    kernels = (identity, single, single, double)
    # Flag labels 00, 01, 10, 11. The last two bins are the first and
    # second projection failures. Hit directions retain the original denominator.
    laws = np.zeros((32, 4, 34))
    for i, kernel in enumerate(kernels):
        laws[:, i, :32] = kernel
        if i >= 2:
            laws[:, i, 32] = 1 - single.sum(axis=1)
        laws[:, i, 33] = 1 - laws[:, i, :33].sum(axis=1)
    assert laws.min() >= -2e-9
    laws = np.maximum(laws, 0.0)
    stats, arrays = {}, {}
    for name, a, b, linked in [
        ("layered", rates[0], rates[1], False),
        ("fresh_fair", rates[2], rates[2], False),
        ("reused_fair", rates[2], rates[2], True),
    ]:
        # A fresh blank dimer supplies the second read on the independent arm.
        # A parked sharp pointer gives the identical bit on its repeat read.
        history = np.array(
            (
                [1 - a, 0, 0, a]
                if linked
                else [(1 - a) * (1 - b), (1 - a) * b, a * (1 - b), a * b]
            )
        )
        joint = history[None, :, None] * laws
        n = np.array(
            [
                rng.multinomial(STARTS, p.ravel() / p.sum()).reshape(4, 34)
                for p in np.round(joint, 12)
            ],
            dtype=np.uint32,
        )
        arrays[name] = n
        k = joint[:, :, :32].sum(axis=1)
        alive = k.sum(axis=1).mean()
        by_history = joint[:, :, :32].sum(axis=2).mean(axis=0)
        last = [by_history[0], by_history[2], by_history[1] + by_history[3]]
        raw = contrast @ k @ h6 / (contrast @ h6)
        counted = n[:, :, :32].sum(axis=1) / STARTS
        observed = contrast @ counted @ h6 / (contrast @ h6)
        two = q @ k @ q / (q @ q)
        stats[name] = dict(
            flag_weights=rounded(history),
            surviving_last=rounded(last),
            survival=round(alive, 10),
            raw_six=round(raw, 10),
            counted_six=round(observed, 10),
            accepted_six=round(raw / alive, 10),
            raw_two=round(two, 10),
            accepted_two=round(two / alive, 10),
        )
    # Joint 99% coverage across calibration, signed clock quadratures and
    # the three held raw contrasts. The individual total allocations sum to .01.
    cal_eps = sqrt(log(2 * len(rates) / 0.003) / (2 * CAL))
    gain_eps = sqrt((contrast @ contrast) * log(2 * len(stats) / 0.003) / (2 * STARTS))
    c = np.array(calibration) / CAL
    pred = (1 - c[0]) * (1 - c[1])
    # Product of two probabilities is 1-Lipschitz in each argument.
    assert abs(stats["layered"]["counted_six"] - pred) <= gain_eps + 2 * cal_eps + 1e-8
    assert (
        abs(stats["reused_fair"]["counted_six"] - (1 - c[2]) ** 2)
        > gain_eps + 2 * cal_eps
    )
    return dict(
        owners=[260, 261],
        messages=messages,
        calibration=calibration,
        calibration_starts=CAL,
        angular_starts=STARTS,
        calibration_halfwidth=cal_eps,
        gain_halfwidth=gain_eps,
        predicted_clear_from_calibration=pred,
        banks=stats,
    ), arrays


def run():
    clock_seed, visibility_seed = np.random.SeedSequence(SEED).spawn(2)
    clock = oscillator(np.random.default_rng(clock_seed))
    m = sum(len(v["reads"]) * 4 for v in clock["modes"])
    eps = sqrt(2 * log(2 * m / 0.004) / clock["shots"])
    for mode in clock["modes"]:
        for row in mode["reads"]:
            outcomes = np.asarray(row["counts"])
            signed = (outcomes[:, :, 0] - outcomes[:, :, 1]) / clock["shots"]
            assert np.max(abs(signed - row["response"])) <= eps
    clock["signed_halfwidth"] = eps
    vis, banks = visibility(np.random.default_rng(visibility_seed))
    return dict(oscillator=clock, visibility=vis), banks


if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--write", action="store_true")
    save = args.parse_args().write
    data, banks = run()
    body = json.dumps(data, separators=(",", ":"), default=lambda x: x.item()) + "\n"
    record = HERE / "records.json"
    archive = HERE / "directions.npz"
    if save:
        record.write_text(body)
        np.savez_compressed(archive, **banks)
    else:
        assert record.read_text() == body
        with np.load(archive, allow_pickle=False) as prior:
            assert set(prior.files) == set(banks)
            for name, array in banks.items():
                np.testing.assert_array_equal(prior[name], array)
    print(data["visibility"]["banks"])
