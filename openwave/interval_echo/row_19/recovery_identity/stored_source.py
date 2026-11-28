"""Move a right head into an auxiliary archive before native delivery or after ranging."""

import argparse
from itertools import product
import json
from math import log

from mpmath import iv
import numpy as np

from . import (
    access,
    bridge_radar as bridge,
    measure,
    paired_relay as pair,
    paired_survey as radar,
    relay,
)

HERE = relay.HERE
RECEIPT = HERE / "stored_source.json"
PROTOCOL = HERE / "stored_source_protocol.md"
SEED = 191208
SOURCE_SHOTS = 2**54
SHOTS = bridge.SHOTS
PAULI = np.array([np.kron(a, b) for a, b in product(pair.AXES, repeat=2)])
EVENTS = 8 * 4 * (4 * 6 * 10 + 2 * 16 * 2) + 12 * 2
L = log(4 * EVENTS / 0.01)


def switch_state(state):
    """Trace an actual 3-register SWAP on each occupied parity sector."""
    source = []
    protected = []
    late = []
    for rho in state:
        w = rho.trace().real
        if w <= 0:
            raise ValueError("Conditional source needs an occupied sector")
        right = np.einsum("iris->rs", rho.reshape(2, 2, 2, 2)) / w
        before = np.kron(rho, right).reshape((2,) * 6)
        after = before.transpose(0, 2, 1, 3, 5, 4)
        source.append(np.trace(after, axis1=2, axis2=5).reshape(4, 4))
        protected.append(np.trace(after, axis1=1, axis2=4).reshape(4, 4))
        # If action is after all interactions, either discarded wire order
        # gives the same marginal on the surviving first register.
        late.append(np.trace(np.trace(after, axis1=2, axis2=5), axis1=1, axis2=3))
    np.testing.assert_allclose(protected, state, atol=1e-14)
    np.testing.assert_allclose(late, measure.head.marginals(state)[0], atol=1e-14)
    np.testing.assert_allclose(source, measure.product_reference(state), atol=1e-14)
    return np.asarray(source)


def fractions():
    rec = json.loads(access.RECEIPT.read_text())
    settings = {"all_accept": 0.0, "all_fail": 1.0}
    for p, cap in rec["panels"].items():
        for i, row in enumerate(cap["settings"], 1):
            settings[f"{p}/{i}"] = row["failed"] / row["cuts"]
    return settings


def source_bank(state, rng):
    angles = np.einsum("aij,bji->ba", PAULI, state).real
    return [
        [pair.draws([(w + m) / 2, (w - m) / 2, 1 - w], SOURCE_SHOTS, rng) for m in row]
        for w, row in zip(state.trace(axis1=1, axis2=2).real, angles, strict=True)
    ]


def acquire():
    board = relay.FORMATION.federation(1)
    route = relay.FORMATION.region(board)
    paths = pair.paths(board, route, relay.tidal.chart(board, route)[0])
    stocks, steps, _, _ = relay.sources(board, route)
    stocks = stocks[:4]
    accepted = np.array([switch_state(s) for s in stocks])
    np.testing.assert_allclose(accepted[2], accepted[3], atol=1e-14)
    rng = np.random.default_rng(SEED)
    rec = json.loads(pair.RECEIPT.read_text())
    paths_serial = [list(p) for p in paths]
    if rec["paths"] != paths_serial:
        raise ValueError("Changed material geometry")
    moved = {}
    for i, courier in enumerate((-1, 1)):
        channel = pair.channel(board, route, paths, courier)
        moved[courier] = pair.deliver(channel, accepted)
        np.testing.assert_allclose(
            pair.deliver(channel, stocks),
            pair.unpack(rec["installed"][i])[:4],
            atol=3e-12,
            rtol=0,
        )
    sites = [*paths[0], paths[1][0]]
    effects = {q: radar.effects(board, sites, paths, q) for q in (-1, 1)}
    clock, _ = relay.RANGING.clocks(0.5, 12)
    hands = [[pair.draws((x, 1 - x), SHOTS, rng) for x in row] for row in clock[1:]]
    panels = {}
    for key, failure in fractions().items():
        blend = failure * stocks + (1 - failure) * accepted
        source = [source_bank(s, rng) for s in blend]
        banks = {}
        for i, courier in enumerate((-1, 1)):
            failed = pair.unpack(rec["installed"][i])[:4]
            received = failure * failed + (1 - failure) * moved[courier]
            for q, panel in effects.items():
                laws = np.einsum("rbij,sji->srb", panel, received).real
                if failure == 0:
                    np.testing.assert_allclose(laws[2], laws[3], atol=1e-14)
                banks[f"{courier},{q}"] = [
                    [pair.draws(p, SHOTS, rng) for p in stock] for stock in laws
                ]
        panels[key] = dict(source=source, radar=banks)
    d = dict(
        protocol=bridge.digest(PROTOCOL),
        access=bridge.digest(access.RECEIPT),
        delivery=bridge.digest(pair.RECEIPT),
        steps=steps,
        seed=SEED,
        shots=SHOTS,
        source_shots=SOURCE_SHOTS,
        failure=fractions(),
        clock=hands,
        panels=panels,
    )
    d["measurements"] = decode(d)
    return d


def moments(bank, coefficients):
    raw = pair.checked(bank, (2, 16, 3), SOURCE_SHOTS)
    means = np.array((raw[:, :, 0] - raw[:, :, 1]) / SOURCE_SHOTS, dtype=float)
    # Each bounded +/-/refused reading is in [-1,1]. Simultaneous Hoeffding.
    eps = np.sqrt(2 * L / SOURCE_SHOTS)
    center = float(np.sum(coefficients * means))
    radius = float(np.sum(np.abs(coefficients)) * eps)
    return iv.mpf(center) + iv.mpf([-radius, radius])


def two_slopes(values, steps):
    return [
        (values[i + 1] - values[i]) / (2 * steps[k])
        for i, k in ((0, "ordinary"), (2, "held"))
    ]


def decode(d):
    rec = json.loads(pair.RECEIPT.read_text())
    failure = fractions()
    if (
        d["protocol"] != bridge.digest(PROTOCOL)
        or d["delivery"] != bridge.digest(pair.RECEIPT)
        or d["access"] != bridge.digest(access.RECEIPT)
        or d["steps"] != rec["half_steps"]
        or d["seed"] != SEED
        or d["source_shots"] != SOURCE_SHOTS
        or d["shots"] != SHOTS
        or d["failure"] != failure
        or set(d["panels"]) != set(failure)
    ):
        raise ValueError("Changed source or acquisition")
    radar.clock_check(d["clock"])
    ref, *_ = measure.head.jacobi.beam.aperture.source_and_changes()
    k = -measure.head.station.log_blocks(ref)
    coefficients = np.einsum("aij,bji->ba", PAULI, k).real / 4
    results = {}
    for key, r in d["panels"].items():
        if len(r["source"]) != 4 or set(r["radar"]) != {
            f"{c},{q}" for c in (-1, 1) for q in (-1, 1)
        }:
            raise ValueError("Incomplete source or radar panels")
        so, sh = two_slopes([moments(s, coefficients) for s in r["source"]], d["steps"])
        rows = {}
        for name, banks in r["radar"].items():
            banks = pair.checked(banks, (4, 6, 11), SHOTS)
            go, gh = two_slopes([bridge.area(b, L) for b in banks], d["steps"])
            det = gh * so - go * sh
            rows[name] = dict(
                ordinary=radar.limits(go),
                held=radar.limits(gh),
                ratio=bridge.ratio(gh, go),
                residual=radar.limits(det),
                compatible=(0 in det) if 0 not in go and 0 not in so else None,
            )
        results[key] = dict(
            ordinary=radar.limits(so),
            held=radar.limits(sh),
            source_ratio=bridge.ratio(sh, so),
            radar=rows,
        )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--counts-only", action="store_true")
    args = parser.parse_args()
    if args.counts_only:
        result = json.loads(RECEIPT.read_text())
        if decode(result) != result["measurements"]:
            raise ValueError("Count decoder disagrees")
    else:
        result = acquire()
        if args.write:
            with RECEIPT.open("x") as output:
                json.dump(result, output, separators=(",", ":"), allow_nan=False)
                output.write("\n")
        else:
            relay.tidal.verify(result, json.loads(RECEIPT.read_text()))
    for key, v in result["measurements"].items():
        arm = v["radar"]["1,1"]
        print(key, v["source_ratio"], arm["ratio"], arm["compatible"])
