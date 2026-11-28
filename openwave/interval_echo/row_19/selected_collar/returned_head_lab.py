"""Read the native heads returned by each two-leg pulse and record the time of arrival."""

import argparse
import json
from hashlib import sha256

import numpy as np
from matched_einstein_limit.source_implication.shared_radar.measure import preparations

from . import reference_read as heads
from . import six_station_read as six

HERE = six.HERE
SEED = 28110673
SHOTS = 2**60
BATCH = 2**50


def sources():
    prep = preparations()
    ref = prep["equilibrium"]
    _, ends, *_ = heads.measure.select_source(ref)
    return {
        "reference": ref,
        "ordinary": prep["ordinary"],
        "turn_plus": prep["plus"],
        "turn_minus": prep["minus"],
        "source_minus": ends[0],
        "source_plus": ends[1],
    }


def pulses(states):
    """Native first arrivals on both legs, retaining the output parity and heads."""
    survey = six.survey
    native = survey.native
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    ports = (*region, 5, 6)
    columns = np.array([4 * p + i for p in ports for i in range(4)])
    coeff = native.prior.weights()
    by_time = np.zeros((len(states), len(survey.ROUTES), 10, 2, 4, 4), complex)
    trace_check = 0.0
    for b in (0, 1):
        for charge in (-1, 1):
            samples = native.kernel_panels(
                native.generator(board, region, charge, bool(b)), columns, coeff
            )
            hits = {
                route: survey.history.first_hits(samples, *route)
                for route in survey.ROUTES
            }
            for k, route in enumerate(survey.ROUTES):
                for t in range(2, 12):
                    raw = [
                        native.prior.multiply(hits[route[::-1]][t - a], hits[route][a])
                        for a in range(1, t)
                    ]
                    amp = np.array(
                        [
                            (x[0].astype(float) + 1j * x[1].astype(float))
                            / float(survey.history.F)
                            for x in raw
                        ]
                    )
                    out = (
                        np.einsum(
                            "aij,pjk,alk->pil",
                            amp,
                            states[:, b],
                            amp.conj(),
                            optimize=True,
                        )
                        / 2
                    )
                    by_time[:, k, t - 2, b] += out
                    effects = sum(
                        native.prior.multiply(native.prior.adjoint(x), x) for x in raw
                    )
                    E = (effects[0].astype(float) + 1j * effects[1].astype(float)) / (
                        2 * float(survey.history.F)
                    )
                    prob = np.einsum("ij,pji->p", E, states[:, b]).real
                    trace_check = max(
                        trace_check,
                        float(abs(np.trace(out, axis1=-2, axis2=-1) - prob).max()),
                    )
    assert trace_check < 2e-15
    return list(ports), by_time, trace_check


def counts(p, rng):
    # Keep all draws below exact binary64 integer precision. Order the rare
    # returns first and put the unique largest category last for stability.
    p = np.asarray(p, dtype=float)
    assert p.min() >= -1e-15 and abs(p.sum() - 1) < 1e-12
    order = np.argsort(p, kind="stable")
    prob = np.maximum(p[order], 0)
    prob[-1] = 1 - prob[:-1].sum()
    samples = rng.multinomial(BATCH, prob, size=SHOTS // BATCH).sum(axis=0)
    assert samples.sum() == SHOTS
    return samples[np.argsort(order)].tolist()


def pairs(state, rng, timeout):
    tape = []
    for a in heads.AXES[1:]:
        row = []
        for b in heads.AXES[1:]:
            E = [
                np.kron((heads.AXES[0] + s * a) / 2, (heads.AXES[0] + t * b) / 2)
                for s in heads.SIGNS
                for t in heads.SIGNS
            ]
            prob = np.einsum("kij,bji->bk", E, state).real.ravel()
            if timeout:
                prob = np.r_[prob, 1 - prob.sum()]
            row.append(counts(prob, rng))
        tape.append(row)
    return tape


def run():
    prep = sources()
    ports, history, err = pulses(np.array(list(prep.values())))
    rng = np.random.default_rng(SEED)
    tapes = {}
    for (label, state), Y in zip(prep.items(), history, strict=True):
        times = np.trace(Y, axis1=-2, axis2=-1).real.sum(axis=2)
        tapes[label] = {
            "time": [counts(np.r_[p, 1 - p.sum()], rng) for p in times],
            "incoming": pairs(state, rng, False),
            "returned": [pairs(s, rng, True) for s in Y.sum(axis=1)],
        }
    return {
        "protocol_sha256": sha256(
            (HERE / "RETURNED_HEADS.md").read_bytes()
        ).hexdigest(),
        "seed": SEED,
        "shots_per_setting": SHOTS,
        "ports": ports,
        "routes": list(map(list, six.survey.ROUTES)),
        "trace_effect_max_difference": err,
        "tapes": tapes,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    data = run()
    path = HERE / "returned_head_lab.json"
    if cli.parse_args().write:
        path.write_text(json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        stored = json.loads(path.read_bytes())
        assert (
            data["tapes"] == stored["tapes"]
            and data["protocol_sha256"] == stored["protocol_sha256"]
        )
        assert data["seed"] == stored["seed"] and data["routes"] == stored["routes"]
        assert (
            data["ports"] == stored["ports"]
            and data["shots_per_setting"] == stored["shots_per_setting"]
        )
    print(len(data["tapes"]), data["ports"], data["trace_effect_max_difference"])
