"""Separate finite-shot readouts of two-port entropy and bridge work on the controlled source pair."""

import argparse
import json

import numpy as np
from mpmath import iv
from scipy.sparse.linalg import expm_multiply

from . import measure


def untagged_read(physical_blocks, trials, rng):
    # An analyzer distinguishes the actual total local load 1, 2, or 3.
    # Preparations have no retained input-charge label. The unselected block,
    # outside mode and a missed local interference are all unconditional events.
    unit = np.sqrt(np.log(2 * 18 / 0.01) / (2 * trials))
    endpoints = []
    for t in (-0.01, 0.01):
        eigen = []
        weights = []
        data = {}
        sym = anti = 0.0
        for q in (-1, 1):
            d, x, y, z = np.array(physical_blocks[q, t]) / 2
            sym += (d + x) / 2
            anti += (d - x) / 2
            parsed = []
            for name, v in zip(measure.NAMES, (d, x, y, z)):
                if name == "occupation":
                    val = int(rng.binomial(trials, d))
                    hist = [val, trials - val]
                    rate = val / trials
                    eps = unit
                else:
                    hist = rng.multinomial(
                        trials, [(d + v) / 2, (d - v) / 2, 1 - d]
                    ).tolist()
                    rate = (hist[0] - hist[1]) / trials
                    eps = 2 * unit
                data[str(q) + "_" + name] = hist
                parsed.append(iv.mpf([rate - eps - 2e-14, rate + eps + 2e-14]))
            weights.append(parsed[0])
            eigen.extend(measure.modular(parsed)[1][1:])
        eigen.append(1 - sum(weights))
        entropy = -sum(v * iv.ln(v) for v in eigen)
        w = rng.multinomial(trials, [sym, anti, 1 - sym - anti]).tolist()
        data["work_8_9_0"] = w
        work = (8 * w[0] + 9 * w[1]) / trials
        err = 9 * unit + 2e-14
        endpoints.append((entropy, iv.mpf([work - err, work + err]), data))
    ds = endpoints[1][0] - endpoints[0][0]
    dw = endpoints[1][1] - endpoints[0][1]
    assert ds.a > 0 and dw.a > 0
    return dict(
        signed_entropy_difference=str(ds),
        signed_work_difference=str(dw),
        entropy_per_work=str(ds / dw),
        source_charge_tag_retained=False,
        counts=[v[-1] for v in endpoints],
    )


def run():
    iv.dps = 30
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    launches = measure.bridge.beam.indices(sum(rings, []))
    near = measure.bridge.beam.indices(route[1:3])
    ref, _, _ = measure.source()
    origin = np.array(
        [np.array(m, dtype=complex) * p for m, p in zip(ref, (1 - 0.1343, 0.1343))]
    )
    arm = json.loads((measure.HERE / "control.json").read_text())
    d = np.einsum(
        "bm,mij->bij",
        np.array(arm["center_coordinates"]).reshape(2, 15),
        measure.c.iso.basis(False),
    )
    trials = 2**59
    rng = np.random.default_rng(512626)
    # 20 settings: two charges, two endpoints, four entropy analyzers and one work read.
    unit = np.sqrt(np.log(40 / 0.01) / (2 * trials))
    rows = {}
    physical_blocks = {}
    for q in (-1, 1):
        orbit = []
        for b in (False, True):
            h = (-2j) * measure.bridge.native.generator(board, route, q, b)
            start = np.zeros((h.shape[0], 40))
            start[launches] = np.eye(40)
            wave = expm_multiply(h, start, traceA=h.diagonal().sum())[near]
            orbit.append(wave.reshape(2, 4, 10, 4))
        endpoints = []
        for t in (-0.01, 0.01):
            state = origin + t * d
            assert min(np.linalg.eigvalsh(v).min() for v in state) > 0
            occupied = (
                np.einsum("brhxi,bij,bshxj->rs", orbit, state, np.conj(orbit)) / 10
            )
            occ = occupied.trace().real
            x = 2 * occupied[0, 1].real
            y = -2 * occupied[0, 1].imag
            z = (occupied[0, 0] - occupied[1, 1]).real
            physical_blocks[q, t] = (occ, x, y, z)
            data = {}
            measured = []
            for key, v in zip(measure.NAMES, (occ, x, y, z)):
                if key == "occupation":
                    n = int(rng.binomial(trials, occ))
                    estimate = n / trials
                    data[key] = [n, trials - n]
                    scale = 1
                else:
                    n = rng.multinomial(trials, [(occ + v) / 2, (occ - v) / 2, 1 - occ])
                    estimate = (int(n[0]) - int(n[1])) / trials
                    data[key] = n.tolist()
                    scale = 2
                err = scale * unit + 2e-14
                measured.append(iv.mpf([estimate - err, estimate + err]))
            _, spectrum = measure.modular(measured)
            entropy = -sum(v * iv.ln(v) for v in spectrum)
            work = rng.multinomial(trials, [(occ + x) / 2, (occ - x) / 2, 1 - occ])
            work_mean = (8 * int(work[0]) + 9 * int(work[1])) / trials
            work_err = 9 * unit + 2e-14
            work_bound = iv.mpf([work_mean - work_err, work_mean + work_err])
            data["work_8_9_0"] = work.tolist()
            endpoints.append((entropy, work_bound, data))
        ds = endpoints[1][0] - endpoints[0][0]
        dw = endpoints[1][1] - endpoints[0][1]
        assert ds.a > 0 and dw.a > 0
        rows[str(q)] = dict(
            signed_entropy_difference=str(ds),
            signed_work_difference=str(dw),
            entropy_per_work=str(ds / dw),
            counts=[v[-1] for v in endpoints],
        )
    return dict(
        trials_per_setting=trials,
        simultaneous_coverage=0.99,
        source_coefficients=[-0.01, 0.01],
        charges=rows,
        fair_untagged=untagged_read(physical_blocks, trials, rng),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "counts.json").write_text(data)
    else:
        print(data, end="")
