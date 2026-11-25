"""Independent four-tick replay on two held native preparation arms."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ... import station
from . import measure


def run(weights):
    survey = measure.survey
    previous = json.loads((survey.HERE / "translation.json").read_text())
    x0 = json.loads((survey.HERE.parent / "results.json").read_text())
    x = np.array(x0["x_real"]) + 1j * np.array(x0["x_imag"])
    d = np.array(previous["tangent_center"], dtype=float).mean(axis=1) / survey.SCALE
    d /= np.linalg.norm(d)
    change = np.einsum("ba,aij->bij", d.reshape(2, 9), station.iso.basis(True))
    arm_step = 5e-5
    plus = x + arm_step * change
    minus = x - arm_step * change
    # Both arms lie in the proven 1e-4 reachable neighborhood of this center.
    assert min(np.linalg.eigvalsh(s).min() for s in (minus, plus)) > 0.02
    spectral = (station.entropy(plus) - station.entropy(minus)) / (2 * arm_step)
    board = measure.native.Atlas(1)
    patch = measure.native.chosen(board)
    board = board.board
    route = station.iso.pairing.formation.region(board)
    pos = measure.regional.bridge.beam.indices(
        sum(measure.regional.bridge.beam.screens(board, route), [])
    )
    edges = measure.native.flux_on_edges(board, patch)
    rows = {}
    for charge in (-1, 1):
        for parity in (False, True):
            h, ops = measure.native.op_table(board, route, patch, charge, parity)
            generator = -2j * h
            starts = np.zeros((h.shape[0], len(pos)))
            starts[pos] = np.eye(len(pos))
            w = expm_multiply(generator, starts, traceA=generator.diagonal().sum())
            channels = {(charge, name): ops[name] for name in measure.regional.NAMES}
            channels.update(
                {
                    key: (edges[key[1:]], 4)
                    for key in measure.CHANNELS
                    if key[0] == charge
                }
            )
            for key, (op, unit) in channels.items():
                idx = np.unique(op.nonzero()[0])
                part = w[idx]
                applied = op[idx, :][:, idx] @ part
                observable = sum(
                    a.conj().T @ b
                    for a, b in zip(np.hsplit(part, 10), np.hsplit(applied, 10))
                ) / (10 * unit)
                # Taking the difference before applying the effect avoids
                # cancellation of the common (much larger) mean.
                diff = (plus - minus)[int(parity)] / (2 * arm_step)
                rows[key] = rows.get(key, 0.0) + float(np.trace(observable @ diff).real)
    keys = [(q, name) for q in (-1, 1) for name in measure.regional.NAMES] + list(
        measure.CHANNELS
    )
    slopes = [rows[k] for k in keys]
    predicted = float(np.array(weights, dtype=float) @ slopes)
    assert abs(predicted - spectral) < 1e-6
    return dict(
        arm_step=arm_step,
        material_slopes=slopes,
        spectrum_entropy_slope=spectral,
        material_entropy_slope=predicted,
    )
