"""Test an isolated neutral source calibration through the already open outlet."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from . import outlet
from .outlet_clock import LABELS
from .measure import HERE
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import Z, circuit


def outlet_effects(groups, sectors, clock):
    dim = sum(len(rows) for _, rows, _ in groups)
    effects = np.zeros((5, dim, dim))
    for count, rows, v in groups:
        ids = np.ix_(rows, rows)
        if count == 1:
            for i in range(3):
                effects[i + 1][ids] = (v * clock[:, i]) @ v.T
        else:
            effects[0 if count == 0 else 4][ids] = np.eye(len(rows))
    np.testing.assert_allclose(effects.sum(axis=0), np.eye(dim), atol=2e-13)
    return effects


def aged(h, initial):
    values, v = eigh(h)
    phase = (values[:, None] - values[None, :]) / outlet.PERIOD
    phase = (phase + 0.5) % 1 - 0.5
    coherence = (
        np.exp(-1j * np.pi * (LABELS - 1) * phase)
        * np.sinc(LABELS * phase)
        / np.sinc(phase)
    )
    stock = v.T @ initial
    rho = np.outer(stock, stock.conj()) * coherence
    return values, v, rho


def run_setting(h, stock, phase_d, exchange_read, effects):
    values, v, idle = aged(h, stock)
    phase = v.T @ phase_d @ v
    pulse = phase @ idle @ phase
    delta = pulse - idle
    work = float(values @ delta.diagonal().real)
    local_work = float(np.trace((v.T @ exchange_read @ v) @ idle).real)
    np.testing.assert_allclose(work, local_work, atol=4e-12)
    readers = np.array([v.T @ effect @ v for effect in effects])
    rows = {}
    for name, state in (("idle", idle), ("pulse", pulse)):
        weighted = readers * state.T
        panels = []
        for t in outlet.WAITS:
            spin = np.exp(1j * t * values)
            p = np.einsum("i,kij,j->k", spin, weighted, spin.conj()).real
            np.testing.assert_allclose(p.sum(), 1, atol=4e-12)
            assert min(p) > -3e-14
            panels.append(np.maximum(0, p).round(13).tolist())
        rows[name] = panels
    np.testing.assert_allclose(rows["pulse"][0], rows["idle"][0], atol=3e-12)
    coefficients = readers[2] * delta.T
    constant = np.trace(coefficients).real
    offdiag = np.abs(coefficients - np.diag(coefficients.diagonal())).sum()
    return dict(
        source_energy=work,
        local_exchange_work=local_work,
        panels=rows,
        all_time_constant=float(constant),
        all_time_oscillating_bound=float(offdiag),
        all_time_upper=float(constant + offdiag),
    )


def acquire():
    board, dimers, path, cell, internal = outlet.outlet_geometry()
    sectors = outlet.receiving_cell(board, cell, internal)
    _, one, vcell = sectors[1]
    response = outlet.clock(
        one,
        one[0],
        np.ptp(eigh(hamiltonians(4)[0], eigvals_only=True)) / 10,
        labels=LABELS,
        low_edge=outlet.PERIOD / LABELS,
    )
    _, wing, _, _ = hamiltonians(4)
    values, v = eigh(wing)
    calibrated = circuit(v, 33)
    source = calibrated[:, 0]
    local_phase = np.kron(Z, np.eye(2))
    lines = abs(calibrated.T @ (local_phase @ source)) ** 2
    out = {}
    for name in ("native", "entrance_cut", "no_reversal"):
        words, h = outlet.material(
            board, dimers, path, cell, internal, name == "entrance_cut"
        )
        n = len(words)
        if name == "no_reversal":
            # Replacing cost 33 by 32 in 33(1-d e) adds d e after removing the constant.
            h += np.diag(
                [int(d == 1 and dimers[-1][0] in w) for d in (0, 1) for w in words]
            )
        groups = outlet.read_groups(words, sectors, cell)
        effects = outlet_effects(groups, sectors, response)
        guide = np.zeros((4, 2, 2, 2))
        guide[0, 0] = source.reshape(2, 2)
        stock = outlet.inject(guide, words, dimers, cell, vcell[:, :1])[:, 0]
        phase = np.kron(Z, np.eye(n))
        exchange_read = np.kron(np.array([[0.0, 1.0], [1.0, 0.0]]), np.eye(n))
        np.testing.assert_allclose(phase @ h @ phase - h, exchange_read, atol=1e-13)
        out[name] = run_setting(h, stock, phase, exchange_read, effects)
    return dict(
        protocol_sha256=sha256((HERE / "TRANSMISSION.md").read_bytes()).hexdigest(),
        isolated_excitation_energies=(values - values[0]).tolist(),
        isolated_phase_probabilities=lines.tolist(),
        isolated_mean_gain=float(lines @ (values - values[0])),
        waits=list(outlet.WAITS),
        settings=out,
    )


def decode(record):
    print(
        "isolated stimulus",
        record["isolated_mean_gain"],
        record["isolated_phase_probabilities"],
    )
    for name, setting in record["settings"].items():
        p = np.array(setting["panels"]["pulse"])
        i = np.array(setting["panels"]["idle"])
        print(
            name,
            "work",
            setting["source_energy"],
            "light excesses",
            (p[:, 2] - i[:, 2]).round(9),
            "all-time excess upper bound",
            setting["all_time_upper"],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    decode(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "transmission.json").read_text()), result)
