"""Read the actual heavy core and external receiver on a continuously connected source."""

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from . import cold_ports, open_ports, outlet, wing
from .measure import HERE
from .outlet_clock import LABELS
from ..joined_decay.measure import core_pulse
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import Z, circuit

SETTINGS = ("native", "source_cut", "exit_cut", "no_reversal")
TICKS = 2 * np.pi * np.array([162, 163, 164]) / outlet.PERIOD
WAITS = np.r_[0, TICKS, 32 + TICKS]


def generator(prefix, tail, words, dimers, name):
    n = len(words)
    if name == "no_reversal":
        tail = tail + np.diag(
            [int(d == 1 and dimers[-1][0] in word) for d in (0, 1) for word in words]
        )
    # The C,D triangle contributes -(1-c)d, relative to the orientation
    # independent source offset. The core and C occupy the eight leading states.
    abc = prefix[::4, ::4]
    ck = np.arange(8) % 2
    d = np.repeat([0.0, 1.0], n)
    bond = np.kron(-(1 - ck), d) if name in ("native", "exit_cut") else np.zeros(16 * n)
    return np.kron(np.eye(8), tail) + np.kron(abc, np.eye(2 * n)) + np.diag(bond)


def prepare(h, initial, k):
    values, v = eigh(h, driver="evd")
    shape = len(k)
    rows = np.flatnonzero((abs(k) ** 2).sum(1) > 0)
    installed = v.reshape(8, shape, -1)
    local = np.vstack([k[rows] @ subsystem for subsystem in installed])
    supported = np.concatenate([i * shape + rows for i in range(8)])
    amp = v.T @ initial
    phase = (values[:, None] - values[None, :]) / outlet.PERIOD
    phase = (phase + 0.5) % 1 - 0.5
    rho = np.exp(-1j * np.pi * (LABELS - 1) * phase)
    rho *= np.sinc(LABELS * phase) / np.sinc(phase)
    rho *= amp[:, None] * amp[None, :]
    post = local @ rho @ local.conj().T
    del rho, local, phase
    weights, states = eigh(post, driver="evd")
    assert weights.min() > -1e-12
    admission = float(np.trace(post).real)
    # Floating-point compression only. The omitted positive mass provides a
    # uniform absolute error bound on any original-launch event probability.
    nonnegative = np.maximum(0, weights)
    start = int(np.searchsorted(np.cumsum(nonnegative), 1e-10))
    discarded = nonnegative[:start].sum()
    kept = np.zeros((len(values), len(weights) - start), complex)
    kept[supported] = states[:, start:] * np.sqrt(nonnegative[start:])
    return (
        values,
        v,
        kept,
        dict(
            admission=admission,
            rank=len(weights) - start,
            discarded_weight_upper=float(discarded + 1e-11),
        ),
    )


@dataclass
class Apparatus:
    dimers: tuple
    words: list
    groups: list
    sectors: dict
    response: np.ndarray
    energy_read: np.ndarray
    collapse: np.ndarray
    h: np.ndarray
    initial: np.ndarray
    cold_filter: np.ndarray


def apparatus(name):
    board, dimers, cell, _, sites, contacts, cut = open_ports.aperture()
    internal = [e for e in contacts if set(e) <= set(cell)]
    sectors = outlet.receiving_cell(board, cell, internal)
    _, one, cv = sectors[1]
    core = hamiltonians(4)[0]
    cev, coreframe = eigh(core)
    energy_read = circuit(coreframe, 32)
    collapse = circuit(core_pulse(coreframe), 32)
    threshold = np.ptp(cev) / 10
    response = outlet.clock(
        one, one[0], threshold, labels=LABELS, low_edge=outlet.PERIOD / LABELS
    )
    heads = wing.installation()[2]
    _, five = eigh(heads["native"])
    original = five[:, 0].reshape(4, 2, 2, 2)
    seams = [e for e in contacts if name != "exit_cut" or e not in cut]
    words, tail = outlet.connected_wing(board, dimers, sites, seams)
    prefix = heads[
        "no_reversal"
        if name == "no_reversal"
        else "boundary_cut"
        if name == "source_cut"
        else "native"
    ]
    h = generator(prefix, tail, words, dimers, name)
    initial = outlet.inject(original, words, dimers, cell, cv[:, :1]).T.reshape(-1)
    groups = outlet.read_groups(words, sectors, cell)
    k = cold_ports.zero_bin(groups, sectors)
    return Apparatus(
        dimers, words, groups, sectors, response, energy_read, collapse, h, initial, k
    )


def acquire_setting(name):
    apparatus_state = apparatus(name)
    h = apparatus_state.h
    n = len(apparatus_state.words)
    energies, v, stock, preparation = prepare(
        h, apparatus_state.initial, apparatus_state.cold_filter
    )
    groups, sectors, response = (
        apparatus_state.groups,
        apparatus_state.sectors,
        apparatus_state.response,
    )
    energy_read, collapse = apparatus_state.energy_read, apparatus_state.collapse
    hamiltonian = sparse.csr_matrix(h)
    del h
    rank = stock.shape[-1]
    base_work = float(np.vdot(stock, hamiltonian @ stock).real)
    rows = {}
    for arm in ("idle", "core_pulse", "D_phase"):
        launched = stock.reshape(4, 2, 2, n, rank)
        if arm == "core_pulse":
            launched = np.einsum("ba,acdmr->bcdmr", collapse, launched)
        elif arm == "D_phase":
            launched = launched * Z.diagonal()[None, None, :, None, None]
        launched = launched.reshape(-1, rank)
        work = float(np.vdot(launched, hamiltonian @ launched).real) - base_work
        spectral = v.T @ launched
        panels = []
        for t in WAITS:
            psi = v @ (np.exp(-1j * t * energies)[:, None] * spectral)
            core_amplitudes = np.einsum(
                "ab,acdr->dbcr", energy_read, psi.reshape(4, 2, 2 * n, rank)
            )
            # flight_read checks unit mass. Keep the intervention's original
            # admission probability outside, and restore it on every output.
            norm = np.linalg.norm(core_amplitudes) ** 2
            rec = outlet.flight_read(
                core_amplitudes.reshape(2 * n, -1) / np.sqrt(norm),
                groups,
                sectors,
                response,
            )
            law = np.asarray(rec["core_by_outlet"]) * norm
            panels.append(law.round(12).tolist())
        rows[arm] = dict(work_per_original_launch=work, core_by_outlet=panels)
    return dict(preparation=preparation, arms=rows)


def collect():
    settings = {}
    for name in SETTINGS:
        settings[name] = acquire_setting(name)
    return dict(
        protocol_sha256=sha256((HERE / "JOINED_OUTLET.md").read_bytes()).hexdigest(),
        waits=WAITS.tolist(),
        clock_false_light=json.loads((HERE / "outlet_clock.json").read_text())[
            "false_light_ceiling"
        ],
        settings=settings,
    )


def decode(data):
    for key, setting in data["settings"].items():
        admission = setting["preparation"]["admission"]
        print(key, "cold admission", admission, "rank", setting["preparation"]["rank"])
        rows = setting["arms"]
        idle = np.array(rows["idle"]["core_by_outlet"])
        for arm, results in rows.items():
            p = np.asarray(results["core_by_outlet"])
            np.testing.assert_allclose(p.sum(axis=(1, 2)), admission, atol=2e-10)
            np.testing.assert_allclose(
                p[0].sum(axis=0), idle[0].sum(axis=0), atol=3e-10
            )
            print(
                arm,
                "work",
                results["work_per_original_launch"],
                "light excess",
                np.round((p - idle)[:, :, 2].sum(1), 8),
                "lower plus light",
                np.round(p[:, 1:3, 2].sum(1), 7),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--setting", choices=SETTINGS)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.setting:
        if not args.acquire:
            parser.error("a single-setting run needs --acquire")
        result = acquire_setting(args.setting)
    else:
        result = collect()
        decode(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "joined_outlet.json").read_text()), result)
