"""Qualify the complete available native approach to the external receiving cell."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from . import outlet, transmission
from .material import route, edge
from .measure import HERE
from .outlet_clock import LABELS
from ..passive_decay.measure import compare_records, hamiltonians
from ..passive_decay.readers import X, Z, circuit

NAMES = ("idle", "X", "Y", "Z")
BASIS = (np.eye(2), -1j * X, -1j * np.array([[0, -1j], [1j, 0]]), -1j * Z)


def aperture():
    board, dimers, _, cell, _ = outlet.outlet_geometry()
    heads = route((3,), crossing=True)[2]
    owner = dimers[-1][0] // 12
    reserved = {p for pair in (*dimers[:-1], *heads.keys()) for p in pair}
    available = (set(range(12 * owner, 12 * (owner + 1))) - reserved) | set(dimers[-1])
    sites = sorted(available | set(cell))
    contacts = sorted(
        edge(a, b) for a in sites for b in board.neighbours[a] if b in sites and a < b
    )
    cut = [link for link in contacts if len(set(link) & set(dimers[-1])) == 1]
    return board, dimers, cell, sorted(available), sites, contacts, cut


def analyze_clock(h, initial, effects, limit, cold_read=None, waits=outlet.WAITS):
    energies, v, idle = transmission.aged(h, initial)
    preparation = None
    if cold_read is not None:
        accepted = v.T @ cold_read @ v
        post = accepted @ idle @ accepted.conj().T
        admission = float(np.trace(post).real)
        assert 0 < admission <= 1 + 1e-12
        idle = post / admission
        cold_ground = v[:, 0]
        preparation = dict(
            admission=admission,
            static_ground_cold_clock=float(
                np.linalg.norm(cold_read @ cold_ground) ** 2
            ),
            static_ground_inventory_one=float(
                np.vdot(cold_ground, effects[1:4].sum(0) @ cold_ground).real
            ),
        )

    gates = [v.T @ np.kron(op, np.eye(len(energies) // 2)) @ v for op in BASIS]
    readers = [v.T @ e @ v for e in effects]
    # Pairs of local amplitudes span every SU(2) intervention on D.
    states = {
        (i, j): gate @ idle @ gates[j].conj().T
        for i, gate in enumerate(gates)
        for j in range(i + 1)
    }
    change = energies[:, None] - energies[None, :]
    resolution = outlet.PERIOD / LABELS
    indicators = [
        abs(change) < resolution,
        change <= -resolution,
        (change >= resolution) & (change < limit),
        change >= limit,
    ]
    records = {}
    extrema = []
    for a, (name, gate) in enumerate(zip(NAMES, gates, strict=True)):
        initial_weights = idle.diagonal().real
        two_readings = abs(gate) ** 2 * initial_weights[None, :]
        classes = [float(two_readings[mask].sum()) for mask in indicators]
        records[name] = dict(
            work=float(energies @ (states[a, a] - idle).diagonal().real),
            unresolved_loss_light_gain_heavy_gain=classes,
            panels=[],
        )
        np.testing.assert_allclose(sum(classes), 1, atol=3e-12)
    for t in waits:
        spin = np.exp(1j * energies * t)
        active = [r * spin[:, None] * spin.conj()[None, :] for r in readers]
        light = active[2]
        form = np.empty((4, 4))
        for (i, j), state in states.items():
            form[i, j] = form[j, i] = np.einsum("ij,ji", light, state).real
            if i == j:
                p = [np.einsum("ij,ji", read, state).real for read in active]
                np.testing.assert_allclose(sum(p), 1, atol=3e-12)
                records[NAMES[i]]["panels"].append(np.maximum(0, p).round(13).tolist())
        extrema.append(float(eigh(form, eigvals_only=True)[-1] - form[0, 0]))
    for name in NAMES[1:]:
        np.testing.assert_allclose(
            np.einsum(
                "kij,ji->k", readers, states[NAMES.index(name), NAMES.index(name)]
            ),
            np.einsum("kij,ji->k", readers, idle),
            atol=3e-12,
        )
    result = dict(fixed=records, all_D_unitaries_max_light_excess=extrema)
    if preparation is not None:
        result["cold_preparation"] = preparation
    return result


def acquire():
    board, dimers, cell, available, sites, contacts, cut = aperture()
    receiving = [e for e in contacts if set(e) <= set(cell)]
    sectors = outlet.receiving_cell(board, cell, receiving)
    _, one, cv = sectors[1]
    core, wing, *_ = hamiltonians(4)
    limit = np.ptp(eigh(core, eigvals_only=True)) / 10
    response = outlet.clock(
        one, one[0], limit, labels=LABELS, low_edge=outlet.PERIOD / LABELS
    )
    soft = (one - one[0] > 1e-8) & (one - one[0] < limit)
    ceiling = float(response[~soft, 1].max())
    values, v = eigh(wing)
    wing_reader = circuit(v, 33)
    guide = np.zeros((4, 2, 2, 2))
    guide[0, 0] = wing_reader[:, 0].reshape(2, 2)
    calibration = {}
    for name, op in zip(NAMES, BASIS, strict=True):
        lines = abs(wing_reader.T @ np.kron(op, np.eye(2)) @ wing_reader[:, 0]) ** 2
        calibration[name] = dict(
            lines=lines.round(13).tolist(), input=float(lines @ (values - values[0]))
        )
    settings = {}
    for name in ("native", "exit_cut", "no_reversal"):
        seams = [e for e in contacts if name != "exit_cut" or e not in cut]
        words, h = outlet.connected_wing(board, dimers, sites, seams)
        if name == "no_reversal":
            h += np.diag(
                [int(d == 1 and dimers[-1][0] in w) for d in (0, 1) for w in words]
            )
        effects = transmission.outlet_effects(
            outlet.read_groups(words, sectors, cell), sectors, response
        )
        initial = outlet.inject(guide, words, dimers, cell, cv[:, :1])[:, 0]
        settings[name] = analyze_clock(h, initial, effects, limit)
    return dict(
        protocol_sha256=sha256((HERE / "OPEN_PORTS.md").read_bytes()).hexdigest(),
        available=available,
        sites=sites,
        blank_contacts=contacts,
        cut=cut,
        threshold=limit,
        clock_false_light=ceiling,
        waits=list(outlet.WAITS),
        isolated_D_sources=calibration,
        settings=settings,
    )


def verdict(data):
    selected = []
    idx = data["waits"].index(32)
    cut = data["settings"]["exit_cut"]["fixed"]
    for name in NAMES[1:]:
        n = data["settings"]["native"]["fixed"]
        pulse, idle = (np.array(n[key]["panels"])[:, 2] for key in (name, "idle"))
        contrast = pulse - idle
        cut_contrast = (
            np.array(cut[name]["panels"])[:, 2] - np.array(cut["idle"]["panels"])[:, 2]
        )
        np.testing.assert_allclose(cut_contrast, 0, atol=3e-12)
        _, _, light, heavy = n[name]["unresolved_loss_light_gain_heavy_gain"]
        if (
            contrast[idx] - data["clock_false_light"] > 0.01
            and n[name]["work"] > 0
            and light > heavy
        ):
            selected.append(name)
        print(
            name,
            "work",
            n[name]["work"],
            "light/heavy gain",
            light,
            heavy,
            "light excess all waits",
            contrast.round(7),
        )
    print(
        "all single-D rotations",
        data["settings"]["native"]["all_D_unitaries_max_light_excess"],
    )
    print("fixed qualifying sources", selected)
    return selected


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    result = acquire()
    verdict(result)
    if args.acquire:
        with args.acquire.open("x") as f:
            json.dump(result, f, allow_nan=False, separators=(",", ":"))
            f.write("\n")
    else:
        compare_records(json.loads((HERE / "open_ports.json").read_text()), result)
