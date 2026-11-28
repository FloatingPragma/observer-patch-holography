"""Compare the complete passive heat records at three native meter settings."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from . import calorimeter
from ..material import edge, route
from ..measure import rebuild
from ...joined_decay.measure import core_pulse
from ...passive_decay.material import columns, commission
from ...passive_decay.measure import hamiltonians, physical
from ...passive_decay.readers import circuit, rotate

HERE = Path(__file__).resolve().parent
MOUNTS = ((32, 32), (33, 32), (32, 33))
TOL = 1e-10


def spectator_edges(board, dimers, heads):
    sites = set(sum(dimers, ()))
    return [
        edge(
            p,
            min(
                q
                for q in board.neighbours[p]
                if q not in sites and edge(p, q) not in heads
            ),
        )
        for p in dimers[-1]
    ]


def meter_clock(h, frame):
    """Interrogate a balanced superposition with two relative-phase analyzers."""
    values, vectors = eigh(h)
    prepared = frame @ np.array([1, 1]) / np.sqrt(2)
    incoming = vectors.T @ prepared
    angle = 0.0
    holds = []
    for t in 4.0 ** np.arange(7) / 64:
        state = frame.T @ (vectors @ (np.exp(-1j * t * values) * incoming))
        pc, ps = (abs(state[0] + phase * state[1]) ** 2 / 2 for phase in (1, 1j))
        phase = np.arctan2(2 * ps - 1, 2 * pc - 1)
        angle = phase + 2 * np.pi * round((4 * angle - phase) / (2 * np.pi))
        holds.append([float(t), float(pc), float(ps)])
    return dict(wait_cosine_sine=holds, gap=float(angle / t))


def acquire():
    board, dimers, heads, emitter, distances = calorimeter.apparatus()
    base = rebuild(board, dimers[:-1], route((3,), crossing=True)[2], False, True)
    settings = calorimeter.generators(board, dimers, heads, emitter, base)
    _, guide = eigh(base)
    _, core = eigh(hamiltonians(4)[0])
    pulse = np.kron(circuit(core_pulse(core), 32), np.eye(8))
    guide_stocks = {"pulse": pulse @ guide[:, 0], "idle": guide[:, 0]}
    core_frame = circuit(core, 32)
    contacts = spectator_edges(board, dimers, heads)
    protocol = (HERE / "DETUNING.md").read_bytes()
    result = dict(
        protocol_sha256=sha256(protocol).hexdigest(),
        distances=distances,
        spectator_edges=contacts,
        waits=list(calorimeter.WAITS),
        column_errors=[columns(), columns(comparison=True)],
        mounts=[],
    )
    base_meter = calorimeter.native_meter()
    for costs in MOUNTS:
        pins = dict(zip(contacts, costs, strict=True))
        bias = np.diag([0.0, float(costs[1] - costs[0])])
        meter = base_meter + bias
        _, eigenvectors = eigh(meter)
        ground = eigenvectors[:, 0]
        frame = rotate(np.arctan2(ground[1], ground[0]))
        clock = meter_clock(meter, frame)
        np.testing.assert_allclose(
            clock["gap"], np.ptp(np.linalg.eigvalsh(meter)), atol=TOL, rtol=0
        )
        reader = np.kron(core_frame, np.kron(np.eye(8), frame))
        mount = dict(
            costs=list(costs),
            clock=clock,
            preparation=commission(board, dimers, heads | pins),
            settings={},
        )
        for name, deposits, reference in settings:
            if name == "no_reversal":
                changed_pins = dict.fromkeys(contacts, 32)
                h = reference
            else:
                changed_pins = pins
                h = reference + np.kron(np.eye(32), bias)
            diagonal = physical(board, dimers, deposits | changed_pins)
            np.testing.assert_array_equal(diagonal - diagonal[0], np.diag(h))
            values, vectors = eigh(h)
            setting = {}
            for arm, stock in guide_stocks.items():
                initial = np.kron(stock, frame[:, 0])
                incoming = vectors.T @ initial
                laws = []
                for wait in calorimeter.WAITS:
                    state = vectors @ (np.exp(-1j * wait * values) * incoming)
                    law = (abs(reader.T @ state) ** 2).reshape(4, 8, 2).sum(axis=1)
                    np.testing.assert_allclose(law.sum(), 1.0, atol=TOL, rtol=0)
                    laws.append(law.round(13).tolist())
                setting[arm] = laws
            mount["settings"][name] = setting
        result["mounts"].append(mount)
    return result


def summary(data):
    """Read the specified contrasts from all eight detector bins."""
    idx = data["waits"].index(32)
    for mount in data["mounts"]:
        panels = mount["settings"]

        def hot(setting, arm):
            return float(np.asarray(panels[setting][arm])[idx, :, 1].sum())

        background = max(
            hot(setting, arm)
            for setting, stocks in panels.items()
            for arm in stocks
            if (setting, arm) != ("native", "pulse")
        )
        joint = np.asarray(panels["native"]["pulse"])[idx, 1:3, 1].sum()
        print(
            "costs",
            mount["costs"],
            "gap",
            round(mount["clock"]["gap"], 7),
            "joint",
            round(joint, 7),
            "max control hot",
            round(background, 7),
            "difference",
            round(joint - background, 7),
        )
        print(
            "pulse excess at each setting",
            {s: round(hot(s, "pulse") - hot(s, "idle"), 7) for s in panels},
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.acquire:
        record = acquire()
        with args.acquire.open("x") as output:
            json.dump(record, output, separators=(",", ":"), allow_nan=False)
            output.write("\n")
    else:
        record = json.loads((HERE / "detuning.json").read_text())
    summary(record)
