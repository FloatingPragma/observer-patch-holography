"""Independently count powers at the triangle feeding the remote calorimeter."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from .calorimeter import HERE, apparatus, generators
from ..material import route
from ..measure import rebuild
from ...joined_decay.measure import core_pulse
from ...passive_decay.measure import hamiltonians
from ...passive_decay.readers import X, Z, rotate

T = 32
BINS = 1024
SHOTS = 2**16
SEED = 26112610


def acquire():
    board, dimers, heads, emitter, _ = apparatus()
    original_heads = route((3,), crossing=True)[2]
    original = rebuild(board, dimers[:-1], original_heads, False, True)
    settings = generators(board, dimers, heads, emitter, original)
    _, core = eigh(hamiltonians(4)[0])
    pulse = np.kron(core_pulse(core), np.eye(8))
    _, guide = eigh(original)
    cold = guide[:, 0]
    meter_ground = (rotate(np.pi / 4) @ Z)[:, 0]
    stocks = (
        ("pulse", np.kron(pulse @ cold, meter_ground)),
        ("idle", np.kron(cold, meter_ground)),
    )
    y = np.array([[0, -1j], [1j, 0]])
    # Native phase on the occupied rail takes X exchange eigenvectors to Y read ports.
    y_frame = np.diag([1, 1j]) @ (rotate(np.pi / 4) @ Z)
    np.testing.assert_allclose(
        y_frame.conj().T @ y @ y_frame, np.diag([1, -1]), atol=1e-13
    )
    assert emitter == dimers[-2][1]
    power = np.kron(np.eye(16), np.kron(np.diag([1, 0]), y / 2))
    energy = np.kron(np.eye(32), (np.eye(2) - X) / 2)
    rng = np.random.default_rng(SEED)
    times = (np.arange(BINS) + 0.5) * T / BINS
    panels = []
    for name, _, h in settings[:2]:
        np.testing.assert_allclose(1j * (h @ energy - energy @ h), power, atol=1e-13)
        comm = h @ power - power @ h
        second = h @ comm - comm @ h
        bound = float(np.linalg.norm(second, 2)) * T * (T / BINS) ** 2 / 24
        values, v = eigh(h)
        for arm, psi in stocks:
            arriving = v @ (
                np.exp(-1j * values[:, None] * times) * (v.T @ psi)[:, None]
            )
            output = np.einsum(
                "fi,defb->deib", y_frame.conj(), arriving.reshape(16, 2, 2, BINS)
            )
            probabilities = (abs(output) ** 2).sum(axis=0)
            # Counts are in the order negative, zero, positive.
            bins = np.array(
                [probabilities[0, 1], probabilities[1].sum(axis=0), probabilities[0, 0]]
            ).T
            counts = [rng.multinomial(SHOTS, p / p.sum()).tolist() for p in bins]
            panels.append(
                dict(apparatus=name, source=arm, quadrature_bound=bound, counts=counts)
            )
    return dict(
        seed=SEED,
        bins=BINS,
        shots_per_bin=SHOTS,
        stop=T,
        sampling_half_width=float(
            T * np.sqrt(np.log(2 * 4 / 0.01) / (2 * BINS * SHOTS))
        ),
        panels=panels,
    )


def decode(data, panel):
    counts = np.array(panel["counts"])
    assert counts.shape == (data["bins"], 3)
    assert np.all(counts.sum(axis=1) == data["shots_per_bin"])
    return (counts[:, 2] - counts[:, 0]).sum() * data["stop"] / (2 * counts.sum())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    args = parser.parse_args()
    if args.acquire:
        data = acquire()
        args.acquire.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    else:
        data = json.loads((HERE / "power_counts.json").read_text())
    for panel in data["panels"]:
        print(
            panel["apparatus"],
            panel["source"],
            decode(data, panel),
            data["sampling_half_width"] + panel["quadrature_bound"],
        )
