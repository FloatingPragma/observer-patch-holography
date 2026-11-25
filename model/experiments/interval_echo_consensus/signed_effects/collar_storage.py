"""Encode a neutral material pointer as a free, conserved seam-memory record."""

import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.experiments.material_inertia.dark_collar import label
from .sequential import require

LAW = candidate.native_law
RECEIPT = Path(__file__).with_suffix(".json")
PAIR = LAW.SRC.index((1, 1))
BITS = [LAW.SRC.index(q) for q in ((2, 0), (0, 2))]
MEMORIES = np.stack([label(2, parity) for parity in (1, -1)], axis=1)
POINTER_OWNERS = (1, 2, 4, 5)


def geometry():
    """Use the actual four pointer seams and two disjoint retrieval contacts apiece."""
    board = federation(0)
    stations = []
    for owner in POINTER_OWNERS:
        a, b, x, y = [12 * owner + j for j in (4, 5, 3, 11)]
        require(board.degree[b] == 5, "pointer reference unavailable")
        require(b in board.neighbours[a] and x in board.neighbours[a]
                and y in board.neighbours[b], "pointer service outside incidence")
        stations.append(dict(pointer=[a, b], scratch=[x, y]))
    return stations


@lru_cache(maxsize=None)
def gates(exchange):
    hopping, contrast = LAW.local_term_parts()
    h = exchange * hopping + candidate.KAPPA * contrast
    reference = np.diag(np.tile([0., 8., 32.], LAW.DIM // 3))
    z = expm(-1j * np.pi / 64 * reference)
    basis = z @ expm(-0.5j * np.pi * h) @ z
    # A calibrated blank contact exchanges a scratch load with the endpoint.
    # Both compared laws have unit contrast phase at their swap time.
    swap = expm(-1j * np.pi / exchange * h)[:, :LAW.N_SRC]
    transposed = [LAW.SRC.index((b, a)) for a, b in LAW.SRC]
    np.testing.assert_allclose(swap, np.eye(LAW.DIM)[:, transposed], atol=3e-13, rtol=0)
    return h, basis, LAW.echo()


def conserved_columns():
    h, basis, echo = gates(1.)
    embedded = np.kron(MEMORIES, np.eye(LAW.N_SRC))
    kinetic = np.tile(np.arange(LAW.N_SRC) != PAIR, 2)
    energies = np.repeat([32., 33.], LAW.N_SRC) * kinetic
    np.testing.assert_allclose(h @ embedded, embedded * energies, atol=2e-13, rtol=0)
    target = np.kron(MEMORIES, np.eye(LAW.N_SRC)[:, [PAIR]]) * np.exp(-1j * np.pi / 4)
    np.testing.assert_allclose((echo @ basis)[:, BITS], target, atol=8e-14, rtol=0)
    return [[int(q), int(even), int(odd)]
            for q, (even, odd) in enumerate(zip(*energies.reshape(2, -1), strict=True))]


def recovered(exchange, encoded, turns):
    """Native encode, swap in an arbitrary ternary pair, dwell, restore, decode.

    The scratch pair carries the stored (1,1) during the dwell. At the inverse
    swaps it recovers the displaced loads. A nonblank final seam is a reject.
    """
    h, basis, echo = gates(exchange)
    encoder = echo @ basis if encoded else echo
    write = encoder[:, BITS].reshape(11, 9, 2)
    np.testing.assert_allclose(write[:, np.arange(9) != PAIR, :], 0, atol=3e-14, rtol=0)
    # Two disjoint, calibrated scratch swaps place the supplied loads on the
    # pointer. Columns below span both bits and all nine scratch words.
    start = np.zeros((11, 9, 9, 2), complex)
    for q in range(9):
        start[:, q, q, :] = write[:, PAIR, :]
    state = (expm(-1j * turns * np.pi * h) @ start.reshape(99, 18)).reshape(11, 9, 9, 2)
    # Check all nonmatching material outputs before undoing the scratch swaps.
    restored = np.zeros_like(state)
    for q in range(9):
        np.testing.assert_allclose(state[:, np.arange(9) != q, q], 0, atol=2e-13, rtol=0)
        restored[:, PAIR, q, :] = state[:, q, q, :]
    out = (encoder.conj().T @ restored.reshape(99, 18)).reshape(99, 9, 2)
    outcomes = np.concatenate([abs(out[:9])**2, np.sum(abs(out[9:])**2, axis=0,
                                                        keepdims=True)])
    np.testing.assert_allclose(outcomes.sum(axis=0), 1, atol=3e-13, rtol=0)
    return outcomes


def report():
    rng = np.random.default_rng(518117)
    shots = 4096
    table = []
    for exchange in (1., .5):
        for encoded in (True, False):
            for turns in (0, 1, 2):
                outcomes = recovered(exchange, encoded, turns)
                banks = []
                for q, scratch in enumerate(LAW.SRC):
                    for bit in range(2):
                        p = outcomes[:, q, bit].round(12)
                        require(np.isclose(p.sum(), 1., rtol=0, atol=1e-11),
                                "incomplete terminal distribution")
                        counts = rng.multinomial(shots, p / p.sum()).tolist()
                        banks.append(dict(scratch_word=list(scratch), bit=bit, counts=counts))
                table.append(dict(exchange=exchange, parity_encoding=encoded,
                                  dwell_in_pi=turns, banks=banks))
    # Both scratch branches have total load two. After a pi dwell, the
    # nonbalanced branch reverses the relative phase of the two memory bits.
    h, _, _ = gates(1.)
    scratch = np.eye(9)[:, [PAIR, LAW.SRC.index((2, 0))]].sum(axis=1) / np.sqrt(2)
    plus = MEMORIES.sum(axis=1) / np.sqrt(2)
    state = (expm(-1j * np.pi * h) @ np.kron(plus, scratch)).reshape(11, 9)
    memory = MEMORIES.T @ state
    reduced = memory @ memory.conj().T
    np.testing.assert_allclose(reduced, np.eye(2) / 2, atol=2e-13, rtol=0)
    return dict(stations=geometry(), native_energies=conserved_columns(),
                equal_total_load_probe_memory_purity=float(np.trace(reduced @ reduced).real.round(12)),
                outcome_labels=[list(q) for q in LAW.SRC] + ["written collar"],
                starts_per_input=shots, seed=518117, banks=table)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as stream:
            stream.write(json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n")
    else:
        require(json.loads(RECEIPT.read_text()) == result, "stored pointer acquisition differs")
        print("Conserved pointer records replay on all nine endpoint disturbances.")


if __name__ == "__main__":
    main()
