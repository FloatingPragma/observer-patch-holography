"""Phase-sensitive and loop-invariant instruments on transported vacancy fibers."""

import argparse
import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from .holonomy import HERE, analyze_bright, move, witness


def exchange_block():
    indices = [law.SRC.index(word) for word in ((2, 1), (1, 2))]
    return law.seam_term(8)[np.ix_(indices, indices)]


def rotation(theta):
    """Two addressed phases remove the pair offset relative to other sites."""
    dwell = 2 * theta
    flip = np.diag([1, 1j])
    return (
        np.exp(8.5j * dwell)
        * flip
        @ expm(-1j * dwell * exchange_block())
        @ flip.conj().T
    )


def analyzer(vector):
    """Fold a three-port target to its first port using its two native bonds."""
    result = np.eye(3, dtype=complex)
    for j in (1, 0):
        a, b = (result @ vector)[j : j + 2]
        phase = np.diag(np.exp(-1j * np.angle([a, b])))
        result[j : j + 2] = (
            rotation(np.arctan2(abs(b), abs(a))) @ phase @ result[j : j + 2]
        )
    np.testing.assert_allclose(result @ vector, [1, 0, 0], atol=4e-14)
    return result


def transport(vector, perm):
    output = np.empty_like(vector)
    output[perm] = vector
    return output


def folded_click(vector, root, rotations):
    work = vector.copy()
    for parent, child, a, b in rotations:
        # The exact trigonometric entries of the implemented pair rotation.
        c, s = np.sqrt(np.array([a, b]) / (a + b))
        p, q = work[parent], work[child]
        work[parent] = c * p + s * q
        work[child] = -s * p + c * q
    return float(abs(work[root]) ** 2)


def records(level=1, shots=2**18, seed=511307):
    board = federation(level)
    info, perm = witness(board)
    root = info["root"]
    source, receiver = info["source_triple"], info["receiving_triple"]
    for triple in (source, receiver):
        assert all(b in board.neighbours[a] for a, b in zip(triple, triple[1:]))
    u = np.ones(3) / np.sqrt(3)
    omega = np.exp(2j * np.pi / 3)
    plus = np.array([1, omega, omega**2]) / np.sqrt(3)
    probe = (u + plus) / np.sqrt(2)
    meter = (u + 1j * plus) / np.sqrt(2)
    prep, read = analyzer(probe).conj().T, analyzer(meter)
    v = np.zeros(board.n, dtype=complex)
    v[source] = prep[:, 0]
    received = transport(v, perm)
    np.testing.assert_allclose(received[receiver], probe, atol=5e-14)
    local = info["wedge_roots"][1]
    approach = info["approaches"][1]
    ab, bc = [
        move(board, approach + [a, b, local] + approach[-2::-1])
        for a, b in (receiver[:2], receiver[1:])
    ]
    rng = np.random.default_rng(seed)

    def clicks(probability):
        yes = int(rng.binomial(shots, probability))
        return [yes, shots - yes]

    orders = []
    for first, second in ((ab, bc), (bc, ab)):
        final = transport(transport(received, first), second)[receiver]
        p = float(abs((read @ final)[0]) ** 2)
        orders.append(dict(probability=p, counts=clicks(p)))

    # Send both columns of the local zero-sum doublet; retain the referee.
    basis = np.stack((plus, plus.conj()), axis=1)
    home = np.zeros((board.n, 2), dtype=complex)
    home[source] = basis
    received_code = transport(home, perm)
    overlap = basis.conj().T @ received_code[receiver]
    fidelity = float(abs(np.trace(overlap) / 2) ** 2)
    dephased = float(sum((np.vdot(row, row).real / 2) ** 2 for row in received_code))
    bright = np.ones(board.n) / np.sqrt(board.n - 1)
    bright[root] = 0
    dark = np.zeros(board.n)
    dark[source[:2]] = [1 / np.sqrt(2), -1 / np.sqrt(2)]
    prepared = np.sqrt(3 / 8) * bright + np.sqrt(5 / 8) * dark
    read_port, rotations = analyze_bright(board, root)
    p = folded_click(prepared, read_port, rotations)
    p_after = folded_click(transport(prepared, perm), read_port, rotations)
    np.testing.assert_allclose([p, p_after], 3 / 8, atol=2e-14)
    # Terminal click followed by dyadic routing with fair native bits.
    bins = p_after * np.array([3 / 4, 1 / 8, 1 / 8]) + (1 - p_after) * np.array(
        [1 / 8, 5 / 8, 1 / 4]
    )
    return dict(
        board_ports=board.n,
        root=root,
        source=source,
        receiver=receiver,
        exchanges=info["contact_pulses"],
        shots=shots,
        seed=seed,
        loop_orders=orders,
        doublet_entanglement_fidelity=fidelity,
        dephased_doublet_entanglement_fidelity=dephased,
        bright_probability=p,
        transported_bright_probability=p_after,
        ternary_probabilities=bins.tolist(),
        ternary_counts=rng.multinomial(shots, bins).tolist(),
        direct_coarse_probability=float(sum(bins[:2])),
        direct_coarse_counts=clicks(sum(bins[:2])),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    text = json.dumps(records(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "holonomy_counts.json").write_text(text)
    else:
        print(text, end="")
