"""Realize coexistent complex effects with a single load on a twelve-port cell."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell


def port_path():
    board = cell()

    def extend(path):
        if len(path) == board.n:
            return path
        for p in sorted(set(board.neighbours[path[-1]]) - set(path)):
            complete = extend(path + [p])
            if complete:
                return complete
        return None

    path = extend([0])
    assert path is not None and all(board.degree[p] == 5 for p in path)
    return list(map(int, path))


def native_block():
    """Check the isolated one-position columns and the empty spectators."""
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    h = exchange + candidate.KAPPA * contrast
    basis = np.eye(law.DIM)[:, [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]]
    block = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ block, atol=2e-13, rtol=0)
    np.testing.assert_allclose(block, [[8.5, -0.5], [-0.5, 8.5]], atol=2e-13, rtol=0)
    np.testing.assert_allclose(h[:, law.SRC.index((0, 0))], 0, atol=2e-13)
    return block


def effect_root(effect):
    np.testing.assert_allclose(effect, effect.conj().T, atol=1e-13, rtol=0)
    values, vectors = np.linalg.eigh(effect)
    if values.min() < -1e-13:
        raise ValueError("Negative effect or complement")
    return (vectors * np.sqrt(np.maximum(values, 0))) @ vectors.conj().T


def compile_isometry(target, block):
    """Adjacent row elimination and inverse playback using actual seam evolution."""
    rails, dimension = target.shape
    assert rails <= len(port_path())
    columns = np.eye(rails, dimension, dtype=complex)
    np.testing.assert_allclose(target.conj().T @ target, np.eye(dimension), atol=2e-13)
    reduced = np.array(target, dtype=complex, copy=True)
    steps = []
    for col in range(dimension):
        for j in range(rails - 1, col, -1):
            a, b = reduced[j - 1 : j + 1, col]
            if abs(b) < 2e-15:
                continue
            phase = np.exp(-1j * np.angle([a, b]))
            theta = np.arctan2(abs(b), abs(a))
            c, s = np.cos(theta), np.sin(theta)
            reduced[j - 1 : j + 1] = np.array([[c, s], [-s, c]]) @ (
                phase[:, None] * reduced[j - 1 : j + 1]
            )
            steps.append((j, theta, phase))

    diagonal = np.diag(reduced[:dimension])
    columns[:dimension] = np.diag(np.exp(1j * np.angle(diagonal)))
    np.testing.assert_allclose(reduced, columns, atol=4e-14, rtol=0)
    for j, theta, phase in reversed(steps):
        # Inverse real Givens from a positive exposure, with addressed degree phases.
        dwell = (-2 * theta) % (4 * np.pi)
        contact = expm(-1j * dwell * block)
        flip = np.diag([1, 1j])
        rotation = np.exp(8.5j * dwell) * flip @ contact @ flip.conj().T
        columns[j - 1 : j + 1] = (phase.conj()[:, None] * rotation) @ columns[
            j - 1 : j + 1
        ]
    np.testing.assert_allclose(columns, target, atol=3e-12, rtol=0)
    return columns, len(steps)


def three_effects(first, second):
    size = len(first)
    return np.vstack(
        [effect_root(e) for e in (first, second, np.eye(size) - first - second)]
    )


def report():
    path = port_path()
    block = native_block()
    identity = np.eye(4)
    u = np.array([1, 1j, 1, -1j]) / 2
    v = np.array([1, 0, 1j, 0]) / np.sqrt(2)
    first = identity / 10 + np.outer(u, u.conj()) / 4
    second = identity / 7 + np.outer(v, v.conj()) / 3
    target = three_effects(first, second)
    joint, joint_count = compile_isometry(target, block)
    # Independent binary circuit and an output basis change on the yes block.
    binary, binary_count = compile_isometry(
        three_effects(first + second, 0 * first), block
    )
    fourier = np.exp(0.5j * np.pi * np.outer(np.arange(4), np.arange(4))) / 2
    rotated_target = target.copy()
    rotated_target[:4] = fourier @ rotated_target[:4]
    rotated, rotated_count = compile_isometry(rotated_target, block)
    expected = (first, second, first + second)
    effects = [rows.conj().T @ rows for rows in (joint[:4], joint[4:8], binary[:4])]
    residual = max(
        np.linalg.norm(a - b, 2) for a, b in zip(expected, effects, strict=True)
    )
    np.testing.assert_allclose(rotated[:4].conj().T @ rotated[:4], first, atol=2e-12)
    rng = np.random.default_rng(47343)
    records = {}
    for name, instrument in (
        ("joint", joint),
        ("binary", binary),
        ("rotated", rotated),
    ):
        probabilities = abs(instrument @ u) ** 2
        counts = rng.multinomial(2**18, probabilities / probabilities.sum())
        records[name] = {
            "counts_in_path_order": counts.tolist(),
            "shots": int(counts.sum()),
        }
    return dict(
        path=path,
        input_ports=path[:4],
        dimension=4,
        joint_output_blocks=[path[i : i + 4] for i in (0, 4, 8)],
        adjacent_contact_counts=[joint_count, binary_count, rotated_count],
        effect_operator_error=float(residual),
        commutator_norm=float(np.linalg.norm(first @ second - second @ first, 2)),
        disjoint_probabilities=[float(np.vdot(u, e @ u).real) for e in expected],
        counted_preparation="(1,i,1,-i)/2 on the four input ports",
        simulated_counts=records,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(report(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_suffix(".json").write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
