"""Two coherent arms take opposite paths through a guarded face's native private exposures."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation, rotations
from ..measure import compile_turn, half_turns


HERE = Path(__file__).resolve().parent
KAPPA = candidate.KAPPA
STEP = 0.2
STRENGTH = 0.35
SLICES = (64, 256, 1024)
SHOTS = 2**22
ANGLES = np.arange(3) * 2 * np.pi / 3
# One path from +x to -x in each arm surrounds the degeneracy once.
SETTINGS = {"enclosing": (0, 1), "displaced": (2, 1), "twice": (0, 2), "parked": (0, 0)}


def private_face():
    board = federation(0)
    return next(face for face in board.faces()
                if all(board.degree[p] == 5 for p in face))


def terms(sign, contrast_strength=KAPPA):
    """Extract native seam terms on all single-defect words; check closure column by column."""
    exchange, contrast = candidate.native_law.local_term_parts()
    local = exchange + contrast_strength * contrast
    words = np.ones((3, 3), dtype=int)
    np.fill_diagonal(words, 1 + sign)
    word_index = {tuple(word): index for index, word in enumerate(words)}
    edges = {}
    for a, b in combinations(range(3), 2):
        matrix = np.zeros((3, 3))
        for column, word in enumerate(words):
            initial = 3 * word[a] + word[b]
            if np.max(np.abs(local[9:, initial])) > 2e-13:
                raise ArithmeticError("A seam wrote outside the blank code")
            for na, nb in product(range(3), repeat=2):
                value = local[3 * na + nb, initial]
                if abs(value) < 2e-13:
                    continue
                output = word.copy()
                output[[a, b]] = na, nb
                matrix[word_index[tuple(output)], column] += value
        edges[a, b] = matrix
    return edges, words


def prepare(edges, words):
    state = expm(-0.5j * np.pi * edges[1, 2])[:, 1]
    # A common private exposure on port 2 sets the same relative phase for either load.
    return np.exp(-0.5j * np.pi * words[:, 2]**2) * state


def propagate(state, step, words, directions, offset, contrast_strength=KAPPA):
    phases = np.cos(directions[:, None] - ANGLES) + offset * np.cos(ANGLES)
    dwell = STEP * STRENGTH * (3 + phases) / KAPPA
    if dwell.min() < -1e-15:
        raise ArithmeticError("Native exposures must be nonnegative")
    # Actual n^2 charges, including the occupied background's common dynamical phase.
    half_kicks = np.exp(-0.5j * contrast_strength * dwell @ (words**2).T)
    for kick in half_kicks:
        state = kick * (step @ (kick * state))
    return state


def probabilities(left, right):
    # Each bin resolves both the interferometer exit and the address of the signed defect.
    joint = np.abs(np.stack((left + right, left - right)))**2 / 4
    np.testing.assert_allclose(joint.sum(), 1, atol=1e-11, rtol=0)
    return joint / joint.sum()


def bank(joint, rng):
    return {
        "probabilities": joint.round(12).tolist(),
        "counts": rng.multinomial(SHOTS, joint.flatten()).reshape(joint.shape).tolist(),
    }


def spatial_reads(face, state, edge, rng):
    """Read every address after two full-cell spatial interventions on the same orbital."""
    cycle = min(r for r in rotations() if r[face[0]] == face[1] and r[face[1]] == face[2])
    turn = compile_turn(cycle)
    a, b = half_turns()
    pulse = expm(-1j * np.pi * edge)[:2, :2]
    initial = np.zeros(12, dtype=complex)
    initial[list(face)] = state

    def arm(word):
        result = initial.copy()
        for endpoints in word:
            ports = list(endpoints)
            result[ports] = pulse @ result[ports]
        return result

    return {"full_turn": bank(probabilities(initial, arm(turn * 3)), rng),
            "half_turn_orders": bank(probabilities(arm(a + b), arm(b + a)), rng),
            "face_turn_permutation": list(cycle)}


def scan(initial, edges, words, contrast, rng):
    step = expm(-1j * STEP * sum(edges.values()))
    results = {}
    for setting, (offset, half_turns) in SETTINGS.items():
        banks = {}
        for slices in SLICES:
            directions = np.pi * half_turns * (np.arange(slices) + 0.5) / slices
            left = propagate(initial, step, words, directions, offset, contrast)
            right = propagate(initial, step, words, -directions, offset, contrast)
            banks[str(slices)] = bank(probabilities(left, right), rng)
        results[setting] = banks
    return results


def acquire():
    rng = np.random.default_rng(20261125)
    face = private_face()
    record = {"carrier": "20 cells", "face": list(face), "contrast": KAPPA,
              "seam_exposure": STEP, "strength": STRENGTH, "shots_per_bank": SHOTS,
              "seed": 20261125, "exits": ["bright", "dark"], "runs": {},
              "changed_drive_contrast_zero": {}}
    for sign in (-1, 1):
        edges, words = terms(sign)
        state = prepare(edges, words)
        results = {"material_spatial_reads": spatial_reads(face, state, edges[0, 1], rng)}
        results.update(scan(state, edges, words, KAPPA, rng))
        record["runs"][str(sign)] = results
    for sign in (-1, 1):
        edges, words = terms(sign)
        changed_edges, _ = terms(sign, 0)
        record["changed_drive_contrast_zero"][str(sign)] = scan(
            prepare(edges, words), changed_edges, words, 0, rng
        )
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(acquire(), indent=2) + "\n")
