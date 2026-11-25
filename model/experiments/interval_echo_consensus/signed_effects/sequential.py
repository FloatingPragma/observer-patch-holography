"""Replay a complete signed native instrument, its poststates and repeat records."""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell
from . import measure

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECEIPT = HERE / "sequential.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


GRID_BITS = 40
GRID = 1 << GRID_BITS
SAMPLING = ("integer weights round(p*2**40) for every outcome except the most "
            "probable, which takes the complement; each draw is a PCG64 integer "
            "below 2**40 compared with the cumulative weights")


def grid_weights(probabilities, tolerance=3e-12):
    """Integer weights on the 2**-40 grid with an explicit complement.

    Propagated probabilities differ in their last bits between linear algebra
    libraries. Rounding moves each probability by at most 2**-41 < 5e-13, and
    the complement by the other roundings plus the norm defect. Both stay
    inside the 3e-12 norm budget that the propagation already checks.
    """
    p = np.asarray(probabilities, dtype=float).ravel()
    weights = np.rint(p * GRID).astype(np.int64)
    top = int(np.argmax(weights))
    weights[top] = GRID - (weights.sum() - weights[top])
    require(np.all(weights >= 0) and np.max(abs(weights / GRID - p)) <= tolerance,
            "sampling weights leave the propagated law")
    return weights


def grid_counts(rng, shots, probabilities, tolerance=3e-12):
    """Exact draws from weights/2**40 using integers only.

    No floating-point operation enters a draw, so a given seed gives the same
    counts on every platform. Zero-weight outcomes stay empty.
    """
    edges = np.cumsum(grid_weights(probabilities, tolerance))
    draws = rng.integers(0, GRID, size=shots, dtype=np.int64)
    return np.bincount(np.searchsorted(edges, draws, side="right"), minlength=edges.size)


def grid_binomial(rng, trials, probability):
    """Successes among trials with success weight on the 2**-40 grid."""
    return int(grid_counts(rng, trials, [probability, 1 - probability])[0])


def sources():
    paths = [
        "model/candidate_laws/contrast_memory_spherical_v1/law.py",
        "model/candidate_laws/interval_echo/candidate.py",
        "model/experiments/interval_echo_consensus/signed_effects/measure.py",
        "model/experiments/interval_echo_consensus/signed_effects/sequential.py",
    ]
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


@lru_cache(maxsize=2)
def native_generators(sign):
    """The full native 99-state column check precedes the twelve-position reduction."""
    require(sign in (-1, 1), "signed one-defect sector required")
    contacts, references = measure.native_controls()
    index = 0 if sign == -1 else 1
    contact = np.zeros((cell().n, cell().n))
    contact[:2, :2] = contacts[index]
    reference = 8 * np.eye(cell().n)
    reference[:2, :2] = references[index]
    require(1 in cell().neighbours[0], "source contact absent")
    return contact, reference


def word(sign, inverse=False):
    """Positive native waits implement U or U*, in chronological order."""
    contact, reference = native_generators(sign)
    pulses = [(contact, np.pi / 2), (reference, np.pi / 20),
              (contact, 3 * np.pi / 2)]
    if inverse:
        # Both actual generators have integer spectra and period 2*pi.
        pulses = [(h, (2 * np.pi - t) % (2 * np.pi))
                  for h, t in reversed(pulses)]
    u = np.eye(cell().n, dtype=complex)
    for h, t in pulses:
        u = expm(-1j * t * h) @ u
    return u


def instrument(sign, swapped=False):
    """One Kraus operator per exhaustive sharp address, including unused rails.

    The actual read is between U and its positive-wait inverse. A retained
    address is a classical label. Postread swapping changes no induced effect.
    """
    u, undo = word(sign), word(sign, inverse=True)
    require(np.linalg.norm(undo @ u - np.eye(cell().n), 2) < 2e-12,
            "positive-wait undo failed")
    contact, _ = native_generators(sign)
    after = expm(-1j * np.pi * contact) if swapped else np.eye(cell().n)
    ops = []
    for address in range(cell().n):
        pointer = np.zeros((cell().n, cell().n))
        pointer[address, address] = 1
        ops.append(after @ undo @ pointer @ u)
    require(np.linalg.norm(sum(k.conj().T @ k for k in ops) - np.eye(cell().n), 2)
            < 3e-12, "incomplete instrument")
    return ops


def joint_probabilities(ops, state):
    """Two uses retain every ordered address pair and their full residual state."""
    return np.array([[np.linalg.norm(second @ first @ state) ** 2
                      for second in ops] for first in ops])


def matrix_data(a):
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def report():
    n = cell().n
    inputs = {
        "port0": np.array([1, 0], complex),
        "port1": np.array([0, 1], complex),
        "plus": np.array([1, 1], complex) / np.sqrt(2),
        "plus_i": np.array([1, 1j], complex) / np.sqrt(2),
    }
    rng = np.random.default_rng(925031)
    shots = 16384
    cases = []
    effect_error = repeat_error = map_error = 0.0
    for sign in (-1, 1):
        plain = instrument(sign)
        changed = instrument(sign, swapped=True)
        effects = [k.conj().T @ k for k in plain]
        effect_error = max(effect_error, *(np.linalg.norm(k.conj().T @ k - e, 2)
                            for k, e in zip(changed, effects, strict=True)))
        # Audit the full twelve-position maps, not just their sampled inputs.
        for k, e in zip(plain, effects, strict=True):
            map_error = max(map_error, np.linalg.norm(k - e, 2))
        for name, small in inputs.items():
            state = np.r_[small, np.zeros(n - 2)]
            for arm, ops in (("repeatable", plain), ("postread_swap", changed)):
                law = joint_probabilities(ops, state)
                require(np.all(np.isfinite(law)) and np.all(law >= 0)
                        and abs(law.sum() - 1) < 5e-12, "invalid complete law")
                mismatch = float(law.sum() - np.trace(law))
                if arm == "repeatable":
                    repeat_error = max(repeat_error, abs(mismatch))
                else:
                    expected = (3 + sign * np.sqrt(5)) / 8
                    require(abs(mismatch - expected) < 3e-12,
                            "same-effects disturbance control missed")
                # Exact rank-one branching has the same conditional repeat
                # error after either first outcome. Sample its binary tree;
                # the full propagated law above is retained without rescaling.
                p0 = float(np.vdot(state, effects[0] @ state).real)
                disagreement = 0 if arm == "repeatable" else (3 + sign * np.sqrt(5)) / 8
                analytic = np.array([[p0 * (1 - disagreement), p0 * disagreement],
                                     [(1 - p0) * disagreement,
                                      (1 - p0) * (1 - disagreement)]])
                require(np.max(abs(law[:2, :2] - analytic)) < 3e-12,
                        "rank-one branching disagrees with full instrument")
                first0 = grid_binomial(rng, shots, p0)
                crossed0 = grid_binomial(rng, first0, disagreement)
                crossed1 = grid_binomial(rng, shots - first0, disagreement)
                counts = np.zeros((n, n), dtype=int)
                counts[:2, :2] = [[first0 - crossed0, crossed0],
                                  [crossed1, shots - first0 - crossed1]]
                cases.append({
                    "sign": sign, "preparation": name, "arm": arm,
                    "joint_address_probabilities": law.tolist(),
                    "joint_address_counts": counts.tolist(),
                    "shots": shots, "repeat_disagreement": mismatch,
                    "first_outcome_poststates_active_rails": [
                        matrix_data(np.outer((k @ state)[:2], (k @ state)[:2].conj()))
                        for k in ops[:2]],
                })
    require(effect_error < 3e-12 and map_error < 3e-12 and repeat_error < 3e-12,
            "instrument or same-effects gate failed")
    return {
        "construction": "Interval Echo signed protected contact/reference menu",
        "execution": "full native block reduction, matrix instrument, synthetic counts",
        "custody": "producer-side software hashes; no external authentication",
        "source_sha256": sources(), "seed": 925031,
        "sampling_grid_bits": GRID_BITS, "sampling": SAMPLING,
        "addresses": list(range(n)), "joint_order": ["first address", "second address"],
        "poststate_coordinates":
            "subnormalized matrices on rails 0,1; all other entries and "
            "first-outcome poststates 2..11 are zero on the four tested inputs",
        "chronological_U": [["contact 0,1", "pi/2"],
                            ["degree reference 0", "pi/20"],
                            ["contact 0,1", "3*pi/2"]],
        "chronological_undo": [["contact 0,1", "pi/2"],
                               ["degree reference 0", "39*pi/20"],
                               ["contact 0,1", "3*pi/2"]],
        "pointer": "complete sharp address between U and undo",
        "same_effects_operator_error": float(effect_error),
        "Lueders_Kraus_operator_error": float(map_error),
        "repeatable_disagreement_error": float(repeat_error),
        "pulse_duration_error_contract":
            "complete two-use cq half-trace error <= epsilon_prep + 218*delta "
            "if each of at most fourteen native waits has absolute error <= delta",
        "cases": cases,
    }


def compare(saved, actual, path="receipt"):
    """Exact source/discrete agreement and bounded float serialization variation."""
    if isinstance(actual, dict):
        require(isinstance(saved, dict) and saved.keys() == actual.keys(),
                f"{path}: field mismatch")
        for key, value in actual.items():
            compare(saved[key], value, f"{path}.{key}")
    elif isinstance(actual, list):
        require(isinstance(saved, list) and len(saved) == len(actual),
                f"{path}: shape mismatch")
        for index, (a, b) in enumerate(zip(saved, actual, strict=True)):
            compare(a, b, f"{path}[{index}]")
    elif isinstance(actual, float):
        require(type(saved) in (float, int) and np.isfinite(saved)
                and abs(saved - actual) <= 1e-11, f"{path}: floating replay mismatch")
    else:
        require(type(saved) is type(actual) and saved == actual,
                f"{path}: exact source or discrete replay mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="explicit first-record generation")
    args = parser.parse_args()
    result = report()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        with RECEIPT.open("x") as output:
            output.write(text)
    else:
        require(RECEIPT.is_file(), "receipt absent; use explicit --write")
        saved = json.loads(RECEIPT.read_text())
        compare(saved, result)
        print("signed sequential instrument: 16 complete banks replay; "
              "same effects, distinct poststates and repeat records")


if __name__ == "__main__":
    main()
