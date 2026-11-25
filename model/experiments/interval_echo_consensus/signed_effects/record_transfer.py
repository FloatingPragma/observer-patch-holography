"""Copy a two-rail address into native material dimers before a final-only read."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell
from . import sequential

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECEIPT = HERE / "record_transfer.json"
PORTS = (0, 1, 6, 10, 4, 5)
POINTERS = ((6, 10), (4, 5))  # park, marked-record rail
STATES = tuple(combinations(PORTS, 3))
INDEX = {frozenset(x): j for j, x in enumerate(STATES)}
DIM = len(STATES)


def pulse(kind, a, b, ratio):
    return (kind, a, b, Fraction(ratio))


def inverse(word):
    # Contact spectra are {0,8,9}; reference spectra are {0,8,32}.
    return [(kind, a, b, (2 - ratio) % 2) for kind, a, b, ratio in reversed(word)]


def basis_rotation(a, b):
    # S* Rx(pi/2) S* is Hadamard up to one common scalar on the dimer.
    return [pulse("reference", b, b, Fraction(3, 16)),
            pulse("contact", a, b, Fraction(1, 2)),
            pulse("reference", b, b, Fraction(3, 16))]


def controlled_phase(a, b, refocus=True, coupled=True):
    # Parity refocusing removes exchange and keeps the unequal-load phase.
    word = []
    if coupled:
        word.append(pulse("contact", a, b, Fraction(1, 34)))
    if refocus:
        word.append(pulse("reference", a, a, Fraction(1, 8)))
    if coupled:
        word.append(pulse("contact", a, b, Fraction(1, 34)))
    if refocus:
        word.append(pulse("reference", a, a, Fraction(1, 8)))
    word += [pulse("reference", a, a, Fraction(1, 16)),
             pulse("reference", b, b, Fraction(1, 16))]
    return word


def copy_word(pointer, refocus=True, coupled=True):
    a, b = POINTERS[pointer]
    basis = basis_rotation(a, b)
    return basis + controlled_phase(0, b, refocus, coupled) + inverse(basis)


def source_word():
    return [pulse("contact", 0, 1, Fraction(1, 2)),
            pulse("reference", 0, 0, Fraction(1, 20)),
            pulse("contact", 0, 1, Fraction(3, 2))]


@lru_cache(maxsize=None)
def generator(sign, kind, a, b):
    sequential.require(sign in (-1, 1) and a in PORTS and b in PORTS,
                       "unsupported signed register or port")
    marks = [set(x) for x in STATES]
    loads = np.array([[1 + sign * int(p in x) for p in PORTS] for x in marks])
    if kind == "reference":
        return np.diag(8 * loads[:, PORTS.index(a)] ** 2)
    sequential.require(kind == "contact" and b in cell().neighbours[a],
                       "non-native contact")
    law = candidate.native_law
    h0, h1 = law.local_term_parts()
    full = h0 + candidate.KAPPA * h1
    result = np.zeros((DIM, DIM))
    for col, marked in enumerate(marks):
        pair = tuple(1 + sign * int(p in marked) for p in (a, b))
        local = law.SRC.index(pair)
        sequential.require(np.max(abs(full[law.N_SRC:, local])) < 3e-13,
                           "blank memory not reducing")
        for row, out in enumerate(law.SRC):
            amplitude = full[row, local]
            if abs(amplitude) < 1e-14:
                continue
            sequential.require(all(x in (1, 1 + sign) for x in out),
                               "native contact exits signed hard-core sector")
            target = marked - {a, b}
            target.update(p for p, value in zip((a, b), out, strict=True) if value != 1)
            sequential.require(len(target) == 3, "mark count changed")
            result[INDEX[frozenset(target)], col] += amplitude
    sequential.require(np.linalg.norm(result - result.T) < 3e-13,
                       "contact is not Hermitian")
    return result


@lru_cache(maxsize=None)
def unitary(sign, word):
    result = np.eye(DIM, dtype=complex)
    for kind, a, b, ratio in word:
        result = expm(-1j * float(ratio) * np.pi * generator(sign, kind, a, b)) @ result
    return result


def execute(sign, word):
    return unitary(sign, tuple(word))


def embedding(pointer_bits=(0, 0)):
    v = np.zeros((DIM, 2), complex)
    for source in (0, 1):
        marked = frozenset((source, *(POINTERS[j][bit] for j, bit in enumerate(pointer_bits))))
        v[INDEX[marked], source] = 1
    return v


def pointer_word(marked):
    words = [int(a in marked) + 2 * int(b in marked) for a, b in POINTERS]
    return 4 * words[0] + words[1]


def history(sign, swapped=False, refocus=True, coupled=True):
    read = source_word()
    out = []
    for j in (0, 1):
        out += read + copy_word(j, refocus=refocus, coupled=coupled) + inverse(read)
        if j == 0 and swapped:
            out.append(pulse("contact", 0, 1, 1))
    return out


def branch_maps(sign, swapped=False, refocus=True, coupled=True):
    evolved = execute(sign, history(sign, swapped, refocus, coupled)) @ embedding()
    maps = []
    for label in range(16):
        branch = evolved.copy()
        branch[[pointer_word(set(x)) != label for x in STATES], :] = 0
        maps.append(branch)
    return maps


def qubit_maps(sign, swapped=False):
    maps = branch_maps(sign, swapped)
    result = {}
    for first in (0, 1):
        for second in (0, 1):
            # A mark at the cross rail records source port 0.
            bits = (1 - first, 1 - second)
            label = 4 * (1 + bits[0]) + (1 + bits[1])
            result[(first, second)] = embedding(bits).conj().T @ maps[label]
    return result


def source_hashes():
    names = [
        "model/candidate_laws/contrast_memory_spherical_v1/law.py",
        "model/candidate_laws/interval_echo/candidate.py",
        "model/executable_model/carrier.py",
        "model/experiments/interval_echo_consensus/signed_effects/sequential.py",
        "model/experiments/interval_echo_consensus/signed_effects/record_transfer.py",
    ]
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}


def serialized_word(word):
    return [[kind, a, b, str(ratio)] for kind, a, b, ratio in word]


def report():
    primitive = candidate.native_law.numeric_checks()
    sequential.require(primitive["U_unitary"] and primitive["reduction_residual"] < 3e-12,
                       "native reduction check failed")
    cases = []
    max_map_error = 0.0
    rng = np.random.default_rng(925041)
    shots = 8192
    probes = {"port0": [1, 0], "port1": [0, 1],
              "plus": np.array([1, 1]) / np.sqrt(2),
              "plus_i": np.array([1, 1j]) / np.sqrt(2)}
    for sign in (-1, 1):
        for swapped in (False, True):
            maps = branch_maps(sign, swapped)
            # One final pointer-occupation read is the only projection.
            sequential.require(np.linalg.norm(sum(k.conj().T @ k for k in maps)
                                               - np.eye(2), 2) < 5e-12,
                               "incomplete final record law")
            direct = qubit_maps(sign, swapped)
            target = sequential.instrument(sign)
            x = np.array([[0, 1], [1, 0]])
            for (first, second), k in direct.items():
                expected = target[second][:2, :2] @ (x if swapped else np.eye(2)) @ target[first][:2, :2]
                # Two -CNOT scalars cancel in the whole history.
                max_map_error = max(max_map_error, float(np.linalg.norm(k - expected, 2)))
            for name, raw in probes.items():
                state = np.asarray(raw, complex)
                p = np.array([np.linalg.norm(k @ state) ** 2 for k in maps])
                allowed = [5, 6, 9, 10]
                sequential.require(p.sum() > 1 - 5e-12 and p.sum() < 1 + 5e-12,
                                   "norm drift")
                sequential.require(p[[j for j in range(16) if j not in allowed]].sum() < 1e-24,
                                   "invalid material pointer outcome")
                # Exact rank-one tree sampling; no rescaling of propagated p.
                first0 = float(np.vdot(state, target[0][:2, :2] @ state).real)
                mismatch = (3 + sign * np.sqrt(5)) / 8 if swapped else 0
                n0 = int(rng.binomial(shots, first0))
                c0, c1 = int(rng.binomial(n0, mismatch)), int(rng.binomial(shots - n0, mismatch))
                counts = np.zeros(16, int)
                counts[[10, 9, 6, 5]] = [n0 - c0, c0, c1, shots - n0 - c1]
                cases.append({"sign": sign, "preparation": name,
                              "arm": "native_swap" if swapped else "nondemolition",
                              "all_16_pointer_word_probabilities": p.tolist(),
                              "all_16_pointer_word_counts": counts.tolist(), "shots": shots})
        # A sharp source probes leakage when the exchange-cancelling pulses are omitted.
        reference = embedding()[:, 0]
        failed = execute(sign, copy_word(0, refocus=False)) @ reference
        good = execute(sign, copy_word(0)) @ reference
        controls = {
            "sign": sign,
            "no_refocus_population_outside_one_mark_per_dimer":
                float(sum(abs(failed[i]) ** 2 for i, x in enumerate(STATES)
                          if not all(len(set(x) & set(pair)) == 1
                                     for pair in ((0, 1), *POINTERS)))),
            "nondemolition_population_outside_one_mark_per_dimer":
                float(sum(abs(good[i]) ** 2 for i, x in enumerate(STATES)
                          if not all(len(set(x) & set(pair)) == 1
                                     for pair in ((0, 1), *POINTERS)))),
        }
        omitted = execute(sign, copy_word(0, coupled=False)) @ embedding()
        actual_copy = execute(sign, copy_word(0)) @ embedding()

        def first_pointer_law(columns):
            return np.array([[sum(abs(columns[i, col]) ** 2 for i, marked in enumerate(STATES)
                                  if int(POINTERS[0][0] in marked)
                                  + 2 * int(POINTERS[0][1] in marked) == label)
                              for label in range(4)] for col in range(2)])

        for name, columns in (("native", actual_copy), ("cross_contact_omitted", omitted)):
            law = first_pointer_law(columns)
            controls[name + "_source_conditioned_pointer_laws"] = law.tolist()
            controls[name + "_source_record_total_variation"] = float(abs(law[0] - law[1]).sum() / 2)
        sequential.require(controls["native_source_record_total_variation"] > 1 - 3e-12
                           and controls["cross_contact_omitted_source_record_total_variation"] < 3e-12,
                           "source-coupling discrimination failed")
        sequential.require(controls["no_refocus_population_outside_one_mark_per_dimer"] > 0.001,
                           "refocusing control lacks sensitivity")
        cases.append({"controls": controls})
    sequential.require(max_map_error < 3e-12, "deferred-read branch map differs")
    return {
        "scope": "two populated source rails only; three prepared matched-sign defects on six ports",
        "source_rails": [0, 1], "pointer_dimers": [list(x) for x in POINTERS],
        "carrier": "native twelve-port cell; every unlisted contact parked",
        "complete_three_mark_states": [list(x) for x in STATES],
        "source_sha256": source_hashes(),
        "execution": "native blank-memory reduction and full twenty-state unitary; synthetic counts",
        "projection": "one final exhaustive four-pointer-port occupation read; no intermediate read",
        "pointer_word": "4*(park1+2*cross1)+(park2+2*cross2); all sixteen words retained",
        "native_copy_word_first_pointer": serialized_word(copy_word(0)),
        "native_two_record_word": serialized_word(history(-1)),
        "same_word_for_both_signs": True,
        "deferred_vs_projective_branch_operator_error": max_map_error,
        "pulse_error_contract": "complete two-record cq half-trace error <= "
            "min(1, epsilon_prep+793*delta) for at most 37 native waits each "
            "with absolute duration error <= delta",
        "custody": "producer-side software hashes; no external authentication",
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as output:
            output.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        sequential.compare(json.loads(RECEIPT.read_text()), result)
        print("native material record transfer: final-only two-pointer read replays")


if __name__ == "__main__":
    main()
