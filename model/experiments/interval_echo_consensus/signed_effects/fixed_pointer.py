"""A fixed neutral material pointer for an unknown vacancy/excess address."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import hashlib
import json
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell
from . import sequential

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RECEIPT = HERE / "fixed_pointer.json"
PORTS = (0, 1, 6, 10, 4, 5, 7, 11)
POINTERS = ((6, 10), (4, 5))  # fixed input loads (2,0)
HOST, REFERENCE = 7, 11


def physical_states():
    # Fixed total eight, with three zeros, two ones and three twos.
    result = []
    for zeros in combinations(range(8), 3):
        for ones in combinations([j for j in range(8) if j not in zeros], 2):
            row = [2] * 8
            for j in zeros:
                row[j] = 0
            for j in ones:
                row[j] = 1
            result.append(tuple(row))
    return tuple(sorted(result))


STATES = physical_states()
INDEX = {state: j for j, state in enumerate(STATES)}
LOADS = np.asarray(STATES, dtype=int)
DIM = len(STATES)


def pulse(kind, a, b, ratio):
    return (kind, a, b, Fraction(ratio))


def inverse(word):
    # All ternary blank-seam eigenvalues {0,8,9,32,33} are integers.
    return [(kind, a, b, (2 - ratio) % 2) for kind, a, b, ratio in reversed(word)]


def diagonal_word(a, b, total, complete_refocus=True):
    """Exactly exp(-i*pi*total*V) on every ternary endpoint word."""
    tau = Fraction(total) / 4
    echo = [pulse("contact", a, b, tau), pulse("reference", a, a, Fraction(1, 8)),
            pulse("contact", a, b, tau), pulse("reference", a, a, Fraction(1, 8))]
    if not complete_refocus:
        return echo + echo
    return (echo + [pulse("reference", a, a, Fraction(1, 32))]
            + echo + [pulse("reference", a, a, Fraction(63, 32))])


def neutral_basis(a, b):
    return [pulse("reference", b, b, Fraction(1, 64)),
            pulse("contact", a, b, Fraction(1, 2)),
            pulse("reference", b, b, Fraction(1, 64))]


def copy_word(pointer, corrected=True, complete_refocus=True):
    a, b = POINTERS[pointer]
    basis = neutral_basis(a, b)
    word = basis + diagonal_word(0, b, Fraction(2, 65), complete_refocus) + inverse(basis)
    if corrected:
        word += [pulse("reference", 0, 0, Fraction(1, 32))]
        word += diagonal_word(0, HOST, Fraction(133, 2210), complete_refocus)
    return word


def source_word():
    return [pulse("contact", 0, 1, Fraction(1, 2)),
            pulse("reference", 0, 0, Fraction(1, 20)),
            pulse("contact", 0, 1, Fraction(3, 2))]


def history(swapped=False):
    word = []
    for pointer in (0, 1):
        read = source_word()
        word += read + copy_word(pointer) + inverse(read)
        if pointer == 0 and swapped:
            word.append(pulse("contact", 0, 1, 1))
    return word


@lru_cache(maxsize=1)
def check_native_reduction():
    law = candidate.native_law
    first, second = law.local_term_parts()
    full = first + candidate.KAPPA * second
    blank = full[:, :law.N_SRC]
    expected = law.seam_term(candidate.KAPPA)
    sequential.require(np.linalg.norm(blank[:law.N_SRC] - expected, 2) < 4e-13
                       and np.linalg.norm(blank[law.N_SRC:], 2) < 4e-13,
                       "native ternary blank-memory reduction failed")
    return expected


@lru_cache(maxsize=None)
def swap_rows(a, b):
    sequential.require(b in cell().neighbours[a], "contact outside native incidence")
    x, y = PORTS.index(a), PORTS.index(b)
    swapped = LOADS.copy()
    swapped[:, [x, y]] = swapped[:, [y, x]]
    return np.asarray([INDEX[tuple(row)] for row in swapped])


def apply(word, columns):
    """Exact two-endpoint block formula on the entire 560-state sector.

    Opposite loads 0 and 2 exchange into all admitted placements. No signed
    hard-core or one-mark-per-dimer projection is made between pulses.
    """
    native = check_native_reduction()
    out = np.array(columns, complex, copy=True)
    for kind, a, b, ratio in word:
        time = float(ratio) * np.pi
        na = LOADS[:, PORTS.index(a)]
        if kind == "reference":
            out *= np.exp(-1j * candidate.KAPPA * time * na**2)[:, None]
            continue
        sequential.require(kind == "contact", "unknown pulse")
        nb = LOADS[:, PORTS.index(b)]
        rows = swap_rows(a, b)
        unequal = na != nb
        diagonal = np.diag(native)[3 * na + nb]
        changed = np.cos(time / 2) * out + 1j * np.sin(time / 2) * out[rows]
        changed *= np.exp(-1j * time * diagonal)[:, None]
        changed[~unequal] = out[~unequal]
        out = changed
    return out


def embedding(bits=(0, 0)):
    """Logical order: vacancy0, vacancy1, excess0, excess1.

    The balancing reference is input data. Every pointer and the host have
    the same sharp preparation in all four logical columns.
    """
    result = np.zeros((DIM, 4), complex)
    for sign_index, sign in enumerate((-1, 1)):
        for address in (0, 1):
            values = {port: 1 for port in PORTS}
            values[address] = 1 + sign
            values[REFERENCE] = 1 - sign
            for (a, b), bit in zip(POINTERS, bits, strict=True):
                values[a], values[b] = (2, 0) if bit == 0 else (0, 2)
            result[INDEX[tuple(values[p] for p in PORTS)], 2 * sign_index + address] = 1
    return result


def pointer_labels():
    words = [LOADS[:, PORTS.index(a)] + 3 * LOADS[:, PORTS.index(b)] for a, b in POINTERS]
    return words[0] + 9 * words[1] + 81 * LOADS[:, PORTS.index(HOST)]


LABELS = pointer_labels()


def branch_maps(swapped=False):
    columns = apply(history(swapped), embedding())
    result = []
    for label in range(243):
        branch = columns.copy()
        branch[LABELS != label] = 0
        result.append(branch)
    return result


def logical_maps(swapped=False):
    maps = branch_maps(swapped)
    result = {}
    for first in (0, 1):
        for second in (0, 1):
            bits = (1 - first, 1 - second)
            label = 81 + (2 + 4 * bits[0]) + 9 * (2 + 4 * bits[1])
            result[(first, second)] = embedding(bits).conj().T @ maps[label]
    return result


def target_maps(swapped=False):
    tables = [sequential.instrument(sign) for sign in (-1, 1)]
    result = {}
    x = np.array([[0, 1], [1, 0]]) if swapped else np.eye(2)
    for first in (0, 1):
        for second in (0, 1):
            target = np.zeros((4, 4), complex)
            for sign_index, table in enumerate(tables):
                target[2*sign_index:2*sign_index+2, 2*sign_index:2*sign_index+2] = (
                    table[second][:2, :2] @ x @ table[first][:2, :2])
            result[(first, second)] = target
    return result


def source_hashes():
    names = ["model/candidate_laws/contrast_memory_spherical_v1/law.py",
             "model/candidate_laws/interval_echo/candidate.py",
             "model/executable_model/carrier.py",
             "model/experiments/interval_echo_consensus/signed_effects/sequential.py",
             "model/experiments/interval_echo_consensus/signed_effects/fixed_pointer.py"]
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}


def coherence_control(corrected):
    # Native neutral-dimer preparation and analyzer on the source/reference
    # seam (0,11). Total charge is fixed; pointer states are the same.
    start = embedding()[:, [2]]  # source excess at 0, reference zero
    basis = neutral_basis(0, REFERENCE)
    state = apply(basis, start)
    state = apply(copy_word(0, corrected=corrected), state)
    state = apply(inverse(basis), state)
    ref = LOADS[:, PORTS.index(REFERENCE)]
    return [float(np.linalg.norm(state[ref == value])**2) for value in range(3)]


def report():
    rng = np.random.default_rng(925051)
    shots = 4096
    preparations = {}
    for sign_index, sign in enumerate((-1, 1)):
        for name, pair in {"port0": [1, 0], "port1": [0, 1],
                           "plus": np.array([1, 1])/np.sqrt(2),
                           "plus_i": np.array([1, 1j])/np.sqrt(2)}.items():
            vec = np.zeros(4, complex)
            vec[2*sign_index:2*sign_index+2] = pair
            preparations[f"{sign}:{name}"] = vec
    preparations["balanced_reference_plus"] = np.array([1, 0, 1, 0])/np.sqrt(2)
    preparations["balanced_reference_plus_i"] = np.array([1, 0, 1j, 0])/np.sqrt(2)
    phase = np.exp(-1j * 133 * np.pi / 130)  # one common phase for two copies
    max_map_error = 0.0
    banks = []
    for swapped in (False, True):
        maps = branch_maps(swapped)
        logical = logical_maps(swapped)
        targets = target_maps(swapped)
        sequential.require(np.linalg.norm(sum(k.conj().T @ k for k in maps)-np.eye(4), 2)
                           < 5e-12, "incomplete fixed-pointer instrument")
        for key, actual in logical.items():
            max_map_error = max(max_map_error, float(np.linalg.norm(actual-phase*targets[key], 2)))
        for name, state in preparations.items():
            law = np.array([np.linalg.norm(k @ state)**2 for k in maps])
            p = {key: float(np.linalg.norm(k @ state)**2) for key, k in targets.items()}
            p0 = p[(0, 0)] + p[(0, 1)]
            sequential.require(0.01 < p0 < 0.99, "sampling probe loses its first branch")
            q0 = p[(0, 1)] / p0 if swapped else 0
            q1 = p[(1, 0)] / (1-p0) if swapped else 0
            n0 = sequential.grid_binomial(rng, shots, p0)
            c0 = sequential.grid_binomial(rng, n0, q0)
            c1 = sequential.grid_binomial(rng, shots-n0, q1)
            counts = np.zeros(243, int)
            counts[[141, 105, 137, 101]] = [n0-c0, c0, c1, shots-n0-c1]
            sequential.require(abs(law.sum()-1) < 5e-12,
                               "complete fixed-pointer law lost norm")
            outside = [j for j in range(243) if j not in (141,105,137,101)]
            sequential.require(law[outside].sum() < 1e-23, "pointer or host escape")
            banks.append({"preparation": name, "arm": "native_swap" if swapped else "repeatable",
                          "all_243_pointer_host_probabilities": law.tolist(),
                          "all_243_pointer_host_counts": counts.tolist(), "shots": shots})
    sequential.require(max_map_error < 4e-12, "coherent unknown-sign branch map failed")
    wrong = apply(copy_word(0, complete_refocus=False), embedding())
    code = np.hstack([embedding(bits) for bits in ((0,0),(1,0))])
    residual = wrong - code @ (code.conj().T @ wrong)
    controls = {
        "incomplete_refocus_escape_by_logical_input": np.sum(abs(residual)**2, axis=0).tolist(),
        "corrected_balancing_analyzer_reference_load_probabilities": coherence_control(True),
        "uncorrected_balancing_analyzer_reference_load_probabilities": coherence_control(False),
    }
    sequential.require(min(controls["incomplete_refocus_escape_by_logical_input"][j] for j in (0,2)) > 1e-5,
                       "mixed-sign refocusing control lacks sensitivity")
    sequential.require(controls["corrected_balancing_analyzer_reference_load_probabilities"][0] > 1-3e-12
                       and controls["uncorrected_balancing_analyzer_reference_load_probabilities"][2] > 1-3e-12,
                       "native balanced-charge phase control failed")
    return {
        "scope": "two populated source rails, unknown sign; fixed neutral pointers and unit host",
        "ports": list(PORTS), "full_fixed_charge_multiset_dimension": DIM,
        "pointer_dimers": [list(pair) for pair in POINTERS],
        "fixed_pointer_initial_loads": [[2,0],[2,0]], "host": HOST, "host_load": 1,
        "balancing_reference": REFERENCE,
        "reference_scope": "untouched input charge reference during acquisition; native final analyzer only in coherence control",
        "controls_same_for_every_input": True,
        "pointer_word": "q1+9*q2+81*host; q=park_load+3*cross_load; all 243 words retained",
        "readout": "one final sharp occupation instrument on four pointer ports and host",
        "execution": "full 560-state native blank-memory pulse propagation; synthetic tree counts",
        "native_copy_word": [[kind,a,b,str(t)] for kind,a,b,t in copy_word(0)],
        "native_history_waits": len(history()),
        "common_two_copy_phase": [float(phase.real),float(phase.imag)],
        "branch_operator_error_including_cross_sign_coherence": max_map_error,
        "source_sha256": source_hashes(), "seed": 925051,
        "sampling_grid_bits": sequential.GRID_BITS, "sampling": sequential.SAMPLING,
        "custody": "producer-side software binding; external authentication absent",
        "controls": controls, "banks": banks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = report()
    if args.write:
        with RECEIPT.open("x") as output:
            output.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    else:
        sequential.compare(json.loads(RECEIPT.read_text()), result)
        print("fixed neutral pointer: unknown-sign and balanced-charge records replay")


if __name__ == "__main__":
    main()
