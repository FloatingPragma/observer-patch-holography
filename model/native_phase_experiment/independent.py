#!/usr/bin/env python3
"""Portable v2 exact independent OPH expectations, including all 36 pointers.

API: compute(repository_root, compact=True).
CLI: python portable_oracle_v2.py --root <OPH repository> --output <JSON file>
     Add --full for matrix entries and complete raw branch matrices.
     --output - writes the JSON to standard output without creating a file.

Only the standard library is used.  No repository module or source function
is imported or executed.  Three fixed source files supply literal carrier
and circuit data, and their byte hashes are a required snapshot contract.
The endpoint law is reconstructed directly as floor/ceiling averaging.

Compact output still includes every (law,seam,coin,raw) slot, its exact ideal
probability and readout contribution, and its accepted-selector counterpart.
It also includes both states' 36 ideal unnormalized pointer probabilities,
in lexicographic (first,second) order on {0,...,5} squared.  Outcomes with
either pointer outside {0,1} are individual zero-probability abort outcomes.
This is a finite software expectation, not a physical preparation claim.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
import sys


SOURCE_COMMIT = "3cf3c14810671a16e61b7915f69079a1773a12a2"
SOURCE_HASHES = {
    "model/observer/native.py":
        "dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py":
        "7e7cda89350be6330c7425df69bf9647549e6e1890a64eaaf28289be6c640b2f",
    "model/observer/source_discrimination/circuit.py":
        "327b537d01901c223e375cfff6dce882cb9a73f6a3b36e5019750e801377dd7e",
}
NATIVE_PATH, FACES_PATH, CIRCUIT_PATH = SOURCE_HASHES
PORTS = 12
SELECTOR_WORD_BITS = 11
SELECTOR_ATTEMPTS = 16
SELECTOR_ACCEPTED_WORDS = 1680


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _canonical_digest(value):
    return _sha(json.dumps(value, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=True).encode("ascii"))


def _literal(tree, name):
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(target, ast.Name) and target.id == name
                for target in node.targets):
            return ast.literal_eval(node.value)
    raise ValueError(f"Required source literal missing: {name}")


def _state(pair):
    loads = [0] * PORTS
    for port in pair:
        loads[port] += 1
    return tuple(loads)


def _endpoint(state, seam, coin):
    """Direct independent endpoint rule, without calling native.native."""
    a, b = seam
    total = state[a] + state[b]
    low, high = total // 2, (total + 1) // 2
    target = list(state)
    target[a], target[b] = (low, high) if coin == 0 else (high, low)
    return tuple(target), state[a] - target[a]


def _raw(state, seam, current, law):
    code = {0: 0, 1: 1, -1: 2}[current]
    return code if law == "C" else code + 4 * ((state[seam[0]] + state[seam[1]]) % 2)


def _clean(matrix):
    return {position: Q(value) for position, value in matrix.items() if value}


def _scale(matrix, factor):
    return _clean({position: value * factor for position, value in matrix.items()})


def _subtract(left, right):
    return _clean({position: left.get(position, Q()) - right.get(position, Q())
                   for position in left.keys() | right.keys()})


def _trace(matrix):
    return sum((value for (i, j), value in matrix.items() if i == j), Q())


def _hs_squared(matrix):
    return sum((value * value for value in matrix.values()), Q())


def _readout(matrix, observable):
    if len(observable) <= len(matrix):
        return sum((value * matrix.get((j, i), Q())
                    for (i, j), value in observable.items()), Q())
    return sum((value * observable.get((j, i), Q())
                for (i, j), value in matrix.items()), Q())


def _matrix(matrix, compact):
    entries = [[i, j, str(value)] for (i, j), value in sorted(matrix.items())]
    result = {"dimension": 66, "nnz": len(entries), "trace": str(_trace(matrix)),
              "hs_squared": str(_hs_squared(matrix)),
              "entries_sha256": _canonical_digest(entries)}
    if not compact:
        result.update(entries=entries, omitted_entries="exactly zero")
    return result


def _rank(matrix, dimension=66):
    rows = [{j: value for (i, j), value in matrix.items() if i == row}
            for row in range(dimension)]
    rank = 0
    for column in range(dimension):
        pivot = next((i for i in range(rank, dimension) if rows[i].get(column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        denominator = rows[rank][column]
        rows[rank] = {j: value / denominator for j, value in rows[rank].items()}
        for row in range(rank + 1, dimension):
            multiplier = rows[row].get(column, Q())
            if not multiplier:
                continue
            for j, value in rows[rank].items():
                updated = rows[row].get(j, Q()) - multiplier * value
                if updated:
                    rows[row][j] = updated
                else:
                    rows[row].pop(j, None)
        rank += 1
    return rank


def _branch(matrix, mapping):
    output = defaultdict(Q)
    for (i, j), value in matrix.items():
        if i in mapping and j in mapping:
            output[mapping[i], mapping[j]] += value / 60
    return _clean(output)


def _add(target, matrix):
    for position, value in matrix.items():
        target[position] += value


def _multiply(left, right):
    rows = defaultdict(dict)
    for (i, j), value in right.items():
        rows[i][j] = value
    result = defaultdict(Q)
    for (i, k), left_value in left.items():
        for j, right_value in rows.get(k, {}).items():
            result[i, j] += left_value * right_value
    return _clean(result)


def _pointer_effects(states, index):
    """Pull back occupation projectors through the defined inverse beam.

    On each single-occupation spectator pair, inverse beam sends the
    symmetric vector to occupation (1,0), and the antisymmetric vector to
    (0,1).  Their pullbacks are |plus><plus| and |minus><minus|, each with
    rational entries +/-1/2.  The (0,0) and (1,1) sectors are unchanged.
    The remaining 32 effects are zero on this ideal binary carrier.
    """
    effects = [{} for _ in range(36)]
    for i, state in enumerate(states):
        first, second = state[:2]
        if first == second:
            effects[6 * first + second][i, i] = Q(1)
    for spectator in range(2, PORTS):
        first = index[_state((0, spectator))]
        second = index[_state((1, spectator))]
        for i in (first, second):
            for j in (first, second):
                effects[6][i, j] = Q(1, 2)  # (1,0): symmetric spectral projector.
                effects[1][i, j] = Q(1 if i == j else -1, 2)  # (0,1).
    return effects


def validate_pointer_probabilities(probabilities, expected_probability, expected_score):
    """Validate all 36 cells, including aborts that sum/score alone miss.

    Accepts exact Fraction values or their rational string representations.
    Returns True, or raises ValueError identifying the violated invariant.
    """
    if len(probabilities) != 36:
        raise ValueError("pointer_distribution_length")
    probabilities = [Q(value) for value in probabilities]
    if any(value < 0 for value in probabilities):
        raise ValueError("negative_pointer_probability")
    if sum(probabilities, Q()) != Q(expected_probability):
        raise ValueError("pointer_probability_sum")
    if probabilities[6] - probabilities[1] != Q(expected_score):
        raise ValueError("pointer_score")
    if any(probabilities[6 * first + second] != 0
           for first in range(6) for second in range(6)
           if first not in (0, 1) or second not in (0, 1)):
        raise ValueError("nonzero_pointer_abort_probability")
    return True


def compute(root, *, compact=True):
    """Recalculate exact expectations from a repository root.

    No paths from the author's machine are used or emitted.  compact=True
    omits matrix entries but retains all 720 raw outcome rows and their
    probabilities/readout contributions, including exact zero outputs.
    Every input and raw row includes each state's 36 pointer probabilities.
    They are ideal unnormalized probabilities; Main may multiply them by q.
    Raises ValueError on a snapshot mismatch and AssertionError on a failed
    mathematical consistency check.  Writes no files and mutates no source.
    """
    root = Path(root)
    check_count = 0

    def check(name, condition):
        nonlocal check_count
        check_count += 1
        if not condition:
            raise AssertionError(name)

    source_bytes = {}
    for relative_path, expected_hash in SOURCE_HASHES.items():
        data = (root / relative_path).read_bytes()
        if _sha(data) != expected_hash:
            raise ValueError(f"Fixed source snapshot hash mismatch: {relative_path}")
        source_bytes[relative_path] = data
        check("fixed_source_hash", True)
    trees = {name: ast.parse(data.decode("utf-8")) for name, data in source_bytes.items()}
    faces = _literal(trees[FACES_PATH], "ORIENTED_BASE_FACES")
    seams = tuple(sorted({tuple(sorted((face[k], face[(k + 1) % 3])))
                          for face in faces for k in range(3)}))
    degrees = Counter(port for seam in seams for port in seam)
    check("20_faces_30_contacts_five_regular", len(faces) == 20 and len(seams) == 30
          and degrees == Counter({port: 5 for port in range(PORTS)}))
    states = tuple(sorted(_state(pair) for pair in combinations(range(PORTS), 2)))
    whole_states = tuple(sorted(_state(pair) for pair in combinations_with_replacement(range(PORTS), 2)))
    index = {state: i for i, state in enumerate(states)}
    check("66_binary_78_whole_states", len(states) == len(index) == 66 and len(whole_states) == 78)

    # Direct comparison with the current equation read in the fixed native.py.
    # This duplicates the equation, not a source function or imported module.
    endpoint_checks = 0
    for state in whole_states:
        for seam in seams:
            a, b = seam
            for coin in (0, 1):
                target, current = _endpoint(state, seam, coin)
                contrast = state[a] - state[b]
                native_current_equation = (contrast - (2 * coin - 1) * (abs(contrast) % 2)) // 2
                check("endpoint_current_equation", current == native_current_equation
                      and target[a] == state[a] - current
                      and target[b] == state[b] + current
                      and all(target[k] == state[k] for k in range(PORTS) if k not in seam))
                endpoint_checks += 1

    circuit_tree = trees[CIRCUIT_PATH]
    local_states = _literal(circuit_tree, "SEAM_STATES")
    swaps = _literal(circuit_tree, "REPAIR_SWAPS")
    raw_width = _literal(circuit_tree, "RECORD_BITS")
    check("source_raw_width_and_current_encoding", raw_width == {"C": 2, "E": 3}
          and _literal(circuit_tree, "CURRENT_CODES") == {0: 0, 1: 1, -1: 2})
    record_decoding = {
        ("C", 0): {0: 0, 1: 1, 2: -1},
        ("C", 1): {0: 0, 1: 1, 2: -1},
        ("E", 0): {0: 0, 1: 2, 2: -2, 4: -1, 5: 1},
        ("E", 1): {0: 0, 1: 2, 2: -2, 4: 1, 6: -1},
    }
    for coin in (0, 1):
        for local_index, local_state in enumerate(local_states):
            sample = local_state + (2 - sum(local_state),) + (0,) * 9
            target, current = _endpoint(sample, (0, 1), coin)
            code = local_index
            for left, right in swaps["C", coin]:
                code = right if code == left else left if code == right else code
            for law in ("C", "E"):
                final_code = code if law == "C" else code ^ (((code >> coin) & 1) << 5)
                raw = final_code // 8
                read_record = current if law == "C" else sample[0] - sample[1]
                check("local_dilation_and_raw_record", local_states[final_code % 8] == target[:2]
                      and raw == _raw(sample, (0, 1), current, law)
                      and record_decoding[law, coin].get(raw) == read_record)

    transition = defaultdict(Q)
    for j, state in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                target, _ = _endpoint(state, seam, coin)
                check("binary_shell_closed", target in index)
                transition[index[target], j] += Q(1, 60)
    transition = _clean(transition)
    check("P_symmetric", all(transition.get((j, i), Q()) == value
                             for (i, j), value in transition.items()))
    check("P_stochastic", all(sum((value for (i, j), value in transition.items() if i == k), Q()) == 1
                              for k in range(66)))
    trace_p = _trace(transition)
    rho = _scale(transition, 1 / trace_p)
    signs = [1 - 2 * state[0] for state in states]
    sigma = {(i, j): value * signs[i] * signs[j] for (i, j), value in rho.items()}
    delta = _subtract(sigma, rho)
    margin = min(transition[i, i] - sum((abs(value) for (r, c), value in transition.items()
                                      if r == i and c != i), Q()) for i in range(66))
    check("rho_strictly_positive", margin > 0)
    check("states_unit_trace_and_distinct", _trace(rho) == _trace(sigma) == 1 and bool(delta))
    check("same_configuration_distribution", all(rho[i, i] == sigma[i, i] for i in range(66)))
    observable = {}
    for k in range(2, PORTS):
        i, j = index[_state((0, k))], index[_state((1, k))]
        observable[i, j] = observable[j, i] = Q(1)
    check("O_ten_exchange_pairs", len(observable) == 20)

    pointer_effects = _pointer_effects(states, index)
    effect_sum = defaultdict(Q)
    for effect in pointer_effects:
        check("pointer_effect_self_adjoint", all(effect.get((j, i), Q()) == value
                                                 for (i, j), value in effect.items()))
        check("pointer_effect_is_projection", _multiply(effect, effect) == effect)
        _add(effect_sum, effect)
    check("36_pointer_effects_sum_to_identity", _clean(effect_sum) == {(i, i): Q(1) for i in range(66)})
    check("pointer_score_effect_equals_O", _subtract(pointer_effects[6], pointer_effects[1]) == observable)

    def pointer_distribution(matrix):
        values = [_readout(matrix, effect) for effect in pointer_effects]
        check("pointer_all_cells_nonnegative", all(value >= 0 for value in values))
        check("pointer_sum_equals_matrix_trace", sum(values, Q()) == _trace(matrix))
        check("pointer_score_equals_exchange_expectation", values[6] - values[1] == _readout(matrix, observable))
        check("each_pointer_abort_is_exact_zero", all(values[6 * first + second] == 0
              for first in range(6) for second in range(6)
              if first not in (0, 1) or second not in (0, 1)))
        check("pointer_public_validator", validate_pointer_probabilities(values, _trace(matrix),
                                                                          _readout(matrix, observable)))
        return [str(value) for value in values]

    rejected_words = (1 << SELECTOR_WORD_BITS) - SELECTOR_ACCEPTED_WORDS
    timeout = Q(rejected_words, 1 << SELECTOR_WORD_BITS) ** SELECTOR_ATTEMPTS
    accepted = 1 - timeout
    check("finite_selector_equal_seam_fibers", SELECTOR_ACCEPTED_WORDS % len(seams) == 0)
    check("finite_selector_full_probability", accepted + timeout == 1 and 0 < timeout < 1)
    check("finite_selector_geometric_sum", accepted == sum(
        (Q(rejected_words, 1 << SELECTOR_WORD_BITS) ** attempt
         * Q(SELECTOR_ACCEPTED_WORDS, 1 << SELECTOR_WORD_BITS)
         for attempt in range(SELECTOR_ATTEMPTS)), Q()))
    laws = {}
    output_deltas = {}
    for law in ("C", "E"):
        raw_rows = []
        out_rho, out_sigma = defaultdict(Q), defaultdict(Q)
        completeness = [Q()] * 66
        recorded_hs = Q()
        for seam in seams:
            for coin in (0, 1):
                by_raw = {raw: {} for raw in range(1 << raw_width[law])}
                for i, state in enumerate(states):
                    target, current = _endpoint(state, seam, coin)
                    by_raw[_raw(state, seam, current, law)][i] = index[target]
                for raw, mapping in by_raw.items():
                    check("raw_partial_injection", len(mapping) == len(set(mapping.values())))
                    for i in mapping:
                        completeness[i] += Q(1, 60)
                    a, b = _branch(rho, mapping), _branch(sigma, mapping)
                    branch_delta = _subtract(b, a)
                    pa, pb = _trace(a), _trace(b)
                    ea, eb = _readout(a, observable), _readout(b, observable)
                    check("raw_equal_probability", pa == pb)
                    check("zero_probability_zero_matrices", bool(pa) or (not a and not b))
                    check("raw_exterior_difference_zero", all(states[i][0] != states[j][0]
                                                               for i, j in branch_delta))
                    recorded_hs += _hs_squared(branch_delta)
                    row = {
                        "law": law, "seam": list(seam), "coin": coin, "raw": raw,
                        "admitted_columns": len(mapping),
                        "rho_probability": str(pa), "sigma_probability": str(pb),
                        "rho_readout_weight": str(ea), "sigma_readout_weight": str(eb),
                        "accepted_rho_probability": str(accepted * pa),
                        "accepted_sigma_probability": str(accepted * pb),
                        "accepted_rho_readout_weight": str(accepted * ea),
                        "accepted_sigma_readout_weight": str(accepted * eb),
                        "rho_conditional_readout": str(ea / pa) if pa else None,
                        "sigma_conditional_readout": str(eb / pb) if pb else None,
                        "delta_hs_squared": str(_hs_squared(branch_delta)),
                        "rho_entries_sha256": _matrix(a, True)["entries_sha256"],
                        "sigma_entries_sha256": _matrix(b, True)["entries_sha256"],
                        "rho_pointer_probabilities": pointer_distribution(a),
                        "sigma_pointer_probabilities": pointer_distribution(b),
                    }
                    if not compact:
                        row.update(column_map=[[i, j] for i, j in sorted(mapping.items())],
                                   rho=_matrix(a, False), sigma=_matrix(b, False),
                                   delta_sigma_minus_rho=_matrix(branch_delta, False))
                    raw_rows.append(row)
                    _add(out_rho, a)
                    _add(out_sigma, b)
        out_rho, out_sigma = _clean(out_rho), _clean(out_sigma)
        out_delta = _subtract(out_sigma, out_rho)
        output_deltas[law] = out_delta
        factor = out_delta.get(next(iter(delta)), Q()) / delta[next(iter(delta))]
        er, es = _readout(out_rho, observable), _readout(out_sigma, observable)
        check("complete_KdaggerK_identity", all(value == 1 for value in completeness))
        check("channel_trace_preserving", _trace(out_rho) == _trace(out_sigma) == 1)
        check("delta_is_exact_scalar_multiple", out_delta == _scale(delta, factor))
        check("all_raw_probability_retained", sum((Q(row["rho_probability"]) for row in raw_rows), Q()) == 1)
        check("all_raw_readout_retained", sum((Q(row["rho_readout_weight"]) for row in raw_rows), Q()) == er
              and sum((Q(row["sigma_readout_weight"]) for row in raw_rows), Q()) == es)
        check("selector_all_raw_accepted_mass", sum((Q(row["accepted_rho_probability"])
                                                     for row in raw_rows), Q()) == accepted)
        check("selector_all_raw_accepted_readout", sum((Q(row["accepted_rho_readout_weight"])
                                                        for row in raw_rows), Q()) == accepted * er)
        laws[law] = {
            "raw_width": raw_width[law], "raw_branch_count": len(raw_rows),
            "zero_probability_branch_count": sum(row["rho_probability"] == "0" for row in raw_rows),
            "rho_readout": str(er), "sigma_readout": str(es),
            "rho_pointer_probabilities": pointer_distribution(out_rho),
            "sigma_pointer_probabilities": pointer_distribution(out_sigma),
            "readout_gap_sigma_minus_rho": str(es - er),
            "delta_factor": str(factor), "delta_rank": _rank(out_delta),
            "rho": _matrix(out_rho, compact), "sigma": _matrix(out_sigma, compact),
            "delta_sigma_minus_rho": _matrix(out_delta, compact),
            "retained_classical_records_delta_hs_squared": str(recorded_hs),
            "finite_selector": {"accepted_probability": str(accepted),
                                "timeout_probability": str(timeout),
                                "accepted_rho_readout_weight": str(accepted * er),
                                "accepted_sigma_readout_weight": str(accepted * es),
                                "accepted_readout_gap_sigma_minus_rho": str(accepted * (es - er)),
                                "accepted_rho_trace": str(accepted), "accepted_sigma_trace": str(accepted)},
            "raw_branches": raw_rows,
        }

    reference = {
        "trace_P": str(trace_p), "P": _matrix(transition, compact),
        "P_diagonal_histogram": dict(sorted(Counter(str(transition[i, i]) for i in range(66)).items())),
        "P_strict_diagonal_dominance_margin": str(margin),
        "rho": _matrix(rho, compact), "sigma": _matrix(sigma, compact),
        "delta_sigma_minus_rho": _matrix(delta, compact), "delta_rank": _rank(delta),
        "rho_readout": str(_readout(rho, observable)),
        "sigma_readout": str(_readout(sigma, observable)),
        "rho_pointer_probabilities": pointer_distribution(rho),
        "sigma_pointer_probabilities": pointer_distribution(sigma),
        "readout_gap_sigma_minus_rho": str(_readout(delta, observable)),
        "observable": _matrix(observable, compact),
    }
    carrier = {"state_count": 66, "whole_native_state_count": 78,
               "seams": [list(seam) for seam in seams], "degrees": [degrees[i] for i in range(PORTS)],
               "source_transition_path_count": 3960,
               "endpoint_current_equation_cross_checks": endpoint_checks}
    if not compact:
        carrier.update(states=[list(state) for state in states], faces=[list(face) for face in faces])
    return {
        "schema": "oph.one_port.portable_independent_oracle.v2", "compact": bool(compact),
        "source_commit": SOURCE_COMMIT, "source_sha256": dict(SOURCE_HASHES),
        "oracle_sha256": _sha(Path(__file__).read_bytes()),
        "conventions": {
            "state_order": "lexicographically sorted twelve-port binary total-two occupation tuples",
            "matrix_order": "row output, column input; P is symmetric",
            "D": "D[x,x]=(-1)^x[0]", "difference": "sigma-rho, sigma=D rho D",
            "observable": "O[(0,k),(1,k)]=O[(1,k),(0,k)]=1 for k=2..11; all other entries zero",
            "instrument": "uniform seam and fair coin; full raw records retained; source reduction sums all raw blocks",
            "selector": "one finite selector; accepted weights are q times ideal unnormalized contributions; timeout is separate",
            "timeout": "flagged selector failure; no repair or readout contribution is assigned to timeout here",
            "scope": "finite software expectations; no changed preparation/repair/readout circuit and no physical theorem",
        },
        "pointer_readout": {
            "outcome_order": [[first, second] for first in range(6) for second in range(6)],
            "normal_outcomes_in_score_order": [[0, 0], [1, 0], [0, 1], [1, 1]],
            "normal_scores": [0, 1, -1, 0],
            "abort_rule": "either first or second outside {0,1}; all 32 individual probabilities exactly zero",
            "construction": "spectral projectors of the exchange observable pulled back from occupation after inverse beam",
            "pointer_probability_normalization": "ideal unnormalized branch probabilities; raw row sum equals rho_probability or sigma_probability; no q applied",
        },
        "carrier": carrier, "preparation": reference, "one_repair": laws,
        "finite_selector": {
            "word_bits": SELECTOR_WORD_BITS, "attempts": SELECTOR_ATTEMPTS,
            "accepted_word_range_inclusive": [0, SELECTOR_ACCEPTED_WORDS - 1],
            "accepted_word_count": SELECTOR_ACCEPTED_WORDS, "rejected_word_count": rejected_words,
            "accepted_words_per_seam": SELECTOR_ACCEPTED_WORDS // len(seams),
            "one_attempt_rejection_probability": str(Q(rejected_words, 1 << SELECTOR_WORD_BITS)),
            "accepted_probability_formula": "q=1-(368/2048)^16",
            "accepted_probability": str(accepted), "timeout_probability": str(timeout),
            "accepted_probability_per_seam": str(accepted / 30),
            "accepted_probability_per_seam_and_coin": str(accepted / 60),
        },
        "checks": {"all_pass": True, "count": check_count, "failed": []},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="Root of the fixed OPH source repository")
    parser.add_argument("--output", required=True, help="JSON output path, or - for standard output")
    parser.add_argument("--full", action="store_true", help="Include matrix entries and all raw branch matrices")
    args = parser.parse_args(argv)
    result = compute(args.root, compact=not args.full)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output == "-":
        sys.stdout.write(encoded)
    else:
        Path(args.output).write_text(encoded, encoding="utf-8", newline="\n")
        print(json.dumps({"all_pass": result["checks"]["all_pass"],
                          "checks": result["checks"]["count"],
                          "compact": result["compact"], "output_sha256": _sha(encoded.encode("utf-8"))},
                         sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
