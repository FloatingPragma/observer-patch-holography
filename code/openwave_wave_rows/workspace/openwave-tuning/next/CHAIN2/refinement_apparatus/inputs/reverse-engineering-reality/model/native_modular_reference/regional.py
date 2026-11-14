"""Native regional Palm algebras and one compatible modular action.

On the fixed two-occupation shell, a region R has the number-preserving
algebra A_R = direct_sum_q M(local q-configurations) tensor I(outside).
The Palm density is diag(degree/600). Its modular action need not preserve
A_R: a local coherence can see different degree ratios at different outside
configurations. Group outside configurations precisely when their positive
degree columns over local configurations are proportional. The resulting
direct-sum matrix algebra C_R is the smallest modular-invariant star algebra
containing A_R. This is a derived boundary-context refinement, not a new
reference or clock. Single-port and classical occupation algebras have
trivial modular action.

Each block has sigma = p mu tensor nu, with faithful diagonal mu and nu.
The explicit state-preserving conditional expectation is id tensor nu on
diagonal blocks, killing cross blocks. Its Kraus form proves complete
positivity; the formula proves unitality, the bimodule law, idempotence and
sigma preservation. All expectations and restrictions here are algebraic
maps. An implementation of these operations by regional observers, including
coherent access and availability of the boundary context, is not inferred.

Spectral projections of the modular generator isolate equal degree ratios
in a local matrix unit. Products of these components and their adjoints
generate the joint outside-context partition, proving minimality. Thus
R subset S implies C_R subset C_S. Intersections are invariant and the
single ambient modular group agrees on every common observable. Arbitrary
disjoint C_R and C_S need not commute, and no spacetime causal net follows.

The full source word already implements this same action on all C_R.
Finite Palm-preparation operator sandwiches survive every regional partial
trace. They bound the modular action of each actual regional marginal by
the ideal action, even when the actual global modular flow does not exactly
preserve C_R. This comparison keeps reference, preparation, coherent-service
and dimensionless-parameter contracts from full_flow.py.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from tools import strict_json
from ..primitive_observer import transfer
from . import preparation, full_flow

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("regional_results.json")
CAP = (0, 1, 5, 7)


def region(ports):
    if (type(ports) not in (tuple, list) or any(type(p) is not int or p not in range(12)
                                              for p in ports)
            or len(set(ports)) != len(ports)):
        raise ValueError("distinct integer native ports required")
    return tuple(sorted(ports))


@lru_cache(maxsize=1)
def source():
    kernel = transfer.native_counts()
    degree, sigma = preparation.anchor()
    if tuple(60-kernel.counts[i][i] for i in range(66)) != tuple(degree):
        raise ArithmeticError("Palm reference and actual native records disagree")
    seams = {frozenset(edge) for edge in kernel.seams}
    occupied = tuple(tuple(i for i, bit in enumerate(s) if bit) for s in kernel.states)
    if any(d != 10-2*int(frozenset(pair) in seams)
           for pair, d in zip(occupied, degree, strict=True)):
        raise ArithmeticError("Palm generator is not the native adjacency projector")
    return kernel, tuple(degree), tuple(sigma), {pair: i for i, pair in enumerate(occupied)}


@dataclass(frozen=True)
class Block:
    charge: int
    local: tuple
    outside: tuple
    indices: tuple
    local_weights: tuple
    outside_weights: tuple
    probability: F

    @property
    def shape(self):
        return len(self.local), len(self.outside)


@lru_cache(maxsize=4096)
def _blocks(ports):
    _, degree, _, index = source()
    complement = tuple(p for p in range(12) if p not in ports)
    out = []
    for q in range(3):
        local = tuple(combinations(ports, q))
        outside = tuple(combinations(complement, 2-q))
        if not local or not outside:
            continue
        groups = {}
        for y in outside:
            weights = tuple(degree[index[tuple(sorted(x+y))]] for x in local)
            key = tuple(F(w, sum(weights)) for w in weights)
            groups.setdefault(key, []).append(y)
        for mu, ys in sorted(groups.items()):
            ys = tuple(ys)
            indices = tuple(tuple(index[tuple(sorted(x+y))] for y in ys) for x in local)
            column_totals = tuple(sum(degree[row[j]] for row in indices) for j in range(len(ys)))
            total = sum(column_totals)
            nu = tuple(F(w, total) for w in column_totals)
            block = Block(q, local, ys, indices, mu, nu, F(total, 600))
            if any(F(degree[i], 600) != block.probability*mu[a]*nu[b]
                   for a, row in enumerate(indices) for b, i in enumerate(row)):
                raise ArithmeticError("native regional reference does not factor in its context")
            out.append(block)
    if sum(block.probability for block in out) != 1:
        raise ArithmeticError("regional blocks lose native reference probability")
    covered = [i for block in out for row in block.indices for i in row]
    if sorted(covered) != list(range(66)):
        raise ArithmeticError("regional blocks are not a partition of native states")
    return tuple(out)


def blocks(ports):
    return _blocks(region(ports))


def basis_supports(ports, *, completed=True):
    """Disjoint supports of regional matrix units, with value one on support."""
    bs = blocks(ports)
    if completed:
        return tuple(tuple((block.indices[a][j], block.indices[b][j])
                           for j in range(len(block.outside)))
                     for block in bs for a in range(len(block.local))
                     for b in range(len(block.local)))
    by_charge = {}
    for block in bs:
        by_charge.setdefault(block.charge, []).append(block)
    return tuple(tuple(pair for block in group for pair in
                       ((block.indices[a][j], block.indices[b][j])
                        for j in range(len(block.outside))))
                 for group in by_charge.values() for a in range(len(group[0].local))
                 for b in range(len(group[0].local)))


def membership(support, ports, *, completed=True):
    """Exact membership for a zero/one sparse matrix using disjoint basis units."""
    remaining = set(support)
    for unit in basis_supports(ports, completed=completed):
        common = remaining.intersection(unit)
        if common and common != set(unit):
            return False
        remaining.difference_update(unit)
    return not remaining


def noninvariance_witness(ports):
    """A bare regional matrix unit with two unequal native modular ratios."""
    _, _, sigma, _ = source()
    for unit in basis_supports(ports, completed=False):
        ratios = {}
        for i, j in unit:
            ratios.setdefault(sigma[i]/sigma[j], (i, j))
        if len(ratios) > 1:
            first, second = sorted(ratios)[:2]
            return {"matrix_unit_support": unit,
                    "first_pair": ratios[first], "second_pair": ratios[second],
                    "first_ratio": first, "second_ratio": second}
    return None


def bare_reference(ports):
    """Unnormalized regional density blocks after forgetting outside context."""
    bs = blocks(ports)
    return tuple((q, tuple(sum((b.probability*b.local_weights[a]
                               for b in bs if b.charge == q), F())
                           for a in range(len(next(b for b in bs if b.charge == q).local))))
                 for q in sorted({b.charge for b in bs}))


def tracial_completion_distance(ports):
    """Exact disturbance of Palm by the bare tracial complement completion."""
    bs = blocks(ports)
    _, _, sigma, _ = source()
    out = [F()] * 66
    for q, weights in bare_reference(ports):
        group = [b for b in bs if b.charge == q]
        multiplicity = sum(len(b.outside) for b in group)
        for a, w in enumerate(weights):
            for block in group:
                for i in block.indices[a]:
                    out[i] = w/multiplicity
    return sum(abs(x-y) for x, y in zip(out, sigma, strict=True))/2


def locality_counterexample():
    """Disjoint bare regions whose invariant completions do not commute."""
    left, right = (0, 1), (2, 8)
    a = ((64, 64), (63, 63), (62, 62), (61, 61), (60, 60), (59, 59))
    b = ((64, 58), (44, 31), (42, 18), (41, 13))
    if not membership(a, left) or not membership(b, right):
        raise ArithmeticError("native locality witness left its regional algebra")
    ab = {(i, l) for i, j in a for k, l in b if j == k}
    ba = {(k, j) for k, l in b for i, j in a if l == i}
    if ab != {(64, 58)} or ba:
        raise ArithmeticError("native completed-algebra commutator witness changed")
    return {"left_ports": left, "right_ports": right,
            "left_matrix_unit": a, "right_matrix_unit": b,
            "commutator_matrix_unit": [64, 58]}


def conditional_expectation(array, ports):
    """Heisenberg state-preserving CE onto C_R; numerical matrix evaluation."""
    array = np.asarray(array)
    if array.shape != (66, 66):
        raise ValueError("a full native 66-by-66 observable required")
    result = np.zeros((66, 66), dtype=np.result_type(array.dtype, complex))
    for block in blocks(ports):
        for a, rows in enumerate(block.indices):
            for b, columns in enumerate(block.indices):
                value = sum(float(weight)*array[i, j] for weight, i, j in
                            zip(block.outside_weights, rows, columns, strict=True))
                for i, j in zip(rows, columns, strict=True):
                    result[i, j] = value
    return result


def restrict_state(array, ports):
    """CPTP restriction to direct-sum local factors, retaining every context."""
    array = np.asarray(array)
    if array.shape != (66, 66):
        raise ValueError("a full native 66-by-66 state required")
    return tuple(np.array([[sum(array[i, j] for i, j in zip(rows, columns, strict=True))
                            for columns in block.indices] for rows in block.indices])
                 for block in blocks(ports))


def restrict_encoded_state(array, ports):
    """Terminal CPTP restriction, including one explicit nonnative failure flag."""
    array = np.asarray(array)
    if array.shape != (512, 512):
        raise ValueError("a complete 512-by-512 encoded state required")
    codes = full_flow.support()[0]
    native = array[np.ix_(codes, codes)]
    outside = tuple(i for i in range(512) if i not in set(codes))
    failure = np.array([[sum(array[i, i] for i in outside)]])
    return (*restrict_state(native, ports), failure)


def restriction_kraus(ports):
    """Direct-sum output encoding; no context is postselected away."""
    bs = blocks(ports)
    dimension = sum(len(block.local) for block in bs)
    out, offset = [], 0
    for block in bs:
        for j in range(len(block.outside)):
            matrix = np.zeros((dimension, 66), complex)
            for a, row in enumerate(block.indices):
                matrix[offset+a, row[j]] = 1
            out.append(matrix)
        offset += len(block.local)
    return tuple(out)


def regional_budget(ports, s=F(16), attempts=369):
    """Uniform comparison to the modular action of the actual regional state."""
    blocks(ports)
    prepared = preparation.budget(attempts, abs(full_flow.parameter(s)))
    if not prepared["successful_work_reference_faithful"]:
        raise ValueError("preparation does not certify faithful regional references")
    bound = full_flow.checked_bound(s)
    return {"warmup_native_attempts": attempts,
            "modular_parameter": F(s),
            "prepared_regional_modular_half_diamond_upper":
                prepared["modular_channel_half_diamond_error_upper"],
            "source_then_restrict_vs_actual_regional_modular_half_diamond_upper":
                bound["half_diamond_upper"]+prepared["modular_channel_half_diamond_error_upper"],
            "all_regional_contexts_and_arbitrary_probe_references_retained": True,
            "compiled_terminal_nonnative_padding_retained_as_failure": True,
            "actual_global_prepared_modular_flow_preserves_regional_algebra_claimed": False}


def _summary(ports):
    bs = blocks(ports)
    witness = noninvariance_witness(ports)
    _, _, sigma, _ = source()
    for unit in basis_supports(ports):
        if len({sigma[i]/sigma[j] for i, j in unit}) != 1:
            raise ArithmeticError("completed regional algebra is not modular invariant")
    census = lambda values: {str(w): values.count(w) for w in sorted(set(values))}
    return {"ports": list(ports), "bare_algebra_dimension": len(basis_supports(ports, completed=False)),
            "completed_algebra_dimension": len(basis_supports(ports)),
            "bare_algebra_modular_invariant": witness is None,
            "bare_reference_blocks": [{"charge": q, "dimension": len(weights),
                                       "weight_census": census(weights)}
                                      for q, weights in bare_reference(ports)],
            "bare_complement_tracial_completion_Palm_half_trace_distance":
                str(tracial_completion_distance(ports)),
            "blocks": [{"local_charge": b.charge, "local_dimension": b.shape[0],
                        "outside_multiplicity": b.shape[1],
                        "probability": str(b.probability),
                        "local_reference_weight_census": census(b.local_weights),
                        "outside_reference_weight_census": census(b.outside_weights)} for b in bs],
            "bare_noninvariance_witness": witness}


def run():
    regions = ((), (0,), (0, 1), CAP, tuple(range(12)))
    rows = [_summary(r) for r in regions]
    inclusions = []
    for small, large in zip(regions, regions[1:], strict=False):
        if not set(small) <= set(large):
            raise ArithmeticError("receipt regions are not nested")
        if not all(membership(unit, large) for unit in basis_supports(small)):
            raise ArithmeticError("native modular completion violates isotony")
        inclusions.append([list(small), list(large)])
    paths = ("model/native_modular_reference/regional.py",
             "model/native_modular_reference/preparation.py",
             "model/native_modular_reference/full_flow.py",
             "model/primitive_observer/transfer.py", "model/observer/native.py")
    hashes = {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}
    result = {"schema": "oph.native_modular_reference.regional.v1", "all_pass": True,
            "source_sha256": hashes,
            "reference": "diag(degree/600)",
            "source_generator_up_to_scalar": "log(5/4) times native adjacency projector",
            "regions": rows, "checked_nested_inclusions": inclusions,
            "disjoint_completed_algebra_noncommutation": locality_counterexample(),
            "cap_modular_budget": preparation._compact(regional_budget(CAP)),
            "scope": {"minimal_modular_invariant_hull_of_bare_regional_algebra": True,
                      "same_ambient_modular_group_agrees_on_intersections": True,
                      "state_preserving_CP_expectation_constructed": True,
                      "regional_coherent_access_and_context_read_implementation_derived": False,
                      "bare_single_bit_record_algebra_has_nontrivial_modular_clock": False,
                      "disjoint_regional_modular_hulls_claimed_to_commute": False,
                      "passive_repair_identified_with_modular_flow": False,
                      "physical_seconds_temperature_or_BW_identified": False}}
    # Portable retained schema: exact rationals are strings and tuples lists.
    return json.loads(preparation.canonical(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        full_flow.corner.compare_receipt(strict_json.load(args.check), result)
        print("Native regional Palm algebras and common modular action verified.")
    elif args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(full_flow.corner.canonical(result)+"\n")
    else:
        print(full_flow.corner.canonical(result))


if __name__ == "__main__":
    main()
