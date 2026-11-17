"""Exact adjacent-pair localization, exit and transport on the native carrier.

The target is fixed by adjacency before applying repair. First exit is a
classical path diagnostic, not an extra measurement or controller service.
All C/E coin/raw branches, including zero maps, enter the certificates.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from . import contract
from ..observer import native, record_recovery as records
from ..observer import adaptive_record_recovery as adaptive
from ..observer.source_discrimination import circuit
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "localized_results.json"
DENOMINATOR = 60


def basis():
    """Native lexicographic binary total-two basis; no supplied defect register."""
    return native.states(binary=True)


def _state(value):
    if type(value) is not tuple or value not in basis() or any(type(x) is not int for x in value):
        raise ValueError("an original binary total-two state tuple required")
    return value


def occupied(state):
    return tuple(i for i, value in enumerate(_state(state)) if value)


def cluster_indices():
    """The two occupied ports share a native edge, independently of dynamics."""
    edges = set(records.seams())
    return tuple(i for i, state in enumerate(basis()) if occupied(state) in edges)


def pair_distances():
    neighbors = [set() for _ in range(12)]
    for a, b in records.seams():
        neighbors[a].add(b)
        neighbors[b].add(a)
    distances = []
    for state in basis():
        a, b = occupied(state)
        pending, distance = deque([a]), {a: 0}
        while pending:
            source = pending.popleft()
            for target in neighbors[source] - distance.keys():
                distance[target] = distance[source] + 1
                pending.append(target)
        distances.append(distance[b])
    return tuple(distances)


def transition_counts(law):
    """Uniform native population kernel: immutable rows / DENOMINATOR (=60).

    Build from the complete quantum instrument, not a replacement walk.
    Each raw partial injection has squared coin amplitude 1/2, and each
    seam weight is 1/30. C and E have different record blocks.
    """
    law = records.candidate(law)
    counts = [[0] * 66 for _ in range(66)]
    for branch in adaptive.branch_maps(binary=True):
        if adaptive.actions()[branch.action][0] == law:
            for source, target in zip(branch.inputs, branch.outputs, strict=True):
                counts[source][target] += 1
    if any(sum(row) != DENOMINATOR for row in counts):
        raise ArithmeticError("incomplete native population instrument")
    return tuple(map(tuple, counts))


def movement_activity(state):
    """Probability of a nonzero current per uniform-seam native attempt."""
    state = _state(state)
    return F(sum(native.native(state, edge, coin)[0] != state
                 for edge in records.seams() for coin in (0, 1)), DENOMINATOR)


def _probability(value, *, positive=True):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or math.gcd(value.numerator, value.denominator) != 1
            or not 0 <= value <= 1 or positive and value == 0):
        raise ValueError("an original canonical probability Fraction required")
    return F(value.numerator, value.denominator)


def policy_bounds(epsilon, repair_probability=F(1)):
    """Uniform bound for arbitrary committed-history classical feedback.

    Each seam has conditional probability at least epsilon/30 before the
    fair coin; C/E selection may depend on the seam and history. There is
    no intermediate source-changing operation beyond native repair/holds.
    Halt is not a continuing attempt and earns no lifetime extension claim.
    """
    epsilon = _probability(epsilon)
    repair_probability = _probability(repair_probability)
    hazard = epsilon * repair_probability / 15
    return {"epsilon": epsilon, "repair_probability_floor": repair_probability,
            "exit_hazard_lower_per_continuing_attempt": hazard,
            "survival_factor_upper": 1 - hazard,
            "mean_first_exit_attempts_upper": 1 / hazard,
            "conditions": "no intervening quantum control; continuing attempts on every compared path",
            "halted_mass_must_remain_separate": True}


def complete_branch_certificate():
    """Exact diagonal exit effects, complete raw maps, and quantum distinction."""
    cluster, states = set(cluster_indices()), basis()
    result = {}
    all_rows = adaptive.branch_maps(binary=True)
    for law in ("C", "E"):
        rows = [b for b in all_rows if adaptive.actions()[b.action][0] == law]
        exits = {edge: [0] * 66 for edge in records.seams()}
        normalization = {edge: [0] * 66 for edge in records.seams()}
        agreement = [[0] * 66 for _ in range(66)]
        canonical = []
        for row in rows:
            _, edge = adaptive.actions()[row.action]
            if len(set(row.outputs)) != len(row.outputs):
                raise ArithmeticError("noninjective recorded native map")
            canonical.append([list(edge), row.coin, row.raw,
                              list(zip(row.inputs, row.outputs, strict=True))])
            for source, target in zip(row.inputs, row.outputs, strict=True):
                normalization[edge][source] += 1
                exits[edge][source] += int(target not in cluster)
            for i in row.inputs:
                for j in row.inputs:
                    agreement[i][j] += 1
        if any(value != 2 for row in normalization.values() for value in row):
            raise ArithmeticError("incomplete C/E quantum instrument")
        if any(sum(row[i] for row in exits.values()) != 4 for i in cluster):
            raise ArithmeticError("uniform cluster exit effect changed")
        # Injection makes every K* P_out K diagonal, including coherent inputs.
        # There are no off-diagonal terms hidden by a population-only check.
        pair = max(((i, j) for i in cluster for j in cluster if i < j),
                   key=lambda ij: (agreement[ij[0]][ij[1]], -ij[0], -ij[1]))
        result[law] = {
            "raw_maps": len(rows), "zero_raw_maps": sum(not row.inputs for row in rows),
            "primitive_input_addresses": sum(len(row.inputs) for row in rows),
            "maps_sha256": hashlib.sha256(json.dumps(canonical, separators=(",", ":")).encode()).hexdigest(),
            "all_seam_normalization_counts": 2,
            "uniform_exit_effect_on_cluster": "P_A/15",
            "cluster_pair_with_maximum_one_event_record_agreement": list(pair),
            "maximum_cluster_pair_agreement": F(agreement[pair[0]][pair[1]], 60),
            "witness_occupied_ports": [list(occupied(states[i])) for i in pair],
        }
    return result


def uniform_certificate():
    """Exact killed chain, motion, three distance classes, and return times."""
    counts = transition_counts("C")
    if counts != transition_counts("E"):
        raise ArithmeticError("C/E diagonal population kernels differ")
    cluster, distances = cluster_indices(), pair_distances()
    quotient = []
    for distance in (1, 2, 3):
        candidates = {tuple(sum(counts[i][j] for j, d in enumerate(distances) if d == target)
                            for target in (1, 2, 3))
                      for i, d in enumerate(distances) if d == distance}
        if len(candidates) != 1:
            raise ArithmeticError("distance projection is not strongly lumpable")
        quotient.append(candidates.pop())
    killed = tuple(tuple(counts[i][j] for j in cluster) for i in cluster)
    if any(sum(row) != 56 for row in killed):
        raise ArithmeticError("first-exit hazard not constant across cluster")
    if any(counts[i][i] != 52 or sum(counts[i][j] for j in cluster if j != i) != 4
           for i in cluster):
        raise ArithmeticError("within-cluster transport census changed")
    # A moved cluster shares one occupied port with its predecessor.
    if any(len(set(occupied(basis()[i])) & set(occupied(basis()[j]))) != 1
           for i in cluster for j in cluster if i != j and counts[i][j]):
        raise ArithmeticError("cluster motion ceased to be a one-port move")
    stationary = (F(5, 11), F(5, 11), F(1, 11))
    if any(sum(stationary[i] * F(quotient[i][j], 60) for i in range(3)) != stationary[j]
           for j in range(3)):
        raise ArithmeticError("distance stationary law changed")
    # h=0 on A; solve/check the exact hitting-time equations on distances2,3.
    hitting = (F(0), F(18), F(24))
    if any(hitting[i] != 1 + sum(F(quotient[i][j], 60) * hitting[j] for j in range(3))
           for i in (1, 2)):
        raise ArithmeticError("return hitting-time equations failed")
    return {
        "distance_class_sizes": dict(sorted(Counter(distances).items())),
        "distance_transition_counts": quotient, "common_denominator": 60,
        "killed_cluster_indices": cluster, "killed_transition_counts": killed,
        "one_attempt": {"unchanged": F(13, 15), "transport_within_cluster": F(1, 15),
                        "first_exit": F(1, 15)},
        "survival_after_n_attempts": "(14/15)^n",
        "first_exit_mass_at_n_ge_1": "(14/15)^(n-1)/15",
        "mean_first_exit_attempts": F(15),
        "mean_within_cluster_moves_before_exit": F(1),
        "exit_distance_class": 2,
        "mean_return_attempts_after_first_exit": hitting[1],
        "mean_hit_cluster_from_antipodal_pair": hitting[2],
        "stationary_cluster_mass": stationary[0],
        "conditioned_survivor_transport_probability_per_attempt": F(1, 14),
        "return_is_new_entry_not_survival_or_protection": True,
    }


def feedback_counts(law, epsilon):
    """Complete source/current-bit kernel, indices 2*source+memory.

    Prefer (0,5) after a positive current and (0,1) otherwise. Fresh uniform
    exploration has weight epsilon. No hidden cluster/exit flag is read.
    Return integer rows and denominator, keeping all raw outcome branches.
    """
    law, epsilon = records.candidate(law), _probability(epsilon)
    p, q = epsilon.numerator, epsilon.denominator
    rows = [[0] * 132 for _ in range(132)]
    states = basis()
    for branch in adaptive.branch_maps(binary=True):
        name, edge = adaptive.actions()[branch.action]
        if name != law:
            continue
        if not branch.inputs:
            continue  # Retained zero map: no population contribution.
        value = circuit.record_values(law, branch.coin)[branch.raw]
        current = records.record_current(law, (*edge, branch.coin, value))
        for source, target in zip(branch.inputs, branch.outputs, strict=True):
            if current != states[source][edge[0]] - states[target][edge[0]]:
                raise ArithmeticError("decoded record current differs from native transport")
            for memory in (0, 1):
                preferred = (0, 5) if memory else (0, 1)
                weight = p + (30 * (q - p) if edge == preferred else 0)
                rows[2 * source + memory][2 * target + int(current > 0)] += weight
    denominator = 60 * q
    if any(sum(row) != denominator for row in rows):
        raise ArithmeticError("feedback kernel lost complete branch mass")
    return tuple(map(tuple, rows)), denominator


def feedback_history(law, epsilon, horizons=(0, 1, 2, 4, 8, 16, 32, 64)):
    """Unconditioned population and offline first-exit marginals, exactly."""
    if (type(horizons) is not tuple or not horizons
            or any(type(n) is not int or not 0 <= n <= 256 for n in horizons)
            or tuple(sorted(set(horizons))) != horizons):
        raise ValueError("increasing distinct integer horizons from0 through256 required")
    counts, denominator = feedback_counts(law, epsilon)
    sparse = [tuple((j, value) for j, value in enumerate(row) if value) for row in counts]
    cluster = {2 * i + m for i in cluster_indices() for m in (0, 1)}
    initial = basis().index(tuple(int(p in (0, 1)) for p in range(12))) * 2
    full, surviving = [0] * 132, [0] * 132
    full[initial] = surviving[initial] = 1
    mass, result = 1, []
    for n in range(max(horizons) + 1):
        if n in horizons:
            cluster_mass, survival = sum(full[i] for i in cluster), sum(surviving)
            if sum(full) != mass or not 0 <= survival <= cluster_mass <= mass:
                raise ArithmeticError("offline first-exit bookkeeping lost probability")
            result.append({"attempts": n, "survival": F(survival, mass),
                           "ever_exited": 1 - F(survival, mass),
                           "currently_clustered": F(cluster_mass, mass),
                           "returned_and_currently_clustered": F(cluster_mass - survival, mass),
                           "currently_separated": 1 - F(cluster_mass, mass)})
        if n == max(horizons):
            break
        next_full, next_survival = [0] * 132, [0] * 132
        for i, row in enumerate(sparse):
            for j, weight in row:
                next_full[j] += full[i] * weight
                if j in cluster:
                    next_survival[j] += surviving[i] * weight
        full, surviving, mass = next_full, next_survival, mass * denominator
    return result


def exterior_certificate():
    """The candidate fails a reference-like exterior; quantify the failure."""
    core = (0, 1)
    outside = tuple(i for i in range(12) if i not in core)
    states, cluster = basis(), set(cluster_indices())
    from ..native_modular_reference import history_transport
    prepared = history_transport.reference("local")
    uniform, mixture, local = Counter(), Counter(), Counter()
    for i, state in enumerate(states):
        key = tuple(state[p] for p in outside)
        uniform[key] += F(1, 66)
        local[key] += prepared[i]
        if i in cluster:
            mixture[key] += F(1, 30)
    keys = uniform.keys() | mixture.keys()
    distance = sum(abs(uniform[key] - mixture[key]) for key in keys) / 2
    local_distance = sum(abs(local[key] - mixture[key]) for key in keys) / 2
    if any(sum(F(1, 30) * states[i][port] for i in cluster) != F(1, 6)
           for port in range(12)):
        raise ArithmeticError("cluster mixture lost its single-port symmetry")
    return {"chosen_core": core, "exterior_ports": outside,
            "reference": "I66/66, the uniform total-two attempt-stationary state",
            "fixed_core_candidate_exterior": "vacuum with probability1",
            "reference_exterior_vacuum_probability": F(1, 66),
            "fixed_core_exterior_trace_distance": F(65, 66),
            "uniform_cluster_mixture_exterior_trace_distance": distance,
            "shared_preparation_reference": "ideal three-charge diagonal target from history_transport.reference('local'); finite preparation error remains separate",
            "shared_reference_exterior_vacuum_probability": local[(0,) * 10],
            "fixed_core_exterior_distance_from_shared_reference": 1 - local[(0,) * 10],
            "cluster_mixture_exterior_distance_from_shared_reference": local_distance,
            "uniform_cluster_mixture_single_port_marginal": F(1, 6),
            "single_port_agreement_does_not_establish_exterior_agreement": True,
            "reference_like_exterior_established": False,
            "new_localization_detector_implemented": False}


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def run():
    epsilon = F(1, 16)
    exact = uniform_certificate()
    feedback = {law: feedback_history(law, epsilon) for law in ("C", "E")}
    if feedback["C"] != feedback["E"]:
        raise ArithmeticError("decoded-current feedback populations differ between C/E")
    rate = policy_bounds(epsilon)["survival_factor_upper"]
    if any(row["survival"] > rate ** row["attempts"] for row in feedback["C"]):
        raise ArithmeticError("feedback example violates uniform hazard theorem")
    freeze = tuple(int(p in (0, 1)) for p in range(12))
    if any(native.native(freeze, (0, 1), c)[0] != freeze for c in (0, 1)):
        raise ArithmeticError("zero-exploration control is not frozen")
    return _json({
        "schema": "oph.clock_defect_energy.localized.v1", "campaign": contract.binding(),
        "definition": {"carrier_dimension": 66, "target_dimension": 30,
                       "target": "two occupied ports adjacent in the supplied native graph",
                       "time_parameter": "actual native repair count, or explicitly bounded continuing attempts",
                       "first_exit_tag": "offline diagonal path diagnostic; not controller state or new measurement"},
        "complete_C_E_instruments": complete_branch_certificate(),
        "uniform": exact,
        "feedback": {"epsilon": epsilon, "initial_occupied_ports": [0, 1], "initial_memory": 0,
                     "policy": "prefer(0,5) iff previous current>0, else(0,1); epsilon uniform exploration",
                     "controller_states": 2, "histories": feedback,
                     "complete_history_not_materialized": True,
                     "population_DP_retains_all_mass": True,
                     "first_exit_tag_never_read_by_policy": True},
        "arbitrary_feedback_bounds": {"repair_only": policy_bounds(epsilon),
                                      "with_repair_floor": policy_bounds(epsilon, F(1, 16))},
        "quantum_memory_comparison": {
            "source": "model/observer/adaptive_record_recovery.py",
            "bounds_at_512_actual_repairs": {
                law: adaptive.bounds(epsilon, 512, law_scope=law, binary=True)
                for law in ("C", "E", "switching")},
            "conditions": "complete committed tape; persistent exploration; no intervening quantum correction or coherent clock controls",
            "population_lifetime_is_quantum_memory_lifetime": False,
            "finite_time_entanglement_breaking_inferred": False},
        "exterior": exterior_certificate(),
        "negative_controls": {
            "zero_exploration_same_occupied_seam": {"survival": "1 for every n", "movement": "0", "outside_positive_exploration": True},
            "halt": "source freezes but earns zero further actual repair events; retain halted mass separately",
            "returned_cluster": "currently_clustered differs from never_exited; all return mass retained"},
        "scope": {"localized_relational_population_candidate": True,
                  "topological_charge_or_autonomous_protection": False,
                  "reference_like_exterior": False, "physical_space_or_time_calibration": False,
                  "mass_or_energy_derived": False, "full_CP03_accepted": False,
                  "history_refinement_or_active_recovery_classes_excluded": False},
        "source_sha256": bound_pins(contract.source_bindings([
            "model/clock_defect_energy/localized.py", "model/clock_defect_energy/contract.py",
            "model/observer/native.py", "model/observer/record_recovery.py",
            "model/observer/adaptive_record_recovery.py", "model/observer/source_discrimination/circuit.py",
            "model/common_source/slow_geometry.py", "model/native_modular_reference/history_transport.py",
            "model/native_modular_reference/local_preparation.py", "tools/strict_json.py",
            "evidence/native_geometric_source_20260925/spec.json",
            "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
            "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py"])),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", nargs="?", const=str(OUTPUT))
    args = parser.parse_args()
    result = run()
    if args.check:
        contract.check_receipt(Path(args.check), result)
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("adjacent-pair localization: exact complete-branch checks passed")


if __name__ == "__main__":
    main()
