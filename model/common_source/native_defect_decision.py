"""Exact decision experiment for a defect in the existing native load law.

The candidate is one double load above the stationary binary total-two
reference. All 78 states, all thirty seams, both coins and both complete C/E
record menus are retained. No defect register, puncture, link phase, energy
scale or new repair rule is supplied. The source support determines the
configuration-sector obstruction; lifetimes and stationary averages use the
separately declared uniform-attempt calibration, not the production E/F policy.

This exhausts classical configuration invariants, not arbitrary quantum,
record-encoded, actively corrected or enlarged-carrier excitations. The
same collision candidate's survival effect is also checked on the complete
quantum instrument, including arbitrary untouched references. The existing
corrected matching memory is audited separately as a control.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path

from tools import strict_json
from ..observer import native
from ..observer.source_discrimination import circuit
from ..repair_scheduler import coins
from . import record_code, slow_geometry
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def _seams(seams):
    allowed = tuple(slow_geometry.source().seams())
    value = allowed if seams is None else tuple(seams)
    if not value or len(set(value)) != len(value) or any(edge not in allowed for edge in value):
        raise ValueError("a nonempty distinct subset of the actual native seams required")
    return value


def transitions(*, binary=False, seams=None):
    """Complete integer row counts; denominator is twice the seam count."""
    if type(binary) is not bool:
        raise ValueError("an explicit binary-sector Boolean required")
    seams = _seams(seams)
    states = native.states(binary=binary)
    index = {state: i for i, state in enumerate(states)}
    counts = [[0]*len(states) for _ in states]
    for i, state in enumerate(states):
        for seam in seams:
            for coin in (0, 1):
                out, _ = native.native(state, seam, coin)
                counts[i][index[out]] += 1
    return tuple(map(tuple, counts))


def strongly_connected_components(counts):
    """Exact directed reachability; no eigenvalue or numerical rank threshold."""
    size = len(counts)
    if (not size or any(len(row) != size for row in counts)
            or any(type(x) is not int or x < 0 for row in counts for x in row)):
        raise ValueError("a nonempty square nonnegative integer transition matrix required")
    reach = []
    for start in range(size):
        seen, queue = {start}, [start]
        while queue:
            current = queue.pop()
            for target, count in enumerate(counts[current]):
                if count and target not in seen:
                    seen.add(target)
                    queue.append(target)
        reach.append(seen)
    remaining, result = set(range(size)), []
    while remaining:
        first = min(remaining)
        component = tuple(i for i in sorted(remaining) if i in reach[first] and first in reach[i])
        result.append(component)
        remaining.difference_update(component)
    return tuple(result)


def complete_record_maps():
    """All raw quantum branch supports, including identically zero outcomes.

    Each listed integer partial injection has Kraus amplitude 1/sqrt(2).
    Seams are interventions here; the uniform schedule adds weight 1/30.
    """
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    rows = []
    for law in ("C", "E"):
        for seam in _seams(None):
            reference = {(coin, label): tuple(mapping)
                         for coin, label, mapping in coins.coin_maps(law, seam, binary=False)}
            for coin in (0, 1):
                raw_maps = {raw: [] for raw in range(1 << circuit.RECORD_BITS[law])}
                for i, state in enumerate(states):
                    out, (_, current) = native.native(state, seam, coin)
                    raw = circuit.CURRENT_CODES[current]
                    if law == "E":
                        raw += 4*((state[seam[0]]+state[seam[1]]) % 2)
                    raw_maps[raw].append((i, index[out]))
                for raw, mapping in raw_maps.items():
                    label = circuit.record_values(law, coin).get(raw)
                    if tuple(mapping) != reference.get((coin, label), ()):
                        raise ArithmeticError("complete raw map differs from existing C/E instrument")
                    if len({target for _, target in mapping}) != len(mapping):
                        raise ArithmeticError("recorded branch is not a partial injection")
                    rows.append({"law": law, "seam": list(seam), "coin": coin,
                                 "raw": raw, "pairs": [list(pair) for pair in mapping]})
    return rows


def _potential(state):
    return F(sum(x*x for x in state), 2)


def quantum_collision_survival(maps):
    """Exact Heisenberg collision effect on all78 matrix entries, per C/E law."""
    collision = {i for i, state in enumerate(native.states()) if max(state) == 2}
    results = {}
    for law in ("C", "E"):
        effect = Counter()
        binary_to_collision = 0
        rows = [row for row in maps if row["law"] == law]
        for row in rows:
            preimages = {}
            for source, target in row["pairs"]:
                if target in collision:
                    preimages.setdefault(target, []).append(source)
                    binary_to_collision += source not in collision
            # Calculate M*P_collision*M explicitly, including cross terms;
            # do not assume classical diagonality of this quantum effect.
            for sources in preimages.values():
                effect.update((i, j) for i in sources for j in sources)
        expected = Counter({(i, i): 50 for i in collision})
        if binary_to_collision or effect != expected:
            raise ArithmeticError("quantum collision survival effect differs from50*P_collision")
        results[law] = {"raw_branches_checked": len(rows),
                       "integer_effect_entries": [[i, j, value] for (i, j), value in sorted(effect.items())],
                       "effect_denominator": 60,
                       "nonzero_off_diagonal_entries": sum(i != j for i, j in effect),
                       "binary_to_collision_branch_addresses": binary_to_collision}
    return {"projector": "P_collision onto all12 double-load configurations",
            "uniform_one_attempt_Heisenberg_identity": "sum_h K_h^dagger P_collision K_h=(5/6)P_collision",
            "normalization": "K_h=M_h/sqrt(60), complete uniform30-seam/fair-coin instrument",
            "C_E": results,
            "arbitrary_quantum_and_untouched_reference_inputs": True,
            "n_attempt_collision_probability": "(5/6)^n * Tr[(P_collision tensor I_reference) rho]",
            "scope": "no persistent quantum subspace supported wholly on the collision configurations under this uniform unread evolution; other quantum/history sectors remain open"}


def _edge_cost(state, seams):
    return sum((state[a]-state[b])**2 for a, b in seams)


def corrected_memory_control():
    """Exact Gram/cost control on record_code's fixed six matching words."""
    seams = _seams(None)
    for support in record_code.SUPPORTS:
        if any(sum(state[p] for state in support) != 1 for p in range(12)):
            raise ArithmeticError("one-port matching compression changed")
        if any(sum((state[a]-state[b])**2 for state in support) != 2 for a, b in seams):
            raise ArithmeticError("native-edge matching compression changed")
    branches, zero = 0, 0
    distributions = {}
    escapes = {}
    encoded_supports = [set(circuit.encode(state) for state in support)
                        for support in record_code.SUPPORTS]
    for law in ("C", "E"):
        weights = Counter()
        escape = F(0)
        for seam in seams:
            for coin in (0, 1):
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    images = record_code.branch_images(law, seam, coin, raw)
                    counts = {len(row) for row in images}
                    if len(counts) != 1:
                        raise ArithmeticError("record probability depends on logical column")
                    probability = F(len(images[0]), 12)
                    weights[str(probability)] += 1
                    branches += 1
                    zero += not probability
                    for image in images:
                        # K|logical> has amplitude1/sqrt(12) on each image.
                        # Each codeword has amplitude1/sqrt(6) on its support.
                        projected_norm = sum(F(len(set(image) & support)**2, 72)
                                             for support in encoded_supports)
                        escape += (F(len(image), 12)-projected_norm)/180
        distributions[law] = dict(sorted(weights.items()))
        if not 0 < escape < 1:
            raise ArithmeticError("unrecovered code escape control changed")
        escapes[law] = str(escape)
    return {"code": "record_code.SUPPORTS; fixed six nonedge matching words",
            "all_raw_branches_checked": branches, "zero_raw_branches_retained": zero,
            "branch_probability_counts": distributions,
            "one_port_compression": "I6/6", "each_native_edge_cost_compression": "I6/3",
            "unweighted_edge_disagreement_compression": "10*I6",
            "quadratic_load_potential_compression": "I6",
            "every_record_cost_expectation_independent_of_logical_state": True,
            "unrecovered_uniform_C_E_code_escape_probability_from_Pcode_over_6": escapes,
            "conditional_persistence_requires_existing_record_selected_recovery": True,
            "autonomous_defect_or_mass_selected": False,
            "cost_scope": "these scalar record/diagonal costs cannot split logical-state masses; a common energy or mass is not excluded and no physical offset is derived",
            "new_source_word_execution_performed": False}


def report():
    seams = _seams(None)
    states = native.states()
    counts = transitions()
    components = strongly_connected_components(counts)
    closed = [c for c in components if all(not counts[i][j]
              for i in c for j in range(78) if j not in c)]
    binary_indices = tuple(i for i, state in enumerate(states) if max(state) == 1)
    if len(closed) != 1 or closed[0] != binary_indices:
        raise ArithmeticError("native source no longer has the unique binary closed class")
    # Every collision reaches that class, so any scalar constant on every
    # allowed branch is constant on the entire78-state support.
    if sorted(map(len, components)) != [1]*12+[66]:
        raise ArithmeticError("native communicating classes changed")
    binary_counts = transitions(binary=True)
    if any(sum(row) != 60 for row in binary_counts) or any(sum(column) != 60 for column in zip(*binary_counts)):
        raise ArithmeticError("uniform binary stationary law changed")
    collision = (2,)+(0,)*11
    collision_index = states.index(collision)
    survive = F(counts[collision_index][collision_index], 60)
    if survive != F(5, 6) or any(counts[collision_index][i] for i in range(78)
                               if i != collision_index and i not in binary_indices):
        raise ArithmeticError("double-load candidate no longer splits on its first incident seam")
    moves, drops, changes = 0, set(), set()
    witnesses = {}
    for i in binary_indices:
        state = states[i]
        for seam in seams:
            for coin in (0, 1):
                out, (_, current) = native.native(state, seam, coin)
                moves += current != 0
                drops.add(str(_potential(state)-_potential(out)))
                change = _edge_cost(out, seams)-_edge_cost(state, seams)
                changes.add(change)
                name = "edge_cost_rise_witness" if change > 0 else "edge_cost_fall_witness"
                if change and name not in witnesses:
                    witnesses[name] = {"before": list(state), "after": list(out),
                                       "seam": list(seam), "coin": coin, "current": current,
                                       "edge_cost_before": _edge_cost(state, seams),
                                       "edge_cost_after": _edge_cost(out, seams),
                                       "quadratic_load_potential_drop": "0"}
    held = tuple(edge for edge in seams if 0 not in edge)
    control_components = strongly_connected_components(transitions(binary=True, seams=held))
    faces = [face for face in combinations(range(12), 3)
             if all(tuple(sorted(edge)) in seams for edge in combinations(face, 2))]
    if len(faces) != 20 or any((state[a]-state[b])+(state[b]-state[c])+(state[c]-state[a])
                              for state in states for a, b, c in faces):
        raise ArithmeticError("native load-difference face circulation changed")
    maps = complete_record_maps()
    if sum(len(row["pairs"]) for row in maps) != 2*78*60:
        raise ArithmeticError("C/E support enumeration omitted source outcomes")
    source_paths = [Path(__file__), Path(native.__file__), Path(slow_geometry.__file__),
                    Path(coins.__file__), Path(circuit.__file__), Path(record_code.__file__),
                    # Owning research note: provenance only; this finite
                    # calculation does not verify every result in the guide.
                    ROOT/"docs/research/OPH_FIELD_RECONSTRUCTION.md"]
    return {"schema": "oph.native_defect_decision.v1",
            "decision": "native double-load defect decays; no nonconstant configuration invariant or maintenance-energy selection",
            "carrier": {"states": 78, "binary_states": 66, "seams": 30,
                        "primitive_branch_addresses": 4680},
            "sectors": {"component_sizes": sorted(map(len, components)),
                        "closed_component_sizes": list(map(len, closed)),
                        "closed_component_indices": [list(c) for c in closed],
                        "diagonal_invariant_dimension": 1, "unique_closed_sector_binary66": True,
                        "invariant_scope": "functions f(x) satisfying f(output)=f(input) for every admitted seam and coin",
                        "uniform_binary_reference": "1/66 on each binary state, zero on collisions"},
            "collision_candidate": {"initial_state": list(collision),
                        "definition": "one double load above the closed binary reference shell",
                        "one_attempt_survival": str(survive),
                        "mean_attempts_to_split": str(1/(1-survive)),
                        "survival_by_attempt": {str(n): str(survive**n) for n in range(13)},
                        "excess_quadratic_potential_before": "1", "after_split": "0",
                        "double_load_hops_before_splitting": False,
                        "all_12_locations_have_same_split_hazard_by_degree_five": all(sum(p in e for e in seams) == 5 for p in range(12))},
            "transport_cost": {"stationary_nonzero_current_probability": str(F(moves, 66*60)),
                        "stationary_mean_squared_current_per_attempt": str(F(moves, 66*60)),
                        "binary_quadratic_drop_values": sorted(drops),
                        "binary_edge_cost_changes": sorted(changes),
                        "binary_edge_cost_values": sorted({_edge_cost(states[i], seams) for i in binary_indices}),
                        "native_attempt_is_not_energy_unit": True, **witnesses},
            "schedule_control": {"held_port": 0, "removed_incident_seams": 5,
                        "binary_component_sizes": sorted(map(len, control_components)),
                        "changed_source_schedule": True,
                        "interpretation": "holding one port creates occupation sectors by disabling its repairs; this is a changed law, not discovered protection"},
            "load_difference_topology": {"edge_variable": "b_ab=n_a-n_b",
                        "states_checked": 78, "faces_checked_per_state": 20,
                        "face_circulations": [0]*20,
                        "every_closed_loop_telescopes": True,
                        "independent_link_connection_or_puncture_supplied": False,
                        "scope": "the actual load-gradient map only; not an exclusion of other bundles or enlarged states"},
            "corrected_memory_control": corrected_memory_control(),
            "quantum_collision_survival": quantum_collision_survival(maps),
            "complete_record_maps": maps,
            "scope": {"classical_configuration_invariants_exhausted": True,
                        "all_quantum_or_record_encoded_sectors_excluded": False,
                        "physical_energy_or_mass_derived": False,
                        "parameter": "native attempt count under uniform independent seams/fair coins",
                        "same_uniform_schedule_as_production_feedback": False,
                        "positive_all_seam_support_required_for_invariant_obstruction": True,
                        "next_obligation": "identify a nontrivial protected state/observable sector under a declared continuing law; an untested quantum/history sector or a changed law remains possible; propagation and conserved energy with clock/kinetic units must then be derived jointly, not fitted to repair counts"},
            "source_sha256": bound_pins({path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in source_paths})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--verify", "--check", dest="verify", type=Path)
    args = parser.parse_args()
    result = report()
    if args.verify:
        if strict_json.load(args.verify) != result:
            raise ValueError("native defect decision or source custody changed")
        print("Verified exact native defect decision, complete raw maps and bounded interpretation")
    else:
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
