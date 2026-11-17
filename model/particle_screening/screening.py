"""Finite cap screening, all-seam code geometry, and bond-exit certificates.

The tested object is the binary total-two native carrier with a classical
two-port core. All probability and distance calculations are rational. The
code is a supplied quantum input with algebraic flagged recovery; this packet
does not compile a new source encoder or recovery word.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from math import comb
from pathlib import Path

from tools import strict_json
from ..exchange_sector import sector
from ..observer import adaptive_record_recovery as adaptive
from ..observer import record_recovery as records

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "screening_results.json"
SEAMS = frozenset(records.seams())
CORES = tuple(combinations(range(12), 2))
CORE = (0, 1)
PREDECESSOR_FILES = (
    "model/exchange_sector/sector.py",
    "model/exchange_sector/README.md",
    "model/common_source/record_code.py",
    "model/common_source/ENCODED_MEMORY.md",
    "model/observer/adaptive_record_recovery.py",
    "model/observer/record_recovery.py",
    "model/clock_defect_energy/defect_history.py",
)


def _edge(a, b):
    return tuple(sorted((a, b))) in SEAMS


def safe_exterior(center):
    """Exterior ports whose two possible q=1 source pairs are nonedges."""
    center = sector.core(center)
    return tuple(j for j in range(12) if j not in center
                 and all(not _edge(i, j) for i in center))


def sharp_code(center):
    """Two orthogonal equal-amplitude columns in a sharp q=1 sector.

    For two safe exterior ports the columns are one rectangle. For four they
    are the direct sum of two rectangles. No code is claimed for zero safe
    ports. The first port of center is retained as the first core label.
    """
    center = sector.core(center)
    safe = safe_exterior(center)
    if len(safe) not in (2, 4):
        return None
    groups = tuple(zip(safe[::2], safe[1::2], strict=True))
    columns = []
    for logical in (0, 1):
        pairs = []
        for u, v in groups:
            pairs.extend((tuple(sorted((center[0], u if logical == 0 else v))),
                          tuple(sorted((center[1], v if logical == 0 else u)))))
        columns.append(tuple(pairs))
    return tuple(columns)


def _reference_sector_masses(center=CORE):
    sigma = sector.reference(center)
    masses = tuple(sum(w for x, w in zip(sector.BASIS, sigma, strict=True)
                       if sector.charge(x, center) == q) for q in range(3))
    if masses != (F(2, 5), F(5, 9), F(2, 45)) or sum(masses) != 1:
        raise ArithmeticError("transported three-charge reference changed")
    return masses


def _hypergeom(outside_count, cap_size, hits):
    if (hits < 0 or hits > outside_count or cap_size-hits < 0
            or cap_size-hits > 10-outside_count):
        return F()
    return F(comb(outside_count, hits)*comb(10-outside_count, cap_size-hits),
             comb(10, cap_size))


def cap_certificate():
    """Exact minimax cap distances over *all* states in each sharp sector.

    Average a cap's exterior-number distribution over all k-subsets. For any
    sharp-q quantum state the result is hypergeometric and state independent.
    Convexity makes its total variation from the reference a lower bound on
    the worst cap's trace distance. The diagonal conditioned reference is
    permutation invariant on the exterior and attains the bound on every cap.
    """
    masses = _reference_sector_masses()
    rows = []
    for k in range(1, 11):
        target = tuple(sum(masses[q]*_hypergeom(2-q, k, r) for q in range(3))
                       for r in range(min(k, 2)+1))
        if sum(target) != 1:
            raise ArithmeticError("reference cap distribution is not normalized")
        distances = tuple(sum(abs(_hypergeom(2-q, k, r)-target[r])
                              for r in range(len(target)))/2 for q in range(3))
        rows.append({"cap_size": k,
                     "reference_exterior_number_distribution": [str(v) for v in target],
                     "sharp_q_minimax_trace_distances": [str(v) for v in distances]})
    if tuple(F(v) for v in rows[-1]["sharp_q_minimax_trace_distances"]) != (F(3, 5), F(4, 9), F(43, 45)):
        raise ArithmeticError("full exterior result disagrees with prior sector certificate")
    return {
        "exterior_ports": 10,
        "reference_sector_probabilities": [str(v) for v in masses],
        "minimax_domain": "all density matrices supported in one fixed-core sharp q sector",
        "cap_menu": "all ten-choose-k exterior subsets for each k=1..10",
        "distance": "half trace norm of the cap-restricted state and transported diagonal reference",
        "lower_witness": "cap particle-number measurement averaged over every cap",
        "attainer": "diagonal reference conditioned on q; uniform within its q sector",
        "single_specified_cap_has_this_lower_bound": False,
        "rows": rows,
    }


def code_geometry_certificate():
    """Specialize the prior all-seam nonedge theorem to q=1 and check maps.

    For a k=2 exact flagged code, each native-edge occupation projector is
    rank one and has zero compression. Each n_j also compresses to a scalar.
    If only one of the core-j source pairs is a nonedge, n_j has rank one on
    the remaining support, hence its scalar compression is zero. Therefore
    only ports nonadjacent to *both* core ports can have exterior occupation.
    """
    if len(SEAMS) != 30:
        raise ArithmeticError("the native seam menu changed")
    classification = Counter()
    all_branches = adaptive.branch_maps(binary=True)
    if len(all_branches) != 720:
        raise ArithmeticError("complete binary raw C/E branch menu changed")
    neighbors = {port: {v if u == port else u for u, v in SEAMS if port in (u, v)}
                 for port in range(12)}
    if {len(value) for value in neighbors.values()} != {5}:
        raise ArithmeticError("native seam graph degree changed")
    by_law = {law: tuple(b for b in all_branches if adaptive.actions()[b.action][0] == law)
              for law in ("C", "E")}
    branch_checks = Counter()
    for center in CORES:
        safe = safe_exterior(center)
        seen, frontier, distance = {center[0]}, {center[0]}, 0
        while center[1] not in seen and frontier:
            frontier = set().union(*(neighbors[p] for p in frontier))-seen
            seen |= frontier
            distance += 1
        if distance != {4: 1, 2: 2, 0: 3}.get(len(safe)):
            raise ArithmeticError("native graph distance does not match safe cap class")
        classification[("edge" if center in SEAMS else "nonedge", len(safe))] += 1
        columns = sharp_code(center)
        if not safe:
            if columns is not None:
                raise ArithmeticError("a q=1 code was asserted with no safe exterior port")
            continue
        if len(safe) not in (2, 4) or columns is None:
            raise ArithmeticError("unexpected native seam/core geometry")
        term_count = len(columns[0])
        if len(columns[1]) != term_count or set(columns[0]) & set(columns[1]):
            raise ArithmeticError("declared q=1 columns are not orthogonal")
        if any(sector.charge(tuple(int(j in pair) for j in range(12)), center) != 1
               or pair in SEAMS for column in columns for pair in column):
            raise ArithmeticError("declared code leaves q=1 or uses a native edge")
        indices = tuple(frozenset(sector.PAIRS.index(pair) for pair in column)
                        for column in columns)
        for law, branches in by_law.items():
            seam_totals = Counter()
            for branch in branches:
                action_law, seam = adaptive.actions()[branch.action]
                if action_law != law:
                    raise ArithmeticError("branch law disagrees with action label")
                inputs = set(branch.inputs)
                counts = tuple(len(row & inputs) for row in indices)
                if counts[0] != counts[1]:
                    raise ArithmeticError("record effect distinguishes the code columns")
                if len(set(branch.outputs)) != len(branch.outputs):
                    raise ArithmeticError("record branch is not injective")
                weight = F(counts[0], 2*term_count)
                seam_totals[seam] += weight
                branch_checks[(law, len(safe), "all")] += 1
                if weight:
                    branch_checks[(law, len(safe), "positive")] += 1
            if len(seam_totals) != 30 or set(seam_totals.values()) != {F(1)}:
                raise ArithmeticError("complete flagged code probabilities do not sum to one")
    expected_classes = {("edge", 4): 30, ("nonedge", 2): 30, ("nonedge", 0): 6}
    if dict(classification) != expected_classes:
        raise ArithmeticError("safe exterior classification changed")
    example = sharp_code(CORE)
    assert example is not None
    ref_occupation = sum(_reference_sector_masses()[q]*F(2-q, 10) for q in range(3))
    if ref_occupation != F(61, 450):
        raise ArithmeticError("one-port transported reference changed")
    old_code_worst = F(1, 2)-ref_occupation
    edge_code_worst = max(ref_occupation, abs(F(1, 4)-ref_occupation))
    nonedge_code_worst = max(ref_occupation, abs(F(1, 2)-ref_occupation))
    if (old_code_worst, edge_code_worst, nonedge_code_worst) != (F(82, 225), F(61, 450), F(82, 225)):
        raise ArithmeticError("one-port screening distances changed")
    return {
        "fixed_core_example": list(CORE),
        "native_seam_degree": 5,
        "core_classes": [
            {"core_type": kind, "safe_exterior_ports": count, "cores": number}
            for (kind, count), number in sorted(classification.items())],
        "edge_core_safe_ports_example": list(safe_exterior(CORE)),
        "edge_core_code_columns": [[list(pair) for pair in column] for column in example],
        "branch_checks": {f"{law}_safe{n}_{scope}": count
                          for (law, n, scope), count in sorted(branch_checks.items())},
        "all_raw_C_E_branch_effects_are_scalar_on_each_witness": True,
        "branch_image_injection": True,
        "arbitrary_reference_exact_recovery_exists_per_retained_flag": True,
        "source_encoder_and_recovery_compiled_for_this_new_code": False,
        "maximum_safe_exterior_ports_any_core": 4,
        "reference_single_port_occupation": str(ref_occupation),
        "unrestricted_sharp_q1_minimax_worst_single_port_distance": "8/225",
        "prior_two_port_code_worst_single_port_distance": str(old_code_worst),
        "edge_core_qubit_optimal_worst_single_port_distance": str(edge_code_worst),
        "nonedge_two_safe_ports_qubit_optimal_worst_single_port_distance": str(nonedge_code_worst),
        "opposite_zero_safe_ports_exact_q1_qubit_exists": False,
        "optimality_scope": "fixed-core q=1 logical dimension two, all native seams, retained C/E flags, transported diagonal reference, maximum over all single exterior ports",
        "physical_particle_identified": False,
    }


def transport_certificate():
    """One-step co-moving charge and graph-adjacency hazards; no clock claim."""
    if CORE not in SEAMS:
        raise ArithmeticError("initial code core is no longer a seam")
    columns = sharp_code(CORE)
    assert columns is not None
    sharp_after_recovery = [center for center in CORES
                            if all(len(set(pair) & set(center)) == 1
                                   for column in columns for pair in column)]
    if sharp_after_recovery != [CORE]:
        raise ArithmeticError("fixed recovered code has another sharp-q1 core label")
    code_rows = {}
    for law in ("C", "E"):
        logical_rows = []
        for column in columns:
            changed = broken = conserved = F()
            for pair in column:
                x = tuple(int(p in pair) for p in range(12))
                for seam in SEAMS:
                    for coin in (0, 1):
                        y, record = records.recorded_step(law, x, seam, coin)
                        moved = sector.moved_core(CORE, law, record)
                        weight = F(1, len(column)*len(SEAMS)*2)
                        changed += weight*int(moved != CORE)
                        broken += weight*int(moved not in SEAMS)
                        conserved += weight*int(sector.charge(y, moved) == 1)
            logical_rows.append((changed, broken, conserved))
        if len(set(logical_rows)) != 1 or logical_rows[0] != (F(1, 15), F(1, 30), F(1)):
            raise ArithmeticError("screened code co-moving transport changed")
        code_rows[law] = {
            "core_label_changes_per_uniform_native_attempt": str(logical_rows[0][0]),
            "core_leaves_native_edge_per_uniform_native_attempt": str(logical_rows[0][1]),
            "sharp_q1_after_co_moving_native_step": str(logical_rows[0][2]),
            "additional_7_bit_core_reset_needed_after_active_source_recovery_when_label_changes": True,
            "core_reset_archive_entry_bits_if_status_old_new_are_retained": 15,
            "new_code_source_pulse_count_known": False,
        }
    edge_rows = {}
    for law in ("C", "E"):
        values = set()
        for center in SEAMS:
            x = tuple(int(p in center) for p in range(12))
            broken = changed = conserved = F()
            for seam in SEAMS:
                for coin in (0, 1):
                    y, record = records.recorded_step(law, x, seam, coin)
                    moved = sector.moved_core(center, law, record)
                    pair = tuple(p for p, n in enumerate(y) if n)
                    weight = F(1, 2*len(SEAMS))
                    broken += weight*int(pair not in SEAMS)
                    changed += weight*int(pair != center)
                    conserved += weight*int(pair == moved)
            values.add((broken, changed, conserved))
        if values != {(F(1, 15), F(2, 15), F(1))}:
            raise ArithmeticError("q=2 edge pair is not uniformly unbound")
        edge_rows[law] = {
            "edge_pair_first_exit_hazard_per_uniform_native_attempt": "1/15",
            "edge_pair_mean_first_exit_attempts": "15",
            "source_pair_equals_co_moving_core_after_every_step": True,
            "q2_relational_charge_preserved_when_edge_bond_breaks": True,
        }
    return {
        "parameter": "uniform native attempted repair count, not a calibrated observer clock",
        "initial_core": list(CORE),
        "only_core_label_for_which_recovered_fixed_code_is_sharp_q1": list(CORE),
        "screened_active_q1_code": code_rows,
        "uncontrolled_q2_edge_pair": edge_rows,
        "q2_first_exit_survival_formula": "(14/15)^t for integer t>=0, until first nonedge pair",
        "q2_first_exit_is_not_particle_decay": True,
    }


def passive_memory_certificate():
    """Reuse the all-history finite theorem to bound passive qubit memory."""
    rows = {}
    for law in ("C", "E"):
        for attempts in range(1, 100):
            result = adaptive.bounds(F(1), attempts, law_scope=law, binary=True,
                                     logical_dimension=2)
            if result["encoded_final_recovery_trace_error_lower"] > 0:
                rows[law] = {
                    "first_attempt_count_with_strictly_positive_final_decoder_error_lower": attempts,
                    "error_lower_at_that_count": str(result["encoded_final_recovery_trace_error_lower"]),
                    "one_step_pair_agreement_upper": str(result["one_step_agreement_upper"]),
                }
                break
        if law not in rows:
            raise ArithmeticError("passive memory bound failed to become nonzero")
    if (rows["C"]["first_attempt_count_with_strictly_positive_final_decoder_error_lower"],
            rows["E"]["first_attempt_count_with_strictly_positive_final_decoder_error_lower"]) != (23, 12):
        raise ArithmeticError("passive information horizon changed")
    return {
        "scope": "any binary66 encoded qubit, complete C or E flags, uniform seam per attempt, no intermediate quantum control, arbitrary final decoder",
        "parameter": "native attempted repair count, not a physical lifetime",
        "source_theorem": "model.observer.adaptive_record_recovery.bounds with eta=1, binary=True, logical_dimension=2",
        "rows": rows,
    }


def science_result():
    return {"schema": "oph.particle_screening.science.v1",
            "scope": "FINITE_COMPONENT_DIAGNOSTICS",
            "cap_minimax": cap_certificate(),
            "all_seam_code_geometry": code_geometry_certificate(),
            "transport": transport_certificate(),
            "passive_memory": passive_memory_certificate(),
            "cp03_accepted": False, "particle_species_identified": False,
            "rest_energy_identified": False, "physical_lifetime_identified": False}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def campaign_binding(science):
    campaign = strict_json.load(HERE/"campaign.json")
    expected = {"schema": "oph.particle_screening.campaign.v1",
                "campaign_id": "particle-screening-20261028-v1",
                "predecessor_git_commit": "7ab50c9fa326b86344ba0dda6ca56804eb094fb7",
                "scope": "FINITE_COMPONENT_DIAGNOSTICS",
                "whole_program_frozen": False,
                "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
                "integrated_gate_accepted": False,
                "native_family_changed": False,
                "physical_comparison_count_added": 0}
    for key, value in expected.items():
        if type(campaign[key]) is not type(value) or campaign[key] != value:
            raise ValueError("particle-screening campaign scope changed")
    source = campaign["predecessor_source_sha256"]
    if type(source) is not dict or set(source) != set(PREDECESSOR_FILES):
        raise ValueError("particle-screening campaign predecessor set changed")
    for path, digest in source.items():
        if type(digest) is not str or _sha(path) != digest:
            raise ValueError("particle-screening campaign predecessor source changed")
    claim_hash = hashlib.sha256(_canonical(science)).hexdigest()
    if campaign["science_result_sha256"] != claim_hash:
        raise ValueError("particle-screening campaign science result changed")
    return {"campaign_id": campaign["campaign_id"],
            "campaign_sha256": _sha("model/particle_screening/campaign.json"),
            "science_result_sha256": claim_hash}


def source_bindings():
    """Reuse the older packet's pinned import closure, then bind new texts."""
    bound = sector.source_bindings()
    for path in ("model/particle_screening/__init__.py",
                 "model/particle_screening/screening.py",
                 "model/particle_screening/test_screening.py",
                 "model/particle_screening/README.md",
                 "model/particle_screening/campaign.json"):
        bound[path] = _sha(path)
    return dict(sorted(bound.items()))


def run():
    science = science_result()
    return {**science, **campaign_binding(science), "source_sha256": source_bindings()}


def check_receipt(path, expected):
    actual = strict_json.load(path)
    if _canonical(actual) != _canonical(expected):
        raise ValueError("particle-screening receipt differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    options.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    elif args.check:
        check_receipt(args.check, result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
