"""Exact particle-candidate boundaries and an executed active-memory control.

The conserved source/core overlap is a relational charge. Its sharp sectors
are not identified with localized particles. A closed odd exchange in an
uninterrupted native history has a rank-one input fiber. Active quantum
recovery is a different, explicitly accounted grammar.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..observer import native, record_recovery as records
from ..observer import adaptive_record_recovery as adaptive
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..native_modular_reference import history_transport
from ..common_source import encoded_memory as memory, native_routing as routing

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "sector_results.json"
BASIS = native.states(binary=True)
PAIRS = tuple(tuple(i for i, n in enumerate(x) if n) for x in BASIS)
CORES = tuple(combinations(range(12), 2))
CORE = (0, 1)
SEAM = (0, 5)
EXCHANGE_WORD = ((0, 5), (0, 11), (5, 11))
WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))
TOLERANCE = 2e-9
PREDECESSOR_SOURCES = frozenset((
    "model/observer/native.py", "model/observer/record_recovery.py",
    "model/observer/adaptive_record_recovery.py",
    "model/native_modular_reference/history_transport.py",
    "model/clock_defect_energy/continuation.py", "model/clock_defect_energy/defect_history.py",
    "model/common_source/encoded_memory.py", "model/common_source/record_code.py",
    "model/common_source/native_routing.py", "model/observer/source_discrimination/circuit.py",
    "model/observer/source_scalar/compiler.py"))


def core(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(p) is not int for p in value) or value not in CORES):
        raise ValueError("an original sorted pair of distinct ports required")
    return value


def charge(state, center=CORE):
    center = core(center)
    if type(state) is not tuple or state not in BASIS or any(type(n) is not int for n in state):
        raise ValueError("an original binary total-two state required")
    return sum(state[p] for p in center)


def reference(center=CORE):
    center = core(center)
    return tuple(WEIGHTS[charge(x, center)] for x in BASIS)


def moved_core(center, law, record):
    center = core(center)
    current = records.record_current(law, record)
    if current not in (-1, 0, 1):
        raise ValueError("binary native current required")
    if not current:
        return center
    a, b = record[:2]
    return tuple(sorted(b if p == a else a if p == b else p for p in center))


def sector_certificate():
    """Exact full-exterior bounds, attained by conditioned diagonal references.

    The algebra contains N_out=sum_{p outside c} n_p. A state supported in
    q has N_out=2-q with certainty. The reference gives that event mass p_q.
    This proves D>=1-p_q for EVERY such quantum state. The conditioned
    reference attains it even on the full source, hence also on the exterior.
    """
    rows = []
    for q in range(3):
        counts = [sum(charge(x, c) == q for x in BASIS) for c in CORES]
        if len(set(counts)) != 1:
            raise ArithmeticError("charge multiplicities depend on core")
        count = counts[0]
        p = count * WEIGHTS[q]
        for c in CORES:
            sigma = reference(c)
            conditioned = tuple(w / p if charge(x, c) == q else F()
                                for x, w in zip(BASIS, sigma))
            distance = sum(abs(a-b) for a, b in zip(conditioned, sigma)) / 2
            if sum(conditioned) != 1 or distance != 1-p:
                raise ArithmeticError("conditioned reference does not attain the lower bound")
        # Equal weight over all cores, conditioned inside each core sector.
        # This is a finite ensemble identity, not a preparation service.
        source_marginal = tuple(sum(F(charge(x, c) == q, 66*count) for c in CORES) for x in BASIS)
        if source_marginal != (F(1, 66),)*66:
            raise ArithmeticError("sector source marginal is not uniform")
        rows.append({"q": q, "fixed_core_dimension": count,
                     "source_core_dimension": 66*count,
                     "reference_sector_mass": str(p),
                     "full_exterior_distance_minimum": str(1-p),
                     "uniform_reference_exterior_distance_minimum": str(1-F(count, 66)),
                     "outside_particle_number": 2-q,
                     "attained_by_conditioned_reference": True,
                     "all_source_marginal_entries": "1/66"})
    return {"sectors": rows, "cores_checked": 66,
            "source_only_channel_distinguishes_sector_ensembles": False,
            "source_only_statement_requires_core_independent_controls": True,
            "core_plus_destructive_source_basis_read_distinguishes_sharp_q": True,
            "exterior_algebra": "all ten complementary ports, including total exterior number",
            "every_smaller_cap_has_this_distance_claimed": False}


@lru_cache(maxsize=1)
def common_record_certificate():
    """Check the local induction premises against every actual raw C/E map."""
    result = {}
    for law in ("C", "E"):
        pair_checks = moves = 0
        branches = [b for b in adaptive.branch_maps(binary=True)
                    if adaptive.actions()[b.action][0] == law]
        for b in branches:
            _, seam = adaptive.actions()[b.action]
            for (i, j), (k, ell) in combinations(zip(b.inputs, b.outputs), 2):
                before = tuple(a-b for a, b in zip(BASIS[i], BASIS[k]))
                after = tuple(a-b for a, b in zip(BASIS[j], BASIS[ell]))
                if before != after:
                    raise ArithmeticError("common native record changed input difference")
                frozen = {p for p, n in enumerate(before) if n}
                if any(BASIS[i][p] != BASIS[j][p] or BASIS[k][p] != BASIS[ell][p] for p in frozen):
                    raise ArithmeticError("distinct spectator moved under a common record")
                if i != j:
                    moves += 1
                    removed = set(PAIRS[i])-set(PAIRS[j])
                    added = set(PAIRS[j])-set(PAIRS[i])
                    if (len(removed) != 1 or len(added) != 1
                            or not removed <= set(PAIRS[i]) & set(PAIRS[k])
                            or added & (set(PAIRS[i]) | set(PAIRS[k]))):
                        raise ArithmeticError("shared movement is not a shared particle into a common vacancy")
                pair_checks += 1
        result[law] = {"complete_raw_maps": len(branches),
                       "distinct_common_input_pair_checks": pair_checks,
                       "shared_movement_pair_checks": moves}
    return result


def odd_closed_inputs(law, branch):
    """Replay the whole native fiber; a caller cannot delete compatible inputs."""
    branch = history_transport._branch(law, branch)
    odd = []
    for i, j in zip(branch.inputs, branch.outputs):
        ordered = PAIRS[i]
        for record in branch.history:
            current = records.record_current(law, record)
            if current:
                a, b = record[:2] if current == 1 else record[:2][::-1]
                if a not in ordered or b in ordered:
                    raise ArithmeticError("native record is inconsistent with ordered particle lift")
                ordered = tuple(b if p == a else p for p in ordered)
        if tuple(sorted(ordered)) != PAIRS[j]:
            raise ArithmeticError("ordered endpoint disagrees with native branch")
        if i == j and ordered == PAIRS[i][::-1]:
            odd.append(i)
    if odd and len(branch.inputs) != 1:
        raise ArithmeticError("a closed odd exchange has a non-singleton native input fiber")
    return tuple(odd)


def exchange_certificate():
    """General local proof premises plus one complete, positive native example."""
    local = common_record_certificate()
    rows = {}
    for law in ("C", "E"):
        branches = history_transport.word_instrument(law, EXCHANGE_WORD)
        odd = [b for b in branches if odd_closed_inputs(law, b)]
        if len(odd) != 1 or odd[0].inputs != odd[0].outputs or odd[0].weight != F(1, 8):
            raise ArithmeticError("the three-attempt odd exchange witness changed")
        b = odd[0]
        rows[law] = {"nonzero_complete_blocks": len(branches),
                     "complete_raw_slots": (2*(1 << circuit.RECORD_BITS[law]))**3,
                     "odd_flag_history": b.history,
                     "odd_flag_input_pair": PAIRS[b.inputs[0]],
                     "odd_flag_domain_rank": len(b.inputs),
                     "odd_flag_probability_on_its_input": str(b.weight),
                     "odd_flag_probability_on_uniform_source": str(b.weight/66),
                     "odd_flag_is_entanglement_breaking": True}
    return {"all_finite_native_histories_theorem":
            "Any positive closed odd-exchange flag on binary total-two has rank-one source domain.",
            "proof_local_premises": local, "three_attempt_complete_example": rows,
            "final_flag_aware_recovery_can_preserve_dimension_gt_one_on_positive_odd_flag": False,
            "intermediate_quantum_recovery_excluded_by_this_theorem": False,
            "scalar_phase_on_committed_flag_changes_instrument": False}


def code_capacity_certificate():
    """Specialize the existing all-seam flagged-KL occupation-scalar proof.

    On q=1 each exterior n_j has rank2. Sum_out n_j=I, so some scalar
    compression V*n_j V=p_j I is nonzero; its rank gives k<=2.
    Existing two-column code attains equality and has scalar complete effects.
    """
    supports = (((0, 2), (1, 3)), ((0, 3), (1, 2)))
    if any(sum(p in CORE for p in pair) != 1 for col in supports for pair in col):
        raise ArithmeticError("existing code is outside sharp q=1")
    exterior_ranks = {str(p): sum(charge(x) == 1 and x[p] for x in BASIS)
                      for p in range(12) if p not in CORE}
    if set(exterior_ranks.values()) != {2}:
        raise ArithmeticError("q=1 exterior occupation rank changed")
    # Instantiate the pre-existing exact code with the actual complete raw
    # C/E maps, not just the source-executed demonstration seam below.
    indices = tuple(tuple(BASIS.index(tuple(int(p in pair) for p in range(12)))
                          for pair in col) for col in supports)
    laws = {}
    for law in ("C", "E"):
        sums = Counter()
        checked = 0
        directed_effects = set()
        moving_effect = [F() for _ in BASIS]
        for b in adaptive.branch_maps(binary=True):
            actual_law, seam = adaptive.actions()[b.action]
            if actual_law != law:
                continue
            overlaps = [len(set(col) & set(b.inputs)) for col in indices]
            if overlaps[0] != overlaps[1]:
                raise ArithmeticError("native flag reveals the logical code value")
            if len(set(b.outputs)) != len(b.outputs):
                raise ArithmeticError("native branch fails coherent injection")
            value = circuit.record_values(law, b.coin).get(b.raw)
            if value is not None:
                current = records.record_current(law, (*seam, b.coin, value))
                if current and b.inputs:
                    a, c = seam if current == 1 else seam[::-1]
                    if set(b.inputs) != {i for i, x in enumerate(BASIS) if x[a] == 1 and x[c] == 0}:
                        raise ArithmeticError("moving raw effect is not the claimed directed occupation projector")
                    directed_effects.add((a, c))
                    for i in b.inputs:
                        moving_effect[i] += F(1, 60)
            sums[seam] += F(overlaps[0], 4)
            checked += 1
        if len(sums) != 30 or set(sums.values()) != {F(1)} or len(directed_effects) != 60:
            raise ArithmeticError("all-seam corrected code dropped native probability")
        laws[law] = {"raw_effects_checked": checked, "seams": len(sums),
                     "directed_occupation_effects": len(directed_effects),
                     "complete_branch_compression": "p_h I_2, p_h in {0,1/4,1/2}"}
        code_activity = [sum(moving_effect[i] for i in col)/2 for col in indices]
        background = sum(moving_effect)/66
        reference_background = sum(w*a for w, a in zip(reference(), moving_effect))
        if code_activity != [F(1, 6)]*2 or background != F(5, 33):
            raise ArithmeticError("complete native activity differs from independent support count")
        laws[law]["uniform_seam_activity"] = {
            "code_any_logical_input": "1/6", "uniform_source_background": str(background),
            "excess_over_uniform_source_per_attempt": str(F(1, 6)-background),
            "ideal_fixed_core_reference_background": str(reference_background),
            "excess_over_fixed_core_reference_per_attempt": str(F(1, 6)-reference_background),
            "event": "nonzero signed native current; all coin and raw flags retained",
            "coherent_recovery_work_included_in_this_activity": False,
            "rest_energy_inferred": False}
    return {"fixed_core": CORE, "sharp_charge": 1,
            "optimal_exact_all_seam_correctable_dimension": 2,
            "exterior_occupation_ranks": exterior_ranks,
            "achieving_existing_code_supports": supports, "complete_C_E_checks": laws,
            "active_quantum_recovery_required": True,
            "unknown_logical_state_given_as_input": True,
            "autonomous_localized_particle_constructed": False}


def terminal_rows(amplitude, *, setting="boundary"):
    """Complete destructive source read; every row keeps its reference vector.

    Boundary setting: existing source computational read, no rotation.
    Logical setting: existing full inverse encoder and one Y(pi/4) source
    pulse, then the SAME 512-output computational read. These are separate
    terminal experiments; a boundary read is not a coherent QND service.
    """
    value = routing._state(amplitude)
    if type(setting) is not str or setting not in ("boundary", "logical"):
        raise ValueError("declared terminal boundary or logical setting required")
    if setting == "logical":
        value = compiler.execute([compiler.pulse("Y", math.pi/4, [0])], 9,
                                 memory.decoding(value))
    return tuple((i, routing._immutable(value[i:i+1])) for i in range(512))


def charge_read_masses(amplitude, center=CORE):
    center = core(center)
    labels = {circuit.encode(x): charge(x, center) for x in BASIS}
    masses = {str(q): 0. for q in range(3)}
    masses["outside_binary_shell"] = 0.
    for address, row in terminal_rows(amplitude):
        label = str(labels[address]) if address in labels else "outside_binary_shell"
        masses[label] += float(np.vdot(row, row).real)
    return masses


def execution_certificate():
    """Actual source encoder, native repair, recovery and two terminal settings.

    The two columns are an arbitrary-reference isometry test: equality on
    both columns proves the ideal map for any logical/reference input. Every
    floating raw branch/row is retained. Tolerance checks are diagnostics,
    not physical service error certificates.
    """
    code = memory.code_matrix()
    ideal = memory.ideal_code_matrix()
    if np.max(np.abs(code-ideal)) > TOLERANCE:
        raise ArithmeticError("source encoder differs from its declared code")
    rows = {}
    plus = np.array([1., 1.])/math.sqrt(2)
    minus = np.array([1., -1.])/math.sqrt(2)
    for law in ("C", "E"):
        native_blocks = memory.repair_branches(law, SEAM, code, recover=False)
        before, after_without_reset, after_reset = 0., 0., 0.
        unreset_charge_effect = np.zeros((2, 2), complex)
        reset_probability = F()
        max_recovery = 0
        logical_plus = {"plus_input": 0., "minus_input": 0.}
        gram = np.zeros((2, 2), complex)
        reference_gram = np.zeros((2, 2), complex)
        residual = 0.
        for b in native_blocks:
            coin, raw = b["coin"], b["record_code"]
            value = circuit.record_values(law, coin).get(raw)
            target_core = CORE if value is None else moved_core(CORE, law, (*SEAM, coin, value))
            raw_map = b["amplitude"]
            before += charge_read_masses(raw_map, target_core)["1"]/2
            recovered = raw_map if b["abort"] else memory.recovery_unitary(law, SEAM, coin, raw)@raw_map
            probability = memory.branch_weight(law, SEAM, coin, raw)
            if not b["abort"] and target_core != CORE:
                reset_probability += probability
            expected = math.sqrt(float(probability))*ideal
            residual = max(residual, float(np.max(np.abs(recovered-expected))))
            gram += recovered.conj().T@recovered
            reference_gram += sum((row.conj().T@row for _, row in terminal_rows(recovered)),
                                  np.zeros((2, 2), complex))
            after_without_reset += charge_read_masses(recovered, target_core)["1"]/2
            q1_rows = [circuit.encode(x) for x in BASIS if charge(x, target_core) == 1]
            unreset_charge_effect += recovered[q1_rows].conj().T@recovered[q1_rows]
            after_reset += charge_read_masses(recovered, CORE)["1"]/2
            if not b["abort"]:
                max_recovery = max(max_recovery, memory.recovery_resources(law, SEAM, coin, raw)["native_pulses"])
                for name, vector in (("plus_input", plus), ("minus_input", minus)):
                    # All other outcomes still exist; only this named event is
                    # summarized, and completeness is checked separately.
                    output = terminal_rows(recovered@vector, setting="logical")
                    logical_plus[name] += float(np.vdot(output[0][1], output[0][1]).real)
        if (residual > TOLERANCE or np.max(np.abs(gram-np.eye(2))) > TOLERANCE
                or np.max(np.abs(reference_gram-np.eye(2))) > TOLERANCE
                or np.max(np.abs(unreset_charge_effect-(7/8)*np.eye(2))) > TOLERANCE
                or max(abs(before-1), abs(after_reset-1), abs(after_without_reset-7/8)) > TOLERANCE
                or abs(logical_plus["plus_input"]-1) > TOLERANCE
                or logical_plus["minus_input"] > TOLERANCE or reset_probability != F(1, 4)):
            raise ArithmeticError("source execution fails the exact sector/recovery/readout control")
        repair_cost = routing.routed_repair_resources(law, SEAM)["native_pulses"]
        maximum = memory.encoding_resources()["native_pulses"] + repair_cost + max_recovery
        rows[law] = {"complete_raw_branches_retained": len(native_blocks),
                     "native_q1_probability": "1",
                     "after_recovery_without_core_reset_q1_probability": "7/8",
                     "after_recovery_and_explicit_core_reset_q1_probability": "1",
                     "probability_core_value_changes_in_reset": "1/4",
                     "logical_plus_read_probabilities": {"plus_input": "1", "minus_input": "0"},
                     "source_execution_matches_exact_maps_within": TOLERANCE,
                     "retained_native_rows_per_raw_branch": 512,
                     "maximum_native_pulses_through_recovery": maximum,
                     "maximum_native_pulses_with_logical_terminal_read":
                         maximum+memory.decoding_resources()["native_pulses"]+1,
                     "maximum_recovery_native_pulses": max_recovery,
                     "terminal_computational_outcomes": 512,
                     "raw_abort_complement_retained": True}
    return {"initial_core": CORE, "source_executed_seam": SEAM, "laws": rows,
            "services": {"unknown_logical_input_qubits": 1, "initial_blank_code_qubits": 8,
                         "source_code_qubits": 9, "fresh_quantum_bits_per_recovery": 0,
                         "coin_initializations_per_attempt": 1,
                         "native_record_initializations_per_attempt": 1,
                         "terminal_read_bits": 9,
                         "core_label_bits": 7,
                         "extra_core_reset_archive_entries_per_nonabort_attempt": 1,
                         "extra_core_reset_archive_bits_per_entry": 15,
                         "reset_entry_format": "status bit, old 7-bit core, new 7-bit core",
                         "native_records_rewritten": 0,
                         "proper_code_load_and_export_events_per_source_pulse": 2,
                         "source_computational_read_is_declared_service": True},
            "arbitrary_reference_preserved_before_destructive_read": True,
            "full_computational_read_preserves_reference_marginal_after_outcomes_summed": True,
            "logical_and_boundary_settings_are_separate_experiments": True,
            "hardware_service_accuracy_certified": False}


def source_bindings():
    """Bind the static repository Python import closure, plus theorem texts.

    Imports inside function bodies are included. External Python packages
    are execution dependencies, not repository source files. Literal data
    receipts are not used as computational inputs by this packet.
    """
    seeds = ("model/exchange_sector/sector.py", "model/exchange_sector/test_sector.py")
    pending, seen = list(seeds), set()
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        path = ROOT/name
        parts = name[:-3].split("/")
        package = parts[:-1]
        for node in ast.walk(ast.parse(path.read_text())):
            candidates = []
            if isinstance(node, ast.Import):
                candidates.extend(a.name.split(".") for a in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = package[:len(package)-node.level+1] if node.level else []
                base = base + (node.module.split(".") if node.module else [])
                candidates.append(base)
                candidates.extend(base+a.name.split(".") for a in node.names if a.name != "*")
            for candidate in candidates:
                if not candidate or candidate[0] not in ("model", "tools"):
                    continue
                for module in ("/".join(candidate)+".py", "/".join(candidate)+"/__init__.py"):
                    if (ROOT/module).is_file() and module not in seen:
                        pending.append(module)
    seen.update(("model/exchange_sector/README.md", "model/exchange_sector/campaign.json",
                 "model/common_source/ENCODED_MEMORY.md",
                 "model/common_source/record_code.py", "model/common_source/record_information.py",
                 "model/clock_defect_energy/defect_history.py", "model/clock_defect_energy/DEFECT_HISTORY.md",
                 "model/clock_defect_energy/continuation.py",
                 "Lean/Geometry/SourceEncodedMemory.lean"))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(seen)}


def campaign_binding():
    value = strict_json.load(HERE/"campaign.json")
    expected = {"campaign_id": "exchange-sector-20261028-v1",
                "schema": "oph.exchange_sector.campaign.v1",
                "predecessor_git_commit": "82a8d90588d819c432276136489121dbf8789072",
                "scope": "FINITE_COMPONENT_DIAGNOSTICS",
                "boundary_status": "FROZEN_DIAGNOSTIC_CONTRACT",
                "whole_program_frozen": False, "integrated_gate_accepted": False,
                "native_family_changed": False, "physical_comparison_count_added": 0}
    for key, item in expected.items():
        if type(value[key]) is not type(item) or value[key] != item:
            raise ValueError("exchange-sector campaign scope changed")
    sources = value["source_sha256"]
    if type(sources) is not dict or set(sources) != PREDECESSOR_SOURCES:
        raise ValueError("source-bound exchange-sector campaign required")
    for name, digest in sources.items():
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("exchange-sector predecessor source changed")
    return {"campaign_id": value["campaign_id"],
            "contract_sha256": hashlib.sha256((HERE/"campaign.json").read_bytes()).hexdigest()}


def run():
    return {"schema": "oph.exchange_sector.v1", "scope": "FINITE_COMPONENT_DIAGNOSTICS",
            **campaign_binding(),
            "predecessor_git_commit": "82a8d90588d819c432276136489121dbf8789072",
            "native_family_changed": False, "physical_comparisons_added": 0,
            "cp03_accepted": False, "cp05_accepted": False,
            "sector_boundaries": sector_certificate(), "closed_exchange": exchange_certificate(),
            "optimal_sharp_sector_code": code_capacity_certificate(),
            "source_execution": execution_certificate(), "source_sha256": source_bindings()}


def check_receipt(path, expected):
    actual = strict_json.load(path)
    canonical = lambda value: json.dumps(value, sort_keys=True, allow_nan=False)
    if canonical(actual) != canonical(expected):
        raise ValueError("exchange-sector receipt differs from exact replay")


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
