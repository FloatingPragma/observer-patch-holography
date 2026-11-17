"""Continue the actual C/E FIFO with a record-updated transported core.

Only binary66 is admitted.  The local prepared reference depends on the
charge in an unordered two-port core.  Each nonzero committed current swaps
the core endpoints before its record can leave the FIFO.  This implements
the reference *metadata* update of history_transport, not a new source gate.
All native records, selectors, law tags, holds and exported records remain
in the complete kernel; the long runner computes an explicitly named
diagonal marginal.  A missing predecessor clock word is a service status,
never an executed clock or a discarded repair outcome.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from math import lcm
from pathlib import Path

from model.completion import observer
from model.completion.family import FamilyParameters
from model.native_modular_reference import history_clock, history_transport
from model.observer import native
from .contract import binding, check_receipt, source_bindings

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "continuation_results.json"
BASIS = native.states(binary=True)
CORES = tuple(combinations(range(12), 2))
CORE_INDEX = {core: i for i, core in enumerate(CORES)}
INITIAL_CORE = (0, 1)
FULL_TO_BINARY = {observer.STATES.index(x): i for i, x in enumerate(BASIS)}
LOCAL_WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))


def _core(core):
    if (type(core) is not tuple or len(core) != 2
            or any(type(p) is not int for p in core) or core not in CORE_INDEX):
        raise ValueError("an ordered encoding of an unordered two-port core required")
    return core


def reference(core=INITIAL_CORE):
    """Faithful extension, not the history posterior or a new preparation."""
    core = _core(core)
    return tuple(LOCAL_WEIGHTS[x[core[0]]+x[core[1]]] for x in BASIS)


def update_core(core, record):
    core = _core(core)
    if record.kind == "timeout":
        if record != observer.Record("timeout"):
            raise ValueError("invalid timeout")
        return core
    if not observer._valid_record(record):
        raise ValueError("complete native record required")
    if not record.current:
        return core
    a, b = record.seam
    return tuple(sorted(b if p == a else a if p == b else p for p in core))


@lru_cache(maxsize=1)
def compiled_cores():
    """Inventory only: the four source words certified by the predecessor."""
    result = {INITIAL_CORE}
    for law in ("C", "E"):
        for branch in history_clock.branches(law):
            permutation = history_transport.port_permutation(law, branch)
            result.add(tuple(sorted(permutation[p] for p in INITIAL_CORE)))
    return frozenset(result)


def clock_service(core):
    """Inventory at predecessor parameter4; identity source action in both cases."""
    return ("predecessor_word_available" if _core(core) in compiled_cores()
            else "unavailable_service")


@dataclass(frozen=True)
class Branch:
    weight: F
    pairs: tuple[tuple[int, int], ...]
    record: observer.Record
    memory: tuple[observer.Record, ...]
    exported: observer.Record | None
    halted: bool
    core: tuple[int, int]
    clock_status: str


def branches(memory, core, params):
    """Complete binary restriction of the existing observer, including zeros.

    The new core is computed from the new current before the caller releases
    the old FIFO slot. No quantum source value controls this metadata update.
    """
    _core(core)
    result = []
    for old in observer.branches(memory, params):
        pairs = tuple((FULL_TO_BINARY[i], FULL_TO_BINARY[j])
                      for i, j in old.pairs if i in FULL_TO_BINARY)
        after = update_core(core, old.record)
        result.append(Branch(old.weight, pairs, old.record, old.memory,
                             old.exported, old.halted, after, clock_service(after)))
    return tuple(result)


def apply_branch(density, branch):
    if any(type(i) is not int or not 0 <= i < 66
           for key in density for i in (key[0], key[2])):
        raise ValueError("binary66 source coordinates required")
    return observer.apply_branch(density, branch)


def replay(tape, params, density=None):
    """Replay a supplied complete tagged tape without resetting the work."""
    if type(tape) is not tuple:
        raise ValueError("complete record tuple required")
    if density is None:
        density = {(i, 0, i, 0): w for i, w in enumerate(reference())}
    memory, exported, core, statuses = (), [], INITIAL_CORE, []
    for position, record in enumerate(tape):
        selected = tuple(b for b in branches(memory, core, params) if b.record == record)
        if len(selected) != 1:
            raise ValueError("record is absent or ambiguous in the complete kernel")
        branch = selected[0]
        density = apply_branch(density, branch)
        core, memory = branch.core, branch.memory
        if branch.exported is not None:
            exported.append(branch.exported)
        statuses.append(branch.clock_status)
        if branch.halted and position != len(tape)-1:
            raise ValueError("a halted replay cannot consume later records")
    return {"density": density, "memory": memory, "exported": tuple(exported),
            "core": core, "clock_statuses": tuple(statuses)}


def advance_marginal(blocks, params):
    """Source/reference/core/FIFO marginal; earlier exports are traced here.

    Input keys are (core, one-slot memory). Grouping is valid because the
    future selector uses only the last current sign and the old slot exports.
    The complete exported-history instrument is still branches()/replay().
    """
    if params.memory_capacity != 1 or params.overflow != "export_oldest":
        raise ValueError("one exported FIFO slot required for this marginal")
    groups = {}
    for (core, memory), density in blocks.items():
        _core(core)
        # Validate memory, including entries otherwise hidden by grouping.
        observer.branches(memory, params)
        key = core, observer.proposal(memory)
        if key not in groups:
            groups[key] = defaultdict(F)
        observer.add_into(groups[key], density)
    result = {}
    for (core, preferred), density in groups.items():
        memory = _representative_memory(int(preferred == native.F))
        for branch in branches(memory, core, params):
            value = apply_branch(density, branch)
            if value:
                key = branch.core, branch.memory
                if key not in result:
                    result[key] = defaultdict(F)
                observer.add_into(result[key], value)
    return {key: dict(value) for key, value in result.items()}


def _representative_memory(bit):
    # Only the selector is consumed; the exported synthetic slot is never
    # presented as an actual record in a complete history or receipt.
    return (observer.Record("repair", "explore", (0, 1), "C", 0, 1, 1),) if bit else ()


def elementary_certificate():
    """Check every binary branch and core; matrix-unit identities follow.

    For each supported i->j, sigma_after[j]=sigma_before[i]. Thus ratios,
    and K sigma^(is)=sigma_after^(is) K, agree on every matrix unit and
    tensor with an arbitrary untouched reference. Raw posterior corners
    equal w P_range sigma_after P_range; normalization is only a scalar.
    """
    if reference() != history_transport.reference():
        raise ArithmeticError("reference differs from the prepared predecessor")
    checked, matrix_units, completeness = 0, 0, 0
    for law in ("C", "E"):
        for seam in observer.SEAMS:
            effects = [F()] * 66
            for coin, label, current, pairs78 in observer._raw_maps(law, seam):
                pairs = tuple((FULL_TO_BINARY[i], FULL_TO_BINARY[j])
                              for i, j in pairs78 if i in FULL_TO_BINARY)
                if len(set(j for _, j in pairs)) != len(pairs):
                    raise ArithmeticError("native branch is not injective")
                record = observer.Record("repair", "explore", seam, law, coin, label, current)
                for i, _ in pairs:
                    effects[i] += F(1, 2)
                for core in CORES:
                    before, after = reference(core), reference(update_core(core, record))
                    for i, j in pairs:
                        if before[i] != after[j]:
                            raise ArithmeticError("core metadata failed reference transport")
                        checked += 1
                    matrix_units += len(pairs)**2
            if effects != [F(1)]*66:
                raise ArithmeticError("complete binary native instrument lost mass")
            completeness += 66
    return {"core_locations": len(CORES), "transported_columns_checked": checked,
            "matrix_units_certified_by_column_identity": matrix_units,
            "complete_native_input_columns_checked": completeness,
            "all_untouched_reference_extensions": True,
            "weights_holds_and_tagged_mixtures_preserve_identity": True}


def _events(params):
    """Exact source/next-bit/swap marginal of the full tagged kernel."""
    table = []
    for bit in (0, 1):
        rows = [defaultdict(F) for _ in BASIS]
        for branch in branches(_representative_memory(bit), INITIAL_CORE, params):
            if not branch.weight:
                continue
            edge = branch.record.seam if branch.record.current else None
            for i, j in branch.pairs:
                rows[i][j, int(branch.record.current > 0), edge] += branch.weight
        if any(sum(row.values(), F()) != 1 for row in rows):
            raise ArithmeticError("compressed transition omitted outcomes")
        table.append(tuple(dict(row) for row in rows))
    return tuple(table)


def _integer_kernel(params):
    if params.memory_capacity != 1 or params.overflow != "export_oldest":
        raise ValueError("one exported slot required")
    events = _events(params)
    denominator = lcm(*(w.denominator for rows in events for row in rows for w in row.values()))
    kernel = []
    for i in range(66):
        for core in CORES:
            for bit in (0, 1):
                row = defaultdict(int)
                for (j, next_bit, edge), weight in events[bit][i].items():
                    if edge is None:
                        after = core
                    else:
                        a, b = edge
                        after = tuple(sorted(b if p == a else a if p == b else p for p in core))
                    target = (j*66+CORE_INDEX[after])*2+next_bit
                    row[target] += int(weight*denominator)
                if sum(row.values()) != denominator:
                    raise ArithmeticError("integer row is not stochastic")
                kernel.append(tuple(row.items()))
    return denominator, tuple(kernel)


def compressed_run(params, horizon=24):
    """Exact diagonal source/core/last-current-bit marginal, no sampled path.

    The source is prepared once. Integer counts avoid repeated rational
    reduction; the transition denominator is carried explicitly. All future
    repair, hold and selector outcomes contribute, including unavailable
    clock-menu statuses. This does not retain the full exported tape in RAM.
    """
    if type(horizon) is not int or not 0 <= horizon <= 64:
        raise ValueError("an integer model-attempt horizon in 0..64 required")
    step_denominator, kernel = _integer_kernel(params)
    denominator = lcm(*(w.denominator for w in reference()))
    counts = [0] * (66*66*2)
    for i, w in enumerate(reference()):
        counts[(i*66+CORE_INDEX[INITIAL_CORE])*2] = int(w*denominator)
    snapshots = []
    menu = compiled_cores()
    checkpoints = {0, 1, 2, 4, 8, 12, 24, horizon}
    for attempt in range(horizon+1):
        if sum(counts) != denominator:
            raise ArithmeticError("complete marginal probability changed")
        if attempt in checkpoints:
            occupied = {CORES[(k//2) % 66] for k, value in enumerate(counts) if value}
            unavailable = sum(value for k, value in enumerate(counts)
                              if CORES[(k//2) % 66] not in menu)
            snapshots.append({"model_attempts": attempt, "total_probability": "1",
                              "positive_joint_states": sum(bool(n) for n in counts),
                              "reachable_cores": len(occupied),
                              "unavailable_clock_service_probability": str(F(unavailable, denominator)),
                              "expected_repairs": str(attempt*params.repair_probability),
                              "expected_holds": str(attempt*(1-params.repair_probability)),
                              "committed_records": attempt, "exported_records": max(0, attempt-1),
                              "retained_records": min(1, attempt)})
        if attempt == horizon:
            break
        following = [0] * len(counts)
        for i, count in enumerate(counts):
            if count:
                for j, multiplicity in kernel[i]:
                    following[j] += count*multiplicity
        counts, denominator = following, denominator*step_denominator
    digest = hashlib.sha256(json.dumps({"denominator": str(denominator),
                                        "counts": list(map(str, counts))},
                                       sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"parameters": params.as_dict(), "horizon": horizon,
            "integer_transition_denominator": step_denominator,
            "final_exact_distribution_sha256": digest, "snapshots": snapshots}


def last_record_counterexample(law="C"):
    """Identical endpoint matrix unit and last record require different cores."""
    if law not in ("C", "E"):
        raise ValueError("C or E required")
    first_label = 0 if law == "C" else -1
    histories = (((0, 5, 0, first_label), (6, 7, 0, 0)), ((0, 5, 0, 1), (6, 7, 0, 0)))
    params = FamilyParameters(c_weight=F(int(law == "C")))
    branches_old = history_transport.word_instrument(law, ((0, 5), (6, 7)))
    selected = [next(b for b in branches_old if b.history == h) for h in histories]
    out = (15, 51)
    reports = []
    for h, old in zip(histories, selected, strict=True):
        tape = tuple(observer.Record("repair", "explore", (a, b), law, coin, label,
                                     history_transport.records.record_current(law, (a, b, coin, label)))
                     for a, b, coin, label in h)
        actual = replay(tape, params)
        if any(i not in old.outputs for i in out):
            raise ArithmeticError("counterexample does not share the output matrix unit")
        ext = history_transport.extended_reference(law, old, history_transport.reference())
        if reference(actual["core"]) != ext or actual["memory"] != (tape[-1],):
            raise ArithmeticError("counterexample failed the actual exported-memory join")
        reports.append({"history": h, "core": actual["core"],
                        "reference_ratio_on_shared_output_unit": str(ext[out[0]]/ext[out[1]]),
                        "actual_history_probability": str(observer.trace(actual["density"])),
                        "exported_records": len(actual["exported"])})
    if reports[0]["reference_ratio_on_shared_output_unit"] == reports[1]["reference_ratio_on_shared_output_unit"]:
        raise ArithmeticError("last-record control no longer distinguishes reference frames")
    return {"law": law, "shared_output_matrix_unit": out, "histories": reports,
            "same_final_source_and_last_record_determine_reference": False}


def unavailable_witness():
    """Positive exploration reaches an uncompiled core in one real attempt."""
    params = FamilyParameters(c_weight=F(1))
    record = observer.Record("repair", "explore", (0, 7), "C", 0, 1, 1)
    actual = replay((record,), params)
    mass = observer.trace(actual["density"])
    if mass <= 0 or actual["core"] in compiled_cores():
        raise ArithmeticError("unavailable-service witness is absent")
    return {"record": {**record.__dict__, "seam": list(record.seam)},
            "core": actual["core"], "probability": str(mass),
            "status": actual["clock_statuses"][-1],
            "inventory_query_source_action": "identity",
            "clock_action_executed": False, "repair_outcome_retained": True}


def run():
    members = {"C": FamilyParameters(c_weight=F(1)),
               "E": FamilyParameters(c_weight=F(0)),
               "mixed_delayed": FamilyParameters(repair_probability=F(3, 4))}
    if _events(members["C"]) != _events(members["E"]):
        raise ArithmeticError("C/E diagonal compressed kernels unexpectedly differ")
    return {"schema": "oph.clock_defect_energy.continuation.v1", "binding": binding(),
            "source_sha256": source_bindings((
                "model/clock_defect_energy/continuation.py",
                "model/clock_defect_energy/test_continuation.py",
                "model/clock_defect_energy/contract.py", "tools/strict_json.py",
                "model/completion/observer.py", "model/completion/family.py",
                "model/native_modular_reference/history_transport.py",
                "model/native_modular_reference/history_clock.py",
                "model/observer/record_recovery.py", "model/repair_scheduler/coins.py")),
            "elementary_certificate": elementary_certificate(),
            "compiled_clock_core_inventory": sorted(compiled_cores()),
            "clock_menu_modular_parameter": str(history_clock.PARAMETER),
            "last_record_counterexamples": {law: last_record_counterexample(law) for law in ("C", "E")},
            "unavailable_service_witness": unavailable_witness(),
            "runs": {name: compressed_run(params) for name, params in members.items()},
            "storage": {"unordered_core_states": 66, "core_bits_upper": 7,
                        "core_times_feedback_states": 132, "controller_bits_upper": 8,
                        "fifo_capacity": 1, "complete_semantic_record_alphabet_upper": 513,
                        "fifo_bits_upper": 10,
                        "environmental_history_records_after_n_attempts": "max(0,n-1)",
                        "bounded_closed_total_storage": False},
            "scope": {"carrier_dimension": 66, "source_preparations": 1,
                      "full_extension_is_posterior": False,
                      "stationary_reference_replaced": False,
                      "new_clock_words_executed": False, "coherent_interventions_admitted": False,
                      "hardware_storage_or_controller_compiled": False,
                      "integrated_cp01_accepted": False,
                      "physical_clock_or_energy_calibration": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    elif args.check:
        check_receipt(OUTPUT, result)
        print("continuation receipt replay: PASS")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
