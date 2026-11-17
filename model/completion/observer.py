"""Continuing finite-memory C/E observer, with an explicit exported history.

The source is H_78; the controller and committed FIFO are classical. Each
branch is sqrt(weight) times a partial injection, tensored with identity on
any reference. Exact sparse density propagation uses rational coefficients.
An environmental output retains overwritten records: bounded controller
memory is not bounded total storage. All parameters count model attempts,
not physical duration. This is an ideal instrument execution, not a newly
compiled source-word implementation of its ancillary services.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from model.observer import native, quantum_mixing
from model.repair_scheduler import coins
from .family import FamilyParameters, binding


HERE = Path(__file__).resolve().parent
STATES = native.states()
SEAMS = quantum_mixing.seams()
IDENTITY = tuple((i, i) for i in range(len(STATES)))


@dataclass(frozen=True)
class Record:
    kind: str
    selector: str = "none"
    seam: tuple[int, int] | None = None
    law: str | None = None
    coin: int | None = None
    label: int | None = None
    current: int = 0


@dataclass(frozen=True)
class Branch:
    weight: F
    pairs: tuple[tuple[int, int], ...]
    record: Record
    memory: tuple[Record, ...]
    exported: Record | None
    halted: bool = False


def proposal(memory):
    """The only controller input is the last retained current, never loads."""
    return native.F if memory and memory[-1].current > 0 else native.E


def _valid_record(record):
    if not isinstance(record, Record):
        return False
    if record.kind == "hold":
        return record == Record("hold")
    if (record.kind != "repair" or record.selector not in ("explore", "feedback")
            or record.seam not in SEAMS or record.law not in ("C", "E")
            or type(record.coin) is not int or record.coin not in (0, 1)
            or type(record.label) is not int or type(record.current) is not int):
        return False
    return any((record.coin, record.label, record.current) == row[:3]
               for row in _raw_maps(record.law, record.seam))


def _memory_update(memory, record, params):
    if len(memory) < params.memory_capacity:
        return (*memory, record), None
    if params.overflow != "export_oldest":
        raise ValueError("halt-on-full is handled before an attempted write")
    return (*memory[1:], record), memory[0]


@lru_cache(maxsize=128)
def _raw_maps(law, seam):
    """Keep every semantic coin/record outcome, including zero operators."""
    observed = {(coin, label): pairs
                for coin, label, pairs in coins.coin_maps(law, seam, binary=False)}
    labels = range(-1, 2) if law == "C" else range(-2, 3)
    result = []
    for coin in (0, 1):
        for label in labels:
            current = label if law == "C" else (label-(2*coin-1)*(abs(label) % 2))//2
            result.append((coin, label, current, observed.get((coin, label), ())))
    return tuple(result)


def branches(memory, params, *, preferred=None):
    """One complete committed instrument; a declared proposal intervention
    changes only the preferred seam and retains the uniform exploration floor.
    Full memory with the halt policy emits a timeout without source evolution.
    A hold is an actual committed record and sets the last-current control to0.
    """
    params.validate()
    if (not isinstance(memory, tuple) or len(memory) > params.memory_capacity
            or any(not _valid_record(record) for record in memory)):
        raise ValueError("a valid bounded classical record memory required")
    if preferred is not None and preferred not in SEAMS:
        raise ValueError("intervention must name an actual native seam")
    if len(memory) == params.memory_capacity and params.overflow == "halt":
        return (Branch(F(1), IDENTITY, Record("timeout"), memory, None, True),)
    preferred = proposal(memory) if preferred is None else preferred
    result = []
    hold = Record("hold")
    after, exported = _memory_update(memory, hold, params)
    result.append(Branch(1-params.repair_probability, IDENTITY, hold, after, exported))
    choices = [("explore", seam, params.exploration_probability/len(SEAMS))
               for seam in SEAMS]
    choices.append(("feedback", preferred, 1-params.exploration_probability))
    for selector, seam, selection_weight in choices:
        for law, law_weight in (("C", params.c_weight), ("E", 1-params.c_weight)):
            for coin, label, current, pairs in _raw_maps(law, seam):
                record = Record("repair", selector, seam, law, coin, label, current)
                after, exported = _memory_update(memory, record, params)
                weight = params.repair_probability*selection_weight*law_weight/2
                result.append(Branch(weight, pairs, record, after, exported))
    return tuple(result)


def apply_branch(density, branch):
    """Density keys are (source_left, reference_left, source_right, reference_right).
    Complex coefficients are allowed; retained experiments use exact fractions.
    Positivity is a preparation obligation. Linear maps also accept matrix units.
    """
    if not branch.weight or not branch.pairs:
        return {}
    mapping = dict(branch.pairs)
    result = {}
    for (i, a, j, b), value in density.items():
        if i in mapping and j in mapping:
            result[mapping[i], a, mapping[j], b] = branch.weight*value
    return result


def trace(density):
    return sum((value for (i, a, j, b), value in density.items()
                if i == j and a == b), F(0))


def add_into(target, source):
    for key, value in source.items():
        target[key] += value


def reference_marginal(density):
    result = defaultdict(F)
    for (i, a, j, b), value in density.items():
        if i == j:
            result[a, b] += value
    return dict(result)


def source_marginal(density):
    result = defaultdict(F)
    for (i, a, j, b), value in density.items():
        if a == b:
            result[i, j] += value
    return dict(result)


def summed_density(blocks):
    result = defaultdict(F)
    for density in blocks.values():
        add_into(result, density)
    return dict(result)


def advance_marginal(blocks, params):
    """Exact marginal on source, reference and the one-slot retained record.

    The old FIFO records and earlier output tape are traced out in THIS output
    only; `branches`/`history` still specify their complete joint instrument.
    Grouping by the controller's selected seam is exact because a one-slot
    export replaces all of the old memory and the old value is exported.
    """
    if params.memory_capacity != 1 or params.overflow != "export_oldest":
        raise ValueError("the efficient marginal runner requires one exported slot")
    groups = {}
    for memory, density in blocks.items():
        key = proposal(memory)
        if key not in groups:
            groups[key] = defaultdict(F)
        add_into(groups[key], density)
    result = {}
    for preferred, density in groups.items():
        for branch in branches((), params, preferred=preferred):
            after = apply_branch(density, branch)
            if not after:
                continue
            if branch.memory not in result:
                result[branch.memory] = defaultdict(F)
            add_into(result[branch.memory], after)
    return {memory: dict(density) for memory, density in result.items()}


def history(density, records, params, *, interventions=None):
    """Replay any complete recorded word, without conditioning away its mass.
    Exports are returned in order. No exponential output tape is materialized.
    An impossible or zero-weight word returns the zero operator.
    """
    memory, exported = (), []
    current = dict(density)
    interventions = {} if interventions is None else interventions
    for position, record in enumerate(records):
        options = branches(memory, params, preferred=interventions.get(position))
        found = [branch for branch in options if branch.record == record]
        if len(found) != 1:
            return {}, memory, tuple(exported)
        branch = found[0]
        current = apply_branch(current, branch)
        memory = branch.memory
        if branch.exported is not None:
            exported.append(branch.exported)
    return current, memory, tuple(exported)


def regional_view(records, ports):
    """Read-only restriction of the shared committed history, with event indices.
    Agreement here assumes the same tape; it does not build remote replicas.
    """
    ports = frozenset(ports)
    return tuple((position, record) for position, record in enumerate(records)
                 if record.kind == "repair" and set(record.seam) <= ports)


def completeness(memory, params):
    """Partial injection plus sum of branch effects I is an exact certificate
    for every matrix unit, hence arbitrary source/reference input and histories.
    """
    effects = [F(0)]*len(STATES)
    rows = branches(memory, params)
    for branch in rows:
        inputs = [i for i, _ in branch.pairs]
        outputs = [j for _, j in branch.pairs]
        if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
            raise ArithmeticError("branch ceased to be a partial injection")
        for i in inputs:
            effects[i] += branch.weight
    if effects != [F(1)]*len(STATES):
        raise ArithmeticError("complete branch effects do not sum to identity")
    return {"semantic_branches_including_zero": len(rows),
            "nonzero_operator_branches": sum(bool(b.weight and b.pairs) for b in rows),
            "zero_operator_or_weight_branches": sum(not (b.weight and b.pairs) for b in rows),
            "source_dimension": len(effects), "sum_effects": "I_78",
            "reference_dimension_limit": None}


def entangled_preparation():
    """Supplied rank-two source/reference state with a current-distinguishing pair."""
    pair = [None, None]
    for i, state in enumerate(STATES):
        if state[4] == 1 and state[0] == 1:
            pair[0] = i
        if state[4] == 1 and state[1] == 1:
            pair[1] = i
    if None in pair:
        raise ArithmeticError("source preparation support changed")
    return {(pair[a], a, pair[b], b): F(1, 2) for a in (0, 1) for b in (0, 1)}


def run(preparation, params, steps):
    if type(steps) is not int or not 0 <= steps <= 32:
        raise ValueError("an exact bounded horizon of0..32 model attempts required")
    blocks = {(): preparation}
    initial_reference = reference_marginal(preparation)
    summaries = []
    for attempt in range(1, steps+1):
        blocks = advance_marginal(blocks, params)
        density = summed_density(blocks)
        if trace(density) != trace(preparation):
            raise ArithmeticError("trace changed")
        if reference_marginal(density) != initial_reference:
            raise ArithmeticError("unconditional untouched reference changed")
        source = source_marginal(density)
        summaries.append({"model_attempt": attempt, "trace": str(trace(density)),
                          "reachable_last_records": len(blocks),
                          "expected_repair_events": str(attempt*params.repair_probability),
                          "expected_hold_events": str(attempt*(1-params.repair_probability)),
                          "committed_records_per_history": attempt,
                          "source_coherence_l1": str(sum((abs(v) for (i, j), v in source.items()
                                                        if i != j), F(0))),
                          "joint_source_reference_offdiagonal_l1": str(sum(
                              (abs(v) for (i, a, j, b), v in density.items() if i != j), F(0))),
                          "exported_fifo_records_per_history": max(0, attempt-1),
                          "retained_fifo_records_per_history": 1})
    return blocks, summaries


def report():
    params = FamilyParameters(c_weight=F(1, 2), repair_probability=F(3, 4),
                              exploration_probability=F(1, 2), memory_capacity=1,
                              overflow="export_oldest")
    certificates = {}
    for c in (F(0), F(1, 2), F(1)):
        p = FamilyParameters(c_weight=c, repair_probability=F(3, 4),
                             exploration_probability=F(1, 2), memory_capacity=1,
                             overflow="export_oldest")
        for current in (-1, 0, 1):
            memory = (Record("repair", "feedback", native.E, "C", 0, current, current),)
            certificates[f"C_weight={c},last_current={current}"] = completeness(memory, p)
    initial = entangled_preparation()
    blocks, summaries = run(initial, params, 6)
    final = summed_density(blocks)
    endpoint_runs = {}
    for name, c in (("C", F(1)), ("E", F(0))):
        endpoint = FamilyParameters(c_weight=c, repair_probability=F(1),
                                    exploration_probability=F(1, 2), memory_capacity=1,
                                    overflow="export_oldest")
        _, endpoint_runs[name] = run(initial, endpoint, 6)
    first = [(b, apply_branch(initial, b)) for b in branches((), params)]
    selected = next((b, rho) for b, rho in first if trace(rho) and b.record.kind == "repair"
                    and b.record.seam == native.E and b.record.current == 1)
    b, rho = selected
    probability = trace(rho)
    conditional_reference = {str(k): str(v/probability)
                             for k, v in reference_marginal(rho).items()}
    # Prefix example beyond the old two-event halt. Holds have nonzero weight;
    # their old contents are exported and never fed back into the controller.
    records = (Record("hold"),)*6
    path, memory, exported = history(initial, records, params)
    if trace(path) != F(1, 4)**6 or exported+memory != records:
        raise ArithmeticError("history probability or export reconstruction failed")
    # Independent two-event reference countercontrol shared with CP-13.
    uniform = {(i, 0, i, 0): F(1, 66) for i, x in enumerate(STATES) if max(x) == 1}
    adaptive = FamilyParameters(c_weight=F(1), repair_probability=F(1),
                                exploration_probability=F(1, 2), memory_capacity=1,
                                overflow="export_oldest")
    adaptive_blocks, _ = run(uniform, adaptive, 2)
    source = source_marginal(summed_density(adaptive_blocks))
    changes = [source.get((i, i), F(0))-F(1, 66)
               for i, state in enumerate(STATES) if max(state) == 1]
    if max(map(abs, changes)) != F(1, 960):
        raise ArithmeticError("adaptive-reference countercontrol changed")
    source_paths = [Path(__file__), HERE/"family.py",
                    HERE.parent/"observer/native.py", HERE.parent/"repair_scheduler/coins.py"]
    return {"schema": "oph.completion.observer.v1", "gate": "CP-01", "campaign_binding": binding(),
            "status": "BOUNDED_COMPONENT_EXECUTED_GATE_OPEN", "parameters": params.as_dict(),
            "source_sha256": {path.relative_to(HERE.parents[1]).as_posix():
                              hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
            "instrument_certificates": certificates, "six_attempt_history": summaries,
            "C_E_six_attempt_histories": endpoint_runs,
            "untouched_reference_final": {str(k): str(v) for k, v in reference_marginal(final).items()},
            "conditional_reference_disturbance_control": {
                "branch_probability": str(probability), "selected_current": 1,
                "conditional_reference": conditional_reference,
                "unconditional_reference_unchanged": True},
            "complete_six_hold_word": {"probability": str(trace(path)),
                "exports": len(exported), "retained": len(memory),
                "reconstructs_complete_word": exported+memory == records},
            "adaptive_reference_control": {"two_attempt_changed_binary_probabilities": sum(bool(x) for x in changes),
                "maximum_deviation_from_uniform": str(max(map(abs, changes))),
                "uniform_source_reference_preserved": False},
            "resource_boundary": {"source_dimension": 78, "demonstrated_reference_dimension": 2,
                "controller_reads": "last retained current only", "memory_slots": 1,
                "history_representation": "complete Kraus transition kernel and exact prefix replay; no exponential history tree materialized",
                "total_history_storage_bounded": False,
                "exported_record_count_after_n_attempts": "max(0,n-1)",
                "supplied_services": ["initial source/reference state", "classical random selection",
                    "committed C/E instrument", "FIFO writes and export", "classical record feedback"],
                "new_native_source_word_execution": False,
                "physical_clock_or_energy_normalization": False,
                "remote_record_replication_constructed": False,
                "q3_common_code_history": "unexecuted, not used by this component"}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = report()
    if args.check:
        from tools import strict_json
        if strict_json.load(args.check) != result:
            raise SystemExit("observer receipt differs")
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    if not args.check and not args.write:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
