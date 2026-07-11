#!/usr/bin/env python3
"""Dynamic-circuit preflight for the primary record-gated Cayley kernel.

The implementation expands each random choice into balanced circuit variants.
This keeps every random source explicit in the receipt and avoids hiding a
software RNG inside the analysis.  No function in this module connects to IBM
or submits a job.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, qpy, transpile
from qiskit_aer import AerSimulator

from generative_repair_kernel import (
    CayleyRepairModel,
    builtin_cayley_models,
    open_loop_cayley_null,
    record_gated_cayley_kernel,
    sha256_json,
)


SCHEMA_VERSION = "oph.record-gated-cayley-circuit.v1"
SECTOR_QUBITS = 3


@dataclass(frozen=True)
class CayleyCircuitRecipe:
    model: str
    protocol: str
    initial_state: int
    disturbance_slot: int
    disturbance_generator: int | None
    second_slot: int
    repair_generator: int
    decision_state: int
    decision_generator: int
    open_loop_second_generator: int | None
    heated_state: int
    decision_record: int
    final_state: int
    state_encoding: tuple[int, ...]


def validate_encoding(model: CayleyRepairModel, encoding: tuple[int, ...]) -> None:
    if len(encoding) != model.size or len(set(encoding)) != model.size:
        raise ValueError("state encoding must map each semantic state to one physical code")
    if min(encoding) < 0 or max(encoding) >= 2**SECTOR_QUBITS:
        raise ValueError("state encoding does not fit the three-qubit register")


def _action(model: CayleyRepairModel, state: int, generator: int) -> int:
    return model.right_action[generator][state]


def build_recipe(
    model: CayleyRepairModel,
    protocol: str,
    initial_state: int,
    disturbance_slot: int,
    second_slot: int,
    state_encoding: tuple[int, ...] | None = None,
) -> CayleyCircuitRecipe:
    """Build one equal-weight variant of L_G R_G or L_G squared.

    There are ``2*|S|`` slots per lazy step.  The first ``|S|`` slots are
    duplicate stays and the second half applies each generator once.  For the
    record-gated second stage both halves propose each generator once, so the
    proposal is uniform while circuit multiplicity remains matched to the null.
    """

    model.validate()
    allowed_protocols = {
        "record_gated",
        "open_loop_heat",
        "delayed_record",
        "shuffled_record",
        "inverted_record",
    }
    if protocol not in allowed_protocols:
        raise ValueError(f"protocol must be one of {sorted(allowed_protocols)}")
    if initial_state not in range(model.size):
        raise ValueError("initial_state is outside the model")
    degree = len(model.generators)
    slots = 2 * degree
    if disturbance_slot not in range(slots) or second_slot not in range(slots):
        raise ValueError("lazy-step slot is outside the balanced range")
    if state_encoding is None:
        state_encoding = tuple(range(model.size))
    validate_encoding(model, state_encoding)

    disturbance_generator = None if disturbance_slot < degree else disturbance_slot - degree
    heated = (
        initial_state
        if disturbance_generator is None
        else _action(model, initial_state, disturbance_generator)
    )
    repair_generator = second_slot % degree
    repair_target = _action(model, heated, repair_generator)
    decision_state = initial_state if protocol == "delayed_record" else heated
    decision_generator = (
        (repair_generator + 1) % degree if protocol == "shuffled_record" else repair_generator
    )
    decision_target = _action(model, decision_state, decision_generator)
    decision = int(model.mismatch[decision_target] < model.mismatch[decision_state])
    if protocol == "inverted_record":
        decision = 1 - decision

    if protocol != "open_loop_heat":
        open_loop_generator = None
        final = repair_target if decision else heated
    else:
        open_loop_generator = None if second_slot < degree else second_slot - degree
        final = (
            heated
            if open_loop_generator is None
            else _action(model, heated, open_loop_generator)
        )

    return CayleyCircuitRecipe(
        model=model.name,
        protocol=protocol,
        initial_state=initial_state,
        disturbance_slot=disturbance_slot,
        disturbance_generator=disturbance_generator,
        second_slot=second_slot,
        repair_generator=repair_generator,
        decision_state=decision_state,
        decision_generator=decision_generator,
        open_loop_second_generator=open_loop_generator,
        heated_state=heated,
        decision_record=decision,
        final_state=final,
        state_encoding=state_encoding,
    )


def _apply_code(circuit: QuantumCircuit, register: QuantumRegister, code: int) -> None:
    for bit in range(SECTOR_QUBITS):
        circuit.x(register[bit]) if (code >> bit) & 1 else circuit.id(register[bit])


def _apply_xor(circuit: QuantumCircuit, register: QuantumRegister, source: int, target: int) -> None:
    mask = source ^ target
    for bit in range(SECTOR_QUBITS):
        circuit.x(register[bit]) if (mask >> bit) & 1 else circuit.id(register[bit])


def build_circuit(recipe: CayleyCircuitRecipe, name: str | None = None) -> QuantumCircuit:
    sector = QuantumRegister(SECTOR_QUBITS, "sector")
    decision_qubit = QuantumRegister(1, "decision_q")
    heated_record = ClassicalRegister(SECTOR_QUBITS, "heated")
    decision_record = ClassicalRegister(1, "decision")
    final_record = ClassicalRegister(SECTOR_QUBITS, "final")
    circuit = QuantumCircuit(
        sector,
        decision_qubit,
        heated_record,
        decision_record,
        final_record,
        name=name
        or (
            f"{recipe.model}__{recipe.protocol}__h{recipe.initial_state}__"
            f"u{recipe.disturbance_slot}__s{recipe.second_slot}"
        ),
    )

    initial_code = recipe.state_encoding[recipe.initial_state]
    heated_code = recipe.state_encoding[recipe.heated_state]
    final_code = recipe.state_encoding[recipe.final_state]
    _apply_code(circuit, sector, initial_code)
    _apply_xor(circuit, sector, initial_code, heated_code)
    circuit.measure(sector, heated_record)

    if recipe.decision_record:
        with circuit.if_test((heated_record, heated_code)):
            circuit.x(decision_qubit[0])
    else:
        with circuit.if_test((heated_record, heated_code)):
            circuit.id(decision_qubit[0])
    circuit.measure(decision_qubit[0], decision_record[0])

    if recipe.protocol != "open_loop_heat":
        with circuit.if_test((decision_record[0], 1)):
            _apply_xor(circuit, sector, heated_code, final_code)
        # Match the open-loop arm's unconditional three-operation envelope.
        for bit in range(SECTOR_QUBITS):
            circuit.id(sector[bit])
    else:
        # Keep the same record-conditioned control-flow envelope while making
        # the actual second lazy-heat move independent of the record.
        with circuit.if_test((decision_record[0], 1)):
            for bit in range(SECTOR_QUBITS):
                circuit.id(sector[bit])
        _apply_xor(circuit, sector, heated_code, final_code)

    circuit.reset(decision_qubit[0])
    circuit.measure(sector, final_record)
    circuit.metadata = {
        "schema_version": SCHEMA_VERSION,
        "recipe": asdict(recipe),
        "claim_boundary": "engineered record-feedback process; standard QM predicts the circuit",
    }
    return circuit


def qpy_hash(circuit: QuantumCircuit) -> str:
    buffer = io.BytesIO()
    qpy.dump(circuit, buffer)
    return hashlib.sha256(buffer.getvalue()).hexdigest()


def parse_counts(counts: Mapping[str, int]) -> list[dict[str, int]]:
    rows = []
    for key, count in counts.items():
        groups = key.split()
        if len(groups) != 3:
            raise ValueError(f"unexpected count key {key!r}")
        final_bits, decision_bits, heated_bits = groups
        rows.append(
            {
                "heated": int(heated_bits, 2),
                "decision": int(decision_bits, 2),
                "final": int(final_bits, 2),
                "count": int(count),
            }
        )
    return rows


def build_catalog(
    model: CayleyRepairModel,
    protocol: str,
    encoding: tuple[int, ...] | None = None,
) -> list[tuple[CayleyCircuitRecipe, QuantumCircuit]]:
    degree = len(model.generators)
    catalog = []
    for initial in range(model.size):
        for disturbance_slot in range(2 * degree):
            for second_slot in range(2 * degree):
                recipe = build_recipe(
                    model,
                    protocol,
                    initial,
                    disturbance_slot,
                    second_slot,
                    encoding,
                )
                catalog.append((recipe, build_circuit(recipe)))
    return catalog


def frozen_catalog_kernel(
    model: CayleyRepairModel,
    protocol: str,
    encoding: tuple[int, ...] | None = None,
) -> np.ndarray:
    """Compute the exact row kernel implied by the balanced public catalog."""

    counts = np.zeros((model.size, model.size), dtype=float)
    for recipe, _ in build_catalog(model, protocol, encoding):
        counts[recipe.initial_state, recipe.final_state] += 1.0
    return counts / counts.sum(axis=1, keepdims=True)


def run_ideal_preflight(
    model: CayleyRepairModel,
    protocol: str,
    shots: int,
    seed: int,
    encoding: tuple[int, ...] | None = None,
) -> dict[str, Any]:
    catalog = build_catalog(model, protocol, encoding)
    circuits = [circuit for _, circuit in catalog]
    backend = AerSimulator(seed_simulator=seed)
    isa = transpile(circuits, backend=backend, optimization_level=0, seed_transpiler=seed)
    result = backend.run(isa, shots=shots).result()

    transition_counts = np.zeros((model.size, model.size), dtype=int)
    inconsistent = 0
    total = 0
    raw_rows = []
    for (recipe, circuit), transpiled_circuit in zip(catalog, isa):
        counts = result.get_counts(transpiled_circuit)
        rows = parse_counts(counts)
        for row in rows:
            count = row["count"]
            total += count
            semantic_final = recipe.state_encoding.index(row["final"])
            transition_counts[recipe.initial_state, semantic_final] += count
            if (
                row["heated"] != recipe.state_encoding[recipe.heated_state]
                or row["decision"] != recipe.decision_record
                or row["final"] != recipe.state_encoding[recipe.final_state]
            ):
                inconsistent += count
        raw_rows.append(
            {
                "name": circuit.name,
                "recipe": asdict(recipe),
                "qpy_sha256": qpy_hash(circuit),
                "logical_depth": circuit.depth(),
                "operation_counts": {
                    str(key): int(value) for key, value in circuit.count_ops().items()
                },
                "raw_counts": dict(counts),
            }
        )

    empirical = transition_counts / transition_counts.sum(axis=1, keepdims=True)
    if protocol == "record_gated":
        expected = record_gated_cayley_kernel(model)
    elif protocol == "open_loop_heat":
        expected = open_loop_cayley_null(model)
    else:
        expected = frozen_catalog_kernel(model, protocol, encoding)
    receipt = {
        "schema_version": SCHEMA_VERSION,
        "receipt_class": "ideal_dynamic_circuit_preflight_not_hardware_evidence",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "model": asdict(model),
        "protocol": protocol,
        "shots_per_balanced_variant": shots,
        "seed": seed,
        "catalog_size": len(catalog),
        "transition_counts": transition_counts.tolist(),
        "empirical_kernel": empirical.tolist(),
        "expected_kernel": expected.tolist(),
        "max_kernel_error": float(np.max(np.abs(empirical - expected))),
        "inconsistent_record_or_feedback_fraction": inconsistent / total,
        "circuits": raw_rows,
        "claim_boundary": (
            "This validates a finite record-conditioned implementation against a frozen open-loop "
            "null. It is not a test of OPH against unrestricted quantum mechanics."
        ),
    }
    receipt["receipt_sha256"] = sha256_json(receipt)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run ideal preflight for the Cayley repair circuit.")
    parser.add_argument("--model", choices=sorted(builtin_cayley_models()), default="s3")
    parser.add_argument(
        "--protocol",
        choices=[
            "record_gated",
            "open_loop_heat",
            "delayed_record",
            "shuffled_record",
            "inverted_record",
        ],
        default="record_gated",
    )
    parser.add_argument("--shots", type=int, default=64)
    parser.add_argument("--seed", type=int, default=509)
    parser.add_argument("--encoding", type=int, nargs="+")
    parser.add_argument("--json-out", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    model = builtin_cayley_models()[args.model]
    encoding = tuple(args.encoding) if args.encoding is not None else None
    receipt = run_ideal_preflight(model, args.protocol, args.shots, args.seed, encoding)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
