from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


pytest.importorskip("qiskit")
pytest.importorskip("qiskit_aer")

PROGRAMS = Path(__file__).resolve().parents[1] / "programs"
sys.path.insert(0, str(PROGRAMS))
MODULE_PATH = PROGRAMS / "record_gated_cayley_circuits.py"
SPEC = importlib.util.spec_from_file_location("record_gated_cayley_circuits", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
circuits = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = circuits
SPEC.loader.exec_module(circuits)

import generative_repair_kernel as kernel  # noqa: E402


def empirical_recipe_kernel(model, protocol: str) -> np.ndarray:
    catalog = circuits.build_catalog(model, protocol)
    counts = np.zeros((model.size, model.size), dtype=float)
    for recipe, _ in catalog:
        counts[recipe.initial_state, recipe.final_state] += 1.0
    return counts / counts.sum(axis=1, keepdims=True)


def test_balanced_recipe_catalog_exactly_realizes_frozen_kernels() -> None:
    for model in kernel.builtin_cayley_models().values():
        primary = empirical_recipe_kernel(model, "record_gated")
        null = empirical_recipe_kernel(model, "open_loop_heat")
        assert np.allclose(primary, kernel.record_gated_cayley_kernel(model), atol=1e-12)
        assert np.allclose(null, kernel.open_loop_cayley_null(model), atol=1e-12)


def test_record_controls_are_frozen_and_not_aliases_of_primary() -> None:
    for model in kernel.builtin_cayley_models().values():
        primary = circuits.frozen_catalog_kernel(model, "record_gated")
        controls = [
            circuits.frozen_catalog_kernel(model, protocol)
            for protocol in ("delayed_record", "shuffled_record", "inverted_record")
        ]
        assert all(np.max(np.abs(control - primary)) > 0.0 for control in controls)


def test_dynamic_circuit_exposes_state_record_decision_and_feedback() -> None:
    model = kernel.builtin_cayley_models()["z5"]
    # Initial 2, stay disturbance, proposal -1 is a descending repair.
    recipe = circuits.build_recipe(
        model,
        protocol="record_gated",
        initial_state=2,
        disturbance_slot=0,
        second_slot=1,
    )
    assert recipe.decision_record == 1
    assert recipe.final_state == 1
    circuit = circuits.build_circuit(recipe)
    assert circuit.num_qubits == 4
    assert circuit.num_clbits == 7
    assert circuit.count_ops()["if_else"] == 2
    assert circuit.count_ops()["measure"] == 7
    assert circuit.metadata["recipe"]["heated_state"] == 2


def test_nontrivial_encoding_changes_physical_codes_not_semantic_kernel() -> None:
    model = kernel.builtin_cayley_models()["s3"]
    encoding = (5, 2, 7, 0, 6, 3)
    catalog = circuits.build_catalog(model, "record_gated", encoding)
    counts = np.zeros((model.size, model.size), dtype=float)
    for recipe, _ in catalog:
        assert recipe.state_encoding == encoding
        counts[recipe.initial_state, recipe.final_state] += 1.0
    empirical = counts / counts.sum(axis=1, keepdims=True)
    assert np.allclose(empirical, kernel.record_gated_cayley_kernel(model), atol=1e-12)


def test_ideal_dynamic_preflight_has_zero_process_error() -> None:
    model = kernel.builtin_cayley_models()["z5"]
    receipt = circuits.run_ideal_preflight(
        model=model,
        protocol="record_gated",
        shots=1,
        seed=509,
    )
    assert receipt["max_kernel_error"] == 0.0
    assert receipt["inconsistent_record_or_feedback_fraction"] == 0.0
