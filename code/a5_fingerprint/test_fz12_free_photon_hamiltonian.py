from __future__ import annotations

import ast
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import fz12_free_photon_hamiltonian as producer
import pytest
import verify_fz12_free_photon_hamiltonian_independent as independent

VERIFIER = producer.HERE / "verify_fz12_free_photon_hamiltonian_independent.py"


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(independent.canonical_json_bytes(value))


def resign(value: dict, field: str = "receipt_sha256") -> None:
    value[field] = independent.self_digest(value, field)


def independent_run(
    *,
    source: Path = independent.DEFAULT_SOURCE,
    remainder: Path = independent.DEFAULT_REMAINDER,
    free_photon_lean: Path = independent.DEFAULT_FREE_PHOTON_LEAN,
    receipt: Path = independent.DEFAULT_RECEIPT,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--source",
            str(source),
            "--remainder",
            str(remainder),
            "--free-photon-lean",
            str(free_photon_lean),
            "--receipt",
            str(receipt),
        ],
        cwd=producer.HERE,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def mutated_receipt(tmp_path: Path, mutation) -> Path:
    value = json.loads(independent.DEFAULT_RECEIPT.read_text())
    mutation(value)
    resign(value)
    path = tmp_path / "fz12_free_photon_hamiltonian_receipt.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_canonical(path, value)
    return path


def test_producer_and_independent_rebuild_are_byte_exact() -> None:
    built = producer.build_receipt()
    assert producer.verify_committed_receipt() == built
    assert independent.DEFAULT_RECEIPT.read_bytes() == producer.canonical_json_bytes(
        built
    )
    parents = [
        independent.load_parent(
            independent.DEFAULT_SOURCE, independent.PARENT_SPECS[0]
        ),
        independent.load_parent(
            independent.DEFAULT_REMAINDER, independent.PARENT_SPECS[1]
        ),
    ]
    lean_raw = independent.load_lean(independent.DEFAULT_FREE_PHOTON_LEAN)
    assert independent.expected_receipt(parents, lean_raw) == built


def test_separate_replay_does_not_import_producer_and_passes() -> None:
    tree = ast.parse(VERIFIER.read_text(encoding="utf-8"))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imports.add(node.module)
    assert "fz12_free_photon_hamiltonian" not in imports
    result = independent_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FZ12_TRANSVERSE_OSCILLATOR_SEPARATE_REPLAY_VALID"


def test_exact_basis_free_transverse_hamiltonian_contract() -> None:
    receipt = producer.build_receipt()
    transverse = receipt["transverse_sector"]
    assert transverse["formula"] == "P_T(k) = I - khat tensor khat"
    assert transverse["global_basis_selected"] is False
    assert all(row["numerical_rank"] == 2 for row in transverse["numerical_rows"])
    assert all(
        row["expected_spectrum"] == ["0", "1", "1"]
        for row in transverse["numerical_rows"]
    )
    assert all(
        row["expected_trace"] == "2"
        and row["idempotence_check_passed"] is True
        and row["transversality_check_passed"] is True
        and row["trace_check_passed"] is True
        and row["spectrum_check_passed"] is True
        for row in transverse["numerical_rows"]
    )

    hamiltonian = receipt["hamiltonian_contract"]
    assert hamiltonian["transverse_dimension_for_nonzero_k"] == 2
    assert hamiltonian["exact_frequency_relation"] == "omega_a(k)^2=Lambda_a(k)"
    assert hamiltonian["zero_mode"] == "Lambda_a(0)=omega_a(0)=0"
    assert hamiltonian["frequency_branch"].endswith(">=0")
    assert hamiltonian["transversality_preserved"] is True
    assert hamiltonian["declared_energy_first_variation_on_generator_zero"] is True
    assert hamiltonian["trajectory_energy_conservation_proved"] is False
    assert hamiltonian["auxiliary_parameter_only"] is True


def test_wave_packet_stays_inside_certified_domain_and_is_degenerate() -> None:
    packet = producer.build_receipt()["wave_packet_check"]
    assert packet["fixture_role"] == (
        "dimensionless auxiliary transverse-oscillator wave-packet check"
    )
    assert packet["q_support"] == "0 < q <= 1"
    assert packet["physical_wave_packet_claimed"] is False
    assert packet["local_transverse_witness_only"] is True
    assert packet["global_polarization_basis_claimed"] is False
    assert packet["center_shift_check_passed"] is True
    assert packet["center_shift_residual_clamped"] == "0"
    assert packet["two_local_transverse_norm_ratio_check_passed"] is True
    assert packet["local_transverse_orthogonality_check_passed"] is True
    assert packet["local_transverse_constraint_check_passed"] is True
    assert packet["ulp_residuals_clamped"] is True
    assert packet["time_reversed_positive_phase_is_also_solution"] is True


def test_lean_source_is_fixed_kernel_checked_and_admission_free() -> None:
    raw = independent.load_lean(independent.DEFAULT_FREE_PHOTON_LEAN)
    source = raw.decode("utf-8")
    assert len(raw) == independent.LEAN_BYTES
    assert independent.sha256(raw) == independent.LEAN_SHA256
    assert "theorem transversePolarization_finrank" in source
    assert "theorem photonModeFrequency_sq" in source
    assert "theorem photonModeEnergy_firstVariation_generator_zero" in source
    assert "#print axioms photonModeEnergy_firstVariation_generator_zero" in source
    assert "sorry" not in source
    assert "admit" not in source


def test_every_physical_attachment_and_comparison_stays_fail_closed() -> None:
    receipt = producer.build_receipt()
    scope = receipt["scope"]
    assert scope["declared_oscillator_relation_proved"] is True
    assert scope["physical_prediction"] is False
    assert scope["comparison_permitted"] is False

    boundary = receipt["physical_boundary"]
    assert all(value is False for key, value in boundary.items() if key != "statement")
    assert "conditional oscillator generator" in boundary["statement"]
    exposure = receipt["exposure_boundary"]
    assert exposure["comparison_inputs"] == []
    assert all(
        value is False for key, value in exposure.items() if key != "comparison_inputs"
    )


def test_semantic_guards_have_matching_fail_closed_mutations(
    tmp_path: Path,
) -> None:
    controls = producer.build_receipt()["semantic_receipt_mutation_guards"]
    assert controls["producer_executes_counterfactuals"] is False
    assert controls["test_suite_exercises_guards"] is True
    assert controls["guards"] == [
        "longitudinal rank-three fiber",
        "non-idempotent projector",
        "transverse-direction-split stiffness",
        "negative square-root branch",
        "mutated FZ-12 symbol sign",
        "q support above one",
        "physical clock promotion",
        "electron or interaction promotion",
    ]

    mutations = (
        lambda value: value["hamiltonian_contract"].__setitem__(
            "transverse_dimension_for_nonzero_k", 3
        ),
        lambda value: value["transverse_sector"]["numerical_rows"][0].__setitem__(
            "idempotence_check_passed", False
        ),
        lambda value: value["hamiltonian_contract"].__setitem__(
            "polarization_degeneracy", "split transverse stiffness"
        ),
        lambda value: value["hamiltonian_contract"].__setitem__(
            "frequency_branch", "omega_a(k)=-sqrt(Lambda_a(k))"
        ),
        lambda value: value["exact_symbol_contract"].__setitem__(
            "lambda_hat", "-(1/5) sum_e [1-cos(q (w_e dot n))]"
        ),
        lambda value: value["wave_packet_check"].__setitem__(
            "q_support", "0 < q <= 2"
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "auxiliary_parameter_identified_with_physical_clock", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "electron_dispersion_selected", True
        ),
    )
    for index, mutation in enumerate(mutations):
        path = mutated_receipt(tmp_path / str(index), mutation)
        result = independent_run(receipt=path)
        assert result.returncode == 1
        assert "semantic drift" in result.stderr


@pytest.mark.parametrize(
    "mutation",
    [
        lambda value: value["scope"].__setitem__("physical_prediction", True),
        lambda value: value["scope"].__setitem__("comparison_permitted", True),
        lambda value: value["transverse_sector"].__setitem__(
            "global_basis_selected", True
        ),
        lambda value: value["wave_packet_check"].__setitem__(
            "global_polarization_basis_claimed", True
        ),
        lambda value: value["wave_packet_check"].__setitem__(
            "observed_center_shift_6_digits", "0"
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "transverse_fiber_identified_with_laboratory_photon", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "source_selected_hamiltonian", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "electromagnetic_interaction_or_pair_kinematics_proved", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "gauge_or_gauss_constraint_quotient_proved", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "physical_masslessness_proved", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "trajectory_energy_conservation_proved", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "reality_pairing_k_and_minus_k_proved", True
        ),
        lambda value: value["physical_boundary"].__setitem__(
            "canonical_symplectic_form_or_kinetic_normalization_selected", True
        ),
        lambda value: value["exposure_boundary"].__setitem__(
            "comparison_permitted", True
        ),
        lambda value: value["exposure_boundary"].__setitem__(
            "comparison_inputs", ["public-data.csv"]
        ),
        lambda value: value["parent_pins"][2].__setitem__("bytes", 1),
        lambda value: value["parent_pins"][0].__setitem__("sha256", "sha256:00"),
        lambda value: value["semantic_receipt_mutation_guards"].__setitem__(
            "test_suite_exercises_guards", False
        ),
        lambda value: value.__setitem__("unregistered_claim", True),
    ],
)
def test_resigned_semantic_promotions_and_drifts_are_rejected(
    tmp_path: Path, mutation
) -> None:
    path = mutated_receipt(tmp_path, mutation)
    result = independent_run(receipt=path)
    assert result.returncode == 1
    assert "semantic drift" in result.stderr


@pytest.mark.parametrize(
    ("parent", "spec_index", "self_field", "argument"),
    [
        (independent.DEFAULT_SOURCE, 0, "receipt_sha256", "source"),
        (independent.DEFAULT_REMAINDER, 1, "receipt_sha256", "remainder"),
    ],
)
def test_resigned_parent_mutation_fails_fixed_raw_pin(
    tmp_path: Path,
    parent: Path,
    spec_index: int,
    self_field: str,
    argument: str,
) -> None:
    value = json.loads(parent.read_text())
    value["audit_mutation"] = True
    resign(value, self_field)
    path = tmp_path / parent.name
    write_canonical(path, value)
    kwargs = {argument: path}
    result = independent_run(**kwargs)
    assert result.returncode == 1
    assert any(
        message in result.stderr
        for message in ("parent byte-count drift", "parent raw hash drift")
    )
    assert independent.PARENT_SPECS[spec_index]["sha256"] != independent.sha256(
        path.read_bytes()
    )


def test_lean_text_mutation_fails_fixed_raw_pin(tmp_path: Path) -> None:
    source = independent.DEFAULT_FREE_PHOTON_LEAN.read_text(encoding="utf-8")
    source = source.replace("= 2 := by", "= 3 := by", 1)
    path = tmp_path / "SeamCurrentFreePhotonLift.lean"
    path.write_bytes(source.encode("utf-8"))
    result = independent_run(free_photon_lean=path)
    assert result.returncode == 1
    assert "Lean raw hash drift" in result.stderr


def test_unresigned_receipt_mutation_fails_self_digest(tmp_path: Path) -> None:
    value = copy.deepcopy(json.loads(independent.DEFAULT_RECEIPT.read_text()))
    value["scope"]["physical_prediction"] = True
    path = tmp_path / "fz12_free_photon_hamiltonian_receipt.json"
    write_canonical(path, value)
    result = independent_run(receipt=path)
    assert result.returncode == 1
    assert "self-digest drift" in result.stderr
