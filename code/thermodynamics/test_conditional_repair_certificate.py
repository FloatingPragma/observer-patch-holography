"""Adversarial tests for the conditional-repair thermodynamics
certificate."""

from __future__ import annotations

import copy
import json
import random
from fractions import Fraction
from pathlib import Path

import pytest

import conditional_repair_certificate as cert

HERE = Path(__file__).resolve().parent


@pytest.fixture(scope="module")
def receipt() -> dict:
    return json.loads(cert.RECEIPT_PATH.read_text())


@pytest.fixture(scope="module")
def rebuilt() -> dict:
    return cert.build_receipt()


def test_receipt_matches_fresh_rebuild(receipt, rebuilt):
    assert receipt == rebuilt


def test_receipt_self_digest(receipt):
    body = {k: v for k, v in receipt.items() if k != "receipt_sha256"}
    assert receipt["receipt_sha256"] == cert.tagged_sha256(
        cert.canonical_json_bytes(body)
    )


def test_mutated_receipt_fails_digest(receipt):
    mutated = copy.deepcopy(receipt)
    mutated["law_map"]["second"] = "entropy always increases"
    body = {k: v for k, v in mutated.items() if k != "receipt_sha256"}
    assert mutated["receipt_sha256"] != cert.tagged_sha256(
        cert.canonical_json_bytes(body)
    )


def test_status_records_open_receipts(receipt):
    assert "COMMON_REFERENCE_REALIZATION_AND_ENERGY_CLOCK_OPEN" in (
        receipt["status"]
    )
    assert set(receipt["open_receipts"]) == {
        "THERMO-GLOBAL",
        "THERMO-COMMON-REFERENCE",
        "THERMO-REALIZATION",
        "THERMO-ENERGY-CLOCK",
        "THERMO-LOW-T-REFINEMENT",
    }


def test_strict_descent_control_is_recorded(receipt):
    control = receipt["strict_descent_control"]
    assert control["normalizer"] == {"a": "b", "b": "b"}
    assert control["entropy_after"] == "0 (deterministic)"


def test_kernel_by_independent_route():
    """Rebuild the kernel as an explicit conditional probability table
    and compare with the certificate's construction."""

    rng = random.Random(7)
    for _ in range(10):
        size = rng.randint(3, 8)
        labels = [rng.randrange(3) for _ in range(size)]
        pi = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        kernel, _ = cert.heat_bath(labels, pi)
        total = sum(pi)
        for x in range(size):
            for y in range(size):
                # conditional probability of y given the fibre of x
                joint = pi[y] / total if labels[y] == labels[x] else (
                    Fraction(0)
                )
                fibre_mass = (
                    sum(
                        pi[z] / total
                        for z in range(size)
                        if labels[z] == labels[x]
                    )
                )
                assert kernel[x][y] == joint / fibre_mass


def test_kernel_contracts_chi_squared():
    """A second contraction functional: the chi-squared divergence to
    the reference also contracts under the kernel, exactly over the
    rationals."""

    rng = random.Random(11)
    for _ in range(10):
        size = rng.randint(3, 7)
        labels = [rng.randrange(2) for _ in range(size)]
        pi = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        total = sum(pi)
        pi_n = [w / total for w in pi]
        kernel, _ = cert.heat_bath(labels, pi)
        raw = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        tot = sum(raw)
        p = [r / tot for r in raw]
        pushed = [
            sum(p[x] * kernel[x][y] for x in range(size))
            for y in range(size)
        ]
        chi_before = sum(
            (p[i] - pi_n[i]) ** 2 / pi_n[i] for i in range(size)
        )
        chi_after = sum(
            (pushed[i] - pi_n[i]) ** 2 / pi_n[i] for i in range(size)
        )
        assert chi_after <= chi_before


def test_lean_binding_paths_exist(receipt):
    root = HERE.parent.parent
    for path in (
        "Lean/Thermodynamics/FiniteConditionalRepair.lean",
        "Lean/Thermodynamics/StationaryRealization.lean",
        "Lean/Thermodynamics/FirstLawIdentity.lean",
        "Lean/Thermodynamics/FluctuationTheorems.lean",
        "Lean/Thermodynamics/CapFirstLaw.lean",
        "Lean/Thermodynamics/EinsteinPremiseLink.lean",
        "Lean/EventAlgebra/PartitionPinchingCP.lean",
    ):
        assert (root / path).exists()
    bindings = receipt["lean_bindings"]
    assert "heatBath_secondLaw" in bindings["conditional_repair"]
    assert "stationary_secondLaw" in bindings["stationary_realization"]
    assert "firstLaw_split" in bindings["first_law"]
    assert "kraus_complete" in bindings["record_channel"]
    assert "integral_fluctuation" in bindings["fluctuation"]
    assert "heatBath_cap_clausius" in bindings["cap_first_law"]
    assert "thermoFirstLawData_passes" in bindings["einstein_premise_link"]


def test_bound_lean_declarations_exist_in_sources(receipt):
    """Every declaration named in lean_bindings exists in its file."""
    root = HERE.parent.parent
    for binding in receipt["lean_bindings"].values():
        path_part, names = binding.split(":", 1)
        source = (root / "Lean" / path_part.strip()).read_text(
            encoding="utf-8"
        )
        for name in names.split(","):
            token = name.strip()
            assert token, "empty declaration token"
            assert token in source, f"{token} missing from {path_part}"


def test_integral_fluctuation_independent_route():
    """The exact integral fluctuation identity, rebuilt with the sum
    grouped by fibre rather than by pair."""

    rng = random.Random(23)
    for _ in range(8):
        size = rng.randint(3, 7)
        labels = [rng.randrange(2) for _ in range(size)]
        pi_raw = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        total = sum(pi_raw)
        pi = [w / total for w in pi_raw]
        kernel, _ = cert.heat_bath(labels, pi)
        raw = [Fraction(rng.randint(1, 9)) for _ in range(size)]
        tot = sum(raw)
        p = [r / tot for r in raw]
        q = [
            sum(p[x] * kernel[x][y] for x in range(size))
            for y in range(size)
        ]
        # group by fibre: sum_x p x K x y exp(-sigma) collapses to
        # sum over y of q y once the x-sum is carried out per fibre
        acc = Fraction(0)
        for lab in set(labels):
            fibre = [i for i in range(size) if labels[i] == lab]
            for y in fibre:
                inner = sum(
                    p[x] * kernel[x][y] * (pi[x] / p[x]) * (q[y] / pi[y])
                    for x in fibre
                )
                acc += inner
        assert acc == 1


def test_realization_probe_receipt_matches_rebuild():
    import collar_matrix_realization_probe as probe

    committed = json.loads(probe.PROBE_PATH.read_text())
    rebuilt = probe.build_probe()
    assert committed == rebuilt


def test_realization_probe_self_digest():
    import collar_matrix_realization_probe as probe

    committed = json.loads(probe.PROBE_PATH.read_text())
    body = {
        k: v for k, v in committed.items() if k != "receipt_sha256"
    }
    assert committed["receipt_sha256"] == cert.tagged_sha256(
        cert.canonical_json_bytes(body)
    )


def test_realization_probe_states_receipt_open():
    import collar_matrix_realization_probe as probe

    committed = json.loads(probe.PROBE_PATH.read_text())
    assert committed["receipt_target"] == "THERMO-REALIZATION"
    assert "RECEIPT_OPEN" in committed["status"]
    assert "RAW_CHAIN_REDUCIBLE" in committed["status"]
    assert committed["measurements"]["off_fibre_mass_max"] == 0.0
    assert committed["measurements"]["protected_datum_cardinality"] == 1
    assert any(
        "reducible" in blocker
        for blocker in committed["inherited_blockers"]
    )


def test_realization_probe_exhausts_committed_field_subset_projections():
    import collar_matrix_realization_probe as probe

    committed = json.loads(probe.PROBE_PATH.read_text())
    audit = committed["raw_coarsening_audit"]
    assert audit["quotient_count"] == 15
    assert len(audit["rows"]) == 15
    assert audit["nontrivial_irreducible_syntactic_quotient_count"] == 4
    assert audit["nontrivial_irreducible_reversible_count"] == 0

    selected = audit["selected_raw_equilibrium_probe"]
    assert selected["packet_fields"] == ["repair_load_bucket"]
    assert selected["state_count"] == 8
    assert selected["irreducible"] is True
    assert selected["aperiodic"] is True
    assert selected["stationary_min"] > 0.0
    assert selected["detailed_balance_max_err"] > probe.DB_TOL
    assert selected["microscopic_reversibility_claimed"] is False
    assert selected["protected_charge_test_nontrivial"] is False
    assert (
        selected["common_reference_with_state_optimizer_identified"]
        is False
    )


def test_realization_probe_recurrent_route_is_singleton_freezeout():
    import collar_matrix_realization_probe as probe

    committed = json.loads(probe.PROBE_PATH.read_text())
    recurrent = committed["recurrent_class_audit"]
    assert recurrent["fine_quotient_closed_class_count"] == 1
    assert recurrent["closed_class_sizes"] == [1]
    assert recurrent["closed_class_state_indices"] == [[12]]
    assert recurrent["nontrivial_recurrent_restriction_available"] is False


def test_support_graph_audit_detects_closed_classes():
    import collar_matrix_realization_probe as probe

    matrix = [
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 1.0, 0.0],
    ]
    assert probe.closed_communicating_classes(matrix) == [[0, 1]]


def test_realization_probe_pins_enforced(monkeypatch, tmp_path):
    import collar_matrix_realization_probe as probe

    tampered = tmp_path / "finite_repair_transition_matrix.npz"
    tampered.write_bytes(
        probe.MATRIX_PATH.read_bytes() + b"\x00"
    )
    monkeypatch.setattr(probe, "MATRIX_PATH", tampered)
    with pytest.raises(cert.ThermoError):
        probe.build_probe()


def test_tampered_kernel_rejected_by_fluctuation(monkeypatch):
    """A perturbed kernel entry breaks the exact fluctuation
    identities."""

    original = cert.heat_bath

    def tampered(labels, pi):
        kernel, mass = original(labels, pi)
        kernel[0][0] += Fraction(1, 997)
        return kernel, mass

    monkeypatch.setattr(cert, "heat_bath", tampered)
    with pytest.raises(cert.ThermoError):
        cert.fluctuation_certificate()


def test_tampered_kernel_rejected(monkeypatch):
    """Perturb one kernel entry and confirm the producer fails closed."""

    original = cert.heat_bath

    def tampered(labels, pi):
        kernel, mass = original(labels, pi)
        kernel[0][0] += Fraction(1, 1000)
        return kernel, mass

    monkeypatch.setattr(cert, "heat_bath", tampered)
    with pytest.raises(cert.ThermoError):
        cert.kernel_algebra_certificate()
