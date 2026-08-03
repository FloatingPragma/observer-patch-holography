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
        "Lean/Thermodynamics/FirstLawIdentity.lean",
        "Lean/EventAlgebra/PartitionPinchingCP.lean",
    ):
        assert (root / path).exists()
    bindings = receipt["lean_bindings"]
    assert "heatBath_secondLaw" in bindings["conditional_repair"]
    assert "firstLaw_split" in bindings["first_law"]
    assert "kraus_complete" in bindings["record_channel"]


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
