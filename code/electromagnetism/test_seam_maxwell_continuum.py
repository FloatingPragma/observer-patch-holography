"""Adversarial controls for the continuum field bridge and its custody."""
import copy
import json
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import seam_maxwell_continuum as producer
import verify_seam_maxwell_continuum as verifier


@pytest.fixture(scope="module")
def geometry():
    return verifier.exact_geometry()


@pytest.fixture(scope="module")
def receipt():
    return json.loads(producer.OUTPUT.read_text())


def test_complete_receipt_replays_independently(receipt, geometry):
    assert verifier.verify(receipt, geometry=geometry) == 120


@pytest.mark.parametrize("mutation", ["omitted_case", "promoted_status", "wrong_bound",
                                      "missing_assumption", "forged_result", "nan",
                                      "changed_source", "missing_custody"])
def test_false_green_receipts_fail(receipt, geometry, mutation):
    value = copy.deepcopy(receipt)
    if mutation == "omitted_case":
        value["cases"].pop()
    elif mutation == "promoted_status":
        value["status"] = "PHYSICAL_MAXWELL_ESTABLISHED"
    elif mutation == "wrong_bound":
        value["universal_bounds"]["frequency"] = "zero error"
    elif mutation == "missing_assumption":
        value["assumptions"].pop(0)
    elif mutation == "forged_result":
        value["cases"][-1]["propagator_error"] = 0.0
    elif mutation == "nan":
        value["cases"][-1]["forcing_error"] = float("nan")
    elif mutation == "changed_source":
        value["source"]["sha256"] = "0" * 64
    else:
        value["implementation"].pop("Lean/Screen/SeamMaxwellContinuum.lean")
    with pytest.raises(ValueError):
        verifier.verify(value, geometry=geometry)


def test_actual_producer_matrices_match_exponentials_and_forcing():
    # Off-axis and infrared data; compare matrices, not only norms that could
    # hide a simultaneous orientation mistake in two propagators.
    for values in [(0, 0, 0), (1e-6, -2e-6, 3e-6), (1, 2, -2), (12, -8, 4)]:
        k = np.array(values, float)
        _, p, j = producer.blocks(k)
        for a in (0.125, 1.0):
            omega = np.sqrt(producer.symbol(a, k))
            g, _ = verifier.generator(k, omega)
            for t in (-2.0, 0.0, 3.0):
                assert np.allclose(producer.propagator(omega, t, p, j), expm(t * g), atol=1e-12)
                aug = np.block([[g, np.eye(6)], [np.zeros((6, 12))]])
                assert np.allclose(producer.constant_forcing(omega, t, p, j),
                                   expm(t * aug)[:6, 6:], atol=1e-12)


def test_wrong_sign_and_wrong_normalization_are_discriminated():
    k = np.array([1.0, 2.0, -2.0])
    g, p = verifier.generator(k, np.sqrt(producer.symbol(0.25, k)))
    wrong_sign = g.copy()
    wrong_sign[3:, :3] *= -1
    assert np.linalg.norm(expm(wrong_sign).conj().T @ expm(wrong_sign) - np.eye(6)) > 1
    assert np.linalg.norm(g @ g + producer.symbol(0.25, k) * p) < 1e-12
    assert np.linalg.norm((2 * g) @ (2 * g) + producer.symbol(0.25, k) * p) > 1


def test_longitudinal_forcing_and_charge_continuity():
    k = np.array([1.0, 2.0, -2.0])
    _, p, j = producer.blocks(k)
    omega = np.sqrt(producer.symbol(0.5, k))
    current = np.array([2.0, 1.0, -1.0], complex)
    initial = np.array([1, 0, 2, 2, -1, 0], complex)
    force = np.concatenate([-current, np.zeros(3)])
    t = 1.25
    field = (producer.propagator(omega, t, p, j) @ initial +
             producer.constant_forcing(omega, t, p, j) @ force)
    rho = 1j * k @ initial[:3] - t * (1j * k @ current)
    assert abs(1j * k @ field[:3] - rho) < 1e-12
    assert abs(1j * k @ field[3:]) < 1e-12
    # A charge history without the continuity term fails the same Gauss test.
    assert abs(1j * k @ field[:3] - 1j * k @ initial[:3]) > 1


def test_real_pairing_is_required():
    k = np.array([1.0, 2.0, -2.0])
    omega = np.sqrt(producer.symbol(0.5, k))
    _, p, j = producer.blocks(k)
    _, pn, jn = producer.blocks(-k)
    u, un = producer.propagator(omega, 0.3, p, j), producer.propagator(omega, 0.3, pn, jn)
    assert np.linalg.norm(un - u.conj()) < 1e-12
    assert np.linalg.norm(u - u.conj()) > 0.1  # wrong even sign for the curl fails


@pytest.mark.parametrize("scale", [0.0, -1.0, float("inf"), float("nan")])
def test_invalid_refinement_scale_rejected(scale):
    with pytest.raises(ValueError):
        producer.symbol(scale, np.array([1.0, 0.0, 0.0]))


@pytest.mark.parametrize("text", ['{"status":1,"status":2}', '{"value":NaN}'])
def test_noncanonical_json_rejected(tmp_path, text):
    path = tmp_path / "receipt.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        verifier.load_receipt(path)
