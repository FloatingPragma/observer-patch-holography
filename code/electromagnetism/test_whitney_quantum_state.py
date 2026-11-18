"""Falsification tests for the selected interacting quantum initial state."""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
import os
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_whitney_quantum_state as verifier
import whitney_quantum_state as producer


@pytest.fixture(scope="module")
def receipt():
    return verifier.load()


def change(packet, keys, value):
    row = packet
    for key in keys[:-1]:
        row = row[key]
    row[keys[-1]] = value


def test_receipt_replays_independently(receipt):
    result = verifier.verify(receipt)
    assert result["accepted"] is True
    assert result["gaussian_sigma"] == "1/2"
    assert result["real_configuration_dimension"] == 56
    assert result["gaussian_norm_squared"] == "1"
    assert result["matter_l2_coefficient"] == "1/5"
    assert result["matter_l4_coefficient"] == "3/35"
    assert result["volume_in_Qsqrt5"] == ["10", "10/3"]
    assert result["initial_state_constructed"] is True
    assert result["initial_observables_provided"] is True
    assert result["quantum_time_history"] is False
    assert result["observer_history"] is False
    assert result["physical_comparison"] is False
    assert result["analytic_operator_domain_proved_by_numeric_replay"] is False
    assert max(row["max_quadrature_replay_difference"] for row in result["numeric_diagnostics"]) < 2e-10


def test_producer_matches_exact_evidence_and_replays(receipt):
    fresh = producer.build()
    assert fresh["source_pins"] == receipt["source_pins"]
    assert fresh["exact_observables"] == receipt["exact_observables"]
    assert verifier.verify(fresh)["accepted"] is True


def test_simplex_cross_terms_and_wick_recurrence():
    # Omitting the six mixed pairs or treating a complex scalar as real gives
    # a different fourth moment, even if a producer were changed coherently.
    assert verifier.standard_gaussian_even_moment(0) == 1
    assert verifier.standard_gaussian_even_moment(6) == 15
    assert verifier.simplex_moment((0, 0, 0, 0)) == 1
    assert verifier.simplex_moment((2, 2, 0, 0)) == Q(1, 210)
    simplex_fourth = 4*Q(1, 35)+12*Q(1, 210)
    assert simplex_fourth == Q(6, 35)
    assert 8*Q(1, 2)**4*simplex_fourth == Q(3, 35)
    assert 3*Q(1, 2)**4*simplex_fourth != Q(3, 35)


@pytest.mark.parametrize("keys,value", [
    (("state", "real_configuration_dimension"), 55),
    (("state", "real_configuration_dimension"), 56.0),
    (("state", "complex_matter_coordinates"), True),
    (("state", "rho_power"), "-1"),
    (("state", "gaussian_prefactor_power"), "-28"),
    (("state", "gaussian_exponent_denominator_factor"), "2"),
    (("state", "quantum_time_history"), True),
    (("state", "quantum_time_history"), 0),
    (("state", "observer_history"), True),
    (("state", "physical_comparison"), True),
    (("state", "empirical_prediction"), True),
    (("state", "initial_state_constructed"), 1),
    (("state", "initial_observables_provided"), False),
    (("parameters", "sigma"), "1"),
    (("parameters", "charge"), .25),
    (("exact_observables", "normalized_gaussian_norm_squared"), "2"),
    (("exact_observables", "complex_gaussian_wick_factors", "fourth"), "3"),
    (("exact_observables", "coefficients_per_volume", "matter_l4"), "2/35"),
    (("exact_observables", "coefficients_per_volume", "quartic_potential"), "3/140"),
    (("exact_observables", "volume_in_Qsqrt5"), ["10", "10/6"]),
    (("numerical_amplitudes", "scope"), "exact amplitudes with certified error"),
    (("numerical_amplitudes", "samples", 1, "coordinates_exact", 0), "1/2"),
    (("numerical_amplitudes", "samples", 1, "coordinates_exact", 30), .5),
    (("numerical_amplitudes", "samples", 0, "evaluations", 0, "quadrature_order"), 4.0),
    (("numerical_amplitudes", "samples", 0, "evaluations", 0, "quadrature_order"), True),
    (("numerical_amplitudes", "samples", 0, "log_g_numeric"), True),
    (("numerical_amplitudes", "samples", 0, "evaluations", 0, "log_rho_numeric"), float("nan")),
    (("source_pins",), []),
    (("numerical_amplitudes",), []),
    (("numerical_amplitudes", "samples", 0), []),
    (("numerical_amplitudes", "samples", 0, "evaluations", 0), []),
])
def test_scientific_and_type_mutations_rejected(receipt, keys, value):
    packet = deepcopy(receipt)
    change(packet, keys, value)
    with pytest.raises(ValueError):
        verifier.verify(packet)


def test_arbitrarily_small_coherent_moment_corruption_rejected(receipt):
    packet = deepcopy(receipt)
    wrong = Q(1, 5)+Q(1, 10**1000)
    evidence = packet["exact_observables"]
    evidence["coefficients_per_volume"]["matter_l2"] = str(wrong)
    evidence["expectations_in_Qsqrt5"]["matter_l2"] = [str(10*wrong), str(Q(10, 3)*wrong)]
    with pytest.raises(ValueError, match="exact Gaussian"):
        verifier.verify(packet)


def test_wrong_half_density_numeric_samples_rejected(receipt):
    packet = deepcopy(receipt)
    # Corrupt all orders coherently; comparing quadrature orders alone would
    # accept this wrong rho^{-1} amplitude.
    for row in packet["numerical_amplitudes"]["samples"]:
        for sample in row["evaluations"]:
            sample["log_f_numeric"] = row["log_g_numeric"]-sample["log_rho_numeric"]
    with pytest.raises(ValueError, match="independent simplex"):
        verifier.verify(packet)


def test_small_density_corruption_rejected(receipt):
    packet = deepcopy(receipt)
    for sample in packet["numerical_amplitudes"]["samples"][2]["evaluations"]:
        sample["log_rho_numeric"] += 1e-6
    with pytest.raises(ValueError, match="independent simplex"):
        verifier.verify(packet)


@pytest.mark.parametrize("mode", ["missing", "wrong_hash", "changed_parent_bytes"])
def test_pins_replayed_fresh(receipt, monkeypatch, mode):
    packet = deepcopy(receipt)
    relative = "code/electromagnetism/whitney_quantum_state.py"
    if mode == "missing":
        del packet["source_pins"][relative]
    elif mode == "wrong_hash":
        packet["source_pins"][relative] = "0"*64
    else:
        original = Path.read_bytes
        def changed(path):
            data = original(path)
            return data+b"\n" if path == verifier.ROOT/relative else data
        monkeypatch.setattr(Path, "read_bytes", changed)
    with pytest.raises(ValueError, match="source pin"):
        verifier.verify(packet)


@pytest.mark.parametrize("fragment", [
    '{"quantum_time_history":false,"quantum_time_history":true}',
    '{"state":{"rho_power":"-1/2","rho_power":"-1"}}',
    '{"log_f_numeric":NaN}', '{"log_f_numeric":Infinity}',
    '{"log_f_numeric":-Infinity}',
    '{"log_f_numeric":1e9999}', '{"log_f_numeric":-1e9999}',
])
def test_ambiguous_or_nonfinite_json_rejected(tmp_path, fragment):
    path = tmp_path/"bad.json"
    path.write_text(fragment, encoding="utf-8")
    with pytest.raises(ValueError):
        verifier.load(path)


@pytest.mark.parametrize("value", [True, 1, .5, "2/4", "0.5", "1/0"])
def test_exact_rationals_require_canonical_strings(value):
    with pytest.raises(ValueError):
        verifier.rational(value)


def test_loader_and_source_reads_use_explicit_utf8(receipt, monkeypatch):
    original = Path.read_text
    def checked(path, *args, **kwargs):
        assert kwargs.get("encoding") == "utf-8"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", checked)
    assert verifier.verify(verifier.load())["accepted"] is True


def test_verifier_filespec_import_from_unrelated_directory(tmp_path):
    # The postdiction builder loads verifiers by absolute file spec; it must
    # work without adding the producer directory to sys.path.
    code = "\n".join((
        "import importlib.util", "import sys",
        "spec=importlib.util.spec_from_file_location('initial_state_replay',sys.argv[1])",
        "module=importlib.util.module_from_spec(spec)", "spec.loader.exec_module(module)",
        "assert module.verify(module.load())['accepted'] is True",
        "assert 'whitney_quantum_state' not in sys.modules",
        "assert 'whitney_interacting_quantum' not in sys.modules",
    ))
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    subprocess.run([sys.executable, "-c", code, str(Path(verifier.__file__).resolve())],
                   cwd=tmp_path, env=env, check=True, capture_output=True, text=True)
