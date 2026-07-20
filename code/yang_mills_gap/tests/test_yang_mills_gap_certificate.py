from __future__ import annotations

import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("collar_gap_certificate", HERE / "collar_gap_certificate.py")
assert SPEC and SPEC.loader
CERTIFICATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CERTIFICATE)


def manifest(name: str) -> dict:
    return json.loads((HERE / "manifests" / name).read_text())


def test_atomic_calibration_has_the_announced_exact_floor() -> None:
    receipt = CERTIFICATE.validate(manifest("atomic_4d_ising_calibration.json"))
    assert receipt["active_type_count"] == 244
    assert receipt["c_floor"] == "1"
    assert receipt["eta_upper"] == "1/2"
    assert receipt["gap_lower"] == "1/2"
    assert receipt["physical_clay_receipt"] is False


def test_binary_float_input_is_rejected() -> None:
    payload = manifest("atomic_4d_ising_calibration.json")
    payload["calibration_type_family"]["template"]["rate_lower"] = 1.0
    try:
        CERTIFICATE.validate(payload)
    except ValueError as exc:
        assert "floating" in str(exc)
    else:
        raise AssertionError("float proof input was accepted")


def test_noncontractive_influence_is_rejected() -> None:
    payload = manifest("atomic_4d_ising_calibration.json")
    payload["calibration_type_family"]["template"]["influences"][0]["upper"] = "1"
    try:
        CERTIFICATE.validate(payload)
    except ValueError as exc:
        assert "< 1" in str(exc)
    else:
        raise AssertionError("noncontractive table was accepted")


def test_physical_placeholder_fails_closed() -> None:
    try:
        CERTIFICATE.validate(manifest("physical_compact_gauge_uninstantiated.json"))
    except ValueError as exc:
        assert "type_table" in str(exc)
    else:
        raise AssertionError("uninstantiated physical manifest was accepted")
