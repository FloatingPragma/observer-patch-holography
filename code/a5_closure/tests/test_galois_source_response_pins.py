"""Fail-closed controls for the audited source pins of the source-response replay."""
import hashlib
import shutil
from pathlib import Path
import sys
import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import galois_source_response_control as control


def test_pins_bind_the_current_bytes():
    assert set(control.SOURCE_PINS) == {
        "DERIVATION.md", "geometry.py", "response.json", "verify_response.py"}
    assert {name: hashlib.sha256((control.MODEL / name).read_bytes()).hexdigest()
            for name in control.SOURCE_PINS} == control.SOURCE_PINS


@pytest.mark.parametrize("name", sorted(control.SOURCE_PINS))
def test_any_stale_pin_rejects_the_replay(monkeypatch, name):
    monkeypatch.setitem(control.SOURCE_PINS, name, "0" * 64)
    with pytest.raises(ValueError, match="audited source drift"):
        control.verify()


def test_derivation_text_is_audited_but_never_a_replay_input(monkeypatch, tmp_path):
    # The prose is bound by digest only. Rewriting it while keeping the pins
    # consistent must leave every executed clause count unchanged.
    for name in control.SOURCE_PINS:
        shutil.copy(control.MODEL / name, tmp_path / name)
    (tmp_path / "DERIVATION.md").write_text("# replacement prose\n")
    monkeypatch.setattr(control, "MODEL", tmp_path)
    monkeypatch.setattr(control, "SOURCE_PINS", {
        name: hashlib.sha256((tmp_path / name).read_bytes()).hexdigest()
        for name in control.SOURCE_PINS})
    result = control.verify()
    assert result["plus"] == result["minus"]
    assert result["plus"]["naturality_pairs"] == 720
    assert result["plus"]["compositions"] == 3600
