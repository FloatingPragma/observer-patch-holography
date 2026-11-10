"""The theorem crosswalk and its declared source cannot silently lose files."""
import pytest
from . import certificate


def test_local_theorem_dependencies_and_toolchain_are_retained():
    hashes=certificate.input_hashes()
    assert set(certificate.LEAN_INPUTS)<=hashes.keys()
    assert all(len(hashes[name])==64 for name in certificate.LEAN_INPUTS)


def test_missing_theorem_input_fails_before_export(monkeypatch):
    monkeypatch.setattr(certificate,'LEAN_INPUTS',certificate.LEAN_INPUTS+('Lean/QFT/absent_dependency.lean',))
    with pytest.raises(ValueError,match='missing two-cell source dependency'):
        certificate.input_hashes()
