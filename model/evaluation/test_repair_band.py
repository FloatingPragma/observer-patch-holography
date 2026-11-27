"""The cell Gram frame must be carried by the actual registered operator."""
import json
from pathlib import Path

import numpy as np
import pytest

from model.executable_model import cell, laws
from . import repair_band


def test_archive_covers_the_registered_response_of_every_candidate():
    measured = repair_band.run()
    saved = json.loads(Path(repair_band.__file__).with_suffix(".json").read_text())
    assert set(measured) == set(saved) == set(laws.names())
    for name, row in measured.items():
        archived = saved[name].copy()
        errors = archived.pop("frobenius_errors")
        assert archived == {key: value for key, value in row.items() if key != "frobenius_errors"}
        assert row["lifted_rank"] == 3
        assert row["frobenius_errors"] == pytest.approx(errors, rel=0, abs=1e-12)


def test_wrong_slow_projector_fails_even_with_the_same_spectrum(monkeypatch):
    candidate = laws.load("null_hopping")
    kinetic = cell().laplacian().toarray()/2
    perm = [1, 0, *range(2, 12)]
    other = kinetic[np.ix_(perm, perm)]
    np.testing.assert_allclose(np.linalg.eigvalsh(other), np.linalg.eigvalsh(kinetic), atol=1e-14)
    monkeypatch.setattr(candidate, "generator", lambda _: other)
    monkeypatch.setattr(candidate, "response", lambda _: other)
    with pytest.raises(ArithmeticError, match="incompatible generator or Gram frame"):
        repair_band.measure(candidate)
    monkeypatch.setattr(candidate, "response", lambda _: kinetic)
    with pytest.raises(ArithmeticError, match="declared response differs"):
        repair_band.measure(candidate)
