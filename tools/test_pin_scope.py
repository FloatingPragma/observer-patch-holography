"""The pin rule: prose, configuration, tests and the rule file are never pinned."""
import hashlib

import pytest

from tools import pin_scope


@pytest.mark.parametrize("path", [
    "docs/research/X.json", "paper/a.tex", "model/x/README.md", "model/x/PROOF.md",
    "model/x/test_packet.py", "conftest.py", ".github/workflows/ci.yml", "tools/pin_scope.py",
    "code/pin_scope.py", "./model/x/README", "model\\x\\notes.rst",
])
def test_exempt_paths(path):
    assert pin_scope.exempt(path)


@pytest.mark.parametrize("path", [
    "model/x/packet.py", "model/x/results.json", "model/x/witness.npz", "tools/replay_model_receipts.py",
    "code/m1_source_realization/receipt.json", "model/requirements.txt", "Lean/X.lean",
])
def test_bound_paths(path):
    assert not pin_scope.exempt(path)


def test_bound_pins_filters_maps_and_leaves_other_values():
    pins = {"model/x/packet.py": "a", "model/x/README.md": "b", "model/x/test_x.py": "c"}
    assert pin_scope.bound_pins(pins) == {"model/x/packet.py": "a"}
    assert pin_scope.bound(list(pins)) == ["model/x/packet.py"]
    assert pin_scope.bound_pins(["model/x/README.md"]) == ["model/x/README.md"]
    assert pin_scope.bound_pins(None) is None


def test_repo_sha256_reads_the_current_file():
    root = pin_scope.__file__ and __import__("pathlib").Path(__file__).resolve().parents[1]
    expected = hashlib.sha256((root / "model/requirements.txt").read_bytes()).hexdigest()
    assert pin_scope.repo_sha256("model/requirements.txt") == expected
