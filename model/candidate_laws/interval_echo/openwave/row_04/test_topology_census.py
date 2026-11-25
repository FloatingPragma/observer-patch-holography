import json
import subprocess
import sys

import numpy as np
import pytest

from model.executable_model import federation
from . import open_pair
from . import topology_census as census


def test_all_cells_and_all_layouts_match_the_archived_exact_census():
    result = census.census()
    assert result == json.loads(census.ARCHIVE.read_text(encoding="utf-8"))
    assert [(b["eligible_cells"], b["cells"], b["certified_layouts"], b["eligible_layouts"])
            for b in result["boards"]] == [(20, 20, 104, 104), (80, 80, 384, 384), (320, 320, 1504, 1504)]
    assert all(not b["ineligible_cells"] for b in result["boards"])
    assert result["sensitivity_controls"] == {"bright_sign": 1992, "single_contact_unseal": 21912}


def test_integer_columns_agree_with_the_full_28680_state_native_sector():
    board = federation(0)
    written = set(open_pair.geometry()[0])
    sector, kinetic, contrast, _, _ = open_pair.full_parts()
    for state in ({(156, 157): 1}, {(156, 166): 1}, {(156, 157): 1, (156, 166): -1},
                  {(156, 157): 1, (156, 166): 1}):
        vector = np.zeros(sector.P, dtype=int)
        for word, amplitude in state.items():
            vector[sector.index(*word)] = amplitude
        integers = census.action(board, written, state)
        for actual, expected in zip(integers, (2 * kinetic @ vector, contrast * vector), strict=True):
            reconstructed = np.zeros(sector.P, dtype=int)
            for word, amplitude in actual.items():
                reconstructed[sector.index(*word)] = amplitude
            np.testing.assert_array_equal(reconstructed, expected)


def test_each_unsealed_contact_reaches_a_missing_full_sector_destination():
    board = federation(0)
    written = set(open_pair.geometry()[0])
    dark = {(156, 157): 1, (156, 166): -1}
    protected = census.action(board, written, dark)
    for contact in written:
        changed = census.action(board, written - {contact}, dark)
        assert changed != protected
        assert set(changed[0]) - dark.keys()


def test_same_certificate_rejects_a_changed_native_binary_rule(monkeypatch):
    law = census.candidate.native_law
    kinetic, contrast = law.local_term_parts()
    changed = kinetic.copy()
    changed[law.SRC.index((0, 1)), law.SRC.index((0, 1))] += .125
    monkeypatch.setattr(law, "local_term_parts", lambda: (changed, contrast))
    with pytest.raises(ArithmeticError, match="Native blank binary column"):
        census.native_blocks()


def test_replay_failure_is_nonzero_and_preserves_the_archive(tmp_path):
    # Isolate custody and CLI failure from the scientific census above.
    archive = tmp_path / "certificate.json"
    original = b'{"wrong":true}\n'
    archive.write_bytes(original)
    code = ("from pathlib import Path; from model.candidate_laws.interval_echo.openwave.row_04 "
            "import topology_census as c; c.ARCHIVE=Path(__import__('sys').argv[1]); "
            "c.census=lambda: {'right': True}; __import__('sys').argv=['census']; c.main()")
    result = subprocess.run([sys.executable, "-c", code, str(archive)], capture_output=True, text=True)
    assert result.returncode != 0 and "differs from its archived certificate" in result.stderr
    assert archive.read_bytes() == original
