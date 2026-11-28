"""Regressions for the local SWAP and the meaning of the spectator records."""

import numpy as np

from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.fixture import (
    FanPath,
    archived_fixture,
)
from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.local_launch.measure import (
    fresh_write,
    terminal,
)


def test_swapping_in_a_fresh_qubit_traces_the_old_code_but_keeps_spectators():
    with archived_fixture() as archived:
        path = FanPath(archived, 3, 1, True)
    source = path.blocks[0, 0]
    stock = np.zeros_like(path.source())
    data = stock[source.rows].reshape(4, -1, 2)
    # Correlate two distinct joint spectator words with orthogonal old
    # code directions even in the same stock input. The other stock is refused.
    data[:, 0, 0] = path.code[:, 0] / np.sqrt(2)
    data[:, -1, 0] = path.code[:, 1] / np.sqrt(2)
    fresh, probability = fresh_write(path, stock)
    np.testing.assert_allclose(probability, 0.5)
    decoded = path.decode(fresh, 0).reshape(-1, 2, 4, 2)
    spectator_density = np.zeros((len(decoded), len(decoded)))
    spectator_density[0, 0] = spectator_density[-1, -1] = 0.25
    for i in range(2):
        np.testing.assert_allclose(decoded[:, 1 - i, :, i], 0, atol=1e-15)
        factor = decoded[:, i, :, i]
        np.testing.assert_allclose(factor @ factor.T, spectator_density, atol=1e-15)
    np.testing.assert_array_equal(path.decode(fresh), 0)


def test_which_single_station_must_report_when_its_record_marks_the_input():
    with archived_fixture() as archived:
        path = FanPath(archived, 3, 1, True)
    block = path.blocks[2, 0]
    wave = np.zeros_like(path.source())
    # Leave one perfectly distinguishing word in the first spectator while
    # the middle observer is already in its vacuum on both input branches.
    for i in range(2):
        one_word = np.eye(block.shape[0])[i]
        quiet = np.ones(block.shape[1]) / np.sqrt(block.shape[1])
        wave[block.rows, i] = np.kron(one_word, np.kron(quiet, path.code[:, i]))
    read = terminal(path, wave, 1.0, 1.0, 0.0)
    source_address, middle_address = [cell[0] for cell in path.cells[:2]]
    assert read["inclusion_minimal_entangled_comparisons"] == [[source_address]]
    for entry in read["comparisons"]:
        if source_address in entry["vacuum_observers"]:
            np.testing.assert_allclose(entry["conditional_negativity"], 0.5)
        elif entry["vacuum_observers"] in ([], [middle_address]):
            np.testing.assert_allclose(entry["conditional_negativity"], 0, atol=1e-15)
