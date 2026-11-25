"""Check a free, overlapping-contact evolution and sensitive readout controls."""

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from . import collar_storage as storage


def test_memory_survives_an_actual_mobile_neighbour():
    law = candidate.native_law
    hop, contrast = law.local_term_parts()
    h = np.kron(hop + 8 * contrast, np.eye(3))
    # Three material ports (a,b,c), with written collar on ab and blank collar on bc.
    # bc acts on the same port b and changes the material at the stored record.
    h += np.kron(np.eye(33), law.seam_term(8))
    signed = np.zeros((11, 2))
    for coin, difference in law.LABELS:
        if abs(difference) == 2:
            row = 1 + law.LABELS.index((coin, difference))
            signed[row] = [0.5, 0.5 * np.sign(difference)]
    loads = np.eye(27)[:, [np.ravel_multi_index((1, 1, 2), (3, 3, 3))]]
    wave = (expm(-2j * h) @ np.kron(signed, loads)).reshape(11, 27, 2)
    for bit in range(2):
        material = signed[:, bit] @ wave[:, :, bit]
        assert np.linalg.norm(wave[:, :, bit] - signed[:, bit, None] * material) < 1e-12
        assert np.sum(abs(material.reshape(3, 3, 3)[:, 2, :])**2) > 1e-6


def test_disturbance_controls_at_each_scratch_word():
    for exchange, encoded, expected in ((1., True, 0.), (1., False, 1.),
                                        (.5, True, .25), (.5, False, .5)):
        data = storage.recovered(exchange, encoded, 1)
        for q in range(9):
            for bit in range(2):
                target_error = 0 if q == storage.PAIR else expected
                assert abs(data[storage.BITS[1-bit], q, bit] - target_error) < 3e-13
                assert abs(data[storage.BITS[bit], q, bit] - (1-target_error)) < 3e-13


def test_station_geometry_and_deposit():
    assert len({tuple(station["pointer"]) for station in storage.geometry()}) == 4
    assert storage.conserved_columns() == [
        [q, 0, 0] if q == 4 else [q, 32, 33] for q in range(9)]
