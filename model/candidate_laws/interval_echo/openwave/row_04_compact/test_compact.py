import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import dark_bands, gateway_table, geometry, port_reader, prepare


def test_four_port_source_and_local_occupation_reader():
    local, axes = geometry()
    car = federation(0)
    tips = gateway_table(car)[0]
    available = [pair for pair in axes if set(tips).isdisjoint(pair)]
    target = np.zeros(12)
    target[list(available[0])] = 0.5
    target[list(available[1])] = -0.5
    source, _ = prepare(local, tips, target)
    for port in range(12):
        _, effect = port_reader(local, port)
        np.testing.assert_allclose(
            effect @ abs(source) ** 2, target[port] ** 2, atol=1e-12
        )
    action = candidate.generator(car)[:, :12] @ source
    np.testing.assert_allclose(action[:12], 51 * source, atol=1e-12)
    np.testing.assert_allclose(action[12:], 0, atol=1e-12)


def test_internal_beat_requires_two_distinct_dark_energies():
    local_h = candidate.generator(cell()).toarray()
    rich = dark_bands(local_h, [0, 1, 2])
    ordinary = dark_bands(local_h, [1, 3, 10])
    assert len(ordinary) == 1
    assert len(rich) == 3
    np.testing.assert_allclose(rich[-1][0] - rich[0][0], np.sqrt(5), atol=1e-12)
