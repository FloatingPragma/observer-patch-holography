"""Check local continuous charge symmetry on written seams, and its blank control."""

import numpy as np
from scipy.linalg import expm

from . import law


def test_written_gauss_generators_on_all_local_columns():
    # 2E = 2(I-WW*) has integer entries. The written blocks are E G E.
    write = law._write_pattern()
    size = len(write)
    twice_complement = 2 * np.eye(size, dtype=int) - write @ write.T
    reversal, contrast = law._label_ops()
    identity = np.eye(law.N_SRC, dtype=int)
    twice_exchange = np.kron(np.eye(law.N_LAB, dtype=int) - reversal, identity)
    field_cost = np.kron(contrast, identity)
    flux = [(1 - 2*c) * (d % 2) for c, d in law.LABELS]
    left = np.array([2*a + f for f in flux for a, b in law.SRC])
    right = np.array([2*b - f for f in flux for a, b in law.SRC])
    counts = np.array([[word.count(load) for _ in law.LABELS for word in law.SRC] for load in range(3)])
    for middle in (twice_exchange, field_cost):
        block = twice_complement @ middle @ twice_complement
        for value in (left, right, *counts):
            assert np.all(block * (value[:, None] - value[None, :]) == 0)


def test_complete_load_and_memory_records_distinguish_uncompensated_kick():
    exchange, contrast = law.local_term_parts()
    frequency = exchange + 8 * contrast
    rng = np.random.default_rng(34921)
    for written in (False, True):
        labels = law.LABELS if written else [(0, 0)]
        sl = slice(law.N_SRC, None) if written else slice(None, law.N_SRC)
        free = expm(-0.7j * frequency[sl, sl])
        state = rng.normal(size=len(free)) + 1j*rng.normal(size=len(free))
        state /= np.linalg.norm(state)
        base = abs(free @ state)**2
        load = np.array([a for label in labels for a, b in law.SRC])
        flux = np.repeat([(0.5-c)*(d % 2) for c, d in labels], law.N_SRC)
        def distance(charge):
            measured = abs(free @ (np.exp(0.5j*np.pi*charge) * state))**2
            return np.abs(measured-base).sum()/2
        assert distance(load) > 1e-3
        if written:
            assert distance(load + flux) < 1e-12


def test_an_echo_reference_can_detect_the_defect_inventory_offset():
    """An addressed controlled echo arms a blank/written Ramsey comparator."""
    kappa = 8
    dwell = np.pi / 2
    exchange, contrast = law.local_term_parts()
    blank = np.zeros(law.DIM)
    blank[[law.SRC.index((0, 2)), law.SRC.index((2, 0))]] = 1/np.sqrt(2)
    echo = law.echo()
    written = echo @ blank
    inventory = np.array([a != 1 for a, b in law.SRC], int) + np.array([b != 1 for a, b in law.SRC], int)
    defects = np.tile(inventory, law.N_LAB + 1)
    total = np.tile([a+b for a,b in law.SRC], law.N_LAB + 1)
    for coupling, probability in [(0, 1), (1, 0)]:
        h = exchange + kappa*contrast + np.diag(coupling*defects + 2.3*total)
        free = expm(-1j*dwell*h)
        bright = np.linalg.norm(free @ blank + echo @ (free @ written))**2 / 4
        np.testing.assert_allclose(bright, probability, atol=1e-12)
