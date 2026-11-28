"""Read the same two exterior ports with an ordinary mobile source present."""

import numpy as np
from scipy.sparse.linalg import expm_multiply

from openwave.interval_echo.row_19.selected_collar import midpoint_clocks as clocks


def compare(adjacency, lap, ports, times):
    psi, _ = clocks.force.compact_mode(lap)
    q = 2  # A distinct load label lets the observer detect just the signal.
    contact = clocks.force.meter.contact_terms()[q]
    h, x, y = clocks.force.meter.two_loads(lap*0.5, adjacency, contact, q)
    initial = np.column_stack([psi[x]*(y == p) for p in ports])
    blank = np.eye(lap.shape[0])[:, ports]
    np.testing.assert_allclose(initial.T @ initial, np.eye(2), rtol=0, atol=1e-14)
    results = {}
    for name, op, first, detector in [
        ("ordinary", h, initial, y), ("empty", lap*0.5, blank, np.arange(lap.shape[0]))
    ]:
        gen = -1j*op
        results[name] = {}
        state = first
        last = 0
        for t in times:
            step = gen*(t-last)
            state = expm_multiply(step, state, traceA=step.diagonal().sum())
            square = abs(state)**2
            np.testing.assert_allclose(square.sum(axis=0), 1., atol=2e-12, rtol=0)
            # Rows = emitter, columns = detector. Sum out every final ordinary position.
            results[name][str(t)] = np.round(np.stack(
                [square[detector == p].sum(axis=0) for p in ports], axis=1), 12).tolist()
            last = t
    return results
