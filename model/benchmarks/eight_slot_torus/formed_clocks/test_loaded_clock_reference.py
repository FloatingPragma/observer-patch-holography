"""Check the complete local exchange instrument against an explicit stationary memory."""
import math

import numpy as np
from scipy import sparse

from . import clock_loader, loaded_clock
from .test_clock_loader_literal import _launch, _local_gates, _native_pair_step, _routing


def _memory_effects():
    # Joint basis (rail,memory) = 00, 01, 10, 11. Only equal-total-charge words mix.
    charge = np.diag([0, 1, 1, 2])
    exchange = np.eye(4)
    exchange[1:3, 1:3] = np.array([[1., 1.], [1., -1.]])/math.sqrt(2)
    np.testing.assert_array_equal(charge @ exchange, exchange @ charge)
    effects = []
    for h in (1, 1j, -1, -1j):
        effect = np.eye(4, dtype=complex)/4
        effect[1, 2], effect[2, 1] = h/4, h.conjugate()/4
        np.testing.assert_array_equal(charge @ effect, effect @ charge)
        assert np.linalg.eigvalsh(effect).min() >= -1e-15
        effects.append(effect)
    np.testing.assert_allclose(np.sum(effects, axis=0), np.eye(4), atol=1e-15, rtol=0)
    return exchange, effects


def test_natal_memory_and_same_root_measurement_on_both_full_coordinates():
    shape, root = (5, 5, 3, 3), (2, 2, 1, 1)
    route = _routing(shape)
    _, pair = _launch(shape, root, route)
    gate, effects = _memory_effects()
    i = 8*np.ravel_multi_index(root, shape)
    mask = np.zeros(route.shape[0])
    mask[i:i+8] = 1
    # 01 initially carries the Born amplitude. 00 and 11 have zero support throughout.
    words = [pair, sparse.csr_matrix(pair.shape, dtype=complex)]
    schedule = clock_loader.schedule(shape, root)
    native, _ = clock_loader.prepare_clock(shape, root, local_reference=True)
    for age in range(29):
        if age == 6:
            local = [sparse.diags(w.diagonal()*mask, format='csr') for w in words]
            words = [w-local[r]+sum(gate[r+1,s+1]*m for s,m in enumerate(local))
                     for r,w in enumerate(words)]
        if age in schedule:
            words = [_local_gates(w, shape, schedule[age], r) for r,w in enumerate(words)]
        if age < 28:
            words = [_native_pair_step(w, route, age, r) for r,w in enumerate(words)]
    native = native.reshape(2,-1)/math.sqrt(2)
    for r,w in enumerate(words):
        np.testing.assert_allclose(w.diagonal(), native[r], atol=2e-12, rtol=0)
    # Each terminal joint effect acts only when the two material endpoints reach the memory.
    for age in (28, 29, 128):
        if age == 29:
            words = [_native_pair_step(w, route, 28, r) for r,w in enumerate(words)]
        if age == 128:
            for j in range(29,128):
                words = [_native_pair_step(w, route, j, r) for r,w in enumerate(words)]
        local = np.zeros((4,8), complex)
        local[1:3] = np.stack([w.diagonal()[i:i+8] for w in words])
        rho = local @ local.conj().T
        measured = np.array([np.trace(e @ rho).real for e in effects])
        np.testing.assert_allclose(measured, loaded_clock.local_ports(local[1:3]), atol=2e-13, rtol=0)
        assert abs(measured.sum()-1/8) < 2e-12
        hand = complex(measured[0]-measured[2], measured[1]-measured[3])
        expected = np.exp(-1j*(age-28)*math.pi/4096)/16
        assert abs(hand-expected) < 2e-12
        reduced = np.trace(rho.reshape(2,2,2,2), axis1=1, axis2=3)
        assert reduced[0,1] == 0 and reduced[1,0] == 0
        accessible = np.array([np.trace(e @ np.diag(np.diag(rho))).real for e in effects])
        np.testing.assert_allclose(accessible, loaded_clock.local_ports(local[1:3], dephased=True),
                                   atol=2e-13, rtol=0)
        np.testing.assert_allclose(accessible, np.full(4,1/32), atol=2e-12, rtol=0)
