"""Audit the loader's exported instructions with two explicit constituent coordinates."""
import math

import numpy as np
import pytest
from scipy import sparse

from . import clock_loader


def _arrival(shape, site, slot):
    other = list(site)
    axis = slot // 2
    other[axis] = (other[axis] + (1 if slot % 2 == 0 else -1)) % shape[axis]
    return np.ravel_multi_index(tuple(other), shape) * 8 + (slot ^ 1)


def _routing(shape):
    rows = [_arrival(shape, site, slot) for site in np.ndindex(shape) for slot in range(8)]
    n = len(rows)
    return sparse.csr_matrix((np.ones(n), (rows, np.arange(n))), shape=(n, n))


def _shift(shape, site, axis, distance):
    target = list(site)
    target[axis] = (target[axis] + distance) % shape[axis]
    return tuple(target)


def _launch(shape, root, route, focus_a=1., focus_b=1., angle=math.pi / 4):
    """Independent inward trains on two sites, then the physical root rotation."""
    size = math.prod(shape)
    equal = np.full(8, math.sqrt(1 / 8))
    coin = 2 * np.outer(equal, equal) - np.eye(8)
    ordinary = sparse.kron(sparse.eye(size), sparse.csr_matrix(coin), format='csr')
    there = _shift(shape, root, 0, 1)
    trains = []
    for slot, f, offset in ((4, focus_a, -1), (5, focus_b, 1)):
        arr = np.zeros((size, 8))
        for site, probability in ((root, f), (there, 1 - f)):
            i = np.ravel_multi_index(_shift(shape, site, 2, offset), shape)
            arr[i] = math.sqrt(probability) * coin[:, slot]
        trains.append(sparse.csr_matrix(arr.reshape(-1, 1)))
    assert trains[0].multiply(trains[1]).nnz == 0
    incoming = trains[0] @ trains[1].T
    advanced = route @ ordinary @ incoming @ ordinary.T @ route.T
    # Each head-on coefficient writes its own root. A record identifies the wire independently of rail.
    vacuum = advanced.tolil()
    marked = {}
    for site in (root, there):
        z = np.ravel_multi_index(site, shape)
        x = advanced[8 * z + 5, 8 * z + 4]
        vacuum[8 * z + 5, 8 * z + 4] = math.cos(angle) * x
        diagonal = np.zeros(8 * size, complex)
        diagonal[8 * z:8 * z + 8] = math.sin(angle) * x * equal
        marked[site] = sparse.diags(diagonal, format='csr')
    refusals = {'other': float(sparse.linalg.norm(marked[there])**2),
                'unmarked': float(sparse.linalg.norm(vacuum)**2)}
    return refusals, marked[root]


def _loop(shape, root):
    loop = np.zeros((math.prod(shape), 8))
    x = _shift(shape, root, 0, 1)
    y = _shift(shape, root, 1, 1)
    xy = _shift(shape, x, 1, 1)
    for start, a, end, b, sign in ((root, 0, x, 1, 1), (x, 2, xy, 3, -1), (xy, 1, y, 0, 1), (y, 3, root, 2, -1)):
        loop[np.ravel_multi_index(start, shape), a] = sign / math.sqrt(8)
        loop[np.ravel_multi_index(end, shape), b] = sign / math.sqrt(8)
    return loop.ravel()


def _local_gates(pair, shape, instructions, rail):
    # The ordinary and unmatched entries of the two-particle matrix are untouched.
    diag = pair.diagonal()
    before = diag.copy().reshape(-1, 8)
    after = before.copy()
    for site, blocks in instructions.items():
        i = np.ravel_multi_index(site, shape)
        gate = blocks[rail]
        np.testing.assert_allclose(gate.conj().T @ gate, np.eye(8), atol=2e-13, rtol=0)
        after[i] = gate @ before[i]
    return pair + sparse.diags(after.ravel() - diag, format='csr')


def _native_pair_step(pair, route, age, rail):
    # After a first mark the pair is supported on the matched diagonal of the full matrix.
    # Retain both endpoint indices when routing the partial-hop amplitude.
    diag = pair.diagonal().reshape(-1, 8)
    assert np.allclose(sparse.linalg.norm(pair)**2, np.vdot(diag, diag).real, rtol=0, atol=3e-13)
    mu = (math.pi / 16 + rail * math.pi / 2048) / 2
    bright = np.mean(diag, axis=1, keepdims=True)
    result = np.exp(-1j * mu) * bright - np.exp(1j * mu) * (diag - bright)
    coined = pair + sparse.diags((result - diag).ravel(), format='csr')
    routed = route @ coined @ route.T
    theta = math.pi * (1 - math.cos(math.pi * (age + 1) / 128)) / 4 if age < 128 else math.pi / 2
    return ((1 + np.exp(-2j * theta)) * coined + (1 - np.exp(-2j * theta)) * routed) * .5


def _run(shape, root, instructions, end_age):
    route = _routing(shape)
    refusal, born = _launch(shape, root, route)
    assert math.isclose(refusal['unmarked'], .5, abs_tol=2e-14) and refusal['other'] == 0
    states = []
    for rail in range(2):
        pair = born
        for age in range(end_age + 1):
            if age in instructions:
                pair = _local_gates(pair, shape, instructions[age], rail)
            if age < end_age:
                pair = _native_pair_step(pair, route, age, rail)
        states.append(pair.diagonal())
    return np.stack(states)


PANELS = (((5, 5, 3, 3), (2, 2, 1, 1)), ((6, 7, 4, 5), (5, 6, 3, 4)))


@pytest.mark.parametrize('shape,root', PANELS)
def test_both_literal_constituents_and_the_original_wire_refusals(shape, root):
    route = _routing(shape)
    for a, b, angle in ((1., 1., math.pi / 4), (.5, .75, math.pi / 4), (0., 1., math.pi / 4), (.5, .75, 0.)):
        raw, born = _launch(shape, root, route, a, b, angle)
        raw['natal'] = float(sparse.linalg.norm(born)**2)
        native = clock_loader.ordinary_launch(shape, root, a, b, angle)
        assert set(native) == set(raw)
        for key in native:
            assert abs(native[key] - raw[key]) < 2e-14
        assert abs(sum(raw.values()) - 1) < 2e-14


@pytest.mark.parametrize('shape,root', PANELS)
@pytest.mark.parametrize('mode', ('exact', 'common', 'phase-cut', 'passive'))
def test_literal_carrier_replays_loaded_and_controlled_rail_fields(shape, root, mode):
    gates = clock_loader.schedule(shape, root, mode)
    assert all(1 <= age <= 28 for age in gates)
    actual = _run(shape, root, gates, 28)
    native, _ = clock_loader.prepare_clock(shape, root, mode)
    np.testing.assert_allclose(actual, native.reshape(2, -1) / math.sqrt(2), atol=2e-12, rtol=0)
    loop = _loop(shape, root)
    overlaps = actual @ loop
    if mode == 'passive':
        assert np.max(abs(overlaps)) < 2e-13
    if mode != 'exact':
        return
    np.testing.assert_allclose(actual, overlaps[:, None] * loop, atol=2e-12, rtol=0)
    mu = np.array([math.pi / 32, math.pi / 32 + math.pi / 4096])
    # Entanglement fidelity of the heralded channel against the known natal rail rotation.
    choi_overlap = np.sum(np.exp(1j * mu) * overlaps) / math.sqrt(2)
    assert abs(abs(choi_overlap)**2 - 1) < 2e-12
    later = _run(shape, root, gates, 133)
    for s in (root, _shift(shape, _shift(shape, root, 0, 1), 1, 1)):
        i = 8 * np.ravel_multi_index(s, shape)
        a, b = later[:, i:i + 8] / math.sqrt(2)
        ports = [np.vdot(a + h * b, a + h * b).real / 4 for h in (1, 1j, -1, -1j)]
        assert abs(sum(ports) - 1 / 8) < 2e-12
        hand = complex(ports[0] - ports[2], ports[1] - ports[3])
        desired = np.exp(1j * (1 - (133 - 28)) * (mu[1] - mu[0]))
        assert abs(hand / abs(hand) - desired) < 2e-11
