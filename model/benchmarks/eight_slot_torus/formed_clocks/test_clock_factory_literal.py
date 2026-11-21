"""Independent multi-origin channel and local pointer checks for the uniform loader."""
import math

import numpy as np
from scipy import sparse

from . import clock_factory
from .test_clock_loader_literal import _native_pair_step, _routing


def _fields(shape, coefficients):
    """Even corners of every translated plaquette, assembled from its signed slots."""
    origins = coefficients.reshape(shape)
    diagonal = np.roll(origins, (1, 1), axis=(0, 1))
    result = np.zeros((*shape, 8), complex)
    result[..., 0] = origins/2
    result[..., 2] = -origins/2
    result[..., 1] = diagonal/2
    result[..., 3] = -diagonal/2
    return result.ravel()


def _angles(age):
    return math.pi*(1-math.cos(math.pi*(age+1)/128))/4 if age < 128 else math.pi/2


def _prediction(shape, route, coefficients, rail, age):
    w = _fields(shape, coefficients)
    elapsed = age-35
    mu = (128+rail)*math.pi/4096
    phase = np.exp(-1j*(sum(map(_angles, range(35)))+mu))
    clock = (-np.exp(1j*mu))**elapsed
    parity = np.exp(-2j*sum(map(_angles, range(35, age))))
    return phase*clock*((1+parity)*w+(1-parity)*(route@w))/2


def test_all_origin_channel_with_untouched_reference():
    # Three reference states, complex spatial/rail entanglement, and a chart
    # including wrap-around plaquettes. Compare full output amplitudes.
    shape = (3, 4, 3, 3)
    volume = math.prod(shape)
    route = _routing(shape)
    random = np.random.default_rng(357401)
    coefficients = random.normal(size=(2, 3, volume))+1j*random.normal(size=(2, 3, volume))
    coefficients /= np.linalg.norm(coefficients)
    gates = clock_factory.uniform_schedule()
    assert 0 not in gates and 35 in gates and set(gates) <= set(range(1, 36))
    for g in gates.values():
        assert g.shape == (2, 8, 8)
        for block in g:
            np.testing.assert_allclose(block.conj().T@block, np.eye(8), atol=2e-13, rtol=0)
    for r in range(2):
        for reference in range(3):
            coeff = coefficients[r, reference]
            pair = sparse.diags(np.repeat(coeff/np.sqrt(8), 8), format="csr")
            for age in range(1, 132):
                pair = _native_pair_step(pair, route, age-1, r)
                if age in gates:
                    old = pair.diagonal().reshape(volume, 8)
                    updated = old@gates[age][r].T
                    pair += sparse.diags((updated-old).ravel(), format="csr")
                if age in (35, 36, 67, 128, 131):
                    actual = pair.diagonal()
                    assert sparse.linalg.norm(pair-sparse.diags(actual)) < 2e-13
                    np.testing.assert_allclose(actual, _prediction(shape,route,coeff,r,age),
                                               atol=2e-12, rtol=0)


def test_two_parities_are_complementary_on_coherent_origins():
    shape = (4, 4, 3, 3)
    volume = math.prod(shape)
    route = _routing(shape)
    basis = np.eye(volume)
    w = np.column_stack([_fields(shape, a) for a in basis])
    for x in ((0,0,0,0), (3,3,2,2)):
        i = np.ravel_multi_index(x, shape)
        np.testing.assert_allclose(clock_factory.packet(shape,x).ravel(), w[:,i], atol=2e-15)
    gram = w.conj().T@w
    np.testing.assert_allclose(gram, np.eye(volume), atol=1e-14)
    # Columns of (W +/- FW)/sqrt(2) are the normalized single-parity loops.
    for parity in (-1,1):
        loop = (w+parity*(route@w))/np.sqrt(2)
        rng = np.random.default_rng(35439)
        a = rng.normal(size=shape)+1j*rng.normal(size=shape)
        expected = a.copy()
        for axis in (0,1):
            expected += parity*(np.roll(a,1,axis)+np.roll(a,-1,axis))/4
        np.testing.assert_allclose(loop.conj().T@(loop@a.ravel()), expected.ravel(), atol=1e-13)
    plane = np.ones(volume)/np.sqrt(volume)
    w_plane = w@plane
    plus = (w+route@w)/2
    # The coherent input and the mixture of orthogonally marked births have
    # identical origin probabilities and different single-parity odds.
    assert abs(np.vdot(plus@plane, plus@plane).real-1) < 2e-13
    assert abs(np.vdot(plus,plus).real/volume-.5) < 2e-13
    checker = (-1.)**np.indices(shape)[:2].sum(axis=0)/np.sqrt(volume)
    assert np.linalg.norm(plus@checker.ravel()) < 1e-13
    assert abs(np.vdot(w_plane,w_plane).real-1) < 2e-13


def test_local_phase_after_forgetting_birthplace_and_time():
    shape = (3, 3, 3, 3)
    volume = math.prod(shape)
    route = _routing(shape)
    local = np.arange(8)
    # Mix every birthplace with its own marker record. At one receiving
    # address only the four attached roots have a possible local signal.
    roots = ((0,0,0,0), (2,0,0,0), (0,2,0,0), (2,2,0,0), (1,1,1,1))
    hands = []
    for age in (150, 4246):
        outcomes = np.zeros(4)
        for origin in roots:
            coeff = np.zeros(volume)
            coeff[np.ravel_multi_index(origin, shape)] = 1
            rails = [_prediction(shape, route, coeff, r, age)[local]/np.sqrt(2) for r in range(2)]
            for i, phase in enumerate((1,1j,-1,-1j)):
                outcomes[i] += np.linalg.norm(rails[0]+phase*rails[1])**2/(4*volume)
        assert abs(sum(outcomes)-1/volume) < 4e-14
        hand = (outcomes[0]-outcomes[2])+1j*(outcomes[1]-outcomes[3])
        assert abs(abs(hand)-1/(2*volume)) < 4e-14
        hands.append(hand)
    # Distinct ages 4096 apart have the same mature spatial effect and
    # opposite coherent rail hands. A wire longer than 4246 services can
    # retain both marker depths, even though this receiver discards them.
    assert abs(sum(hands)) < 5e-13


def test_late_native_birth_requires_the_age_condition():
    # A second distinct marker word is born after the first 35-step table
    # has run. Apply either its own age clock or the original wall clock.
    shape = (3, 3, 3, 3)
    route = _routing(shape)
    origin = np.zeros(math.prod(shape))
    origin[np.ravel_multi_index((2, 1, 2, 0), shape)] = 1
    gates = clock_factory.uniform_schedule()
    target = _fields(shape, origin)
    fidelities = []
    for offset in (0, 35):
        pair = sparse.diags(np.repeat(origin/np.sqrt(8), 8), format="csr")
        for age in range(1, 36):
            pair = _native_pair_step(pair, route, age-1, 0)
            if age+offset in gates:
                before = pair.diagonal().reshape(-1, 8)
                updated = before@gates[age+offset][0].T
                pair += sparse.diags((updated-before).ravel(), format="csr")
        assert abs(sparse.linalg.norm(pair)**2-1) < 3e-13
        fidelities.append(abs(np.vdot(target, pair.diagonal()))**2)
    assert abs(fidelities[0]-1) < 2e-12
    assert fidelities[1] < 2e-25
