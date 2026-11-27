"""Independent small-spin checks for static contact bins and the spectral bound."""

from itertools import combinations, product

import numpy as np

from .apparatus import Guide
from .emitter import reflected
from .release import LEVELS, calibrate, local_bins, low_weight_certificate
from .test_emission import spins, vacuum


def test_every_physical_rail_pair_including_nonneighbour_chords():
    n = 4
    g = Guide(n,'native',.5)
    _, _, gamma, _ = spins(n,.5)
    host = vacuum(g,gamma)
    parent = np.einsum('i,ijk->jk',g.modes[:,0].real,gamma) @ host
    for state, cov in [(host,g.c),(parent,reflected(g.c,g.modes[:,0].real))]:
        physical = np.zeros(2**(2*n),complex)
        for idx, bits in enumerate(product((0,1),repeat=n)):
            occupied = [v for b in bits for v in (1-b,b)]
            physical[int(''.join(map(str,occupied)),2)] = state[idx]
        for u, v in combinations(range(2*n),2):
            pair = np.moveaxis(physical.reshape([2]*(2*n)),(u,v),(0,1)).reshape(4,-1)
            bins = [np.linalg.norm(pair[[0,3]])**2,0,0,
                    np.linalg.norm(pair[1]+pair[2])**2/2, np.linalg.norm(pair[1]-pair[2])**2/2]
            np.testing.assert_allclose(local_bins(cov,u,v),bins,atol=4e-13)
    instrument = calibrate()
    np.testing.assert_array_equal(instrument['(0, 1)'],[0,.5,.5,0,0])
    np.testing.assert_array_equal(instrument['(0, 2)'],[0,0,0,.5,.5])
    assert LEVELS @ instrument['balanced +1'] == 32
    assert LEVELS @ instrument['balanced -1'] == 33


def test_true_static_ternary_spectrum_with_shared_exits_and_unreleased_host_control():
    # Two dimers, one external unit; the two inter-dimer seams reach the same wrong pair.
    n = 2
    edges = [(0,1),(2,3),(1,2),(0,2),(0,4),(1,4)]
    basis = [tuple(w) for w in product(range(3),repeat=5) if sorted(w)==[0,0,1,2,2]]
    row = {w:i for i,w in enumerate(basis)}
    data = [i for i,w in enumerate(basis) if w[4]==1 and all(w[j]+w[j+1]==2 for j in (0,2))]
    h = np.zeros((len(row),len(row)))
    for i,w in enumerate(basis):
        for u,v in edges:
            if w[u] == w[v]:
                continue
            h[i,i] += 8*(w[u]-w[v])**2 + .5
            rev = list(w)
            rev[u],rev[v] = rev[v],rev[u]
            h[row[tuple(rev)],i] -= .5
    other = sorted(set(range(len(row)))-set(data))
    block = h[np.ix_(other,data)]
    np.testing.assert_allclose(np.linalg.norm(block,2)**2,1.,atol=1e-13)
    assert np.linalg.norm(block,2)**2 <= (2**2+2)/4
    # For this deliberately short, double-crossed sample the fixed cut is uninformative.
    b = dict(contacts=[dict(path_positions=[u,v]) for u,v in edges if v<4],
             unchanged_data_host_contacts=2)
    with np.testing.assert_raises(ValueError):
        low_weight_certificate(n,b)
    # Re-test the complete separation inequality on an exactly diagonalizable pure path.
    n = 4
    words = [w for w in product((0,2),repeat=2*n) if sum(w)==2*n]
    row = {w:i for i,w in enumerate(words)}
    data = [i for i,w in enumerate(words) if all(w[j]!=w[j+1] for j in range(0,2*n,2))]
    h = np.zeros((len(row),len(row)))
    for i,w in enumerate(words):
        for j in range(2*n-1):
            if w[j] != w[j+1]:
                h[i,i] += 32.5
                v = list(w)
                v[j],v[j+1] = v[j+1],v[j]
                h[row[tuple(v)],i] -= .5
    values, vec = np.linalg.eigh(h)
    b = dict(contacts=[dict(path_positions=[j,j+1]) for j in range(2*n-1)],
             unchanged_data_host_contacts=0)
    certified = low_weight_certificate(n,b)
    weight = np.linalg.norm(vec[data][:,values<=values[0]+.1],2)**2
    assert 0 < weight < certified['low_weight_upper']


def test_physical_mode_reader_on_true_spin_amplitudes_and_restore():
    from scipy.linalg import expm

    from .mode_reader import analyzer

    n = 4
    g = Guide(n,'native',.5)
    xs, zs, gamma, _ = spins(n,.5)
    vac = vacuum(g,gamma)
    _, bits, holds = analyzer(g)
    circuit = np.eye(2**n, dtype=complex)
    for axis, dwell in holds:
        j = axis//2
        generator = (65/4)*zs[j] @ zs[j+1] if axis%2 else -xs[j]/2
        circuit = expm(-1j*generator*dwell) @ circuit
    for target_mode in range(-1,n):
        state = vac if target_mode < 0 else np.einsum('i,ijk->jk',g.modes[:,target_mode].real,gamma) @ vac
        after = circuit @ state
        selected = after
        for j, x in enumerate(xs):
            sign = bits[j]*(-1 if j==target_mode else 1)
            selected = (selected+sign*x @ selected)/2
        np.testing.assert_allclose(circuit.conj().T @ selected,state,atol=2e-13)
