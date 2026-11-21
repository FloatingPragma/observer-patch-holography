"""Independent physical counter, first-wire census and position-space terminal readout."""
import math
from decimal import Decimal, localcontext
from fractions import Fraction

import numpy as np
import pytest
from scipy import sparse

from . import clock_comparison as comparison
from . import clock_factory, clock_inference
from .test_clock_loader_literal import _routing, _shift


COIN = np.ones((8, 8))/4 - np.eye(8)


@pytest.mark.parametrize("period", (2, 4))
def test_coupled_optical_pointer_and_which_private_position(period):
    shape = (3, 3, 3, period)
    station = (1, 0, 2, period-1)
    route = _routing(shape)
    step = route @ sparse.kron(sparse.eye(math.prod(shape)), COIN, format="csr")
    cursor = np.zeros((*shape, 8))
    cursor[(*station[:3], slice(0,None,2), 6)] = 1/math.sqrt(period)
    cursor[(*station[:3], slice(0,None,2), 7)] = -1/math.sqrt(period)
    mask = np.zeros_like(cursor, dtype=bool)
    mask[(*station[:3], slice(0,None,2))] = True
    mask = mask.ravel()
    # The second coordinate is a quantum pointer. Increasing it on the
    # even subspace uses the same null record for every private site.
    wave = np.zeros((cursor.size, 8))
    wave[:,0] = cursor.ravel()
    expected = {}
    for t in range(1, 11):
        wave = step @ wave
        wave[mask] = np.roll(wave[mask], 1, axis=1)
        power = (wave**2).sum(axis=0)
        np.testing.assert_allclose(power, np.eye(8)[t//2], atol=2e-13, rtol=0)
        if t % 2 == 0:
            expected[t//2] = t
    triggers, _ = comparison.optical_counter(shape,station,tuple(expected))
    assert triggers == expected
    missed, _ = comparison.optical_counter(shape,station,tuple(expected),tag="displaced")
    assert not missed
    at_return = (step @ (step @ cursor.ravel())).reshape(cursor.shape)
    # Write which private point on this return and keep its orthogonal
    # record. The next even incidence sums over these distinct histories.
    next_return = 0.
    for q in range(0, period, 2):
        branch = np.zeros_like(cursor)
        branch[(*station[:3], q)] = at_return[(*station[:3], q)]
        two = step @ (step @ branch.ravel())
        next_return += np.vdot(two[mask],two[mask]).real
    assert abs(next_return - (1 if period == 2 else 5/8)) < 2e-13


def _first_wires(shape, focus, stop):
    """Evolve the two independent feeds, attenuating each ordered (5,4) entrance."""
    volume = math.prod(shape)
    root = (0,)*4
    route = _routing(shape)
    ordinary = route @ sparse.kron(sparse.eye(volume), COIN, format="csr")
    feeds = []
    for weight, slot, distance in ((focus[0],4,-1), (focus[1],5,1)):
        source = np.zeros((volume,8))
        for mouth, p in ((root,weight),(_shift(shape,root,0,1),1-weight)):
            inlet = np.ravel_multi_index(_shift(shape,mouth,2,distance),shape)
            source[inlet] = math.sqrt(p)*COIN[:,slot]
        feeds.append(source.ravel())
    pair = np.outer(*feeds)
    starts = 8*np.arange(volume)
    blocks = starts[:,None]+np.arange(8)
    row, col = blocks[:,:,None], blocks[:,None,:]
    writes = []
    for _ in range(stop):
        collision = pair[starts+5, starts+4].copy()
        writes.append(collision**2/2)
        pair[starts+5, starts+4] = collision/math.sqrt(2)
        onsite = pair[row,col].copy()
        # The matched subspace is invariant and empty in the unmarked
        # branch. Unmatched coincident slots keep their input before routing.
        assert np.max(abs(onsite.diagonal(axis1=1,axis2=2))) < 2e-14
        pair[row,col] = 0
        far = ordinary @ pair @ ordinary.T
        pair.fill(0)
        pair[row,col] = onsite
        pair = far + route @ pair @ route.T
    return np.stack(writes), np.vdot(pair,pair).real


@pytest.mark.parametrize("fi", (0,1))
def test_all_original_cohort_and_four_receiver_words_by_direct_propagation(fi):
    shape, focus = comparison.FIXTURES[fi]
    sites = comparison.stations(shape)
    appointments = np.array(comparison.RETURNS)*2
    history, refused = _first_wires(shape,focus,16)
    q = history.sum()
    assert abs(q+refused-1) < 3e-13
    assert abs(history[1,0]-np.prod(focus)/2) < 3e-14
    assert history[2].sum() < 1e-28 and history[3:].sum() > .01
    origins = np.array(list(np.ndindex(shape)))
    displacements = [np.ravel_multi_index(((site-origins) % shape).T, shape) for site in sites]
    bookings = {}
    for j,t in enumerate(appointments):
        for b in range(len(history)):
            bookings.setdefault(t-b,[]).append((j,b))
    perm = _routing(shape).indices
    phases = np.exp(1j*np.array([128,129])*math.pi/4096)[:,None,None]
    expected = np.zeros((3,4,len(appointments),7))
    for menu,gates in enumerate((clock_factory.uniform_schedule(),{})):
        wave = np.zeros((2,math.prod(shape),8),complex)
        wave[:,0,:] = .25
        for age in range(max(bookings)+1):
            if age in gates:
                for r in (0,1):
                    wave[r] = wave[r] @ gates[age][r].T
            if age in bookings:
                local = np.stack([np.sum(abs(wave[0]+h*wave[1])**2,axis=1)/4
                                  for h in (1,1j,-1,-1j)],axis=1)
                for j,b in bookings[age]:
                    for s,offset in enumerate(displacements):
                        expected[menu,s,j,:4] += history[b] @ local[offset]
            if age == max(bookings):
                break
            bright = wave.mean(axis=-1,keepdims=True)
            coined = phases.conj()*bright-phases*(wave-bright)
            routed = coined.reshape(2,-1)[:,perm].reshape(wave.shape)
            if age < 127:
                theta = math.pi*(1-math.cos(math.pi*(age+1)/128))/4
                e = np.exp(-2j*theta)
                wave = ((1+e)*coined+(1-e)*routed)/2
            else:
                wave = routed
    expected[2,...,:4] = expected[0,...,:4].mean(axis=-1,keepdims=True)
    expected[...,4] = q-expected[...,:4].sum(axis=-1)
    expected[...,5] = refused
    with np.load(comparison.ARCHIVE,allow_pickle=False) as data:
        np.testing.assert_array_equal(data["ticks"][fi],np.broadcast_to(appointments,(4,len(appointments))))
        np.testing.assert_allclose(data["weights"][:,fi]/(1<<38),expected,atol=2e-11,rtol=0)


def test_registered_elapsed_decodes_also_follow_count_only_hoeffding_disks():
    # Independent conservative disk bound conditional on each original
    # four-port total. Original-start cell coverage supplies a lower bound
    # for incidence when translating the source/grid error to this disk.
    with np.load(comparison.ARCHIVE,allow_pickle=False) as data:
        cells = data["counts"][0].reshape(8,9,7)
    n = cells[...,:4].sum(axis=-1)
    K = cells.sum(axis=-1)
    family = 216
    cell_bound = np.sqrt(np.log(4*family*7/1e-6)/(2*K))
    bias = 4*np.sqrt(2)*2e-11/(n/K-4*cell_bound)
    r = np.sqrt(16*np.log(8*family/1e-6)/n) + bias
    centers = 2*((cells[...,0]-cells[...,2])+1j*(cells[...,1]-cells[...,3]))/n
    arcs = np.arcsin(np.minimum(1.,r/abs(centers)))
    arcs[r>=abs(centers)] = np.pi
    candidates = np.arange(0,8192,2)
    for j in range(1,9):
        relative = np.angle(centers[:,j]*centers[:,0].conj())
        error = abs(np.angle(np.exp(1j*(relative[:,None]+candidates*np.pi/4096))))
        kept = candidates[np.all(error <= (arcs[:,0]+arcs[:,j])[:,None],axis=0)]
        expected = 2*(comparison.RETURNS[j]-comparison.RETURNS[0])
        np.testing.assert_array_equal(kept,[expected])
        # A common permutation of the receiver labels costs no phase prior.
        oracle = clock_inference.decode_elapsed(cells[::-1,0],cells[::-1,j],rate=-math.pi/4096,
                                                family_size=family,cell_error=2e-11)
        assert oracle["even_runs_step_two"] == [[expected,expected]]


def test_extreme_binomial_tails_and_unresolved_rare_receivers():
    tail = Fraction(1, 6*216*10**6)
    with localcontext() as ctx:
        ctx.prec = 85
        for n in (2**36, 2**63-1, 2**64+1):
            for k in (1, n-1, n):
                p = Decimal.from_float(clock_inference._lower(k,n,tail))
                if k == 1:
                    exact = 1-(1-p)**n
                elif k == n-1:
                    exact = p**(n-1)*(p+n*(1-p))
                else:
                    exact = p**n
                assert 0 <= exact <= Decimal(tail.numerator)/tail.denominator
    # The native original-start bias can dominate a seemingly precise
    # conditional direction when almost every launch missed the receiver.
    sparse_word = np.array([[10**6, 5*10**5, 0, 5*10**5, 2**62]], dtype=np.int64)
    decoded = clock_inference.decode_elapsed(sparse_word,sparse_word,rate=math.pi/4096,cell_error=2e-11)
    assert decoded["even_runs_step_two"] == [[0,8190]] and decoded["informative_stations"] == 0
    assert decoded["original_start_count_total"] == 2*(2**62+2*10**6)
