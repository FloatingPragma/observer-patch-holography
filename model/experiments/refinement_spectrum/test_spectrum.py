"""Independent graph, many-body instrument, deep-root and negative receipt checks."""

import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from model.experiments.cage import construction as cell
from model.experiments.cage.instruments import fock_operators, lift
from . import flight, instruments as ins, limit, source
from tools.strict_json import load
from .__main__ import HERE, ROOT, compare_spectrum, score, verify_saved


@pytest.mark.parametrize('depth', [0, 1, 2, 3])
def test_inventory_includes_dark_copies_and_both_signs(depth):
    eig = np.linalg.eigvalsh(source.tree(depth).toarray()) / cell.HOP
    bands = source.inventory(depth)
    positive = np.concatenate([np.tile(b['positive'], b['copies']) for b in bands])
    np.testing.assert_allclose(eig, np.sort(np.r_[-positive, positive]), atol=2e-13)


@pytest.mark.parametrize('length', [2, 3, 9, 17, 33])
def test_deep_root_by_high_precision_determinant(length):
    # All other squared roots exceed 1, so the first sign change brackets the light root.
    with mp.workdps(85):
        low, high = mp.mpf(0), mp.mpf(1)
        for _ in range(245):
            x = (low+high)/2
            prev, curr = 1, 1-x
            for _ in range(2, length+1):
                prev, curr = curr, (5-x)*curr-4*prev
            if curr > 0:
                low = x
            else:
                high = x
        assert source.light_rate(length) == pytest.approx(float(mp.sqrt((low+high)/2)), rel=5e-14, abs=0.)
        assert source.inventory(length-1)[-1]['positive'][0] == pytest.approx(
            float(mp.sqrt((low+high)/2)), rel=5e-14, abs=0.)
    assert source.light_rate(length-1)/source.light_rate(length) < 8/3


@pytest.mark.parametrize('branching', [2, 3, 4])
@pytest.mark.parametrize('relative_hop', [.08, .5, 1., 3., 12.])
def test_free_hopping_ratio_bound_against_full_tree(branching, relative_hop):
    depth = 3
    nodes = sum(branching**k for k in range(depth+1))
    T = np.eye(nodes)
    for child in range(1,nodes):
        T[child,(child-1)//branching] = -relative_hop
    all_lines = np.sort(np.linalg.svd(T,compute_uv=False))
    g = relative_hop*np.sqrt(branching)
    minima = []
    complete = []
    for n in range(1,depth+2):
        rates = np.linalg.svd(np.eye(n)-g*np.eye(n,k=-1),compute_uv=False)
        copies = 1 if n == depth+1 else (branching-1)*branching**(depth-n)
        complete.extend(np.tile(rates,copies))
        minima.append(rates[-1])
    np.testing.assert_allclose(all_lines, np.sort(complete), atol=4e-14)
    ratios = np.array(minima[:-1])/minima[1:]
    assert np.all((g < ratios) & (ratios <= g+1))
    distinct = all_lines[np.r_[True,np.diff(all_lines)>1e-9]]
    if g >= 2:
        np.testing.assert_allclose(distinct[:3],minima[-3:][::-1],atol=4e-14)
    else:
        assert distinct[1]/distinct[0] < 3


def test_recorded_ratios_obey_the_bounds_of_the_source():
    light = load(HERE/'records'/'result.json')['result']['light']
    first, second = light[1]/light[0], light[2]/light[1]
    # Successive chain minima at g = 2: above g, at most g + 1, and below 8/3.
    assert 2 < first < 8/3 and 2 < second < 8/3 and abs(first-second) < 1
    assert light[2]/light[0] < 64/9


def test_shallow_family_including_large_hopping():
    for g in (0., .125, 2., 250.):
        # Two children of the root, in the original node basis.
        rate = -g/np.sqrt(2)
        T = np.array([[1.,0.,0.],[rate,1.,0.],[rate,0.,1.]])
        u, middle, heavy = np.sort(np.linalg.svd(T,compute_uv=False))
        assert middle == pytest.approx(1.)
        assert heavy*u == pytest.approx(1.)
        assert middle/u == pytest.approx(heavy/middle)


@pytest.mark.parametrize('law', ['graded', 'no_children', 'unsigned'])
def test_full_labeled_cell_evolves_as_prepared_blocks(law):
    car = cell.cell()
    full = source.microscopic(2, car, law)
    assert full.shape == (6300, 6300)
    assert full.nnz//2 == (21 if law == 'no_children' else 41)*150 + 42*240
    for profile in ('uniform', 'ring', 'port'):
        u, _, p = source.outer_profile(car, profile)
        for n in (1, 2, 3):
            v = np.kron(source.seed(2,n), np.r_[u, np.zeros(120)])
            t = 1.21
            direct = np.vdot(v, expm_multiply(-1j*t*full, v, traceA=0.))
            block = source.reduced(2, n, p, law)
            assert direct == pytest.approx(expm(-1j*t*block)[0,0], abs=2e-13)


def test_many_body_interferometer_reads_both_fills_and_signed_quadratures():
    a, b, fuel = ops = fock_operators(3)
    number = fuel.conj().T @ fuel
    v = (a+1j*b)/np.sqrt(2)
    split = expm(-.25j*np.pi*(v.conj().T@fuel + fuel.conj().T@v))
    h = np.array([[.7, .2j], [-.2j, -.4]])
    t = 2.3
    free = expm(-1j*t*lift(h,ops[:2]))
    vec = np.array([1.,-1j])/np.sqrt(2)
    amp = np.vdot(vec, expm(-1j*t*h) @ vec)
    law = ins.fuel_law(np.array([amp]), (0,2,1,3))[0]
    result = []
    for q in (0,2,1,3):
        program = split.conj().T @ expm(.5j*np.pi*q*number) @ free @ split
        readout = np.diag(program.conj().T @ number @ program).real
        result.append([np.mean(readout[:4]), np.mean(readout[4:])])
    np.testing.assert_allclose(result, law, atol=4e-15)
    delta = np.diff(result,axis=1).ravel()
    assert delta[0]-delta[1]+1j*(delta[2]-delta[3]) == pytest.approx(amp, abs=2e-15)


def test_clock_absolute_counts_distinguish_occupied_reference():
    t = .9
    clock = ins.fuel_law(ins.clock_amplitude(np.array([t])),(0,2,1,3))[0]
    # In the dimer eigenbasis the positive orbital is a spectator. Evolve the
    # negative orbital and the fuel, before averaging the two source occupations.
    orbital, fuel = fock_operators(2)
    number = fuel.conj().T @ fuel
    split = expm(-.25j*np.pi*(orbital.conj().T@fuel + fuel.conj().T@orbital))
    free = expm(1j*t*orbital.conj().T@orbital)
    for i,q in enumerate((0,2,1,3)):
        U = split.conj().T @ expm(.5j*q*np.pi*number) @ free @ split
        by_fill = np.diag(U.conj().T @ number @ U).real.reshape(2,2)
        np.testing.assert_allclose(by_fill.mean(axis=1), clock[i], atol=1e-14)
        assert not np.allclose(by_fill[:,1], clock[i])


def test_full_federation_filter_and_receiver_partitions():
    spec = load(HERE/'protocol.json')
    menu = spec['flight']
    saved = load(HERE/'records'/'flight.json')
    with pytest.raises(ValueError, match='histories'):
        flight.assess(spec, dict(saved,histories=[]))
    with pytest.raises(ValueError, match='histories'):
        flight.assess(spec, dict(saved,histories=[saved['histories'][0]]*3))
    car = cell.federation(menu['federation'])
    sender, neighbor = saved['sender_cell'], saved['receiver_cell']
    assert any({i//12,j//12} == {sender, neighbor} for i,j in car['edges'])
    seam = tuple(saved['seam_seed'])
    assert seam[0]//12 == seam[1]//12 == sender
    mode = car['edges'].index(seam)
    owners = np.array([i//12 if i//12 == j//12 else -1 for i,j in car['edges']]
                      + [center//12 for _,center,_ in car['wedges']])
    masks = [owners == sender, owners == neighbor, owners == -1,
             ~np.isin(owners,[sender,neighbor,-1])]
    full = source.microscopic(2,car,'graded') / cell.HOP
    assert full.shape == (saved['source_modes'],)*2 == (139860,139860)
    assert full.nnz == 2*363330
    internal = source.tree(2).toarray() / cell.HOP
    evolution = expm(-1j*menu['filter_clock_interval']*internal)
    for n,history in enumerate(saved['histories'],1):
        # Literal A-mode addresses on the breadth-first tree, before any reduction.
        v = np.zeros(42,complex)
        if n == 3:
            v[0] = 1
        else:
            first = 5 if n == 1 else 1
            v[first:first+2] = [1/np.sqrt(2),-1/np.sqrt(2)]
        failures = []
        for _ in range(menu['filter_rounds']):
            moved = evolution @ v
            failures.append(np.linalg.norm((v-moved)/2)**2)
            v = (v+moved)/2
        state = np.zeros((42,len(owners)),complex)
        state[:,mode] = v
        wave = state.ravel()
        previous = 0
        for i,t in enumerate(menu['waits']):
            if t > previous:
                wave = expm_multiply((-1j*(t-previous))*full,wave,traceA=0.)
            weights = np.sum(abs(wave.reshape(state.shape))**2,axis=0)
            law = [*failures, *(weights[mask].sum() for mask in masks)]
            np.testing.assert_allclose(law, history['laws'][i], atol=2e-12, rtol=0)
            previous = t


def test_limit_short_time_bound_and_failed_controls():
    H = source.chain(4,5)/cell.HOP
    for t in (.1, 1., 2.):
        err = abs(expm(-1j*t*H)[0,0]-limit.limiting_return(t))
        assert err <= limit.moment_tail(5,t)+2e-14
    with pytest.raises(ValueError,match='reference contrast'):
        ins.clock_read(np.full((4,4,2),100),200)
    with np.load(HERE/'records'/'counts.npz',allow_pickle=False) as z:
        clock = z['clock'].copy()
    clock[15] = clock[14]
    with pytest.raises(ValueError,match='lost a turn'):
        ins.clock_read(clock,2**20)
    with pytest.raises(ValueError,match='visible spectrum'):
        compare_spectrum(dict(lines=[.19394], dc=0), [.19394,1.90321,2.70927],
                         dict(maximum_frequency_error=.003))
    with pytest.raises(ValueError,match='wrong value'):
        verify_saved(dict(low=0.),dict(low=1.74623e-10))
    with pytest.raises(ValueError,match='wrong value'):
        verify_saved(dict(positive=[0.]),dict(positive=[1.74623e-10]))


def test_short_window_uses_the_observers_clock():
    spec = load(HERE/'limit_protocol.json')
    saved = load(HERE/'records'/'limit'/'result.json')
    # Choose a boundary just below the first nonzero *recorded* tick. The first
    # simulator step lies inside it; the observer has only the zero-tick sample.
    with np.load(HERE/'records'/'limit'/'counts.npz',allow_pickle=False) as z:
        ticks = ins.clock_read(z['clock'],spec['shots'])
        assert ticks[1] > spec['clock_step']
        end = (ticks[1]+spec['clock_step'])/2
        limiting = np.array([limit.limiting_return(t) for t in ticks])
        expected = {}
        for n in spec['chain_lengths']:
            observed = ins.quadrature(z[str(n)],spec['shots'])
            err = abs(observed-limiting)
            expected[n] = max(err[ticks <= end])
            published = next(r for r in saved if r['length'] == n)
            assert max(err[ticks <= spec['finite_window_end']]) == pytest.approx(
                published['atoms_versus_limit_over_window'],abs=2e-12)
    _, rows = limit.run(dict(spec,finite_window_end=end))
    for row in rows:
        assert row['atoms_versus_limit_over_window'] == pytest.approx(expected[row['length']],abs=2e-12)


def test_heldout_counts_and_weak_mode_weights_are_checked():
    spec = load(HERE/'protocol.json')
    with np.load(HERE/'records'/'counts.npz', allow_pickle=False) as archive:
        ticks = ins.clock_read(archive['clock'], spec['shots'])
        counts = archive['graded_uniform_3'].copy()
    counts[spec['fit_samples']:] = 0
    with pytest.raises(ValueError, match='held-out'):
        ins.read_lines(ticks, ins.quadrature(counts,spec['shots']), spec)
    results = load(HERE/'records'/'result.json')
    weights = results['reads']['graded_uniform_3']['weights']
    weights[0], weights[-1] = weights[-1], weights[0]
    with pytest.raises(ValueError, match='weights'):
        score(spec,results['carriers'],results['reads'])


@pytest.mark.parametrize('suffix', ['', '.limit'])
def test_public_replay_with_optimization(suffix):
    subprocess.run([sys.executable,'-B','-O','-m','model.experiments.refinement_spectrum'+suffix],
                   check=True, cwd=ROOT)
