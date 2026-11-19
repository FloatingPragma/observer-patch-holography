"""Independent graph, many-body instrument, deep-root and negative receipt checks."""

import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from ..cage import construction as cell
from ..cage.instruments import fock_operators, lift
from . import instruments as ins, limit, source
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
    assert source.light_rate(length-1)/source.light_rate(length) < 8/3


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
    subprocess.run([sys.executable,'-B','-O','-m','model.scalar_channel.refinement_spectrum'+suffix],
                   check=True, cwd=ROOT)
