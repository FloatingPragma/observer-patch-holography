import numpy as np
import pytest

from . import clock_ballistic as run


def test_saddles_reject_three_competitors_on_both_rails():
    for i,name in enumerate(("low","high")):
        q = run.aperture(name)
        assert abs(np.linalg.norm(q)-1) < 1e-13
        for mu in run.MUS:
            k = run.saddles(mu)
            group = np.cos(mu)*np.sin(k)/(4*np.sin(run.energy(k,mu)))
            assert np.max(abs(group-run.V)) < 2e-13
            response = lambda x: abs(np.vdot(q,np.exp(1j*run.OFFSETS*x)))
            assert response(k[i]) > 1
            assert max(response(x) for x in [-k[0],-k[1],k[1-i]]) < 1e-11


def test_one_address_cannot_annihilate_all_unwanted_occupied_slot_directions():
    # Use the actual eight slot coordinates of each native mature mode.
    signs = np.array([0,0,0,0,1,-1,0,0])
    for mu in run.MUS:
        lo,hi = run.saddles(mu)
        def eig(k,s):
            return np.exp(-1j*mu+1j*signs*k)-np.exp(1j*s*run.energy(k,mu))
        nuisance = np.column_stack([eig(-lo,-1),eig(hi,1),eig(-hi,-1)])
        for k,s in [(-lo,-1),(hi,1),(-hi,-1)]:
            assert abs(run.opening(np.array([k]),mu)[(1-s)//2][0]) > 1e-6
        assert np.linalg.matrix_rank(nuisance,tol=1e-10) == 3
        assert np.linalg.norm(nuisance@np.linalg.lstsq(nuisance,eig(lo,1),rcond=None)[0]
                              -eig(lo,1)) < 2e-13


def test_phase_intervals_handle_wrap_and_uninformative_words():
    assert run._increment(None,(-1.,0.)) is None
    arc = run._increment((2.5,2.6),(-2.7,-2.6))
    assert 0 < arc[0] < arc[1] < np.pi
    assert run._distance_bound(arc,np.mean(arc)-2*np.pi) < .11
    assert run._distance_bound(None,0.) == np.pi
    with pytest.raises(ValueError):
        run.amplitudes(24,[1],[0])


def test_retained_moving_observer_counts():
    result = run.experiment()
    assert result["counted_launches"] == 28*2**34
    assert set(result["records"]) == set(run.MENUS)
    for menu in run.MENUS:
        assert len(result["records"][menu]["advance_intervals_radians"]) == 6
    # Preserve the registered negative result in the authenticated acquisition.
    assert result["failed_gates"] == [["high","conditional_amplitude"]]
    assert result["records"]["low"]["mature_waits_outside_amplitude_gate"] == []
    assert result["records"]["high"]["mature_waits_outside_amplitude_gate"] == run.WAITS.tolist()
