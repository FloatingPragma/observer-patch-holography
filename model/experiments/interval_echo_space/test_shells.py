"""Independent source, measurement, alias and degeneracy checks for the clocked flight."""
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from .shell_source import Board, local_source, prepare
from .shells import hands, protocol, read_hand
from .shell_spectroscopy import analyzer, first_line
from .shell_beats import response, window_read


@pytest.mark.parametrize("amount",[1,2])
def test_overlapping_clock_and_exchange(amount):
    # Build the entire image on three physical load ports and two independent seam memories.
    # Both terms touch the middle port; all escaping basis coefficients are retained.
    h0,h1 = law.local_term_parts()
    h = h0+8*h1
    sign = [1+law.LABELS.index((0,s)) for s in (2,-2)]
    basis = [(a,b,c,m,0) for m in sign for a,b,c in
             [(amount,0,0),(0,amount,0),(0,0,amount)]]
    joint = np.zeros((6,6))
    for col,(a,b,c,m,n) in enumerate(basis):
        coeff = {}
        for clock,load,mem in [(True,(a,b),m),(False,(b,c),n)]:
            column = h[:,9*mem+law.SRC.index(load)]
            for idx in np.flatnonzero(abs(column)>1e-13):
                new,label = law.SRC[idx%9],idx//9
                key = (*new,c,label,n) if clock else (a,*new,m,label)
                coeff[key] = coeff.get(key,0.)+column[idx]
        assert all(k in basis or abs(x)<1e-12 for k,x in coeff.items())
        joint[:,col] = [coeff.get(k,0.) for k in basis]
    sign_term, mobile = local_source(1.,8.,amount)
    data = np.zeros((3,3))
    data[1:,1:] = mobile
    expected = np.kron(sign_term,np.eye(3))+np.kron(np.eye(2),data)
    np.testing.assert_allclose(joint,expected,atol=1e-13)
    np.testing.assert_allclose(expm(-.7j*joint),np.kron(expm(-.7j*sign_term),expm(-.7j*data)),atol=3e-13)
    # Replacing the sharp +/-2 clock by a +/-1 label opens a leakage channel.
    dirty = 9*(1+law.LABELS.index((0,1)))+law.SRC.index((0,1))
    assert np.linalg.norm(law.echo()[:9,dirty])>.5


def test_native_pointer_and_degree_phase():
    for q in (0,1,2):
        raw = law.echo()[:,law.SRC.index((q,0))]
        written = raw[9:].reshape(10,9)
        support = [d for row,(_,d) in zip(written,law.LABELS,strict=True) if np.linalg.norm(row)>1e-12]
        assert support==[q,q]
        np.testing.assert_allclose(law.echo()@raw,np.eye(99)[:,law.SRC.index((q,0))],atol=1e-15)
    for q in (1,2):
        edges = [law.SRC.index(a) for a in ((q,0),(0,q))]
        swap = expm(-1j*np.pi*law.seam_term(8))[np.ix_(edges,edges)]
        np.testing.assert_allclose(swap,[[0,1],[1,0]],atol=6e-14)


def test_refinement_and_arbitrary_phase_source():
    matching = protocol("shell_protocol.json")["matching"]
    with pytest.raises(ValueError,match="perfect matching"):
        Board(0,matching[:-1]+[matching[0]])
    board = Board(1,matching)
    rng = np.random.default_rng(9721)
    target = rng.normal(size=board.car.n)+1j*rng.normal(size=board.car.n)
    target[0]=0
    target /= np.linalg.norm(target)
    wave,_ = prepare(board,target,2,root=0)
    np.testing.assert_allclose(wave,target,atol=1e-13)
    assert np.all(np.bincount(board.bins)==24)


def test_joint_hand_and_analyzer_complements():
    angles=np.arange(10)*.7
    joint=hands(local_source(.5,-8.)[0],angles/.5)
    np.testing.assert_allclose(read_hand(joint),angles,atol=1e-13)
    with pytest.raises(ValueError,match="unresolved"):
        read_hand(np.ones((3,4)))
    with pytest.raises(ValueError,match="missing"):
        read_hand(np.zeros((3,4)))
    probes=np.array([[1.,1.],[1.,-1.],[1.,1j]])
    actual=analyzer(probes,np.full(2,.5),np.ones(3),0)
    np.testing.assert_allclose(actual,[[1,0,0],[.25,.25,.5],[.625,.125,.25]],atol=1e-14)
    with pytest.raises(ValueError,match="incomplete"):
        response(np.zeros((2,3,3,4),int),32)


def test_projector_and_beat_controls():
    rotation=np.eye(4)
    rotation[1:3,1:3]=np.array([[3,4],[-4,3]])/5
    a,r,_=first_line(np.array([0,1.,1.,2.]),rotation)
    assert r==2
    np.testing.assert_allclose(a,[0,1,0,0],atol=1e-15)
    t=np.arange(513)/128
    for rest,inc in [(40.,0.),(-80.,-.03),(32.,.04)]:
        reading,band=window_read(np.exp(1j*t*rest),np.exp(1j*t*(rest+inc)),t,128,.015,.0015)
        assert reading["status"]=="resolved"
        assert band[0]<=(rest+inc)**2-rest**2<=band[1]
    with pytest.raises(ValueError,match="unresolved"):
        window_read(np.zeros_like(t),np.ones_like(t),t,128,.015,.0015)
