"""Independent source, measurement, alias and degeneracy checks for the clocked flight."""
import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from .shell_source import Board, local_source, prepare
from .shells import HERE, draws, hands, protocol, read_hand
from .shell_spectroscopy import analyzer, extract, first_line, response
from .shell_beats import window_read


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
    assert window_read(np.zeros_like(t),np.ones_like(t),t,128,.015,.0015)[1] is None


def test_return_counts_before_marginalizing():
    counts = np.zeros((2,5,3,4),dtype=np.int32)
    counts[:,:,2,0] = 32  # All misses are a valid, uninformative cohort.
    hand = np.arange(5,dtype=float)
    np.testing.assert_array_equal(response(counts,32,5),np.zeros(5))
    assert extract([counts,counts],hand,.02,32,5)[1] is None
    damaged = counts.copy()
    damaged[0,0,1,:2] = [-1,1]  # Identical marginals conceal impossible raw bins.
    for wrong in (damaged,counts.astype(float),counts[:,:,:,:1],counts[:,:4],counts*0):
        with pytest.raises(ValueError,match="incomplete"):
            extract([counts,wrong],hand,.02,32,5)
    for shots in (True,0,3.5,2**31):
        with pytest.raises(ValueError,match="shot count"):
            draws(np.random.default_rng(1),np.ones(1),shots)


def test_common_relativistic_shift_and_blind_clock():
    t = np.arange(513)/128
    for rest in (40.,160.):
        arms = np.exp(1j*np.outer([rest,np.hypot(rest,2.)],t))
        _,band = window_read(*arms,t,128,.01582,.001502)
        assert band[0]<4.<band[1]
        assert window_read(*arms,np.zeros_like(t),128,.01582,.001502)[1] is None
        assert window_read(*arms,None,128,.01582,.001502)[1] is None


def test_joint_return_histograms_against_native_contact_spectrum():
    # Assemble a separate position generator from native matrix elements, with the
    # energy of the identical vacuum hand subtracted on each individual seam.
    from model.executable_model import federation
    plan = protocol("shell_spectroscopy_protocol.json")
    matching = {tuple(e) for e in protocol("shell_protocol.json")["matching"]}
    car = federation(plan["level"])
    h0,h1 = law.local_term_parts()
    shots = plan["shots"]
    times = np.arange(plan["clock_steps"]+1)/128
    cos,sin = np.cos(times),np.sin(times)
    clock = np.stack([(1+x*cos)*(1+y*sin)/4 for x in (1,-1) for y in (1,-1)],axis=-1)
    # A joint 99% Hoeffding budget for every individual bin of both acquisitions.
    bins = 2*len(plan["law_members"])*2*2*2*2*len(times)*3*4
    tolerance = np.sqrt(np.log(2*bins/.01)/(2*shots))+1e-8
    with np.load(HERE/"shell_spectroscopy_counts.npz") as endpoint, np.load(HERE/"shell_beat_counts.npz") as beats:
        for member in plan["law_members"]:
            b,k = (1.,8.) if member=="native" else (.5,-8.)
            h = b*h0+k*h1
            for q in plan["loads"]:
                matrix = np.diag(k*q*q*(6-car.degree))
                for i,j in car.seams:
                    signed = i//12==j//12 and (i%12,j%12) in matching
                    memory = (1+law.LABELS.index((0,2))) if signed else 0
                    zero,left,right = [9*memory+law.SRC.index(s) for s in [(0,0),(q,0),(0,q)]]
                    matrix[np.ix_([i,j],[i,j])] += h[np.ix_([left,right],[left,right])]-h[zero,zero]*np.eye(2)
                spectrum = np.linalg.eigvalsh(matrix)
                for reference in plan["reference_terms"]:
                    for mode,energy in zip(plan["modes"],spectrum[:2],strict=True):
                        angle = (energy-k*q*q*reference)*times/b
                        for a,quadrature in enumerate((np.cos(angle),np.sin(angle))):
                            ports = np.stack([(1+quadrature)/2,(1-quadrature)/2,np.zeros_like(angle)],axis=-1)
                            predicted = ports[:,:,None]*clock[:,None,:]
                            for stored in (endpoint,beats):
                                counts = stored[f"{member}_{reference}_{q}_{mode}"][a]
                                assert counts.dtype==np.int32 and counts.min()>=0
                                np.testing.assert_array_equal(counts.sum(axis=(1,2)),shots)
                                assert np.max(abs(counts/shots-predicted))<tolerance


def test_all_low_modes_and_optimal_coherent_depth_reader():
    board = Board(0,protocol("shell_protocol.json")["matching"])
    kinetic = (board.inside+board.outside).toarray()
    coarse = np.repeat(np.eye(board.car.cells),12,axis=0)/np.sqrt(12)
    np.testing.assert_allclose(coarse.T@kinetic@coarse,board.dual.toarray()/24,atol=1e-14)
    angular = np.linalg.eigvalsh(board.dual.toarray()/24)
    full = np.linalg.eigvalsh(kinetic)
    gap = (3-np.sqrt(5))/2
    # At each discontinuity of either counting function, use interior thresholds.
    edges = np.unique(np.round(np.r_[0,angular,full,angular/(1+1/gap+angular/gap),gap],12))
    for energy in (edges[:-1]+edges[1:])/2:
        if 0 < energy < gap:
            count = sum(full<energy)
            assert sum(angular<energy) <= count <= sum(angular<(1+1/gap)*energy/(1-energy/gap))
    average = coarse@np.full(board.car.cells,1/np.sqrt(board.car.cells))
    opposite = np.tile(np.r_[np.ones(6),-np.ones(6)],board.car.cells)/np.sqrt(len(average))
    # Opposite imaginary depth phases have identical port probabilities. A coherent
    # instrument can distinguish them. Test its *optimal* separation as the energy falls.
    for eta in (1e-2,1e-4,1e-6):
        incident = np.sqrt(1-eta)*average+1j*np.sqrt(eta)*opposite
        first = np.outer(incident,incident.conj())
        second = first.conj()
        optimum = np.linalg.eigvalsh(first-second).clip(min=0).sum()
        np.testing.assert_allclose(optimum,2*np.sqrt(eta*(1-eta)),atol=1e-13)
        assert optimum <= 2*np.sqrt(np.vdot(incident,kinetic@incident).real/gap)
