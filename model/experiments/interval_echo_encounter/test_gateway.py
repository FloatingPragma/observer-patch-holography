"""Independent source, detector, spectral and failure checks for the remote phase experiment."""
from fractions import Fraction
from math import factorial, log

import numpy as np
import pytest
from scipy import sparse

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector, PairSector
from model.executable_model.carrier import CELL_SEAMS, Carrier, federation
from . import gateway as g, gateway_complement as complement, gateway_message as message
from .gateway_source import band_budget, geometry, high_projection, reference, source


def test_dense_federation_launch_complement_and_clock():
    board = Carrier("two complete cells",24,[*CELL_SEAMS,*[(a+12,b+12) for a,b in CELL_SEAMS],(0,12)])
    pair = PairSector(board)
    layout = geometry(board)
    s,u,v,_,edges = layout
    h,admit,prepared = source(pair,8,layout)
    dense = BinarySector(board,2)
    H = dense.hamiltonian(8)
    np.testing.assert_array_equal(h.toarray(),H)
    energies,basis = np.linalg.eigh(H)
    def evolve(state,t):
        return basis @ (np.exp(-1j*t*energies)[:,None]*(basis.T@state))
    start = np.eye(pair.P)[:,[dense.index[(0,2)]]]
    home = dense.index[tuple(sorted((s,u)))]
    second = dense.index[tuple(sorted((s,v)))]
    assert abs(abs(evolve(start,np.pi)[home,0])**2-admit)<1e-12
    for j,sgn in enumerate((1,-1)):
        psi = np.zeros(pair.P,complex)
        psi[[home,second]]=[1,sgn*1j]
        np.testing.assert_allclose(prepared[:,j]/prepared[home,j],psi,atol=2e-13)
    high = basis[:,energies>=-5]
    exact = high@(high.T@prepared)
    projected,err = high_projection(h,prepared)
    assert max(np.linalg.norm(projected-exact,axis=0))<=err
    receiver = [pair.index(*e) for e in edges]
    budget = band_budget(h,prepared,receiver)
    x = (h+2*sparse.eye(pair.P,format="csr"))/14
    state = projected.copy()
    weights = complement.coefficients()
    for n in range(4):
        evolved = evolve(exact,n*np.pi)
        np.testing.assert_allclose(state,evolved,atol=err+n*1e-8,rtol=0)
        assert np.all(complement.apertures(evolved,receiver)<=budget)
        state = complement.step(x,state,weights)
    assert reference([0,1,2,7,127]) == {i:i for i in (0,1,2,7,127)}
    for bad in ([],[1,1],[2,1],[128],[True]):
        with pytest.raises(ValueError):
            reference(bad)


def test_short_time_coefficients_from_actual_generator():
    sector = PairSector(federation(0))
    s,u,v,b,remote = geometry(sector.car)
    target = [sector.index(*e) for e in remote]
    incident = [i for i,e in enumerate(remote) if b in e]
    for k in (0,8):
        H = (2*sector.hamiltonian(k)).astype(np.int64)
        columns = np.zeros((sector.P,2),np.int64)
        columns[[sector.index(s,u),sector.index(s,v)],[0,1]] = 1
        for degree in range(6):
            if degree<=3:
                assert not np.any(columns[target])
            if degree==4:
                expected = np.zeros((30,2),np.int64)
                expected[incident,0] = 2
                np.testing.assert_array_equal(columns[target],expected)
            if degree==5:
                expected = np.zeros(30,np.int64)
                expected[incident] = -4
                np.testing.assert_array_equal(columns[target,1],expected)
            columns = H@columns
    a = Fraction(2,2**4*factorial(4))
    b_im = Fraction(4,2**5*factorial(5))
    assert -2*len(incident)*a*b_im == Fraction(-1,18432)


def test_receiver_only_adversarial_records_and_controls():
    cfg = g.settings()
    receipt = g.load_json(g.RECEIPT)
    key = "native_plus_32"
    counts = receipt["counts"][key]
    p = np.asarray(receipt["weights"][key])/10**13
    L = log(2*63*81/cfg["alpha"])
    for bad in (counts[:-1], [True]+counts[1:], [-1]+counts[1:],
                [0]*31+[cfg["shots"]],receipt["counts"]["native_minus_32"]):
        with pytest.raises(ValueError):
            g.read(bad,cfg["shots"],p,L)
    blocked = list(receipt["counts"]["blocked_plus_32"])
    blocked[0],blocked[-2] = 1,blocked[-2]-1
    with pytest.raises(ValueError,match="impossible"):
        g.read(blocked,cfg["shots"],np.asarray(receipt["weights"]["blocked_plus_32"])/10**13,L)
    for word in cfg["appointments"]:
        for law in cfg["laws"]:
            weights = [np.asarray(receipt["weights"][f"{law}_{menu}_{word}"]) for menu in cfg["menus"]]
            np.testing.assert_allclose((weights[0]+weights[1])/2,weights[2],atol=32,rtol=0)
            assert all(z[-1]==weights[0][-1] or abs(z[-1]-weights[0][-1])<=32 for z in weights)
    with pytest.raises(ValueError):
        g.summarize(cfg,{key:p},np.zeros((31,2)),{key:counts})
    with pytest.raises(ValueError):
        g.close({"quiet":1},{"quiet":True})


def test_capped_message_reader_uses_counts_only():
    cfg = g.load_json(message.PROTOCOL)
    cap,live,cut = cfg["maximum_starts"],cfg["admissions_per_frame"],cfg["threshold_minus"]
    assert message.decode([cap,live-1,cut+1],cfg)=="erasure"
    assert message.decode([cap,live,cut-1],cfg)=="plus"
    assert message.decode([cap,live,cut],cfg)=="minus"
    for bad in ([live-1,live,cut],[live-1,live-1,cut],[cap+1,live,cut],[cap,live,live+1]):
        with pytest.raises(ValueError):
            message.decode(bad,cfg)
    rng = np.random.default_rng(128)
    assert message.frame(rng,5,8,0,1)==[8,0,0]
    assert message.frame(rng,5,8,1,1)==[5,5,5]
    # A one-success frame capped at two fair attempts has three terminal cases.
    observed = [message.frame(rng,1,2,0.5,1) for _ in range(8192)]
    for row,p in (([1,1,1],0.5),([2,1,1],0.25),([2,0,0],0.25)):
        assert abs(observed.count(row)/len(observed)-p)<0.025
    message.replay()
