import numpy as np
import pytest
from scipy.linalg import expm

from model.contrast_memory_spherical_v1 import carrier, law, sectors
from .compiler import Compiler
from .preparation import execute_source, make_source, raw_event
from .receiver import home_read, ideal_vector, receive


@pytest.mark.parametrize("kind", ["complex","sparse","basis"])
def test_arbitrary_phase_and_all_failure_branches(kind):
    car = carrier.Carrier("test path",19,[(j,j+1) for j in range(18)])
    comp = Compiler(car,0,[1,1,1,1])
    neg = list(range(9,17))
    ev = sectors.SeamEvents(sectors.PairSector(car),8)
    if kind=="complex":
        desired = np.exp(1j*np.arange(8)**2)*np.arange(1,9)
    elif kind=="sparse":
        desired = np.array([0,0,1,0,1j,0,0,0],complex)
    else:
        desired = np.eye(8,dtype=complex)[7]
    desired /= np.linalg.norm(desired)
    inp = ev.start((0,1))
    state,losses,_ = execute_source(comp,ev,neg,desired,"native",inp)
    expected = np.zeros_like(inp)
    expected[[ev.sector.index(0,p) for p in comp.slots]] = -desired/2
    np.testing.assert_allclose(state,expected,atol=2e-13)
    np.testing.assert_allclose(losses,np.r_[.25,.25,.25*abs(desired)**2],atol=2e-13)
    with pytest.raises(ValueError):
        list(make_source(comp,comp.slots,desired,8,"native"))
    with pytest.raises(ValueError,match="contrast coupling"):
        list(make_source(comp,neg,desired,7,"native"))


def test_arbitrary_duration_from_unshifted_blank_law():
    # Put both, one or neither load on edge (0,1), without restricting to marked inputs.
    car = carrier.cell()
    sector = sectors.PairSector(car)
    ev = sectors.SeamEvents(sector,8)
    edge = car.seams.index((0,1))
    angle = .3715
    rng = np.random.default_rng(311)
    before = rng.normal(size=sector.P)+1j*rng.normal(size=sector.P)
    actual = before.copy()
    raw_event(ev,actual,edge,angle)
    ex,con = law.local_term_parts()
    local = expm(-1j*angle*(ex+8*con))[:9,:9]
    expected = np.zeros_like(actual)
    for i,(a,b) in enumerate(zip(sector.aa,sector.bb,strict=True)):
        pair = {int(a),int(b)}
        source = (int(0 in pair),int(1 in pair))
        for final in [(0,0),(0,1),(1,0),(1,1)]:
            if sum(final)!=sum(source):
                continue
            result = sorted((pair-{0,1})|{j for j,n in enumerate(final) if n})
            expected[sector.index(*result)] += local[law.SRC.index(final),law.SRC.index(source)]*before[i]
    np.testing.assert_allclose(actual,expected,atol=2e-13)


@pytest.mark.parametrize("member", ["native","omit_phase"])
def test_receiver_reads_complex_overlap_on_arbitrary_incident_wave(member):
    car = carrier.Carrier("test path",19,[(j,j+1) for j in range(18)])
    comp = Compiler(car,0,[1,1,1,1])
    neg = list(range(9,17))
    sector = sectors.PairSector(car)
    ev = sectors.SeamEvents(sector,8)
    rng = np.random.default_rng(470)
    target = rng.normal(size=8)+1j*rng.normal(size=8)
    incident = rng.normal(size=8)+1j*rng.normal(size=8)
    target /= np.linalg.norm(target)
    incident /= np.linalg.norm(incident)
    start = np.zeros(sector.P,complex)
    start[[sector.index(0,j) for j in comp.slots]] = incident
    state,losses,resources = receive(comp,ev,neg,target,member,start)
    observed = state[sector.index(0,1)]
    expected = -np.vdot(ideal_vector(target,member),incident)/2
    np.testing.assert_allclose(observed,expected,atol=2e-13)
    assert abs(np.vdot(state,state).real+sum(losses)-1)<1e-12
    assert resources["reader_queries"]==10


def test_uniform_obstruction_depends_on_the_event_law():
    car = carrier.cell()
    sector = sectors.PairSector(car)
    events = sectors.SeamEvents(sector,8)
    sharp = events.start((0,1))
    scheduled = sharp.copy()
    unshifted = sharp.copy()
    for seam,angle in [((1,5),np.pi/2),((4,5),3*np.pi/16)]:
        edge = car.seams.index(seam)
        events.apply(scheduled,edge,angle)
        raw_event(events,unshifted,edge,angle)
    marked = [sector.index(0,j) for j in range(1,car.n)]
    np.testing.assert_allclose(scheduled[marked].sum(),1,atol=2e-14)
    np.testing.assert_allclose(unshifted[marked].sum(),0,atol=2e-14)
    np.testing.assert_allclose(np.linalg.norm(unshifted[marked]),1,atol=2e-14)


def test_home_frequency_uses_the_whole_record_in_seam_order():
    seams = [(0,3),(0,5),(1,2)]
    # Two early stops, then no pair and three distinct seam bins.
    counts = np.array([2,3,4,1,9,1],dtype=np.int32)
    assert home_read(counts,2,seams,(0,5),20)==(9,9/20)
    counts[4] -= 3
    counts[0] += 3
    assert home_read(counts,2,seams,(0,5),20)==(6,6/20)
    for corrupt in (counts[2:],-counts,counts.astype(float),counts+1):
        with pytest.raises(ValueError):
            home_read(corrupt,2,seams,(0,5),20)
