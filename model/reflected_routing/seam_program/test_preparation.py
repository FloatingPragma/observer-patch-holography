import numpy as np
import pytest
from scipy.linalg import expm

from model.polynomial_sv import carrier, law, sectors
from .compiler import Compiler
from .preparation import execute_source, make_source, raw_event


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
