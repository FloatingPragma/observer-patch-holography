import numpy as np
import pytest
from scipy.linalg import expm

from model.contrast_memory_spherical_v1 import carrier, law, sectors
from model.seam_program.compiler import Compiler


@pytest.mark.parametrize("quarters", [1,4])
def test_unshifted_99_state_seam_on_arbitrary_binary_input(quarters):
    # Rebuild directly from the original write/echo and record generator.
    tau = quarters*np.pi/4
    exchange, contrast = law.generator_parts()
    echo = law.echo()
    native = echo @ expm(-1j*tau*(exchange+8*contrast)) @ echo.conj().T
    # Four load states tensored with blank memory, including both particles on this seam.
    inputs = [law.SRC.index(pair) for pair in [(0,0),(0,1),(1,0),(1,1)]]
    embedded = np.eye(law.DIM,dtype=complex)[:,inputs]
    swap = np.eye(4)[:,[0,2,1,3]]
    generator = (np.eye(4)-swap)/2-np.diag([0,0,0,16])
    np.testing.assert_allclose(native@embedded, embedded@expm(-1j*tau*generator),atol=3e-14)
    if quarters==1:
        wrong = echo @ expm(-1j*tau*(exchange+7*contrast)) @ echo.conj().T
        assert np.linalg.norm((wrong-native)@embedded)>0.5


def test_nonedge_compilation_on_every_pair_state_with_marker():
    car = carrier.cell()
    comp = Compiler(car,0,[1,1,1,1])
    sector = sectors.PairSector(car)
    events = sectors.SeamEvents(sector,8)
    marked = [sector.index(0,i) for i in range(1,car.n)]
    state = np.eye(sector.P,dtype=complex)[:,marked]
    source = state.copy()
    assert (1,2) not in car.seams
    program = comp.exchange(1,2,1)
    for k,q in program:
        assert 0 not in car.seams[k]
        events.apply(state,k,q*np.pi/4)
    mixed = np.eye(len(marked),dtype=complex)
    phase = np.exp(-1j*np.pi/4)
    mixed[:2,:2] = [[(1+phase)/2,(1-phase)/2],[(1-phase)/2,(1+phase)/2]]
    np.testing.assert_allclose(state,source@mixed,atol=3e-15)


def test_coin_fixes_uniform_mode_including_unused_ports():
    comp = Compiler(carrier.cell(),0,[1,1,1,1])
    gate = comp.columns(comp.build("nominal"))
    np.testing.assert_allclose(gate.sum(axis=0),np.ones(8),atol=1e-14)
    np.testing.assert_allclose(gate.conj().T@gate,np.eye(8),atol=1e-14)
    np.testing.assert_allclose(gate[[0,9,10,11]],0,atol=2e-15)
    with pytest.raises(ValueError):
        comp.exchange(0,3,4)
    with pytest.raises(ValueError):
        Compiler(carrier.cell(),0,[2,1,1,1])


def test_marker_cut_and_changed_words():
    car = carrier.Carrier("line",12,[(j,j+1) for j in range(11)])
    comp = Compiler(car,4,[1,1,1,1])
    with pytest.raises(ValueError,match="disconnects"):
        comp.exchange(1,10,4)
    with pytest.raises(ValueError):
        comp.build("unregistered")
