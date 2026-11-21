from itertools import product
from math import comb

import numpy as np
import pytest

from tools.strict_json import load
from .clocked_apparatus import clock_laws, patches
from .clocked_field import eigenphase, mobile_weights, refinement, scalar_kernel, stopping_law
from .clocked_pins import verify_sources
from .clocked_run import finite_source, receive_columns, translated_word
from .compiler import Compiler
from .gate import real_step
from .preparation import HERE
from .receiver import receive
from .free_marker import leakage, routed_read
from model.executable_model import carrier
from model.candidate_laws.contrast_memory_spherical_v1 import sectors


def test_all_coin_directions_and_endpoints():
    uniform = np.ones(8)/np.sqrt(8)
    reflection = 2*np.outer(uniform,uniform)-np.eye(8)
    for k in (np.zeros(4), np.full(4,np.pi), np.full(4,np.pi/3), np.array([0,.3,-1.4,2.5])):
        shift = np.zeros((8,8),complex)
        shift[np.arange(8)^1,np.arange(8)] = np.exp(1j*np.repeat(k,2)*np.tile([1,-1],4))
        operator = shift@reflection
        phase = eigenphase(k)
        rank, weights = mobile_weights(k)
        basis = np.column_stack([uniform, shift@uniform])
        projector = basis@np.linalg.pinv(basis)
        np.testing.assert_allclose(weights,np.diag(projector).real,atol=2e-14)
        np.testing.assert_allclose(operator@(np.eye(8)-projector),-shift@(np.eye(8)-projector),atol=2e-14)
        assert np.isclose(weights.sum(),rank)
        for n in (0,1,2,5,12):
            amp = uniform@np.linalg.matrix_power(operator,n)@uniform
            assert abs(amp-np.cos(n*phase)) < 3e-14
    assert mobile_weights(np.zeros(4))[0] == mobile_weights(np.full(4,np.pi))[0] == 1
    np.testing.assert_allclose(mobile_weights(np.full(4,np.pi/3))[1],.25)


def test_literal_translated_signed_packets():
    shape = (3,4,2,5)
    for order in range(shape[-1]):
        wave = np.repeat(np.exp(2j*np.pi*order*np.arange(shape[-1])/shape[-1]),8).reshape(-1,8)
        wave /= np.linalg.norm(wave)
        grid = np.zeros((*shape,8),complex)
        grid[0,0,0] = wave
        for n in range(6):
            for pos in np.ndindex(shape[:3]):
                amplitude = np.vdot(wave,grid[pos])
                assert abs(amplitude-scalar_kernel(shape,order,pos,[n])[0]) < 2e-14
            grid = real_step(grid)


def test_clock_tree_and_timeouts():
    verify_sources()
    spec = load(HERE/"clocked_protocol.json")["fixture"]
    host,data,ports,front,bank = patches(0)
    assert len(set(ports)) == data.n and all(p//12 not in {front,bank} for p in ports)
    assert {(ports[a],ports[b]) for a,b in data.seams} == {(a,b) for a,b in host.seams if a in ports and b in ports}
    material,_ = clock_laws(host,front,bank,spec)
    cap = spec["alarm_after_periods"]
    for r in spec["thresholds"]:
        w,timeout = stopping_law(r,cap)
        starts = np.zeros(cap+1)
        for bits in product((0,1),repeat=cap):
            if not r:
                starts[0] += 1
            elif sum(bits)>=r:
                starts[np.flatnonzero(np.cumsum(bits)==r)[0]+1] += 1
        np.testing.assert_array_equal(starts/2**cap,material[r])
        np.testing.assert_array_equal(w,material[r])
        assert timeout == sum(comb(cap,j) for j in range(r))/2**cap
    with pytest.raises(ValueError):
        stopping_law(-1,cap)


def test_clock_long_window_and_small_tail():
    pmf, timeout = stopping_law(3,2)
    assert not pmf.any() and timeout==1
    assert stopping_law(1,64)[1] == 2**-64
    pmf, timeout = stopping_law(1024,2048)
    # By binomial symmetry exactly half the noncentral outcomes time out.
    centre = comb(2048,1024) / 2**2048
    assert pmf[-1] == centre/2
    assert timeout == (2**2048-comb(2048,1024)) / 2**2049
    assert abs(pmf.sum()+timeout-1) < 3e-16


def test_phase_keeps_resolved_modes_at_both_endpoints():
    for angle in (1e-180, 1e-8, np.pi-1e-8):
        # For four identical angles on this interval, E is exactly the angle.
        assert abs(eigenphase(np.full(4,angle))-angle) <= 2*np.spacing(angle)
    rank, weights = mobile_weights(np.full(4,1e-180))
    assert rank == 2
    np.testing.assert_allclose(weights,.25)
    for function in (eigenphase,mobile_weights):
        with pytest.raises(ValueError):
            function([0,0,0,1j])
    assert scalar_kernel((2,2,2,3),0,(0,0,0),0)[0]==1


@pytest.mark.parametrize("shape,order,position,periods", [
    ((2,2,2,0),0,(0,0,0),[1]), ((2,2,2,3),.5,(0,0,0),[1]),
    ((2,2,2,3),True,(0,0,0),[1]), ((2,2,2,3),0,(.5,0,0),[1]),
    ((2,2,2,3),0,(0,0,0),[-1]), ((2,2,2,3),0,(0,0,0),[.5]),
])
def test_field_requires_executable_discrete_settings(shape,order,position,periods):
    with pytest.raises(ValueError):
        scalar_kernel(shape,order,position,periods)


def test_aperture_and_batched_adjoint_on_arbitrary_pair_input():
    host,data,ports,_,_ = patches(0)
    sector,_,_ = finite_source(data,host,ports,64)
    small = Compiler(data,0,[1,1,1,3])
    full = Compiler(data,0,[2,2,2,3])
    events = sectors.SeamEvents(sector,8)
    ports_grid = np.array(full.slots).reshape((*full.shape,8))
    for pos in ((0,0,0),(1,0,0),(0,1,0),(0,0,1)):
        swap = full.columns(translated_word(full,small,pos))
        indices = np.array(full.slots)
        other = ports_grid[pos].ravel()
        for a,b in zip(small.slots,other,strict=True):
            indices[[a-1,b-1]] = indices[[b-1,a-1]]
        np.testing.assert_allclose(swap,np.eye(data.n)[:,indices],atol=2e-13)
    rng = np.random.default_rng(442)
    states = rng.normal(size=(sector.P,2))+1j*rng.normal(size=(sector.P,2))
    states /= np.linalg.norm(states,axis=0)
    neg = list(range(data.n-len(small.slots),data.n))
    wave = np.exp(2j*np.pi*np.arange(24)/24)/np.sqrt(24)
    result,losses,cost = receive_columns(small,events,neg,wave,states)
    for col in range(2):
        ref,stops,single = receive(small,events,neg,wave,"native",states[:,col])
        np.testing.assert_allclose(result[:,col],ref,atol=2e-14)
        np.testing.assert_allclose(losses[:,col],stops,atol=2e-14)
        assert cost["events"]==single["events"] and cost["queries"]==single["reader_queries"]


def test_clock_channel_error_and_fixed_compact_cutoff():
    rows = refinement()
    for row in rows:
        assert row["fixed_error"] <= row["fixed_bound"]
        assert row["clock_error"] <= row["clock_bound"]
    # At fixed private site count an order-one rest frequency grows with spatial refinement.
    phase = eigenphase(np.array([0,0,0,2*np.pi/3]))
    assert 64*phase > 32*phase > 16*phase > 1
    # Independently sum the first-passage law, including a nontrivial coherent input.
    r = 6
    w,timeout = stopping_law(r,100)
    z = np.exp(.71j)
    assert timeout <= 1e-14
    assert abs(w@z**np.arange(101)-(z/(2-z))**r) < 1e-14


def test_complete_private_tower_spacing_and_sign():
    for p in range(2,100):
        s = np.arange(p)
        phase = eigenphase(np.column_stack((np.zeros((p,3)),2*np.pi*s/p)))
        np.testing.assert_allclose(phase[1:],phase[:0:-1],atol=1e-14)
        steps = np.diff(phase[:p//2+1])
        assert np.all(steps > 0) and np.all(np.diff(steps) <= 1e-14)
        assert phase[0]==0


def test_static_channel_rest_and_refinement_curvature():
    from model.candidate_laws.interval_echo.openwave.row_30_autonomous.channel import K, bloch
    a = v = K/2
    root = np.hypot(a,v)
    rest = root-a
    speed_sq = a**3/(4*root)
    slopes = []
    for q in (0,.01,.02,.2,1.,np.pi):
        values,_ = bloch(q)
        middle = (values[2]-values[1])/2
        expected_sq = v*v+2*a*a-2*a*np.hypot(v,a*np.cos(q/2))
        assert abs(middle**2-expected_sq)<1e-13
        if q==0:
            assert abs(middle-rest)<1e-14
        if 0<q<=.02:
            slopes.append((middle-rest)/q**2)
    # Cancel the first finite-spacing correction before testing the Taylor coefficient.
    assert abs((4*slopes[0]-slopes[1])/3-speed_sq/(2*rest))<1e-8


def test_partner_independent_exchange_kernel_threshold():
    for d in range(5):
        hops = np.zeros(6,complex)
        hops[:d] = np.arange(1,d+1)*np.exp(1j*np.arange(d))
        terms = []
        for b in range(6):
            for c in range(b):
                row = np.zeros(6,complex)
                row[b],row[c] = hops[c],hops[b]
                terms.append(row)
        assert 6-np.linalg.matrix_rank(terms)==(6 if d==0 else int(d<3))


def test_free_leakage_against_full_native_generator():
    from model.experiments.interval_echo_encounter.conjugation import binary_generator
    car = carrier.cell()
    sec = sectors.PairSector(car)
    for marker in (0,3,11):
        on = (sec.aa==marker)|(sec.bb==marker)
        for kappa in (-8,0,8):
            block = sec.hamiltonian(kappa)[~on][:,on].toarray()
            info = leakage(car,marker)
            expect = np.repeat([x[0] for x in info["squared_singular_values"]],
                               [x[1] for x in info["squared_singular_values"]])
            np.testing.assert_allclose(np.linalg.eigvalsh(block.T@block),np.sort(expect),atol=2e-14)
        reference = sec.hamiltonian(8)[~on][:,on].toarray()
        for profile in (0,1,"singleton"):
            other = binary_generator(sec,8,profile)[~on][:,on].toarray()
            np.testing.assert_allclose(other,reference,atol=2e-14)


def test_native_marker_reader_on_complete_active_carrier():
    _,car,_,_,_ = patches(0)
    comp = Compiler(car,0,[1,1,1,3])
    sec = sectors.PairSector(car)
    ev = sectors.SeamEvents(sec,8)
    rng = np.random.default_rng(705)
    initial = rng.normal(size=sec.P)+1j*rng.normal(size=sec.P)
    initial /= np.linalg.norm(initial)
    found,final = routed_read(comp,ev,initial)
    marked = [sec.index(0,i) for i in range(1,car.n)]
    np.testing.assert_allclose(found,abs(initial[marked])**2,atol=3e-16)
    initial[marked] = 0
    np.testing.assert_allclose(final,initial,atol=1e-13)


def test_triangulation_hypothesis_on_every_released_federation():
    import networkx as nx
    faces = [tuple(face) for face in nx.enumerate_all_cliques(nx.icosahedral_graph()) if len(face)==3]
    for level in range(3):
        dual = nx.Graph()
        owners = {}
        for i,(a,b,c) in enumerate(faces):
            for edge in ((a,b),(b,c),(c,a)):
                key = tuple(sorted(edge))
                if key in owners:
                    dual.add_edge(i,owners[key])
                else:
                    owners[key] = i
        actual = nx.Graph((a//12,b//12) for a,b in carrier.federation(level).external)
        assert nx.vf2pp_is_isomorphic(actual,dual)
        mid = {edge:max(max(f) for f in faces)+1+i for i,edge in enumerate(owners)}
        refined = []
        for a,b,c in faces:
            u,v,w = (mid[tuple(sorted(edge))] for edge in ((a,b),(b,c),(c,a)))
            refined.extend(((a,u,w),(u,b,v),(w,v,c),(u,v,w)))
        faces = refined
