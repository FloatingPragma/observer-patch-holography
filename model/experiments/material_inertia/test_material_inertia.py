"""Independent full-space and operator tests for the material instruments and the mass claim."""
import os
import subprocess
import sys
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import expm

from tools.strict_json import load
from model.executable_model import carrier
from model.candidate_laws.contrast_memory_spherical_v1 import law, sectors
from . import challenge,flight,instrument,memory,prediction
from .detector import position_bins,query
from .__main__ import HERE,ROOT,check_json
from .source import Ring,pair_hop,spectral_measure


def test_complete_translation_reduction_and_direct_return():
    size,n,kappa = 8,3,2.
    ring = Ring.build(size,n)
    graph = carrier.Carrier('cycle',size,[(i,(i+1)%size) for i in range(size)])
    reference = sectors.BinarySector(graph,n)
    native = reference.hamiltonian(kappa)+6*n*kappa*np.eye(len(ring.states))
    base = ring.hamiltonian('contrast',kappa)
    np.testing.assert_allclose(native,base.toarray(),atol=1e-14)
    assisted = ring.hamiltonian('assisted',kappa)
    # Independently install each three-port projector in the full bitstring basis.
    for i in range(size):
        a,b,c = (1<<((i+j)%size) for j in range(3))
        for s in ring.states:
            if s&a and s&b and not s&c:
                v = np.zeros(len(ring.states))
                v[ring.index[s]],v[ring.index[s^a^c]] = 1,-1
                native += pair_hop(kappa)*np.outer(v,v)
    np.testing.assert_allclose(assisted.toarray(),native,atol=2e-14)
    full = np.linalg.eigvalsh(native)
    values,returns = [],[]
    for ell in range(size):
        E,w = spectral_measure(assisted,ring.frame(ell))
        values.extend(E)
        returns.append(np.exp(-.89j*E)@w)
    np.testing.assert_allclose(sorted(values),full,atol=2e-13)
    start = 2
    masks = [sum(1<<((i+j)%size) for j in range(n)) for i in range(size)]
    positions = [ring.index[m] for m in masks]
    wave = expm(-.89j*native)[:,positions[start]]
    expected = abs(wave[positions])**2
    np.testing.assert_allclose(flight.probabilities(returns,start),np.r_[expected,1-expected.sum()],atol=1e-13)
    # Even-N antipodal patterns may have short translation orbits.
    even = Ring.build(10,4)
    assert any(len(orbit)<even.size for orbit in even.orbits)
    assert sum(even.frame(k).shape[1] for k in range(even.size))==len(even.states)


def test_native_local_memory_and_full_qutrit_embedding():
    V = memory.local_conversion()
    F = np.zeros(9)
    F[4] = 1
    m = np.array([((d==0)-int(abs(d)==2)/np.sqrt(2))/np.sqrt(2) for _,d in law.LABELS])
    G = np.zeros((99,99))
    G[9:,9:] = np.kron(np.outer(m,m),np.diag(F))
    U = law.echo()
    evolved = U@expm(-.71j*G)@U
    np.testing.assert_allclose(evolved[:9,:9],expm(-.71j*V),atol=2e-14)
    assert np.linalg.norm(evolved[9:,:9]) < 2e-14
    car = carrier.cell()
    sector,H = memory.joined(car, 8, True)
    basis = sectors.MultisetSector(car,2)
    # The independent 78-state qutrit Hamiltonian is never reduced to the pair subspace.
    dense = basis.hamiltonian(8)
    for a,b in car.seams:
        v = np.zeros(len(basis.states))
        for x,y,weight in [(1,1,1.),(2,0,-1/np.sqrt(2)),(0,2,-1/np.sqrt(2))]:
            state = [0]*car.n
            state[a],state[b] = x,y
            v[basis.index[tuple(state)]] = weight
        dense += np.outer(v,v)
    order = []
    for s in basis.states:
        sites = np.flatnonzero(s)
        order.append(sector.P+sites[0] if len(sites)==1 else sector.index(*sites))
    np.testing.assert_allclose(dense,H.toarray()[np.ix_(order,order)]+96*np.eye(len(order)),atol=1e-13)
    # Check the independent two-bitplane implementation for every three-load state including doubles.
    size = 6
    ring = memory.ternary_ring(size,3)
    bare = memory.ternary_law(ring,law.seam_term(8),8).toarray()
    qutrit = sectors.MultisetSector(carrier.Carrier('ring',size,[(i,(i+1)%size) for i in range(size)]),3)
    order = [ring.index[sum((1<<i) for i,v in enumerate(s) if v)+
                        sum((1<<(i+size)) for i,v in enumerate(s) if v==2)] for s in qutrit.states]
    np.testing.assert_allclose(bare[np.ix_(order,order)],qutrit.hamiltonian(8))


def test_spectral_measure_requires_a_hermitian_generator_and_isometric_frame():
    from scipy.sparse import eye
    ring = Ring.build(7,3)
    frame = ring.frame(1)
    H = ring.hamiltonian('contrast',8)
    spectral_measure(H,frame)
    # An imaginary scalar commutes with every translation and leaks no sector.
    with pytest.raises(ValueError,match='Hermitian'):
        spectral_measure(H+1j*eye(len(ring.states)),frame)
    # Even an invariant zero generator cannot certify a scaled or duplicated frame.
    with pytest.raises(ValueError,match='isometric'):
        spectral_measure(0*H,2*frame)


def test_registered_position_counts_from_full_coordinate_space():
    # Reconstruct the two- and three-load flights without the source's orbit
    # basis, momentum blocks, Fourier readout or phase-interval reader.
    spec = load(HERE/'flight_protocol.json')
    base = load(HERE/'protocol.json')
    results = load(HERE/'records/flight/result.json')
    size = spec['size']
    kappa,_ = challenge.challenge_inputs()
    delta = 1+2*kappa
    pair_rest = (5+1/delta)/3
    assisted_hop = (1/pair_rest-1/(2*delta))/(2*(1-1/delta**2))
    with np.load(HERE/'records/flight/counts.npz',allow_pickle=False) as counts:
        for n in (2,3):
            states = [frozenset(s) for s in combinations(range(size),n)]
            index = {s:i for i,s in enumerate(states)}
            ports = [index[frozenset((j+d)%size for d in range(n))] for j in range(size)]
            native = np.zeros((len(states),len(states)))
            assisted = np.zeros_like(native)
            for col,s in enumerate(states):
                native[col,col] = 4*kappa*n
                for x in s:
                    for step in (-1,1):
                        y,z = (x+step)%size,(x+2*step)%size
                        if y not in s:
                            native[col,col] += kappa+0.5
                            native[index[s-{x}|{y}],col] -= 0.5
                        elif z not in s:
                            assisted[col,col] += assisted_hop
                            assisted[index[s-{x}|{z}],col] -= assisted_hop
            for member,H in [('contrast',native),('assisted',native+assisted)]:
                key = f'{member}_{n}'
                energies,vectors = np.linalg.eigh(H)
                overlaps = vectors[ports]*vectors[ports[spec['starting_port']]]
                for tick,recorded in zip(results[key]['ticks'],counts[key],strict=True):
                    phases = np.exp(-1j*(energies-energies[0])*(tick/base['tick_divisor']))
                    p = abs(overlaps@phases)**2
                    edges = np.r_[0,np.floor(spec['shots']*np.cumsum(p)+0.5),spec['shots']]
                    np.testing.assert_array_equal(np.diff(edges),recorded)


def test_travelling_predicate_is_reversible_and_retains_breakup():
    for word in range(32):
        for scratch in range(32):
            for pointer in (0,1):
                after,read = query(word,5,scratch,pointer)
                assert after==scratch and read==(pointer^int((word^scratch)==31))
                assert query(word,5,after,read)==(scratch,pointer)
    ring = Ring.build(11,4)
    bins = position_bins(ring)
    expected = np.full(len(ring.states),11)
    for j in range(11):
        s = sum(1<<((j+i)%11) for i in range(4))
        expected[ring.index[s]] = j
    np.testing.assert_array_equal(bins,expected)
    for control in [{'flip_message':True},{'read_b':False}]:
        with pytest.raises(ValueError,match='disputed'):
            position_bins(ring,control)


def test_equal_binding_control_nonzero_offset_and_unseen_motion():
    spec = load(HERE/'challenge_protocol.json')
    base = load(HERE/'protocol.json')
    # The hash challenge happens to give a zero common shift. Exercise a nonzero one independently.
    a = challenge.binding_inbox(spec,base,'contrast',6,0)
    for member in ('contrast','assisted'):
        b = challenge.binding_inbox(spec,base,member,6,7)
        np.testing.assert_array_equal(a,b)
    # Independently derive the split-arm signal for ideal return paths on a small real-space ring.
    ring = Ring.build(7,2)
    H = ring.hamiltonian('contrast',6).toarray()
    tim = .63
    waves = [ring.frame(k).toarray()[:,0] for k in (0,1)]
    U = expm(-1j*tim*H)
    a,b = waves
    overlap = np.vdot(np.kron(U@a,b),np.kron(a,U@b))
    expected = []
    for k in (0,1):
        E,w = spectral_measure(ring.hamiltonian('contrast',6),ring.frame(k))
        expected.append(np.exp(-1j*tim*E)@w)
    assert abs(overlap-expected[0].conj()*expected[1])<2e-13


def test_analytic_pair_and_wrap_certificate():
    for size,kappa in [(9,.6),(12,2.),(15,8.)]:
        ring = Ring.build(size,2)
        H = ring.hamiltonian('contrast',kappa)
        e0 = np.linalg.eigvalsh((ring.frame(0).conj().T@H@ring.frame(0)).toarray())[0]
        delta = 1+2*kappa
        z = 1/delta
        exact = 10*kappa+1-z
        err = (1-z*z)*z**(size-3)/(1+z**(size-2))
        assert -2e-13 <= exact-e0 <= err+2e-13
    # Compare the periodized infinite-chain droplet with the exact finite generator.
    size,n,kappa,ell = 11,3,1.,3
    K = 2*np.pi*ell/size
    delta = 1+2*kappa
    eta = np.arccosh(delta)
    ring = Ring.build(size,n)
    H = ring.hamiltonian('contrast',kappa)
    coeff = [(np.exp(1j*K)*np.sinh((n-j)*eta)+np.sinh(j*eta))/np.sinh(n*eta) for j in range(1,n)]
    vector = []
    for occupied in combinations(range(size),n):
        val = 0j
        for cut in range(n):
            loc = occupied[cut:]+tuple(x+size for x in occupied[:cut])
            val += np.exp(1j*K*loc[0])*np.prod(np.array(coeff)**(np.diff(loc)-1))
        vector.append(val)
    v = np.array(vector)
    E = 4*kappa*n + np.sinh(eta)*(np.cosh(n*eta)-np.cos(K))/np.sinh(n*eta)
    rho = max(abs(x) for x in coeff)
    from math import comb
    cap = np.sqrt(comb(size,n)/size)*n*(delta+rho)*rho**(size-n)/(1-(n-1)*rho**(size-n))
    assert np.linalg.norm(H@v-E*v)/np.linalg.norm(v) <= cap


def test_counter_clock_line_and_receipt_rejections(tmp_path):
    # Simulate equal exposures carrying large rest phases. The two kinetic
    # shifts are exactly 1/2 and 1/64 radians in the final phase records.
    phase = instrument.PhaseLine
    for common_turns in (0,2**27,2**48):
        observed = instrument.band_ratios([phase(common_turns,1.,0.),phase(common_turns,1.5,0.),
                                          phase(2*common_turns,2.,0.),phase(2*common_turns,2+1/64,0.)])
        assert observed['inertia']==[32.,32.,32.]
        assert observed['rest']==[2.,2.,2.]
    spec = load(HERE/'protocol.json')
    shots = spec['shots']
    clock = instrument.clock(spec)
    assert abs(instrument.clock_rate(clock,spec)-1/spec['tick_divisor']) <= prediction.clock_radius(spec)
    invalid = clock.copy()
    invalid[:] = [0,shots]
    with pytest.raises(ValueError):
        instrument.clock_rate(invalid,spec)
    ticks = instrument.tick_menu(spec)
    rate = 1/spec['tick_divisor']
    pointer = instrument.quadratures(np.exp(-1j*ticks*rate*.6),shots)
    estimate = instrument.read_line(pointer,rate,spec)
    assert abs(estimate[0]-.6) < estimate[1]
    pointer[-1] = pointer[-1,::-1].copy()
    assert instrument.read_line(pointer,rate,spec) is None
    with pytest.raises(ValueError):
        instrument.signal(pointer.astype(float),len(pointer),shots)
    with pytest.raises(ValueError):
        instrument.signal(pointer[:-1],len(pointer),shots)
    with pytest.raises(ValueError):
        instrument.counts(np.array([1.1]),shots)
    assert instrument.quotient((1.,.1),(.01,.02)) is None
    assert challenge.score(None,[1,2],.08)=='unresolved'
    assert challenge.score([3.,2.9,3.1],[1.,2.],.08)=='excluded'
    for bad in ([1.,float('nan'),2.],[1.,1.5,1.25],[2.,1.,1.5],[1.,0.,2.],[True,1.,1.]):
        with pytest.raises(ValueError):
            challenge.score(bad,[1.,2.],.08)
    for forecast in ([2.,1.],[0.,1.],[1.,float('inf')],[True,2.]):
        with pytest.raises(ValueError):
            challenge.score([1.,1.,1.],forecast,.08)
    for radii in ((1.,-1.),(float('nan'),0.)):
        with pytest.raises(ValueError):
            instrument.quotient(radii,(1.,0.))
    for size,n,tick in ((19,1,0),(10,5,0),(19,2,-1)):
        with pytest.raises(ValueError,match='bound domain'):
            prediction.native_transport_bound({},size,n,tick,1.,spec,1024)
    target = tmp_path/'test.json'
    target.write_text('{"answer":1}',encoding='utf8')
    with pytest.raises(ValueError):
        check_json(target,{'answer':True})
    archive = tmp_path/'counts.npz'
    good = {'p':clock}
    with archive.open('wb') as f:
        np.savez_compressed(f,p=clock.astype(float))
    with pytest.raises(ValueError):
        challenge.custody(archive,good,False)
    with archive.open('wb') as f:
        np.savez_compressed(f,p=clock[:-1])
    with pytest.raises(ValueError):
        challenge.custody(archive,good,False)


@pytest.mark.parametrize('suffix,args', [('',[]),('.memory',[]),('.challenge',['--stage','predict']),
                                         ('.challenge',['--stage','test']),('.flight',[])])
def test_public_optimized_replay(suffix,args):
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    subprocess.run([sys.executable,'-B','-O','-W','error','-m',
                    'model.experiments.material_inertia'+suffix,*args],cwd=ROOT,env=env,check=True)


def test_scientific_verdicts_include_failures():
    base = load(HERE/'records/result.json')
    assert all(v['same_counts'] and max(v['changed_region_between_reads'])>0 for v in base['federations'].values())
    verdicts = base['audit']['rest_inertia']
    assert verdicts['contrast_16_2_ratios']=='outside tolerance'
    assert verdicts['assisted_16_2_ratios']=='within tolerance'
    assert verdicts['assisted_16_3_ratios']=='outside tolerance'
    assert verdicts['exchange_16_4_ratios']=='unresolved'
    result = load(HERE/'records/challenge/test.json')['readings']
    assert result['contrast_5_1']['verdict']=='unresolved'
    assert result['contrast_5_9']['verdict']=='compatible'
    assert all(r['verdict']=='excluded' for k,r in result.items() if k.startswith('assisted'))
    memory_result = load(HERE/'records/memory/result.json')
    assert memory_result['rest_inertia']=={f'{size}_{n}':'outside tolerance' for size in (12,16) for n in (2,3)}
    measured = load(HERE/'records/flight/result.json')
    for key,row in measured.items():
        assert row['passes']==key.startswith('contrast')
        if key.startswith('contrast'):
            assert max(row['total_variation']) < max(row['binding_only_tv_bound']) < .15
