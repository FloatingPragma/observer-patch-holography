"""Independent Fock-space, graph, observer-record and full replay checks."""
import os
import math
import subprocess
import sys

import numpy as np
import pytest
from scipy import sparse as sp
from scipy.linalg import expm
from scipy.special import jv

from tools.strict_json import load
from ...contrast_memory_spherical_v1 import carrier
from . import apparatus, relay, rest, source
from .__main__ import HERE, ROOT, check_counts, compare, gates


def annihilators(n):
    """Literal 2n fermion bits and n three-level neutral oscillators."""
    dim = 4**n*3**n
    basis = np.arange(dim)
    mats = []
    for j in range(2*n):
        step = 1<<j
        full = basis[(basis & step)!=0]
        signs = np.array([(-1)**int(x & (step-1)).bit_count() for x in full])
        mats.append(sp.csr_matrix((signs,(full-step,full)),shape=(dim,dim)))
    for j in range(n):
        step = 4**n * 3**j
        occ = basis//step%3
        full = basis[occ>0]
        mats.append(sp.csr_matrix((np.sqrt(occ[full]),(full-step,full)),shape=(dim,dim)))
    return mats[:n],mats[n:2*n],mats[2*n:]


def embedding(n):
    dim = 4**n * 3**n
    P = np.zeros((dim,2*n*n))
    for i in range(n):
        for j in range(n):
            P[(1<<i)+(1<<(n+j)),i*n+j] = 1
            P[4**n*(3**i+3**j),n*n+i*n+j] = 1 if i==j else 1/np.sqrt(2)
    return P


@pytest.mark.parametrize('law',['exchange','dephased_contact','parked_neutral'])
def test_local_fock_law_analyzers_and_continuity(law):
    n = 3
    car = carrier.Carrier('triangle',n,[(0,1),(1,2),(0,2)])
    plus,minus,field = annihilators(n)
    P = embedding(n)
    active = np.array([True,False,True])
    H = sp.csr_matrix((len(P),len(P)))
    probes = [[],[],[],[]]
    charge = sum(a.T@a-b.T@b for a,b in zip(plus,minus))
    K = sum(a.T@a+b.T@b+c.T@c for a,b,c in zip(plus,minus,field))
    for i,j in car.seams:
        for channel,ops in enumerate((plus,minus,field)):
            if channel==2 and law=='parked_neutral':
                H += ops[i].T@ops[i]+ops[j].T@ops[j]
            else:
                op = ops[i]+ops[j] if channel<2 else ops[i]-ops[j]
                H += op.T@op
                probes[channel].append(op.T@op/2)
    if law=='parked_neutral':
        probes[2] = [c.T@c for c in field]
    for i in range(n):
        pair = minus[i]@plus[i]
        bosons = field[i]@field[i]/np.sqrt(2)
        if active[i] and law!='dephased_contact':
            op = pair-bosons
            H += op.T@op
            probes[3].append(op.T@op/2)
        else:
            term = pair.T@pair + bosons.T@bosons
            H += term
            probes[3].append(term)
    reduced = source.generator(car,law,active)
    np.testing.assert_allclose(H@P,P@reduced.toarray(),atol=8e-15)
    assert abs(H@charge-charge@H).max()<1e-13
    assert abs(H@K-K@H).max()<1e-13
    assert H[:,0].nnz==0
    rng = np.random.default_rng(632)
    X,Y = rng.normal(size=(2,n,n))+1j*rng.normal(size=(2,n,n))
    Y = (Y+Y.T)/2
    norm = np.linalg.norm(np.stack([X,Y]))
    X,Y = X/norm,Y/norm
    psi = P @ np.r_[X.ravel(),Y.ravel()]
    rows,scale = source.energy_instruments(car,X,Y,law,active)
    direct = []
    for group in probes:
        for probe in group:
            applied = probe@psi
            mean = np.vdot(psi,applied).real
            two = (np.vdot(applied,applied).real-mean)/2
            direct.append([1-mean+two,mean-2*two,two])
    np.testing.assert_allclose(rows,direct,atol=3e-15)
    assert scale@(rows[:,1]+2*rows[:,2]) == pytest.approx(np.vdot(psi,H@psi).real,abs=1e-13)


def test_symmetric_neutral_measurement_and_luders_instrument():
    car = carrier.Carrier('three addresses',3,[(0,1),(1,2)])
    n = car.n
    P = embedding(n)
    packed = P[np.flatnonzero(np.linalg.norm(P,axis=1))]
    H = source.generator(car,'exchange')
    U = expm(-.83j*(packed@H@packed.T))
    initial = np.zeros(2*n*n)
    initial[1] = 1
    wave = next(source.evolve(H,initial,[.83]))
    parts = relay.sectors(wave,np.arange(n),1)
    # A literal occupation read in the packed number basis.
    b = np.array([0]* (4**n*3**n))
    for index in range(len(b)):
        b[index]=index//(4**n*3)%3
    b = b[np.flatnonzero(np.linalg.norm(P,axis=1))]
    first = U@(packed@initial)
    for h,v in enumerate(parts):
        collapsed = first*(b==h)
        np.testing.assert_allclose(packed@v,collapsed,atol=2e-14)
        after = next(source.evolve(H,v,[.83]))
        np.testing.assert_allclose(packed@after,U@collapsed,atol=2e-14)
    assert np.linalg.norm(sum(parts)-wave)<2e-14


def test_exact_dark_sector_and_dense_cell_propagation():
    car = carrier.cell()
    n = car.n
    H = source.generator(car,'exchange')
    U = expm(-.63j*H.toarray())
    ordered,even,odd = (source.initial(car,'cell',k) for k in ('ordered','even','odd'))
    for v in (ordered,even,odd):
        np.testing.assert_allclose(next(source.evolve(H,v,[.63])),U@v,atol=3e-13)
    p,q,r = [np.linalg.norm((U@v)[n*n:])**2 for v in (ordered,even,odd)]
    assert q==pytest.approx(2*p,abs=1e-14) and r < 1e-25
    assert 0 < p <= 29/120
    assert np.max(abs(U-U.T))<2e-13


def test_rest_multiplicity_and_disjoint_packets():
    car = carrier.federation(0)
    Q,L = source.one_body(car)
    values = np.linalg.eigvalsh(Q.toarray())
    assert min(values)==pytest.approx(5-np.sqrt(5),abs=1e-13)
    assert sum(abs(values-min(values))<1e-11)==30
    assert np.linalg.norm(L@np.ones(car.n))==0
    psi,homes = rest.packets(car,'disjoint')
    assert len(set(homes))==4
    assert np.linalg.norm(source.generator(car,'exchange')@psi-2*min(values)*psi)<2e-13


def test_receiver_meter_erases_switch_work_for_each_outcome():
    car = carrier.Carrier('three ports', 3, [(0,1),(1,2)])
    labels = np.array([0,1,1])
    before = source.generator(car, 'exchange', labels==0)
    after = source.generator(car, 'exchange')
    relay.check_switch(before, after, labels, 1)
    delta = after-before
    # A coherent charged/neutral coincidence does carry switch work before the meter.
    wave = np.zeros(2*car.n**2, complex)
    pair = car.n+1
    wave[[pair, car.n**2+pair]] = 1/np.sqrt(2)
    assert np.vdot(wave, delta@wave).real == pytest.approx(-1)
    for part in relay.sectors(wave, labels, 1):
        assert np.vdot(part, delta@part) == 0
    # An additional onsite cost does work even after reading the meter.
    extra = sp.csr_matrix(([1.], ([pair], [pair])), shape=delta.shape)
    with pytest.raises(ValueError, match='receiver contact'):
        relay.check_switch(before, after+extra, labels, 1)


def test_detector_support_and_energy_confidence_budget():
    car = carrier.cell()
    shots = 10
    edges = len(car.seams)
    raw = np.zeros((3*edges+car.n, 3), dtype=np.int32)
    raw[:,0] = shots
    for row in (0, edges, 3*edges):
        corrupt = raw.copy()
        corrupt[row] = [shots-1, 0, 1]
        with pytest.raises(ValueError, match='binary energy detector'):
            apparatus.energies(car, corrupt, shots, 'exchange')
    raw[2*edges] = [shots-1, 0, 1]
    assert apparatus.energies(car, raw, shots, 'exchange')['neutral_energy'] == pytest.approx(4/shots)
    spec = load(HERE/'rest_protocol.json')
    for key, invalid in [('energy_cap',25), ('max_bin_energy',3)]:
        with pytest.raises(ValueError, match='understates'):
            rest.halfwidth(dict(spec, **{key:invalid}), carrier.federation(0))


def test_receiver_rejects_inconsistent_meter_and_unresolved_contrast():
    spec = load(HERE/'relay_protocol.json')
    with np.load(HERE/'records'/'relay'/'counts.npz', allow_pickle=False) as archive:
        records = dict(archive)
    car = carrier.federation(spec['level'])
    receiver = relay.receiver_cell(car)
    immediate = records['receiver_on_joint'][0]
    i = int(np.argmax(immediate[0, :car.cells**2]))
    immediate[0,i] -= 1
    immediate[1,i] += 1
    with pytest.raises(ValueError, match='simultaneous receiver'):
        relay.analyze(spec, records)
    immediate[0,i] += 1
    immediate[1,i] -= 1
    # Preserve the observed herald totals. A raw gain above sixteen counts can
    # be consistent with equal probabilities at the simultaneous precision.
    local = receiver*(car.cells+1)
    for arm, target in [('receiver_on',110), ('receiver_off',90)]:
        for branch in records[arm+'_joint'][1:,2]:
            donor = car.cells**2 + int(np.argmax(branch[car.cells**2:]))
            change = target - int(branch[local])
            assert branch[donor] >= change
            branch[donor] -= change
            branch[local] = target
    checks = relay.analyze(spec, records)['gates']
    assert checks['receiver_contrast']
    assert not checks['resolved_receiver_effect']


def chebyshev(car, state, t, active):
    """Independent matrix action and Bessel series on the certified energy interval [0,26]."""
    i,j = np.asarray(car.seams).T
    A = sp.coo_matrix((np.ones(2*len(i)),(np.r_[i,j],np.r_[j,i])),shape=(car.n,car.n)).tocsr()
    D = sp.diags(np.asarray(A.sum(axis=1)).ravel())
    Q,L = D+A,D-A
    def centered(pair):
        x,y = pair
        charged = Q@x + (Q@x.T).T
        neutral = L@y + (L@y.T).T
        ix = np.diag_indices(car.n)
        charged[ix] += x[ix] - active*y[ix]
        neutral[ix] += y[ix] - active*x[ix]
        return np.stack([charged,neutral])/13-pair
    z = 13*t
    degree = math.ceil(z)
    # |J_k(z)| <= (z/2)^k/k! by the Poisson integral. Bound the rest geometrically.
    while (2*math.exp((degree+1)*math.log(z/2)-math.lgamma(degree+2)) /
           (1-z/(2*(degree+2)))) > 1e-13:
        degree += 1
    prev = state.reshape(2,car.n,car.n)
    curr = centered(prev)
    result = jv(0,z)*prev-2j*jv(1,z)*curr
    for k in range(2,degree+1):
        prev,curr = curr,2*centered(curr)-prev
        result += 2*(-1j)**k*jv(k,z)*curr
    return np.exp(-1j*z)*result


def test_independent_federation_relay_and_rest_energy():
    car = carrier.federation(0)
    home = np.arange(car.n)//12==0
    far = np.arange(car.n)//12==relay.receiver_cell(car)
    wave = chebyshev(car,source.initial(car,0,'ordered'),8,home)
    incoming = np.zeros_like(wave)
    incoming[1] = wave[1]*far[:,None]*far[None,:]
    data = load(HERE/'records'/'relay'/'result.json')['readings']
    shots = load(HERE/'relay_protocol.json')['shots']
    for name,active in [('receiver_on',home|far),('receiver_off',home)]:
        final = chebyshev(car,incoming,1,active)
        p = np.linalg.norm(final[0][far][:,far])**2
        measured = data[name]['later_receiver_pair_counts'][1]/shots
        assert abs(p-measured) < 7*np.sqrt(p*(1-p)/shots)+8/shots
    psi,homes = rest.packets(car,'overlap')
    result = chebyshev(car,psi,32,np.ones(car.n,bool))
    i,j = np.array(car.seams).T
    # Number-conserving neutral bond energy uses both indistinguishable coordinates.
    energy = 2*np.sum(abs(result[1,i]-result[1,j])**2,axis=1)
    outside = (i//12==j//12)&~np.isin(i//12,homes)
    read = load(HERE/'records'/'rest'/'result.json')['histories']['0_overlap']
    assert abs(energy.sum()-read['neutral_energy'][-1]) < read['energy_halfwidth']
    assert abs(energy[outside].sum()-read['remote_neutral_energy'][-1]) < read['energy_halfwidth']


def test_corrupted_counts_clocks_and_heldout_spectrum(tmp_path):
    original = {'reading':np.array([[200,57]],dtype=np.int32)}
    check_counts(tmp_path/'original',original,record=True)
    with pytest.raises(FileExistsError):
        check_counts(tmp_path/'original',original,record=True)
    for artifacts in ({}, {'reading':np.array([[199,58]],np.int32)}, {'reading':original['reading'].astype(float)}):
        np.savez(tmp_path/'counts.npz',**artifacts)
        with pytest.raises(ValueError):
            check_counts(tmp_path,original)
    with np.load(HERE/'records'/'counts.npz',allow_pickle=False) as z:
        clock=z['clock']
        orbital=z['plus']
    spec=load(HERE/'protocol.json')
    ticks=apparatus.read_clock(clock,spec['shots'])
    clock[30]=clock[29]
    with pytest.raises(ValueError,match='clock'):
        apparatus.read_clock(clock,spec['shots'])
    orbital[spec['fit_samples']:] = orbital[0]
    reading=apparatus.lines(orbital,ticks,spec)
    assert reading['heldout_rms']>spec['gates']['heldout_rms']
    result=load(HERE/'records'/'result.json')
    result['spectra']['plus']=reading
    assert not gates(spec,result)['spectrum_plus']
    with pytest.raises(ValueError):
        compare({'field':[False]},{'field':[0]})


def test_coarse_position_retains_every_charge_and_coincidence():
    X=np.zeros((4,4),complex)
    Y=X.copy()
    X[2,0]=np.sqrt(.2)
    Y[0,2]=Y[2,0]=np.sqrt(.25)
    Y[3,3]=np.sqrt(.3)
    np.testing.assert_allclose(source.position(X,Y,np.array([0,0,1,1])),[0,0,.2,0,0,.5,.3])


@pytest.mark.parametrize('suffix',['','.relay','.rest'])
def test_public_replays(suffix):
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    subprocess.run([sys.executable,'-B','-O','-W','error','-m','model.scalar_channel.radiative_repair'+suffix],
                   cwd=ROOT,env=env,check=True)
    folder = HERE/'records'/suffix.lstrip('.')
    assert all(load(folder/'result.json')['gates'].values())
