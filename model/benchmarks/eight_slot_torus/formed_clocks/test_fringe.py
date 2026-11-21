"""Independent literal operators, high-precision eigenphases and detector controls."""
import mpmath as mp
import numpy as np
import pytest

from .. import formation, law
from . import fringe, fringe_resonance, fringe_walk as study


@pytest.mark.parametrize("k", [(0.,0.,0.), (.93,-1.2,.74), (np.pi,0.,-np.pi)])
@pytest.mark.parametrize("mu", [study.MU, fringe_resonance.RES])
def test_two_dimensional_sector_in_full_slot_operator(k, mu):
    eye = np.eye(8)
    routing = law.walk_matrix([*k,0]) @ law.coin(eye)
    contact = law.contact_coin(eye, 2*mu).T
    mature = routing @ contact
    prefix = [eye]
    for age in range(study.OPEN):
        angle = law.hop_angle(age, study.MEMBER.maturation)
        hop = np.exp(-1j*angle)*(np.cos(angle)*eye + 1j*np.sin(angle)*routing)
        prefix.append(hop @ contact @ prefix[-1])

    def block(services):
        if services <= study.OPEN:
            return prefix[services]
        return np.linalg.matrix_power(mature, services-study.OPEN) @ prefix[-1]

    period = block(study.PERIOD)
    walk = study.BrightWalk(k, mu=mu)
    for services in (0,1,127,128,8191,8192,8193,16397):
        cycles, tail = divmod(services, study.PERIOD)
        expected = block(tail) @ np.linalg.matrix_power(period, cycles) @ law.EQUAL
        np.testing.assert_allclose(walk.wave(services), expected, rtol=0, atol=2e-11)


@pytest.mark.parametrize("boost", [0,1])
def test_literal_preparations_and_local_aperture(boost):
    size = 5
    shape = (1,1,size,1,8)
    x = np.arange(size)
    k = 2*np.pi/size
    phase = .48
    weights = np.exp(1j*boost*k*x)*(1+np.exp(1j*(k*x+phase)))/np.sqrt(2*size)
    incident = np.zeros(shape*2, complex)
    packets = np.outer(law.coin(np.eye(8)[4]), law.coin(np.eye(8)[5]))
    for z, amplitude in enumerate(weights):
        incident[0,0,(z-1)%size,0,:,0,0,(z+1)%size,0,:] = amplitude*packets
    free, early = formation.first_writes(incident)
    np.testing.assert_array_equal(early, 0)
    free = law.step_contact(free, study.MEMBER.rail_contact_phases[0])
    free, writes = formation.first_writes(free)
    writes = writes.ravel()
    np.testing.assert_allclose(writes, weights/np.sqrt(2), rtol=0, atol=4e-15)
    assert np.vdot(free,free).real == pytest.approx(.5)

    # Keep a separate material vector per orthogonal outgoing marker.
    branches = np.zeros((size,*shape), complex)
    for z in range(size):
        branches[z,0,0,z,0] = writes[z]*law.EQUAL
    coherent = weights.reshape(1,1,size,1,1)*law.EQUAL
    erased = branches.sum(axis=0)/np.sqrt(size)
    np.testing.assert_allclose(erased, coherent/np.sqrt(2*size), atol=2e-15)
    marked_density = branches.reshape(size,-1).T @ branches.reshape(size,-1).conj()
    for z in range(size):
        for y in range(size):
            if y != z:
                np.testing.assert_array_equal(marked_density[z*8:(z+1)*8,y*8:(y+1)*8],0)

    n = 135
    for age in range(n):
        coherent = law.step_maturing(coherent, study.MEMBER, age)
        branches = np.stack([law.step_maturing(v,study.MEMBER,age) for v in branches])
    walk = study.BrightWalk(study.angles((1,1,size)))
    modes = walk.wave(n)[0,0]
    coherent_click = np.linalg.norm(coherent[0,0,0,0])**2
    marked_click = np.linalg.norm(branches[:,0,0,0,0])**2
    contrast = np.vdot(modes[boost], modes[boost+1])
    point = np.fft.ifft(modes, axis=0)
    marked = np.sum(np.sum(abs(point)**2, axis=1)*np.exp(-1j*k*x))
    assert size*coherent_click == pytest.approx(1+(np.exp(1j*phase)*contrast).real, abs=5e-13)
    assert 2*size*marked_click == pytest.approx(1+(np.exp(1j*phase)*marked).real, abs=5e-13)


def test_continuum_coefficients_with_independent_high_precision_characteristic_polynomial():
    # In the nonorthogonal basis (S, FS), F swaps the columns exactly.
    with mp.workdps(58):
        k = mp.mpf('1e-9')
        a = (3+mp.cos(k))/4
        swap = mp.matrix([[0,1],[1,0]])
        length = study.PERIOD
        for resonant in (False, True):
            mu = mp.pi/(4*length) if resonant else mp.pi/32
            d, e = mp.exp(-1j*mu), -mp.exp(1j*mu)
            coin = mp.matrix([[d,(d-e)*a],[0,e]])
            matrix = mp.eye(2)
            for j in range(1,study.OPEN+1):
                angle = mp.pi*(1-mp.cos(mp.pi*j/study.OPEN))/4
                hop = mp.exp(-1j*angle)*(mp.cos(angle)*mp.eye(2)+1j*mp.sin(angle)*swap)
                matrix = hop*coin*matrix
            matrix = (swap*coin)**(length-study.OPEN)*matrix
            trace = matrix[0,0]+matrix[1,1]
            determinant = mp.det(matrix)
            assert abs(determinant+1j)<mp.mpf('1e-38')
            root = mp.sqrt(trace**2-4*determinant)
            eig = [(trace+r)/2 for r in (root,-root)]
            if resonant:
                velocity = abs(mp.arg(eig[0]/eig[1]))/(2*k*length)
                assert float(velocity) == pytest.approx(fringe_resonance.speed()[0], abs=3e-10)
            else:
                bright = min(eig,key=lambda z:abs(z-1))
                curvature = -mp.arg(bright)/(k*k*length)
                assert float(curvature) == pytest.approx(study.envelope_coefficient(), abs=3e-10)


def test_local_counter_and_invalid_records():
    assert study.triggered_service(0)==0
    assert study.triggered_service(2**24)==2**25
    for booking in (-1,1.5,True):
        with pytest.raises(ValueError):
            study.triggered_service(booking)
    shots = 2**25
    counts = np.tile([0,shots//2,shots//2], (4,1)).astype(np.int64)
    p = counts/shots
    assert fringe.check_panel(counts,shots,20,p,"dark") is None
    false = counts.copy()
    false[0] += [1,-1,0]
    for data in (false, counts[:3], counts[:,::-1], counts.astype(float), counts-1):
        with pytest.raises(ValueError):
            fringe.check_panel(data,shots,20,p,"dark")
    for forecast in (p+1, p*.5, np.full((4,3),np.nan)):
        with pytest.raises(ValueError):
            fringe.check_panel(counts,shots,20,forecast,"bad")
