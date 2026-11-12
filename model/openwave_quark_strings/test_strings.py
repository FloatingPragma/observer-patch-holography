import sympy as s
from .strings import full_string,full_gauss,shortest_paths,hop,inner,pointer_probabilities


def test_actual_CAR_signs_and_reference_are_necessary():
    state=full_string('uc',shortest_paths(1)[1])
    assert all(not full_gauss(state,cell,current) for cell in (0,1) for current in range(12))
    assert full_gauss(state,0,11,omit_reference=True)
    wrong={k:abs(v) for k,v in state.items()}
    assert any(full_gauss(wrong,cell,current) for cell in (0,1) for current in range(8))


def test_native_hopping_lands_on_both_spin_channels_by_axis():
    source=full_string('dc',());moved=hop(source)
    for dest,spin in [(1,0),(2,1),(4,1)]:
        path=shortest_paths(dest)[0]
        amplitude=inner(full_string('dc',path,spin),moved)
        assert s.conjugate(amplitude)*amplitude==s.Rational(1,16)
        assert inner(full_string('dc',path,1-spin),moved)==0
    assert inner(source,moved)==0
    assert inner(moved,moved)==s.Rational(3,8)


def test_phase_readout_is_complete_and_detects_conjugation():
    p=pointer_probabilities(s.Rational(2,3))
    anti=pointer_probabilities(-s.Rational(2,3))
    assert p[:2]==anti[:2] and p[2:]==tuple(reversed(anti[2:]))
    assert p[0]==s.Rational(1,8)
    assert pointer_probabilities(s.Rational(2,3),3)==(s.Rational(1,2),0,s.Rational(1,4),s.Rational(1,4))


def test_Duhamel_bound_on_independent_random_complete_carriers():
    import numpy as np
    from scipy.linalg import expm
    rng=np.random.default_rng(711)
    for size in (4,8,12):
        h0=np.diag(rng.normal(size=size));a=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
        tmat=(a+a.conj().T)/2;lam=.03;t=.01;u=np.eye(size,dtype=complex)[:,0]
        projector=np.eye(size);projector[0,0]=0
        actual=projector@expm(-1j*t*(h0+lam*tmat))@u
        first=-1j*lam*t*projector@tmat@u
        tau=np.linalg.norm(tmat,2);b=np.linalg.norm(h0,2)
        assert np.linalg.norm(actual-first)<=lam*t*t*tau*(2*b+lam*tau)/2
