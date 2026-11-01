from itertools import product

import numpy as np
from scipy.linalg import expm

from ..field_symmetry import symmetry
from . import algebra,valleys


def test_symbolic_contact_projection_uses_general_complex_spinor():
    value=algebra.contact_certificate()
    assert value['complex_spinor_parameters']==4
    assert value['contact_commutator_exact']
    assert value['state_or_particle_number_projection_used'] is False


def test_two_site_contact_norm_and_full_CAR_Duhamel_bound():
    a=[symmetry._annihilator(4,j).matrix().toarray() for j in range(4)]
    n=[x.conj().T@x for x in a]
    u=np.array([1+2j,-2+1j,3-1j,1j])/4
    field=sum(u[j]*a[j] for j in range(4))
    composite=sum(u[j]*n[j^1]@a[j] for j in range(4))
    assert np.linalg.norm(composite,2)<=np.linalg.norm(u)+1e-14
    # Ordinary annihilators substituted for n_opposite*c have a different
    # one-site anticommutator rank (four rather than two).
    q=a[0]@a[0].conj().T+a[0].conj().T@a[0]
    b=n[1]@a[0]
    qb=b@b.conj().T+b.conj().T@b
    assert np.trace(q)!=np.trace(qb)
    hopping=a[0].conj().T@a[2]+a[2].conj().T@a[0]
    interaction=n[0]@n[1]+n[2]@n[3]
    U=.7
    for t in (.03,.2,-.1):
        full=expm(1j*t*(hopping+U*interaction))
        free=expm(1j*t*hopping)
        error=np.linalg.norm(full@field@full.conj().T-free@field@free.conj().T,2)
        assert error<=abs(t)*U*np.linalg.norm(u)+1e-14
        assert error>0


def test_all_eight_valleys_and_even_translation_phase():
    value=valleys.certificate()
    assert value['valley_count']==8
    assert sum(row['kinetic_orientation']==1 for row in value['valleys'])==4
    assert sum(row['kinetic_orientation']==-1 for row in value['valleys'])==4
    assert value['nearest_even_translation_phase_one']
    assert any(valleys.phase(b,(1,0,0))==-1 for b in valleys.VALLEYS)
    for b in valleys.VALLEYS:
        assert valleys.phase(b,(2,-2,4))==1
