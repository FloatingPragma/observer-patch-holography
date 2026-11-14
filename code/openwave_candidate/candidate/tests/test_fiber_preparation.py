import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.geometry import flat_ring
from oph_openwave.fiber_preparation import FiberPreparation


def test_composite_modulus_gauss_and_all_native_intertwiners():
    p=FiberPreparation(Construction(3,((0,1),(1,2),(2,0)),mass=.23),6)
    full,v=p.tiny_isometry()
    assert np.allclose(v.conj().T@v,np.eye(p.n))
    assert np.allclose(full.hamiltonian@v,v@p.reduced.hamiltonian)
    for seam in range(3):
        for (r,k),(s,l) in zip(full.operators[seam],p.reduced.kraus_for(seam)):
            assert r==s and np.allclose(k@v,v@l)
    for pos in range(p.n):
        for coordinate in range(6):
            word=p.record_word(pos,(coordinate,))
            assert full.gauss_residue((pos,word))==(1,0,0,0,0,0)


def test_record_exit_retains_complete_continuation():
    p=FiberPreparation(Construction(3,((0,1),(1,2),(2,0))),6)
    rho=p.reduced.initial_density(); full,outputs=p.record_measurement(rho,0)
    assert sum(x['probability'] for x in outputs)==pytest.approx(1)
    assert sum(x['complement_weight'] for x in outputs)==pytest.approx(5/6)
    for x in outputs:
        assert x['probability']==pytest.approx(1/6)
        assert x['code_weight']==pytest.approx(1/36)
        state=x['poststate_subnormalized']/x['probability']
        evolved,ledger=full.step(state,.12,.2)
        assert np.trace(evolved)==pytest.approx(1)
        assert sum(z['probability'] for z in ledger)==pytest.approx(1)


def test_cut_preserves_code_and_preparation_does_not_enumerate_carrier():
    p=FiberPreparation(flat_ring(1,1),6)
    assert p.resources()['fiber_cardinality']==6**40
    assert not p.resources()['native_only_preparation']
    word=p.record_word(86,(0,)*40)
    g=[int(i==86) for i in range(87)]
    for (a,b),f in zip(p.edges,word):g[a]+=f;g[b]-=f
    assert tuple(x%6 for x in g)==(1,)+(0,)*86
    with pytest.raises(ValueError):p.tiny_isometry()
    rho=p.reduced.initial_density(); branches=p.closed_cut_measurement(rho,(0,1))
    assert sum(x['probability'] for x in branches)==pytest.approx(1)
    assert all(x['complement_weight']==0 for x in branches)
    with pytest.raises(ValueError):p.record_word(True,(0,)*40)
    with pytest.raises(ValueError):p.record_word(0,(True,)*40)
