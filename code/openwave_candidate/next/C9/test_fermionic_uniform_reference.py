import numpy as np
from scipy.linalg import expm
from fermionic_uniform_reference import run,reduction,FermionicSharedRecords,Construction


def test_full_native_reference_and_every_intertwiner():
    r=run();assert all(r['checks'].values());assert r['parameters']['fiber_cardinality']==3
    assert r['full_H_minimum']<r['ground_energy']-1e-4
    for law,out in r['outputs'].items():
        assert abs(sum(x['probability'] for x in out['complete_records'])-1)<1e-12
        assert abs(sum(x['unnormalized_excitation_energy'] for x in out['complete_records'])-out['instrument_work'])<1e-12
        assert out['ground_survival']<1


def test_composite_modulus_and_reference_coherent_stationarity():
    m=FermionicSharedRecords(Construction(3,((0,1),(1,2),(2,0)),.23),7,4)
    masks,v,h,ops,fiber=reduction(m);assert fiber==4 and m.n==80
    for law in ops:
        for full,small in zip(m.operators[law],ops[law]):
            for (_,k),(_,l) in zip(full,small):assert np.max(abs(k@v-v@l))<1e-12
    e,w=np.linalg.eigh(h);psi=v@w[:,0]
    evolved=expm(-.37j*m.hamiltonian)@psi
    assert np.max(abs(np.outer(evolved,evolved.conj())-np.outer(psi,psi.conj())))<1e-12
