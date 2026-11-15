import numpy as np
import pytest
from oph_openwave.construction import Construction
from oph_openwave.modular_native import LocalCoarseModular
from oph_openwave.many_load_native import ManyLoadNative


def test_exact_one_load_recovery_including_cyclic_records():
    c=Construction(3,((0,1),(1,2),(2,0)),.23)
    original=LocalCoarseModular(c,3);many=ManyLoadNative(c,(1,0,0,0,0,0),3)
    permutation=[many.index[(tuple(int(v==pos) for v in range(6)),f)] for pos,f in original.states]
    assert np.array_equal(many.hamiltonian[np.ix_(permutation,permutation)],original.hamiltonian)
    for e,ops in enumerate(original.operators):
        actual={tuple(record.values()):k[np.ix_(permutation,permutation)] for record,k in many.operators['MIX'][e]}
        for record,k in ops:assert np.array_equal(actual[tuple(record.values())],k)


def test_three_unit_tree_complete_cp_gauss_and_actual_evolution():
    model=ManyLoadNative(Construction(3,((0,1),(0,2)),.23),(3,0,0,0,0),7)
    assert model.n==35
    assert np.allclose(model.hamiltonian,model.hamiltonian.T)
    g=model.gauss(model.states[0])
    for state in model.states:assert model.gauss(state)==g
    for law in ('C','E','MIX'):
        for ops in model.operators[law]:
            assert np.allclose(sum(k.T@k for _,k in ops),np.eye(model.n))
        final,ledger=model.step(model.initial_density(),.17,.29,law)
        assert np.isclose(np.trace(final),1)
        assert np.linalg.eigvalsh(final).min()>-1e-12
        assert np.isclose(sum(r['probability'] for r in ledger),1)
    # Moving three units in one Bose hop has sqrt(3), not unit amplitude.
    assert np.isclose(max(abs(model.hamiltonian[:,0])),np.sqrt(3))


def test_complete_cut_decoder_and_three_unit_alias():
    for modulus in (3,7):
        model=ManyLoadNative(Construction(2,((0,1),),.23),(3,0,0),modulus)
        g=model.gauss(model.states[0])
        for state in model.states:
            for region in ((0,),(1,),(2,),(0,2),(0,1,2)):
                q=model.regional_charge(state,region);cut=model.closed_cut_record(state,region)
                assert (q+cut-sum(g[v] for v in region))%modulus==0
                if modulus==7:assert model.decode_bounded((3-q)%7)==3-q
        if modulus==3:
            assert 3%modulus==0
            with pytest.raises(ValueError):model.decode_bounded(0)
        else:
            assert model.decode_bounded(3)==3
            assert model.decode_bounded(4)==-3


def test_two_unit_native_balancing_is_not_independent_unit_sampling():
    model=ManyLoadNative(Construction(2,((0,1),)),(2,0,0),7)
    for coin in (0,1):
        target,j,d=model.target(model.states[0],0,coin)
        assert target[0]==(1,1,0) and j==1 and d==2
