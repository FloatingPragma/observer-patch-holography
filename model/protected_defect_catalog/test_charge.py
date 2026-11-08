"""Complete source charge, exact selection rules and negative controls."""
import numpy as np
import pytest

from . import charge


def test_complete_charge_and_creation_on_actual_source():
    op = charge.operators()
    q2 = op['twice_charge']
    for q, twice_a in op['twice_creation'].items():
        assert np.array_equal(q2 @ twice_a-twice_a @ q2, 2*q*twice_a)
    report = charge.report()
    assert report['charge_spectral_multiplicities'] == dict(zip(map(str, range(-4,5)), [1,8,28,56,70,56,28,8,1]))
    assert report['ideal_created_charge_support'] == [-1,1]
    assert report['unique_selected_ground_charge'] == 0
    assert all(row['physical_mass'] is None for row in report['creation_channels'])


def test_local_Z_is_not_itself_a_charge_eigenoperator():
    op = charge.operators(); q = op['twice_charge']; z = op['local_Z']
    commutator = q @ z-z @ q
    assert np.any(commutator)
    assert not np.array_equal(commutator, 2*z)
    assert not np.array_equal(commutator, -2*z)


def test_wrong_uniform_charge_is_rejected_by_source_H():
    op = charge.operators()
    indices = np.arange(256)
    wrong = sum((np.eye(256, dtype=np.int64)[indices ^ (1 << x)] for x in range(8)),
                np.zeros((256,256), dtype=np.int64))
    assert np.any(wrong @ op['H']-op['H'] @ wrong)


def test_nonconserving_source_mutation_detected(monkeypatch):
    original = charge.source.matrix
    def changed():
        h, scale = original(); h=h.copy()
        h += np.diag([1-2*(word & 1) for word in range(256)])
        return h, scale
    monkeypatch.setattr(charge.source, 'matrix', changed)
    with pytest.raises(ArithmeticError, match='does not conserve'):
        charge.operators()


def test_full_source_including_empty_and_doubly_occupied_sites():
    assert charge.spin_flips(0,2) == ()
    assert charge.spin_flips(15,2) == ()
    assert charge.spin_flips(1,2) == (2,)
    result = charge.full_source_charge()
    assert result['complete_all_seam_full_Fock_columns'] == 6976
    assert result['source_charge_pulled_through_spin_and_Marshall_maps']


def test_local_two_outcome_instrument_prepares_coherent_charge_branches():
    op=charge.operators(); identity=np.eye(256,dtype=np.int64)
    total=sum((a.T@a for a in op['twice_creation'].values()), np.zeros_like(identity))
    assert np.array_equal(total,4*identity)
    x=identity[np.arange(256)^1]
    basis=op['walsh']; charges=np.array(op['charges'])
    p0=basis[:,charges==0]@basis[:,charges==0].T
    for q, twice_a in op['twice_creation'].items():
        assert np.array_equal(twice_a,(identity+q*op['signs'][0]*x)@op['local_Z'])
        pq=basis[:,charges==q]@basis[:,charges==q].T
        # Spectral projectors here have the shared exact denominator256.
        assert np.array_equal(2*pq@op['local_Z']@p0,256*twice_a@p0)


def test_fine_record_grouping_loses_within_charge_coherence():
    op=charge.operators(); basis=op['walsh']/16
    labels=[i for i,q in enumerate(op['charges']) if q==1]
    v=(basis[:,labels[0]]+basis[:,labels[1]])/np.sqrt(2)
    coherent=np.outer(v,v)
    fine=(np.outer(basis[:,labels[0]],basis[:,labels[0]])
          +np.outer(basis[:,labels[1]],basis[:,labels[1]]))/2
    assert np.linalg.norm(coherent-fine) > .7
    assert np.isclose(np.trace(coherent),np.trace(fine))
    assert np.isclose(np.trace(coherent@coherent),1)
    assert np.isclose(np.trace(fine@fine),.5)
