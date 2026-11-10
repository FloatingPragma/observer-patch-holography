from fractions import Fraction
from itertools import product
import pytest
from .model import column,exact_ground,preparation_records
from ..periodic_family.model import local_yukawa


def test_compiler_matches_full_native_columns_and_occupied_controls():
    masks=[0,2**30-1]+[1<<i for i in range(30)]+[(1<<i)|(1<<(i+15)) for i in range(15)]
    for mask,scalar in product(masks,range(5)):
        target=dict(local_yukawa(mask,scalar))
        if scalar:target[mask,scalar]=target.get((mask,scalar),0)+1
        assert column((mask,scalar))==target
    with pytest.raises(ValueError):column((2**30,0))


def test_exact_complete_graph_ground_and_preparation_scope():
    r=exact_ground();lo,hi=map(Fraction,r['ground_energy_interval'])
    assert r['preparation_sector_dimension']==97878
    assert r['Hamiltonian_nonzero_entries']==933620
    assert Fraction(-1114103,100000)<lo<hi<Fraction(-1114102,100000)
    assert Fraction(r['ground_preparation_probability_lower_bound'])>0
    assert not r['global_ground_across_other_conserved_sectors_proved']


def test_full_source_remainder_and_complete_gauge_number_records():
    r=preparation_records();rows=r['complete_grouped_records']
    assert r['full_source_H5_input_squared_norm']==409247068
    assert set(rows)=={'up','down','charged_lepton','all_other_outputs'}
    assert [rows[k]['full_detector_rank'] for k in ('up','down','charged_lepton')]==[288,288,32]
    assert all(Fraction(row['probability_interval'][0])>0 for row in rows.values())
    assert sum(Fraction(v['probability_interval'][0]) for v in rows.values())<=1
    assert sum(Fraction(v['probability_interval'][1]) for v in rows.values())>=1


def test_exact_cyclic_spectrum_keeps_every_original_weight():
    from .quotient import spectrum
    result=spectrum();rows=result['complete_energy_records']
    assert result['equitable_dimension']==240
    assert result['full_physical_cyclic_dimension']==len(rows)==36
    assert all(Fraction(row['probability_interval'][0])>0 for row in rows)
    assert sum(Fraction(row['probability_interval'][0]) for row in rows)<=1
    assert sum(Fraction(row['probability_interval'][1]) for row in rows)>=1
    assert result['outside_cyclic_preparation_sector_record']['probability_interval']==['0','0']
    assert not result['numerical_eigensolver_used']
    assert Fraction(result['ground_filter_success_interval'][0])>Fraction(4772,10000000)


def test_proposed_partition_cannot_hide_a_source_row(tmp_path,monkeypatch):
    import numpy as np
    from . import quotient as q
    with np.load(q.PARTITION,allow_pickle=False) as archive:
        colors=archive['colors'].copy();digest=archive['source_digest'].copy()
    colors[1]=colors[0]
    wrong=tmp_path/'wrong.npz';np.savez_compressed(wrong,colors=colors,source_digest=digest)
    q.quotient.cache_clear()
    with monkeypatch.context() as patch:
        patch.setattr(q,'PARTITION',wrong)
        with pytest.raises(ValueError,match='intertwiner|contiguous'):
            q.quotient()
    q.quotient.cache_clear()


def test_interval_division_cannot_cross_zero():
    from .quotient import divide
    with pytest.raises(ValueError,match='zero-crossing'):
        divide((1,2),(-1,1))
