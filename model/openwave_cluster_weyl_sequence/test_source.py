from fractions import Fraction as F
import pytest
from .source import cluster_identity,tail,sequence_member,uniform_data


def test_full_separated_CAR_intertwiner_retains_each_source_component():
    rows=cluster_identity();assert len(rows)==2
    assert all(row['terms'][kind]['exact_intertwiner']for row in rows for kind in('base','triquark','meson'))
    assert all(row['terms']['triquark']['all_output_occupations']>500 for row in rows)


def test_full_tail_bounds_and_normalizable_sequence_converge():
    d=uniform_data()
    for kind in('proton','pion'):
        assert d[kind]['floor']>d[kind]['energy_upper']
        assert 0<tail(kind,16)<tail(kind,8)<tail(kind,4)<F(1,10**10)
    rows=[sequence_member(r,2**(4*r))for r in(4,8,16)]
    assert all(F(b['full_H5_energy_residual_upper_GeV'])<F(a['full_H5_energy_residual_upper_GeV'])for a,b in zip(rows,rows[1:]))
    assert all(r['minimum_cluster_gap_cells']>3 and r['exact_normalized_product_norm_squared']=='1'for r in rows)
    assert F(rows[-1]['full_H5_energy_residual_upper_GeV'])<F(1,10**8)


def test_invalid_sequence_domains_cannot_silently_remove_a_tail():
    for radius in(0,-1,F(3,2),True):
        with pytest.raises(ValueError):tail('proton',radius)
    for args in((1,4),(4,1),(2,F(5,2))):
        with pytest.raises(ValueError):sequence_member(*args)
