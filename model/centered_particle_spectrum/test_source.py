import numpy as np
import pytest
from . import source,packet


def test_full_source_symmetry_and_particle_hole():
    value=source.symmetry_certificate()
    assert value['signed_cube_group_order']==96
    assert value['origin_stabilizer_order']==6
    assert value['all_eight_particle_hole_CAR_intertwiners_checked'] is True


@pytest.mark.parametrize('kind,expected',[('ground',(4900,120,32)),('local',(3920,656,228))])
def test_complete_projector_ranges(kind,expected):
    witness=packet.load_witness();r=source.reduced(kind,tuple(map(int,witness[kind+'_columns'])))
    assert (len(source.masks(r['up'],r['down'])),r['orbit_dimension'],r['dimension'])==expected
    assert np.array_equal(r['basis'].T@(r['sizes'][:,None]*r['basis']),r['gram'])
    assert np.array_equal(r['projector']@r['basis'],r['denominator']*r['basis'])
    for op in ['spin_plus','eta_minus']:
        assert not np.any(source.ladder(r['up'],r['down'],op)@r['full_basis'])


def test_actual_car_image_has_no_discarded_component():
    witness=packet.load_witness();g=tuple(map(int,witness['ground_columns']));c=tuple(map(int,witness['local_columns']))
    assert source.transition(g,c).shape==(228,32)


@pytest.mark.parametrize('up,down',[(True,4),(4,False),(-1,4),(9,4),(4,2.0)])
def test_sector_admission_rejects_noninteger_or_outside(up,down):
    with pytest.raises(ValueError):source.masks(up,down)


def test_integer_overflow_is_rejected():
    with pytest.raises(OverflowError):source.int_product(np.array([[2**62]],dtype=np.int64),np.array([[4]],dtype=np.int64))


def test_duplicate_projector_columns_are_rejected():
    with pytest.raises(ValueError):source.reduced('ground',(0,)*32)
