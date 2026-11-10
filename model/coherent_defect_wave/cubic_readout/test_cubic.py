from fractions import Fraction as F
from itertools import product

import pytest
import sympy as s

from . import lattice as l
from .fourier import compress,expand
from .readout import experiment,higher_momentum_experiment
from .certificate import geometric_certificate


def test_three_dimensional_complete_records_and_nonflat_phase_support():
    result=experiment()
    assert result['full_tick_amplitude_support']==9344
    assert len(result['complete_joint_records'])==32
    assert sum(F(r['probability_given_positive_Gauss']) for r in result['complete_joint_records'])==1
    assert [r['complex_first_moment'] for r in result['momenta']]==[
        ['-29/30','0'],['0','29/30'],['0','29/30'],['29/30','0'],
        ['0','29/30'],['29/30','0'],['29/30','0'],['0','-29/30']]
    prep=result['preparation']
    assert sum(map(F,prep['original_trial_joint_record_weights']))+F(prep['retained_other_flags_total_probability'])==1
    assert result['physical_scope']['rest_mass_extracted'] is False


def test_fourier_orbit_encoding_equals_full_translation_projector():
    source=l.packet(); evolved=l.tick(source,4)
    for k in ((0,0,0),(1,0,1)):
        compact=compress(evolved,k)
        full=l.projector(evolved,k)
        assert expand(compact)==full
        assert compact.norm_squared()==full.norm_squared()==F(1,8)


def test_finer_cube_has_all64_momenta_and_complete_phase_records():
    result=higher_momentum_experiment()
    assert len(result['momenta'])==64 and len(result['complete_joint_records'])==256
    assert sum(F(row['probability_given_positive_Gauss']) for row in result['complete_joint_records'])==1
    assert all(F(row['unitary_variance'])>=F(59,900) for row in result['momenta'])
    assert result['probability_any_defect_changes_face']=='1/450'
    assert sum(F(row['probability_per_original_ambient_trial']) for row in result['complete_joint_records'])+F(result['Gauss_preparation']['retained_other_flags_total_probability'])==1


def test_local_encoding_has_the_literal_source_inverse():
    for word in l.representatives():
        for index in (l.mode(5,3,1,7),l.mode(0,0,0,0)):
            changed=() if word==1 else (((0,0,0),word),)
            source=l.State({((0,0,0),index,changed):(1,0)})
            assert l.onsite(l.onsite(source),True)==source


def test_encoded_hopping_and_onsite_equal_two_cell_parent_execution():
    from ..species.spatial_fock import PairFockState,positive_basis,exchange,onsite
    from ..species.fock import MODES
    def embed(state):
        out={}
        for (position,index,field),(real,imag) in state.amplitudes:
            assert position in ((0,0,0),(1,0,0))
            cell=position[0];words=tuple(dict(field).get((i,0,0),1) for i in range(2))
            for flips in product(range(2),repeat=2):
                links=tuple(g^(l.STAR if flip else 0) for g,flip in zip(words,flips,strict=True))
                phase=l.chi(index)**flips[cell]
                out[links,(cell*MODES+index,)]=phase*(s.Rational(real.numerator,real.denominator)+s.I*s.Rational(imag.numerator,imag.denominator))
        return PairFockState(out,4*state.norm_denominator)
    index=l.mode(5,3,0,0)
    source=l.State({((0,0,0),index,(((1,0,0),3),)):(1,0)})
    parent=positive_basis((1,3),(index,))
    assert embed(source)==parent
    assert embed(l.hopping_layer(source,4,0,0))==exchange(parent)
    assert embed(l.onsite(source))==onsite(onsite(parent,0),1)


def test_volume_and_joint_translation_controls():
    result=geometric_certificate()
    assert result['bare_wave_vs_complete_translation_squared_difference']=='1/225'
    assert result['fixed_spacing_one_tick_volume_comparison']['translated_complete_sparse_vectors_equal']


def test_invalid_domains_do_not_silently_change_the_program():
    with pytest.raises(ValueError):l.tick(l.packet(),3)
    with pytest.raises(ValueError):l.tick(l.packet(),True)
    with pytest.raises(ValueError):compress(l.packet(),(0,0,2))
    with pytest.raises(ValueError):l.hopping_layer(l.packet(),4,3,0)
    with pytest.raises(ValueError):l.State({((0,0,0),0,()):(0.5,0)})
    with pytest.raises(ValueError):l.State({((0,0,0),0,(((1,0,0),1),)):(1,0)})
