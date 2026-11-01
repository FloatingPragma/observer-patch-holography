from fractions import Fraction as F

import numpy as np

from . import spin as s
from ..retarded_matter_response import response


def test_entire_N5_and_N3_source_spin_actions_and_multiplicities():
    five=s.representation_certificate(5);three=s.representation_certificate(3)
    assert five['full_dimension']==4368
    assert five['spin_multiplicities']=={'5/2':56,'3/2':504,'1/2':1008}
    assert three['spin_multiplicities']=={'3/2':56,'1/2':168}
    assert five['all_source_spin_intertwiners_exact']


def test_high_spin_projection_needs_exact_one_fifth_coefficient():
    assert s.projection_certificate()['D_high_Gram']=='P_plus/5'
    t=s.raising(4,1);d=response.dilation('ground','down').vectors
    u=tuple(x[:56] for x in response.dilation('ground','up').vectors)
    for down,up in zip(d,u,strict=True):
        assert not np.any(t@(5*down-t.T@up))
    assert any(np.any(t@(4*down-t.T@up)) for down,up in zip(d,u,strict=True))


def test_all_ground_creator_and_annihilator_columns_are_exact_ladders():
    value=s.family_certificate()
    assert value['creator_columns_checked']==value['annihilator_columns_checked']==80
    assert value['source_ground_norms_squared_relative_to_G']==[1,4,24,144,576]
    for row in value['normalized_squared_ladder_factors']:
        assert F(row['up_high'])+F(row['down_high'])==F(6,5)
        assert F(row['up_low'])+F(row['down_low'])==1


def test_coherent_ground_superposition_is_not_replaced_by_a_classical_mixture():
    states=s.ground_family()
    # These different M components must coexist as source amplitudes.
    coherent=s.combine((1,states[0]),(2,states[2]),(-1,states[4]))
    for mode in range(16):
        for creation in (False,True):
            assert s.car(coherent,mode,creation)==s.combine(
                (1,s.car(states[0],mode,creation)),
                (2,s.car(states[2],mode,creation)),
                (-1,s.car(states[4],mode,creation)))
    assert s.inner(coherent,coherent)==(1296*(1+4*24+576),0)


def test_direct_full_source_twirl_keeps_spatial_coherences_and_first_moment():
    value=s.twirled_moments()
    assert value['all_16_by_16_spatial_and_spin_entries_checked']
    assert value['retarded_first']=='(A tensor I_spin)/16+I_16/4'
    assert value['local_addition_probability']=='3/4'
