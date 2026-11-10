import pytest
from .lane_collision import certificate,countershift


def test_complete_coin_and_two_incoming_objects_survive_collision():
    result=certificate()
    assert result['complete_local_coin_basis_columns_checked']==121
    assert result['two_object_local_basis_states_retained']==100
    assert result['coincident_two_object_lane_labels']==[0,0,1,4,0,0]


def test_countershift_rejects_truncated_or_outside_lane_domain():
    for invalid in ((0,0,0,0),(0,0,0,0,0),(11,0,0,0,0,0)):
        with pytest.raises(ValueError):countershift(invalid)
