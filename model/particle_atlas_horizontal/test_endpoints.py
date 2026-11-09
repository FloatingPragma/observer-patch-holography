from . import endpoints


def test_primitive_history_and_memory_creation_boundary():
    value = endpoints.certificate()
    assert value['primitive_branches'] == [
        {'weight': '1/2', 'input_seam_loads': [1, 0], 'output_seam_loads': [0, 1],
         'record': [0, 1], 'next_seam': [0, 5]},
        {'weight': '1/2', 'input_seam_loads': [1, 0], 'output_seam_loads': [1, 0],
         'record': [1, 0], 'next_seam': [0, 1]}]
    assert value['loop_C']['probability_on_initial_pair'] == '1/8'
    assert value['loop_E']['orientation_action'] == 'X'
    assert value['charge_commutation_columns_checked'] > 0
    assert not value['loop_is_the_literal_two_event_policy']
    assert not value['spatial_topological_defect_creation_derived']
