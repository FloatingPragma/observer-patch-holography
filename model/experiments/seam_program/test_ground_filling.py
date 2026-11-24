from .ground_filling import check_boards, local_bound


def test_partial_cell_cost_and_minimum_bridge_cuts():
    assert local_bound() == (5, 17)
    assert check_boards() == [3, 3, 3]
