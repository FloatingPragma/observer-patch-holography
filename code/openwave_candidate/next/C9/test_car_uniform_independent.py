from review_car_uniform_fiber import run


def test_composite_modulus_full_fiber_and_positive_reference():
    result=run(4)
    assert result['full_dimension']==80
    assert result['fiber_cardinality']==4
    assert all(result['checks'].values())
    assert result['native_reference_work']>1
    assert result['full_criterion_met'] is False
