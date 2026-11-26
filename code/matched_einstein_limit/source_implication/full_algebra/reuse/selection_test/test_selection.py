from fractions import Fraction

from .measure import record


def test_fixed_volume_deficit_and_observer_counts():
    data = record()
    assert Fraction(data["exact_population_deficit_lower_nats"]) > Fraction("0.015834")
    assert data["counts_objective_deficit_lower_nats"] > 0.0101
    assert abs(data["numerical_volume_difference"]) < 1e-10
    assert 0.0806 < data["numerical_objective_deficit"] < 0.0807
    assert data["product_plus_counts_optimum_then_reuse"] == [11983, 9336]
