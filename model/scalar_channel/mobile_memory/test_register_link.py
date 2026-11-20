from .register_link import compile_link


def test_addressed_cage_link_uses_native_register_routes():
    result = compile_link()
    assert result["oriented_seam"] == result["central_support_ports"] == [10, 20]
    assert result["cycle_checks"] > 0
    assert all(result["routing_gate_counts"])
    assert result["literal_even_basis_columns"] == 32
    assert result["route_residual"] < 2e-11
    assert result["literal_generator_residual"] < 2e-11
    assert result["literal_occupation_residual"] < 2e-11
    assert result["literal_unitary_residual"] < 2e-11
