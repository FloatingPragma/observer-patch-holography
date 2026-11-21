from .asynchronous import exhaustive_patch, verify_native_service


def test_sharp_quotient_and_unguarded_fork():
    verify_native_service()
    result = exhaustive_patch()
    assert result["states"] == 3**9
    assert result["public_quotient_states"] == 3**4
    assert result["joined_one_step_peaks"] == 284
    assert result["guarded_endpoint"] in result["unguarded_endpoints"]
    assert result["unguarded_endpoints"][0][1] != result["unguarded_endpoints"][1][1]
