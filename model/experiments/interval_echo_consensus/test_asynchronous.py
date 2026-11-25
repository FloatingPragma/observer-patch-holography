from .asynchronous import exhaustive_patch, verify_native_service


def test_sharp_quotient_and_unguarded_fork():
    verify_native_service()
    result = exhaustive_patch()
    assert result["states"] == 3**9
    assert result["public_quotient_states"] == 3**4
    assert result["joined_one_step_peaks"] == 284
    assert result["guarded_endpoint"] in result["unguarded_endpoints"]
    assert result["unguarded_endpoints"][0][1] != result["unguarded_endpoints"][1][1]


def test_complete_coherent_sweep_and_joint_phase():
    from .coherent_sweep import run

    result = run()
    assert result["legal_orders"] == 5
    assert result["changed_by_unguarded_reverse"] == 192
    assert {b["visit_layers"] for b in result["board_depths"]} == {12}
    assert result["phase"]["plus_joint"] == [0., .5, .5, 0.]
    assert result["phase"]["minus_joint"] == [.5, 0., 0., .5]
    assert result["phase"]["premature_probe_joint"] == [.25] * 4


def test_unordered_guarded_repairs_preserve_quantum_flag_subspaces():
    from .asynchronous import compare_flag_schedules, entangled_flag_read

    path = compare_flag_schedules()
    assert path["flag_patterns"] == 81 and path["fair_periodic_orders"] == 24
    assert path["cases"] == 1944 and path["cases_with_multiple_commits"] > 0
    assert path["largest_mean_uniform_attempts"] == "6"
    assert path["retry_example"]["reverse_order_writes_by_pass"] == [[0], [2]]
    correlations = entangled_flag_read()["joint_rail_correlations"]
    assert correlations["one_initial_flag_space"] == [1, 1]
    assert correlations["different_initial_flags"] == [1, 0]
