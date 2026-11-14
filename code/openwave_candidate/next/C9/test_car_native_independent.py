from audit_car_native_algebra import run


def test_exact_injective_outcomes_and_gauge_domain_rejection():
    result=run()
    assert all(result['checks'].values())
    assert result['checked_outcome_arrows']>1000
    assert result['wrap_control']['continuous_representative_phase_difference']>.1
    assert result['wrap_control']['discrete_phase_difference']<1e-12
