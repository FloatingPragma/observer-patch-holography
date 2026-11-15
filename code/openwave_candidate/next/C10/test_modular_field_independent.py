from review_modular_field import run


def test_full_sources_keep_feedback_work_and_reject_code_compression():
    result=run()
    assert all(result['checks'].values())
    for experiment in result['experiments']:
        assert all(experiment['checks'].values())
        assert experiment['feedback_position_difference']>1e-5
        assert experiment['uniform_code_electric_leak']>1e-3
        assert experiment['state_error']<1e-12
        assert experiment['ledger_error']<1e-12
        assert abs(experiment['native_work']['MIX']['work'])>.05
    assert result['full_criterion_met'] is False
