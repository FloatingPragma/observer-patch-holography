from positive_reference_attachment import run

def test_current_ring_reference_and_complete_native_lift_rejection():
    r=run();assert all(r['checks'].values())
    assert r['attachment']['negative_filled_modes']==63
    assert abs(r['attachment']['charged_excitation_gap']-.23)<1e-12
    assert abs(r['attachment']['fixed_number_neutral_gap']-.46)<1e-12
    assert r['exact_rejecting_control']['completeness_defect_max']==1
    assert r['exact_rejecting_control']['fair_coin_weight_corrected_completeness_defect_max']==1
    assert not r['full_criterion_met']
