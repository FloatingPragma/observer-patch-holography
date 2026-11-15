from .source import phase_boundary,coherent_extension

def test_all_native_history_effects_and_phase_obstruction():
    x=phase_boundary()
    assert x['raw_history_blocks']==625 and x['phase_twin_total_variation_exact']=='0'

def test_actual_coherent_source_with_complete_export():
    x=coherent_extension()
    assert x['matrix_units_checked']==4356 and x['phase_twin_total_variation']>1-2e-12
    assert all(len(r['all_four_processor_abort_probabilities'])==4 for r in x['records'])
