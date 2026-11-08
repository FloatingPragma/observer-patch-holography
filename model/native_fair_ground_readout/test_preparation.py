from fractions import Fraction as F
import numpy as np
import pytest
from . import preparation as p


def test_complete_budget_preserves_resources_and_failures():
    b=p.budget()
    assert b['formal_raw_words']==25**11 and b['live_raw_words']==66*2**11
    assert b['native_coherent_pulses_upper']==12309
    assert not b['success_postselection'] and not b['growing_coherent_service_embedding']


@pytest.mark.parametrize('terminal',[False,-1,66,F(1)])
def test_invalid_decoded_terminal_is_rejected(terminal):
    with pytest.raises(ValueError):p.ground_word(terminal)


def test_zero_raw_outcome_preserved_on_entangled_columns():
    # Select an actually zero first slot, not a fabricated erasure label.
    edge=p.c.detector.tree()[0]
    zero=next(i for i,(_,pairs) in enumerate(p.m.fair_kernel(edge)) if not pairs)
    columns=np.eye(66,dtype=complex)[:,:2]
    result,records=p.ideal_branch((zero,)+(0,)*10,columns)
    assert result.shape==(66,2) and not np.any(result)
    assert len(records)==11


def test_actual_source_ground_replay():
    result=p.run()
    assert result['exact_original_fair_transfer_ground_checked'] is True
    assert result['all_terminal_routes_exactly_checked']==66
    assert result['ground_trace_distance_from_P_over_56']=='55/56'
    assert result['generic_host_status']=='RESOURCE_STOP'
    assert result['energy_readout_join']['complete_unnormalized_score_contrast_strict_lower']=='3/50'


def test_phase_targets_one_native_occupation_without_a_particle_label():
    join=p.scored_energy_join()
    assert join['exact_exchange_energy_of_ideal_excitation']=='50/33'
    assert join['source_words']==68 and join['native_coherent_pulses_upper']==3877451
    assert join['same_fair_source_ground_phase_and_readout_dictionary'] is True
    assert join['stable_particle_pole_or_mass_identified'] is False
