from fractions import Fraction as F
import sympy as sy
from .source import packet,run,KAPPA,POSITION_RADIUS
from model.openwave_normalizable_particle_records.source import moment


def test_normalizable_trial_has_uniform_isolated_band_weight():
    d=packet()
    assert d['gap']-2*d['perturbation']>F(17,100)
    assert d['band_weight']>F(9997,10000)
    assert d['full_Q_upper']<F(1,2000000)
    assert moment(0,KAPPA)==1 and moment(2,KAPPA)==KAPPA*KAPPA/7
    assert d['position_variance']==655360
    assert d['position_tail']==F(5,512)
    assert d['energy'][1]-d['energy'][0]<F(1,10**30)


def test_native_current_records_keep_every_failure_and_no_filter_is_required():
    value=run();r=value['actual_trial_current_records']
    for key in ('conditional_on_prepared_trial_exact_probabilities','original_P7_over_7_input_exact_probabilities'):
        p={k:F(v)for k,v in r[key].items()}
        assert sum(p.values())==1 and all(v>=0 for v in p.values())
        assert p['native_proton_pion_fields']>0 and p['other_compact_N5']>0
        assert p['probe_failure']>0 and p['noncompact_N5']==0
    assert not r['exact_ground_band_filter_required_for_these_trial_records']
    assert r['current_channel_conditioning_preserves_spatial_packet_exactly']
    assert not value['optional_ground_band_conditioning']['finite_local_hardware_compilation_of_exact_band_filter_proved']
    assert not value['autonomous_on_shell_decay_or_lifetime_proved']


def test_covariant_inverse_rule_and_pole_channel_complement_are_positive():
    d=packet();G=d['G'];J=d['J']
    assert J==J.T
    # Positive inverse-rest action follows from the exact rest spectrum; the
    # actual source is the Gram-covariant matrix, not a supplied mass vector.
    assert all(F(x[0])>0 for x in d['rest']['all_seven_compact_source_intervals'])
    assert d['pole_channel_probability'][0]>0
    assert 0<d['pole_other_probability'][0]<d['pole_other_probability'][1]
    assert 0<d['pole_failure_probability'][0]<d['pole_failure_probability'][1]<1
    assert d['conditioned_position_tail']<F(1,40)
