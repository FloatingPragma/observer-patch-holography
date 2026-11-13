from fractions import Fraction as F
import math
import pytest
from .source import log_decay_interval,complete_phase_blind_control,source_band

def test_complete_history_effects_and_phase_analyzer_positive_control():
    row=complete_phase_blind_control()
    assert row['complete_raw_histories']==15625
    assert row['total_effect_diagonal_exact']==['1']*8
    assert row['every_native_record_probability_difference_exact']=='0'
    control=row['external_analyzer_positive_control']
    assert F(control['plus_probability_exact'])!=F(control['phase_probability_exact'])

def test_outward_decay_bound_and_invalid_mass_substitution():
    for lam in (F(69,70),F(89,90),F(1)):
        lo,hi=log_decay_interval(lam)
        # Independent high-precision check does not define the exact enclosure.
        import decimal
        with decimal.localcontext() as ctx:
            ctx.prec=90
            value=-(decimal.Decimal(lam.numerator)/decimal.Decimal(lam.denominator)).ln()
            assert decimal.Decimal(lo.numerator)/decimal.Decimal(lo.denominator)<=value
            assert value<=decimal.Decimal(hi.numerator)/decimal.Decimal(hi.denominator)
    with pytest.raises(ValueError):log_decay_interval(0)

def test_same_source_band_cannot_hide_a_large_lepton_ratio():
    row=source_band(0)
    assert F(row['all_coherence_decay_rate_ratio_upper_exact'])==F(30,23)
    assert F(row['source_to_any_unitary_channel_trace_distance_lower_exact'])==F(1,252)
    assert all(x['phase_per_attempt_interval_exact']==['0','0'] for x in row['tested_modes'])
    assert all(F(x['four_event_multiplier_exact'])<F(x['two_event_multiplier_exact'])<1 for x in row['tested_modes'])
