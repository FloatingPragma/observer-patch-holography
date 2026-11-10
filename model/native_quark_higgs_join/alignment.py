"""Normalized admitted states and common histories with actual preparation."""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from pathlib import Path

from model.native_higgs_filtered_state.model import filtered_bounds
from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.model import certificate as preparation
from model.native_charge_resolution.bounds import root, round_up
from .bounds import response, atomic, field_interval
from model.native_quark_response.full import encoded

ROOT = Path(__file__).resolve().parents[2]
PURE_LEMMA_INPUTS = {
    'model/native_higgs_filtered_state/model.py':
        '3de6a3df51dcd09d7cbff10f5b9dd0f63a40e4674092789110cdfe912b07a095',
    'model/native_higgs_filtered_state/PROOF.md':
        'f6c1dd98cbfb0e7fa3fb393ebdccb3cc318af854bcac3dc7e5118cca9f9838e5',
}


@lru_cache(None)
def certificate():
    for name, expected in PURE_LEMMA_INPUTS.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected
    prep = preparation()['finite_energy_agreement_filter']
    pair_fidelity = F(prep['conditional_ground_pair_fidelity_lower'])
    epsilon = root(2*(1-pair_fidelity), True)
    delta = F(prep['conditional_one_source_trace_distance_upper'])
    rows = {}
    for family in ('uc', 'dc', 'Q'):
        source = response(family)
        result = filtered_bounds(source['source_excitation_relative_actual_ground'],
                                 source['exact_ground_elementary_bin_probability'], epsilon,
                                 source_remainder=2*F(source['full_hopping_norm']))
        # This pure inequality is reused with this package's own source,
        # full reference and original two-source preparation. No lepton
        # state, narrow-bin selector, depth-three tree or port is transferred.
        result.pop('actual_dressed_W_expectation_difference_upper')
        result['diagnostic_binary64'].pop('actual_dressed_W_expectation_difference_upper')
        qlo = F(source['preparation_success_lower'])
        for channel in ('elementary', 'composite_B_over_2'):
            lower = F(result['actual_'+channel+'_in_bin_probability_lower'])
            result['original_attempt_'+channel+'_in_bin_probability_lower'] = str(-round_up(-qlo*lower))
            result['fair_setting_joint_'+channel+'_in_bin_probability_lower'] = str(-round_up(-qlo*lower/2))
        qa = F(result['ideal_elementary_in_bin_probability_lower'])
        qb = F(result['ideal_composite_B_over_2_in_bin_probability_lower'])
        energy = F(result['energy_center'])
        residual = F(result['ideal_relative_vector_residual_upper'])
        whi = F(source['exact_ground_elementary_bin_probability'][1])
        upper_amplitude = root(whi, True)
        hit_intervals = {
            'elementary': (F(result['actual_elementary_in_bin_probability_lower']),
                           min(F(1), (upper_amplitude+epsilon)**2)),
            'composite_B_over_2': (F(result['actual_composite_B_over_2_in_bin_probability_lower']),
                                   min(F(1), ((energy+residual)*upper_amplitude+2*epsilon)**2/4)),
        }
        local = atomic(family)
        composite_total = field_interval(local['composite_norm_squared']/4, local['source_ground'])
        eta = F(source['exact_ground_trial_vector_error'])
        totals = {'elementary': (F(1, 2)-delta, F(1, 2)+delta),
                  'composite_B_over_2': (max(F(0), composite_total[0]-eta-delta),
                                         min(F(1), composite_total[1]+eta+delta))}
        records = {}
        for channel in totals:
            lo, hi = totals[channel]
            bl, bu = hit_intervals[channel]
            laws = {'pump_failure': (1-hi, 1-lo), 'bin_hit': (bl, bu),
                    'outside_bin': (max(F(0), lo-bu), min(F(1), hi-bl))}
            assert all(0 <= a <= b <= 1 for a, b in laws.values())
            records[channel] = {'prepared_conditional': {key: encoded(val) for key, val in laws.items()},
                                'original_attempt_joint_fair_setting': {
                                    key: encoded((qlo*val[0]/2, val[1]/2)) for key, val in laws.items()},
                                'failed_pump_energy_hit': '0'}
        assert totals['composite_B_over_2'][1] < totals['elementary'][0]
        assert hit_intervals['composite_B_over_2'][1] < hit_intervals['elementary'][0]
        independent = min(F(1), F(result['ideal_normalized_state_trace_distance_upper'])+delta/qa+delta/qb)
        result.update(absolute_active_source_plus_reference_bin=source['absolute_active_source_plus_reference_bin'],
                      independent_source_reference_flagged_normalization_trace_bound=str(round_up(independent)),
                      reference_success_cost='4', reference_failure_cost='0',
                      source_and_common_retained_preparation_environment_included=True,
                      different_classical_setting_label_in_compared_state=False,
                      same_history_policy_may_read_pump_setting=False,
                      common_history_total_variation_upper=result['actual_normalized_state_trace_distance_upper'],
                      complete_setting_records=records,
                      distinct_pump_and_bin_success_probabilities_certified=True,
                      original_pump_success_rates_equated=False,
                      complete_fair_setting_pump_energy_labels=8,
                      all_preparation_failure_records_retained=True)
        rows[family] = result
    return dict(pure_inequality_inputs=PURE_LEMMA_INPUTS,
                pure_inequality_source_state_and_reference_transferred=False,
                original_two_source_preparation=True, additional_P0=False,
                conditional_pair_ground_fidelity_lower=str(pair_fidelity),
                common_retained_environment_vector_error_upper=str(epsilon),
                actual_source_trace_distance_upper=str(delta),
                settings=['A', 'C/2'], source_clock='h_active+h_partner+4N_ref',
                complete_law=['1-q, every fine physical preparation failure',
                              'q/2 Tr(K_X rho K_X^dagger), energy flag zero',
                              'q/2 Tr(P X rho X^dagger), energy flag one',
                              'q/2 Tr((I-P) X rho X^dagger), energy flag zero',
                              'zero for failed pump and energy flag one, each X=A,C/2'],
                conditional_history_domain='same source/reference/preparation-environment instrument policy after bin acceptance; setting label excluded',
                unconditional_complete_laws_claimed_close=False,
                rows=rows)
