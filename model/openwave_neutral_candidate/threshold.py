"""Finite variational thresholds relative to the actual interacting ground."""
import sympy as s
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.lattice import gap_data
from model.coherent_defect_wave.chiral_joint.collective_spin_response.transfer import data as ground_data
from model.coherent_defect_wave.chiral_joint.atomic_spin_geometry.cube import spectrum as cube_spectrum
from model.coherent_defect_wave.chiral_joint.ground_composite_response.spectrum import response as local_response
from model.coherent_defect_wave.chiral_joint.local_reference.reference import certificate as local_reference
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.transfer import data as charged_data
from model.openwave_particle_readouts.readouts import ANCHOR_GEV
from . import source


def certificate():
    atomic=gap_data();ground=ground_data();reference=local_reference()
    assert reference['raw_L']=='1' and reference['raw_B3']=='0' and reference['fermion_parity']=='odd'
    coupling=ground['coupling'];tau=ground['tau'];v=coupling*tau;gap=atomic['gap']
    assert coupling==s.Rational(1,2**34) and tau==360 and gap==ground['gap']
    assert v<gap/5
    r4=25*v**4/gap**3
    cube_gap=s.Rational(cube_spectrum()['gap_interval'][0])
    spin_gap=2*coupling**2*ground['jlo']*cube_gap-2*r4
    assert spin_gap>0
    rows=[]
    for kind,raw_L in (('remove',7),('add',9)):
        local=local_response(kind,0)
        record=local['complete_target_energy_records'][0]
        trial_gap=tuple(map(s.Rational,record['excitation_energy_interval']))
        weight=tuple(4*s.Rational(x) for x in record['joint_probability_interval'])
        assert trial_gap[0]>0 and weight[0]>0
        # A blank-link atomic ground trial gives E_ground<=8E0, without a +v loss.
        low=gap-v;high=trial_gap[1]+v
        assert 0<low<high
        rows.append({'source_neutral_operator':kind,'target_raw_L':raw_L,'target_raw_B3':0,
            'target_fermion_parity':'odd','native_gauge_charge_change':'zero for all twelve actual currents',
            'local_ground_relative_trial_gap_interval':list(map(str,trial_gap)),
            'local_lowest_line_probability_for_fixed_setting_interval':list(map(str,weight)),
            'unchanged_antiperiodic_source_sector_threshold_interval':list(map(str,(low,high))),
            'periodic_variant_zero_momentum_sector_threshold_interval':list(map(str,(low,high))),
            'threshold_decimal_interval':[float(low),float(high)],
            'periodic_source_ground_momentum':[0,0,0],
            'periodic_zero_momentum_sector_nonempty':True,
            'lowest_sector_projector_stationary_under_its_actual_full_H':True,
            'energy_decay_to_lower_state_in_same_sector_possible':False,
            'decay_to_L8_ground_forbidden_by_source_L_conservation':True,
            'arbitrary_charge_carrying_environment_stability_claimed':False,
            'gauge_L_B3_and_translation_preserving_perturbation_norm_budget':str(low/4),
            'remaining_threshold_under_that_perturbation_lower':str(low/2)})
    # A common illustrative energy conversion has explicit electron-anchor ancestry.
    ec=charged_data('ec')['local']['rows'][0]['source_energy_difference_interval']
    anchor=(ec[0]+ec[1])/2;scale=ANCHOR_GEV/anchor
    for row in rows:
        interval=tuple(s.Rational(x)*scale for x in row['periodic_variant_zero_momentum_sector_threshold_interval'])
        row['conditional_GeV_energy_threshold_interval']=list(map(str,interval))
        row['conditional_GeV_decimal_interval']=[float(x) for x in interval]
    return {'unchanged_source':'complete N2 PW119 H0+2^-34 T with original antiperiodic fermions',
            'periodic_variant':'same full H0 and link matrices; only twelve fermion seam signs reversed',
            'actual_ground_reference_for_each_source':'unique full physical ground from the same K2 and quartic-isolation bounds',
            'ground_raw_L':8,'ground_raw_B3':0,'ground_parity':'even',
            'complete_atomic_gap_lower':str(gap),'full_hopping_norm_upper':str(tau),
            'source_hopping_coefficient':str(coupling),'perturbation_norm_upper':str(v),
            'variant_quartic_effective_H_remainder_upper':str(r4),
            'inherited_ground_isolation_positive_check':str(spin_gap),
            'ground_representation_transfer':'L/B3 and translation commute with H0,T,P0; canonical graph uniqueness preserves their P0 actions',
            'ground_energy_interval':'[8E0-v,8E0], with all energy constants retained',
            'sectors':rows,'zero_momentum_trial':source.momentum_trial_certificate(),
            'calibration':{'electron_anchor_GeV':str(ANCHOR_GEV),'native_local_ec_gap_anchor':str(anchor),
                'GeV_per_native_energy':str(scale),'scale_is_supplied_empirical_calibration':True,
                'neutral_threshold_was_used_to_choose_scale':False,
                'GeV_threshold_is_observed_relativistic_rest_mass':False},
            'complete_formal_spectral_preparation':{
                'input':'maximally mixed complete source with neutral apparatus references retained',
                'success':'rank(P_sector_lowest)/D; rank>=1',
                'failure':'1-rank(P_sector_lowest)/D',
                'D':str((5*2**30)**8*119**24),
                'outcomes_exhaust_complete_carrier':True,
                'native_efficient_minimum_filter_executed':False},
            'scope':{'existing_current_neutrality_positive_gap_and_closed_source_protection_joined':True,
                'periodic_variant_has_nonempty_zero_momentum_sector':True,
                'full_sector_minimum_numerically_diagonalized':False,
                'relativistic_mass_shell_or_asymptotic_particle_proved':False,
                'active_or_sterile_neutrino_identified':False,
                'physical_dark_matter_or_cosmological_abundance_identified':False,
                'neutral_sector_proved_electromagnetically_inert':False}}
