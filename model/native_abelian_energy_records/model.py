"""Same-H whole-bin energy and ordered arrival/finite-clock bounds."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import json
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.affine_gauss_charged_clock.model import tail_bound,acceptance_indices,rational_clock_identities
M=2**20
SUBSTEPS=2**10
SEPARATE_RAW_SUBSTEPS=2**19
HALF_WINDOW=F(3,100)

def sqrt_bound(value,upper=False,bits=96):
    value=F(value);assert value>=0
    scale=1<<bits
    n=isqrt(value.numerator*scale*scale//value.denominator)
    if upper and n*n*value.denominator!=value.numerator*scale*scale:n+=1
    return F(n,scale)

def load(relative):return json.loads((ROOT/relative).read_text())

@lru_cache(None)
def inputs():
    return (load('model/native_abelian_field_join/results.json'),
      load('model/affine_gauss_full_spectrum/results.json'),
      load('model/affine_gauss_spectrum/results.json'),
      load('model/affine_gauss_mobile/results.json'),
      load('model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json'))

@lru_cache(None)
def ground_bounds():
    join,full,coarse,mobile,prep=inputs()
    L,U=map(F,join['intervals']['native_selected_E_interval'])
    center=(L+U)/2;epsilon=(U-L)/2
    aff_lo,aff_hi=map(F,full['actual_charged_ground_energy_interval'])
    f0=F(full['complete_ground_preparation']['conditional_ground_window_probability_interval'][0])
    gap0=2-aff_hi;gap=gap0-2*epsilon
    angle=epsilon/(gap0-epsilon)
    # Equal finite ranks turn the one-sided ground-space residual estimate
    # into the full projector norm estimate. Finite G1 covariance prevents
    # splitting of the two source states inside the isolated cluster.
    assert gap>F('.124') and 0<angle<F('.008')
    infidelity=(sqrt_bound(1-f0,True)+angle)**2
    fidelity=1-infidelity
    assert fidelity>F('.9997')
    return {'native_center':str(center),'native_half_width':str(epsilon),
      'comparison_to_immobile_bin_gap_lower':str(gap0),
      'actual_admitted_bin_ground_gap_lower':str(gap),
      'actual_vs_comparison_ground_projector_norm_upper':str(angle),
      'actual_bin_ground_energy_relative_full_neutral_interval':list(map(str,(L+2+aff_lo,U+2+aff_hi))),
      'source_G1_ground_rank':2,'static_reference_spin_multiplicity_retained':2,
      'coherent_native_CAR_spatial_pump_conditional_bin_ground_fidelity_lower':str(fidelity),
      'reference_energy_retained':'4','outside_J_full_charged_ground_identified':False}

@lru_cache(None)
def clock_bounds():
    join,full,coarse,mobile,prep=inputs();ground=ground_bounds()
    L,U=map(F,join['intervals']['native_selected_E_interval'])
    center=(L+U)/2;epsilon=(U-L)/2
    e0lo,e0hi=map(F,prep['preparation']['actual_local_preparation']['ground_energy_interval'])
    hopping=F(45,2147483648)
    Eg_center=4*(e0lo+e0hi)
    Eg_radius=4*(e0hi-e0lo)+hopping
    assert -90<Eg_center-Eg_radius<Eg_center+Eg_radius<0
    baseline_radius=26*Eg_radius
    aff_center=F(full['energy_center_exact_binary_value'])
    aff_radius=F(full['arithmetic_admission']['maximum_admitted_actual_H_residual'])
    target_center=center+2+aff_center
    radius=epsilon+aff_radius+baseline_radius
    next_energy=L+4-baseline_radius
    # EC<=6/31, photon energy<9, hopping norm<3/10 on this sector.
    assert U+4+F(6,31)+9+F(3,10)+baseline_radius<15
    assert L+4-F(3,10)-baseline_radius>0
    low,high=acceptance_indices(target_center,HALF_WINDOW,M)
    dg=(HALF_WINDOW-radius)/16
    de=min((next_energy-target_center-HALF_WINDOW)/16,
           1-(15-target_center+HALF_WINDOW)/16)
    assert dg>0 and de>0
    alpha,beta=tail_bound(M,dg),tail_bound(M,de)
    # Raw free H includes all26 native copies and the field vacuum.
    # |Eg|<90, Evac<45, 0<K<15 give |Hraw|<2400; the explicit
    # scalar phase B=26 Eg_center+Evac has |B|<2400 too.
    tB=F(6600,7) # (pi/8)*2400 < (22/56)*2400
    separate_gate=(M-1)*2*tB**7/(100800*SEPARATE_RAW_SUBSTEPS**6)
    # The explicit scalar B is part of the controlled extra interaction.
    # Common autonomous Hraw is retained; only its controlled function is
    # approximated through the bounded admitted operator Hraw-B.
    gate=(M-1)*F(165,28)**7/(100800*SUBSTEPS**6)
    assert gate<F('3e-12') and separate_gate<F('1e-11')
    a=alpha+gate;b=beta+gate
    fidelity=F(ground['coherent_native_CAR_spatial_pump_conditional_bin_ground_fidelity_lower'])
    good=fidelity*(1-a);bad=(1-fidelity)*b
    source=join['whole_bin_coupling']
    Q=F(source['native_physical_preparation_filter_success_lower'])
    p=F(source['additional_P0_success_lower'])
    w=F(source['bin_weight_unnormalized_interval'][0])
    return {'pointer_dimension':M,'source_step':'pi/8','accepted_pointer_first':low,
      'accepted_pointer_last':high,'accepted_pointer_count':high-low+1,
      'target_center_above_explicit_pointer_phase':str(target_center),
      'ground_radius_including_native_raw_phase_uncertainty':str(radius),
      'native_ground_center_for_pointer_phase':str(Eg_center),
      'native_ground_radius_from_full_hopping_bound':str(Eg_radius),
      'pointer_phase':'B=26 Eg_center+E_vac; exp(+i j pi B/8)',
      'actual_autonomous_clock':'exp(-it Hraw), including all26 native copies and every reference',
      'admitted_shifted_spectrum_interval':['0','15'],
      'actual_excited_energy_lower_above_pointer_phase':str(next_energy),
      'ground_rejection_upper':str(a),'excited_acceptance_upper':str(b),
      'R33_implemented_base_gate':'R33(t0 (Hraw-B)/n)^n on admitted sector; identity extra gate elsewhere',
      'R33_substeps_per_shifted_gate':SUBSTEPS,
      'R33_controlled_unitary_error_upper':str(gate),
      'rational_factor_count_upper':SUBSTEPS*(M-1),
      'separate_raw_and_scalar_alternative':{'substeps_per_factor':SEPARATE_RAW_SUBSTEPS,
        'operator_error_upper':str(separate_gate),'factor_count_upper':2*SEPARATE_RAW_SUBSTEPS*(M-1)},
      'coherent_spatial_pump_clock_true_acceptance_conditional_bin_lower':str(good),
      'coherent_spatial_pump_clock_unconditional_true_acceptance_lower':str(Q*p*w*good),
      'coherent_spatial_pump_clock_accepted_infidelity_upper':str(bad/(good+bad)),
      'long_clock_transferred_by_epsilon_times_duration':False,
      'clock_applied_directly_to_actual_H':True,
      'ground_energy_is_only_in_admitted_native_bin':True}

@lru_cache(None)
def ordered_arrival_bounds():
    join,full,coarse,mobile,prep=inputs();ground=ground_bounds();clock=clock_bounds()
    epsilon=F(ground['native_half_width']);angle=F(ground['actual_vs_comparison_ground_projector_norm_upper'])
    source=join['whole_bin_coupling'];Q=F(source['native_physical_preparation_filter_success_lower']);p=F(source['additional_P0_success_lower'])
    wlo,whi=map(F,source['bin_weight_unnormalized_interval'])
    vector=np.load(ROOT/'model/affine_gauss_full_spectrum/trial_state.npy',allow_pickle=False)
    norm2=sum((F.from_float(float(v))**2 for v in vector),F(0))
    affine_angle=F(full['normalized_trial_sine_angle_to_actual_ground_upper'])
    ph_hi=F(full['actual_ground_nonvacuum_transverse_probability_interval'][1])
    rejection=F(clock['ground_rejection_upper']);false=F(clock['excited_acceptance_upper'])
    rows=[]
    for row in mobile['finite_Peierls_first_integral_controls']:
        x,y,edge=int(row['arrival_vertex']),int(row['source_vertex']),int(row['actual_oriented_edge'])
        time=F(row['time']);eta=F(1,24)
        pvac=F.from_float(float(vector[x*5456]))**2/norm2-affine_angle
        assert pvac>0
        px_hi=F(full['actual_ground_position_records'][x]['actual_ground_position_probability'][1])
        overlap=F(coarse['spectral']['operator_bounds']['radial_and_surface_vacuum_overlaps'][int(edge>=12)][0])
        # These controls start at the reference vertex y=1, EC_y=0.
        assert y==1 and x in (0,2)
        delta=F(26 if x==0 else 33,217)
        vacuum_first=eta*time*overlap*(1-delta**2*time**2/24)
        photon_first=eta*time*sqrt_bound(1-overlap**2,True)
        remainder=F(row['all_orders_remainder_interval'][1])
        affine_ground_arrival=(sqrt_bound(pvac)*vacuum_first
          -sqrt_bound(ph_hi,True)*photon_first-sqrt_bound(px_hi,True)*remainder)
        arrival_upper=eta*time+remainder
        native_ground_arrival=affine_ground_arrival-angle*arrival_upper-epsilon*time
        assert native_ground_arrival>0
        ground_weight=native_ground_arrival**2
        actual_arrival_upper=arrival_upper+epsilon*time
        excited_weight_upper=actual_arrival_upper**2-ground_weight
        assert excited_weight_upper>0
        good=ground_weight*(1-rejection)
        bad=excited_weight_upper*false
        lower=Q*p*wlo*good
        rows.append({'source_vertex':y,'arrival_vertex':x,'arrival_time':str(time),
          'affine_ground_component_after_arrival_amplitude_lower':str(affine_ground_arrival),
          'native_ground_component_after_arrival_amplitude_lower':str(native_ground_arrival),
          'comparison_arrival_amplitude_upper':str(arrival_upper),
          'actual_arrival_amplitude_upper':str(actual_arrival_upper),
          'unconditional_arrival_then_true_ground_clock_acceptance_lower':str(lower),
          'unconditional_arrival_then_true_ground_clock_acceptance_lower_decimal_diagnostic':format(float(lower),'.17g'),
          'clock_accepted_state_infidelity_upper_after_arrival':str(bad/(good+bad)),
          'clock_accepted_state_infidelity_upper_after_arrival_decimal_diagnostic':format(float(bad/(good+bad)),'.17g'),
          'ordered_measurements':['position at stated short source time','finite actual-H energy clock after that position record'],
          'same_time_joint_position_energy_projector_claimed':False})
    return rows


def calculate():
    return {'schema':'oph.native_abelian_energy_response.v1',
      'whole_bin_source_unchanged':True,'source_Hamiltonian_choice':'whole-bin, no complement penalty',
      'ground_subspace':ground_bounds(),'finite_clock':clock_bounds(),
      'ordered_arrival_energy_records':ordered_arrival_bounds(),
      'rational_clock_operator_identities':rational_clock_identities(),
      'complete_grouped_records':['native preparation failure','P0 failure','native CAR pump failure',
        'outside native bin','bin nonarrival','arrival with each of M finite energy-pointer outcomes'],
      'record_law':['1-Q','Q(1-p)','Qp/2','Qp(1/2-w_J)','Qp w_J(1-alpha_arrival)',
        'Qp w_J ||K_k P_x U(t) Phi_J||^2 for every k'],
      'no_op_branches_keep_common_autonomous_Hraw_evolution':True,
      'physical_spin_Dirac_or_mass_identification':False}
