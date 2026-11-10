"""Actual native charge centroid and a complete finite phase estimator."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import factorial
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
from model.native_charge_resolution import bounds as resolution
from model.coherent_defect_wave.chiral_joint.atomic_reference_control.lattice import hopping_certificate
from model.affine_gauss_charged_clock.model import rational_clock_identities
M=2**20
SUBSTEPS=512

def mean_grid_error(size):
    """Uniform linear phase error for true phase in[1/4,3/4]."""
    if not isinstance(size,int) or size<4 or size&(size-1):raise ValueError('power-of-two pointer with at least4 labels required')
    # ln(size/2)=(log2(size)-1)ln2; ln2<7/10 follows from e^(7/10)>2.
    bits=size.bit_length()-1
    assert sum((F(7,10)**j/F(factorial(j)) for j in range(5)),F(0))>2
    return (F(7,2)+F(7*(bits-1),20))/size+F(7,size*size)

@lru_cache(None)
def calculate():
    d=resolution.data();row=resolution.window(128)
    controls=hopping_certificate(2)
    count=controls['all_native_coordinates_blank_to_blank_matrix_elements_checked']
    assert count==59 and controls['every_actual_path_nonempty_and_distinct_edges']
    assert F(controls['full_hopping_norm_upper_bound'])==360
    v=d['v'];r=128*v;alo,ahi=d['energies'][0];center=(alo+ahi)/2
    distance=F(row['selected_atom_distance_to_bin_complement_lower'])
    maximum_offset=r+(ahi-alo)/2
    weight_upper=d['weights'][0][1]
    rho=d['eta']+F(row['complete_other_atom_leakage_norm_upper'])
    weight_lower=F(row['unnormalized_bin_weight_interval'][0])
    selected_moment=v*v*weight_upper/distance
    contamination=maximum_offset*(2*resolution.root(weight_upper,True)*rho+rho*rho)
    error=resolution.round_up((selected_moment+contamination)/weight_lower)
    elo,ehi=d['e'];neutral=d['neutral_shift']
    centroid=(alo-4-8*ehi-error,ahi-4-8*elo+neutral+error)
    prep=json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json').read_text())['preparation']
    assert F(prep['source']['hopping_coefficient'])==F(1,2**34)
    physical=prep['finite_energy_agreement_filter']
    trace=F(physical['conditional_one_source_trace_distance_upper'])
    qlower=F(physical['success_probability_lower'])
    preparation_error=2*r*trace/weight_lower
    prepared=(centroid[0]-preparation_error,centroid[1]+preparation_error)
    grid=mean_grid_error(M)
    gate=(M-1)*F(33,7)**7/(100800*SUBSTEPS**6)
    estimator_error=4*r*(grid+gate)
    finite=(prepared[0]-estimator_error,prepared[1]+estimator_error)
    reported_ground=4*(elo+ehi)-neutral/2
    ground_phase_uncertainty=4*(ehi-elo)+neutral/2
    raw_mean=(alo-error-preparation_error-estimator_error,
              ahi+error+preparation_error+estimator_error)
    decoded_mean=(raw_mean[0]-4-reported_ground,raw_mean[1]-4-reported_ground)
    assert decoded_mean[0]-ground_phase_uncertainty==finite[0]
    assert decoded_mean[1]+ground_phase_uncertainty==finite[1]
    assert gate<F('3e-11')
    assert centroid[1]-centroid[0]<F('3.2e-9')
    assert finite[1]-finite[0]<F('3.4e-9')
    return {'schema':'oph.actual_native_charge_centroid.v1',
      'source_h_and_native_spectral_selection_unchanged':True,
      'absolute_total_source_plus_reference_bin':row['absolute_source_plus_reference_bin'],
      'unnormalized_bin_weight_interval':row['unnormalized_bin_weight_interval'],
      'first_moment_identity':{'actual_blank_link_coordinates_replayed':count,
        'every_actual_hopping_path_nonempty_with_distinct_edges':True,
        'P_allblank_T_P_allblank':'0 on every matter/scalar input',
        'selected_atomic_first_centered_moment':'<psi0,(h+Href-a0)psi0>=0'},
      'complete_atomic_tail_first_moment_absolute_upper':str(selected_moment),
      'actual_ground_and_other_atom_contamination_norm_upper':str(rho),
      'centroid_error_from_exact_selected_atomic_energy':str(error),
      'actual_bin_conditional_source_excitation_centroid_interval':list(map(str,centroid)),
      'diagnostic_actual_centroid_interval':list(map(float,centroid)),
      'actual_source_ground_shift_interval':list(map(str,(-neutral,F(0)))),
      'actual_physical_preparation':{'pair_filter_success_lower':str(qlower),
        'conditional_one_source_trace_distance_to_actual_ground_upper':str(trace),
        'extra_P0_ground_projection_used':False,
        'actual_pumped_bin_weight_interval':[str(weight_lower-trace),str(F(row['unnormalized_bin_weight_interval'][1])+trace)],
        'original_attempt_bin_success_lower':str(qlower*(weight_lower-trace)),
        'conditional_CAR_pump_failure_probability_interval':[str(F(1,2)-trace),str(F(1,2)+trace)],
        'conditional_complete_outside_bin_probability_interval':[
          str(F(1,2)-F(row['unnormalized_bin_weight_interval'][1])-trace),
          str(F(1,2)-weight_lower+trace)],
        'post_bin_trace_distance_to_exact_ground_comparison_upper':str(trace/weight_lower),
        'conditional_centroid_preparation_error_upper':str(preparation_error),
        'source_plus_partner_and_all_filter_failures_retained':True},
      'actual_prepared_bin_conditional_centroid_interval':list(map(str,prepared)),
      'clock':{'pointer_dimension':M,'all_pointer_labels_retained':True,
        'active_generator':'A=h_active+4 N_ref; partner h and preparation references continue autonomously',
        'pointer_scalar_phase':str(center-2*r),'base_time':'pi/(2r), r=128*360*2^-34',
        'true_phase_interval':['1/4','3/4'],'decoded_total_energy':'B+4r k/M, k=0,...,M-1',
        'decoded_source_energy':'B+4r k/M-4-Eg_report, with explicit rational Eg_report',
        'reported_neutral_energy_baseline':str(reported_ground),
        'actual_Eg_distance_from_reported_baseline_upper':str(ground_phase_uncertainty),
        'rational_base_gate':'R33(t0(A-B)/n)^n on admitted one-reference bin; identity extra gate elsewhere',
        'R33_substeps':SUBSTEPS,'R33_factor_count_upper':SUBSTEPS*(M-1),
        'R33_controlled_error_upper':str(gate),
        'uniform_mean_absolute_grid_phase_error_upper':str(grid),
        'energy_estimator_bias_upper_including_rational_error':str(estimator_error),
        'reference_cost_retained':'4','exact_native_Eg_used_as_free_pointer_phase':False,
        'common_autonomous_Htotal_on_all_branches':True,
        'large_source_pointer_matrix_numerically_executed':False,
        'prior_bin_admission_remains_explicit_spectral_service':True},
      'finite_clock_conditional_mean_source_excitation_interval':list(map(str,finite)),
      'finite_pointer_mean_absolute_source_plus_reference_energy_interval':list(map(str,raw_mean)),
      'finite_pointer_mean_decoded_source_energy_interval':list(map(str,decoded_mean)),
      'diagnostic_finite_clock_mean_interval':list(map(float,finite)),
      'complete_law':['1-q','q tr(K_fail rho K_fail^dagger)',
        'q tr((I-P_I) B rho B^dagger (I-P_I))',
        'q tr(K_k P_I B rho B^dagger P_I K_k^dagger) for every k'],
      'post_bin_density_in_law_unnormalized':True,
      'source_partner_and_all_preparation_reference_outputs_retained':True,
      'source_dimension':'D^2, D=2^240 5^8 119^24',
      'retained_reference_pointer_dimension':'2^30 2^32 Mprep 2^3 2^20',
      'Mprep':1197721483502453176152084669,
      'gate_synthesis_and_controller_workspace_counted':False,
      'single_charged_eigenvalue_or_original_Estar_relocalized':False,
      'relativistic_mass_spin_or_SM_identification':False,
      'rational_clock_identities':rational_clock_identities()}
