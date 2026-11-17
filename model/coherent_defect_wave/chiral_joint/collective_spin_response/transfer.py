"""Complete actual-ground spin records with controlled full-source residue."""
from pathlib import Path
import json
import sympy as s

from ..higgs_composites.certificate import validate_pins
from ..family_cell.quotient import outward
from .model import MAX_COUPLING
from .residue import exact_data
from .bounds import residue_errors, probability_interval

ROOT=Path(__file__).resolve().parents[4]
GEOMETRY_RECEIPT='model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/results.json'


def geometry_data():
    receipt=json.loads((ROOT/GEOMETRY_RECEIPT).read_text())
    validate_pins(receipt['input_sha256'])
    transfer=receipt['full_interacting_transfer'];spin=receipt['full_spin_spectrum']
    if (not transfer['unique_full_physical_interacting_ground_proved']
        or transfer['first_excited_cluster_rank'] != 3
        or not transfer['constant_ranks_follow_from_strict_isolation_for_every_positive_admitted_lambda']
        or not transfer['first_three_levels_equal_by_actual_cubic_rotation_and_translation']):
        raise ValueError('the complete rank-one ground and rank-three isolated cluster theorem is required')
    assert transfer['declared_source'] == 'H_lambda=H0+lambda*T, N=2, kappa=1, all three Yukawa coefficients=1'
    assert transfer['exact_ground_momentum_in_pi_units'] == [0,0,0]
    assert sorted(transfer['exact_first_excited_source_momenta_in_pi_units']) == [[0,0,1],[0,1,0],[1,0,0]]
    assert spin['full_spin_dimension'] == 256 and spin['ground_rank'] == 1 and spin['first_excited_rank'] == 3
    assert transfer['all_full_Fock_scalar_link_and_physical_Gauss_sectors_retained']
    return transfer,spin


def data():
    geometry,spin=geometry_data();atomic=exact_data()
    gap=s.Rational(geometry['atomic_separation_lower'])
    tau=s.Rational(geometry['full_hopping_norm_upper']);assert tau == 360
    jlo,jhi=map(s.Rational,geometry['actual_exchange_interval'])
    gl,gh=map(s.Rational,spin['gap_interval'])
    above=s.Rational(spin['complete_spectrum'][2]['energy_interval'][0])+4
    dg=2*jlo*gl;ds=2*jlo*min(gl,above)
    admitted=s.Rational(geometry['admitted_positive_hopping_interval'][1])
    if not (gap>0 and tau>0 and 0<jlo<=jhi and 0<gl<=gh and above>0 and admitted>0):
        raise ValueError('positive ordered full-source gap, exchange and coupling controls required')
    coupling=min(MAX_COUPLING,admitted)
    while True:
        try:
            errors=residue_errors(gap,tau,coupling,dg,ds)
        except ValueError:
            coupling/=2
            continue
        amplitude=errors['line_and_outside_record_amplitude_error']
        if amplitude**2 <= min(atomic['line_interval'][0],atomic['outside_interval'][0])/16:
            break
        coupling/=2
    intervals={name:probability_interval(atomic[name+'_interval'],
        errors['line_and_outside_record_amplitude_error'] if name in ('line','outside')
        else errors['total_success_and_failure_amplitude_error']) for name in ('line','outside','success')}
    intervals['failure']=(1-intervals['success'][1],1-intervals['success'][0])
    assert all(0 < lo <= hi < 1 for lo,hi in intervals.values())
    remainder=errors['full_Q_effective_H_remainder']
    energy=(s.cancel(2*coupling**2*jlo*gl-2*remainder),
            s.cancel(2*coupling**2*jhi*gh+2*remainder))
    separation=s.cancel(2*coupling**2*jlo*above-2*remainder)
    assert 0 < coupling <= admitted and energy[0]>0 and separation>0 and gap-2*coupling*tau>0
    return {'geometry':geometry,'atomic':atomic,'gap':gap,'tau':tau,'jlo':jlo,'jhi':jhi,
        'ground_gap_coefficient':dg,'cluster_separation':ds,'coupling':coupling,
        'admitted':admitted,'errors':errors,'intervals':intervals,'energy_interval':energy,
        'next_cluster_gap':separation}


def certificate():
    controls=data();geometry=controls['geometry']
    intervals={name:tuple(map(s.Rational,outward(bounds,64)))
        for name,bounds in controls['intervals'].items() if name!='failure'}
    intervals['failure']=(1-intervals['success'][1],1-intervals['success'][0])
    dimension=(5*2**30)**8*119**24
    prep=geometry['formal_full_carrier_ground_filter']
    assert int(prep['full_dimension']) == dimension
    assert s.Rational(prep['success_probability']) == s.Rational(1,dimension)
    preparation_probability=s.Rational(1,dimension)
    rows=[{'record':'initial_ground_preparation_failure','probability':str(1-preparation_probability),
        'all_source_and_apparatus_outputs_retained':True}]
    rows.append({'record':'ground_prepared_pump_failure','source_momentum_in_pi_units':[0,0,0],
        'reference':'vacuum','exact_unconditional_probability':'(1-p_line-p_outside)/D',
        'probability_interval':[str(v/dimension) for v in intervals['failure']]})
    for axis in range(3):
        phase=[int(j==axis) for j in range(3)]
        for name in ('line','outside'):
            rows.append({'record':'ground_prepared_pump_'+name,'reference_vector_axis':axis,
                'source_momentum_in_pi_units':phase,'reference_momentum_in_pi_units':phase,
                'joint_momentum_in_pi_units':[0,0,0],
                'exact_unconditional_probability':f'p_{name}/(3D)',
                'probability_interval':[str(v/(3*dimension)) for v in intervals[name]],
                'source_energy_record':('first isolated excitation energy' if name=='line'
                    else 'complete complement of the first excitation cluster')})
    return {'source':{'Hamiltonian':'H0+lambda*T on the complete N2 PW119 source',
            'all_onsite_Yukawa_scalar_and_link_coefficients':'1',
            'hopping_coefficient':str(controls['coupling']),
            'source_clock':'exp(-it H_lambda)',
            'clock_between_pump_and_readout':'all energy and momentum record probabilities hold at every real source time',
            'root_geometry_admitted_upper_coefficient':str(controls['admitted']),
            'detector_coefficient_selected_without_empirical_inputs':True,
            'reference':'unique rank-one full physical ground of this H_lambda',
            'reference_stationary_under_this_same_source_clock':True,
            'unit_hopping_or_atomic_product_reference_substituted':False},
        'complete_full_source_operator_controls':{
            'atomic_gap_lower':str(controls['gap']),'full_hopping_norm_upper':str(controls['tau']),
            'actual_exchange_interval':[str(controls['jlo']),str(controls['jhi'])],
            'effective_ground_separation_lower':str(controls['ground_gap_coefficient']),
            'effective_triplet_separation_lower':str(controls['cluster_separation']),
            'effective_and_continued_ground_ranks':[1,1],
            'effective_and_continued_triplet_ranks':[3,3],
            'same_rank_projection_norm_identity_required':True,
            'full_complement_states_discarded':False,
            'exact_rational_error_bounds':{key:str(value) for key,value in controls['errors'].items()}},
        'actual_interacting_energy_line':{
            'source_momenta_in_pi_units':[[1,0,0],[0,1,0],[0,0,1]],
            'exact_full_source_cluster_rank':3,
            'three_energies_equal_by_cubic_symmetry':True,
            'excitation_energy_interval_in_source_units':outward(controls['energy_interval'],128),
            'next_low_band_energy_separation_lower':str(controls['next_cluster_gap']),
            'first_line_probability_conditioned_on_ground_preparation':list(map(str,intervals['line'])),
            'outside_line_probability_conditioned_on_ground_preparation':list(map(str,intervals['outside'])),
            'total_pump_success_probability_conditioned_on_ground_preparation':list(map(str,intervals['success'])),
            'pump_failure_probability_conditioned_on_ground_preparation':list(map(str,intervals['failure'])),
            'each_reference_axis_line_weight_is_one_third':True,
            'complete_outside_line_weight_strictly_positive':True,
            'spectral_measurement':'joint physical Gauss and actual H_lambda spectral projectors, first isolated rank-three cluster versus its complete complement',
            'ideal_finite_spectral_service_is_declared':True,
            'efficient_native_energy_filter_schedule_executed':False,
            'full_source_eigenvectors_or_full_energy_histogram_numerically_diagonalized':False},
        'complete_preparation_and_detector_records':{
            'full_source_dimension':str(dimension),
            'initial_state':'I/D on the full kinematic carrier, with reference vacuum',
            'initial_filter':'the joint physical Gauss and unique ground-energy projector',
            'ground_preparation_success_probability':str(preparation_probability),
            'initial_failure_branch_control':'identity on retained source/reference; no pump or phase readout; terminal failure record',
            'records':rows,
            'positive_record_count':8,
            'all_other_flag_phase_and_reference_outcomes':{'probability':'0','retained':True},
            'exact_conditional_completeness':'p_line+p_outside+p_failure=1',
            'exact_unconditional_completeness':'(1-1/D)+(p_line+p_outside+p_failure)/D=1',
            'conditioned_initial_filter_failures_omitted':False,
            'all_phase_readout_failures_and_outside_cluster_states_retained':True},
        'scope':{'finite_stationary_interacting_ground_spectral_record':True,
            'native_all_family_spin_operator_and_full_H_error_transfer':True,
            'three_X_modes_identified_with_empirical_generations':False,
            'zero_momentum_rest_mass_or_relativistic_spin_inferred':False,
            'continuous_rotated_spin_symmetry_at_finite_hopping_assumed':False,
            'continuum_thermodynamic_limit_or_GeV_units':False}}
