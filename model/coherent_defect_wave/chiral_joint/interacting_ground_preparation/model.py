"""Actual local ground recipe and finite two-copy energy-agreement controls."""
import json
from pathlib import Path
import sympy as s
from ..higgs_composites.certificate import validate_pins
from ..atomic_reference_control.lattice import local_norm
from ..atomic_momentum_response.model import reference_data
from ..collective_spin_response.transfer import data as response_data
from ..family_cell.quotient import outward, sqrt_interval

ROOT = Path(__file__).resolve().parents[4]
LOCAL_RECEIPT = 'model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/results.json'


def local_recipe():
    receipt = json.loads((ROOT/LOCAL_RECEIPT).read_text())
    validate_pins(receipt['input_sha256'])
    local = receipt['complete_physical_spectrum']
    reference, ground = reference_data()
    if not (local['equitable_dimension']==403 and local['full_physical_cyclic_dimension']==49
        and local['exact_full_graph_and_cyclic_intertwiners']
        and local['positive_Perron_overlap_of_physical_seed']):
        raise ValueError('the complete physical 403-class ground construction is required')
    e = s.Symbol('E')
    p = s.Poly(s.sympify(local['full_H_minimal_polynomial'],locals={'E':e}),e)
    f = s.Poly((e*(e-1))**3-313*(e*(e-1))**2+31708*e*(e-1)-1041468,e)
    assert p.degree()==49 and p.rem(f).is_zero and s.gcd(p,p.diff()).degree()==0
    assert ground['global_physical_ground_rank']==2 and reference['raw_L']=='1'
    assert local['odd_component_ground_energy_interval']==ground['global_physical_ground_energy_interval']
    assert s.Rational(local['ground_filter_success_interval'][0])>0
    assert local['full_seed_component_dimension']==113948 and local['full_local_dimension']==5*2**30
    return {'full_local_dimension':5*2**30,'complete_graph_dimension':113948,
        'equitable_dimension':403,'full_H_cyclic_degree':49,
        'exact_ground_energy_polynomial':str(f.as_expr()),
        'selected_root':'least real root',
        'ground_energy_interval':ground['global_physical_ground_energy_interval'],
        'native_seed':'psi_up=((2048,1)-(1024,2))/sqrt(2)',
        'ground_vector':'g_up = q_E(H_local) psi_up / ||q_E(H_local) psi_up||, q_E(z)=p_H(z)/(z-E)',
        'ground_vector_denominator_strictly_positive':True,
        'opposite_spin_vector':'g_down=S_minus g_up',
        'full_local_ground_and_403_class_vector_joined_by_same_simple_root':True,
        'reference':'two spin-half gauge-neutral fermion modes, each raw L=-1; complete four-state Fock carrier',
        'retained_local_singlet':'(g_up tensor ref_down - g_down tensor ref_up)/sqrt(2)',
        'local_preparation_unitary':'swap full vacuum with this singlet; identity on the complete orthogonal complement',
        'preserved_local_quantities':['Gauss','total raw L','total fermion parity','total spin'],
        'preparation_success_probability':'1',
        'algebraic_state_synthesis_service_declared':True,
        'native_efficient_gate_synthesis_executed':False}


def filter_bounds(trial_probability, energy_norm, gap, leakage=s.Rational(1,1024)):
    p,b,g,e = map(s.sympify,(trial_probability,energy_norm,gap,leakage))
    if any(v.is_Rational is not True for v in (p,b,g,e)) or not (s.Rational(1,2)<p<=1 and b>0 and 0<g<=2*b and 0<e<1):
        raise ValueError('rational p>1/2, full energy norm, positive gap and leakage in (0,1) required')
    length = int(s.ceiling(2*b/(e*g)))
    actual_leak = s.cancel(2*b/(length*g))
    denominator = p**2+(1-p)**2+2*actual_leak**2*p*(1-p)
    fidelity = s.cancel(p**2/denominator)
    trace_distance = sqrt_interval((1-fidelity,1-fidelity),100)[1]
    return {'length':length,'leakage':actual_leak,'success_lower':p**2,
        'ground_pair_fidelity_lower':fidelity,'single_copy_trace_distance_upper':trace_distance}


def certificate():
    local = local_recipe(); response = response_data()
    eta = response['errors']['full_ground_vector_error_after_phase_alignment']
    assert 0<eta<1
    overlap = (1-eta**2/2)**2
    # This bound applies on the full source, without centering by an unknown energy.
    norm = 8*local_norm()+24+response['coupling']*response['tau']
    control = filter_bounds(overlap,norm,response['energy_interval'][0])
    p0 = control['success_lower']; delta = control['single_copy_trace_distance_upper']
    intervals = {key:(max(0,value[0]-delta),min(1,value[1]+delta))
        for key,value in response['intervals'].items()}
    assert all(0<lo<=hi<1 for lo,hi in intervals.values())
    return {'source':{'Hamiltonian':'the same full N2 PW119 H0+lambda*T as collective_spin_response',
            'hopping_coefficient':str(response['coupling']),
            'onsite_Yukawa_scalar_and_link_coefficients':'1','clock':'exp(-it H_lambda)',
            'physical_source_copies':2,'source_copy_dimension':str((5*2**30)**8*119**24),
            'every_source_reference_and_filter_register_retained':True},
        'actual_local_preparation':local,
        'atomic_cube_trial':{
            'source_trial':'u = U times the unique normalized cube singlet at the least root of x^3+20x^2-16x-576',
            'reference_trial':'J u; J is the tensor product of spin-half time reversal on the retained one-reference-particle-per-cell space',
            'initial_joint_state':'Psi = product of eight local source/reference singlets',
            'target_joint_state':'u tensor J u, links blank',
            'exact_initial_target_overlap':'1/16',
            'Householder_vector':'w=(Psi-u tensor J u)/sqrt(15/8)',
            'Householder_unitary':'I-2|w><w| on the complete carrier',
            'preserved_quantities':['Gauss','total raw L','total parity','joint cubic rotations','joint translations'],
            'trial_source_momentum_in_pi_units':[0,0,0],
            'exact_interacting_ground_trial_overlap_lower':str(overlap),
            'full_H_ground_vector_distance_upper':str(eta),
            'full_H_ground_eigenvector_used_as_preparation_input':False},
        'finite_energy_agreement_filter':{
            'full_absolute_H_norm_upper':str(norm),'spectral_gap_lower':str(response['energy_interval'][0]),
            'rational_step':str(1/norm),
            'unitary':'C=(I-i H_lambda/B)(I+i H_lambda/B)^(-1); W=C tensor C^dagger',
            'filter':'F_M(W)=(I+W+...+W^(M-1))/M',
            'explicit_integer_M':str(control['length']),
            'ground_excited_chord_gap_lower':str(response['energy_interval'][0]/norm),
            'mixed_ground_excited_amplitude_upper':str(control['leakage']),
            'success_probability_lower':str(p0),
            'failure_probability_interval':['0',str(1-p0)],
            'conditional_ground_pair_fidelity_lower':str(control['ground_pair_fidelity_lower']),
            'conditional_one_source_trace_distance_upper':str(delta),
            'complete_ancilla_measurement':'uniform M-state control, controlled W^k, inverse uniform preparation; retain every output k=0,...,M-1',
            'success_flag':'ancilla zero','failure_flag':'all other ancilla outcomes, retained separately or coarse-grained',
            'H_energy_shift_or_ground_energy_required':False,
            'finite_algebraic_Cayley_and_state_synthesis_services_declared':True,
            'native_efficient_gate_schedule_executed':False},
        'same_source_readout_after_approximate_preparation':{
            'conditional_probability_intervals':{key:outward(value,64) for key,value in intervals.items()},
            'conditional_to_exact_ground_readout_total_variation_upper':str(delta),
            'bound_holds_for_every_real_clock_time':True,
            'preparation_failure_record_retained':True,
            'second_source_and_spin_references_retained':True,
            'exact_completeness':'q_fail + q_success*(p_line+p_outside+p_pump_fail)=1',
            'prepared_state_asserted_exactly_stationary':False},
        'scope':{'constructive_finite_algebraic_service_schedule':True,
            'coherent_full_source_evolution_numerically_executed':False,
            'physical_native_efficiency_or_mass_identification':False}}
