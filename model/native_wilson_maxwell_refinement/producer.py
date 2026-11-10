"""Replay the fixed-source photon-cutoff theorem and exact local algebra."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
RER=HERE.parents[1]
sys.path.insert(0,str(RER))
from model.native_wilson_maxwell_refinement import bounds,algebra
from model.particle_common_contract.provenance import import_closure
PARENT='model/native_wilson_maxwell/results.json'
PARENT_SHA='eba4963148b62f6841a038119ed5b9baf4a6f142009915b31e49d72dacd12b72'
RESEARCH=HERE/'research_reference.json'
RESEARCH_SHA='e8ae519a5d5247de12cfaa90d8b63f870f0bd72110fe88c39e892fbc429394d9'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def source():
    path=RER/PARENT
    if digest(path)!=PARENT_SHA:raise ValueError('frozen finite source receipt differs')
    saved=json.loads(path.read_text())
    for relative,expected in saved['input_sha256'].items():
        if digest(RER/relative)!=expected:raise ValueError('source input differs: '+relative)
    for relative,expected in saved['owned_source_sha256'].items():
        if digest(path.parent/relative)!=expected:raise ValueError('parent source differs: '+relative)
    return saved


def record_transfer(saved,cut):
    p=F(saved['actual_preparation']['global_code_probability_lower'])
    finite=saved['complete_interaction_records']['exact']
    a=F(finite['ideal_complete_target_amplitude_lower'])
    b=F(finite['dimensionless_mT'])*F(finite['whole_perturbation_norm_over_m_upper'])
    delta=F(cut['state_vector_error_upper'])
    assert 0<delta<a
    target=p*(a-delta)**2
    connected=p*((a-delta)**2-(b+delta)**4)
    assert connected>0
    values={'same_native_trial_code_probability_lower':p,
        'comparison_ideal_target_amplitude_lower':a-delta,
        'comparison_actual_target_probability_lower':target,
        'comparison_actual_target_complement_probability_upper':1-target,
        'comparison_ideal_each_excitation_marginal_upper':(b+delta)**2,
        'comparison_actual_binary_connected_lower':connected,
        'finite_actual_target_probability_lower':F(finite['actual_unconditional_target_probability_lower']),
        'comparison_photon_overflow_probability_upper':delta**2}
    exact={key:str(bounds.rounded(value,'upper' in key)) for key,value in values.items()}
    return {'exact':exact,'diagnostic_binary64':{k:float(F(v)) for k,v in exact.items()},
        'actual_native_trial_and_full_references_retained':True,
        'same_p_code_exact_mixture_in_both_sources':True,
        'conditioned_on_rare_record':False,'target':[1,1,1],
        'finite_executed_pointer_labels':1161,
        'common_overflow_padded_labels':1451,
        'padding_changes_parent_executed_apparatus':False,
        'nonstationary_free_vacuum_dressing':True,
        'coupled_ground_radiation_or_asymptotic_scattering':False}


def calculate():
    saved=source();cut=bounds.cutoff_bound(3)
    assert F(cut['state_vector_error_upper'])<F(1863,10**12)
    algebra_data=algebra.calculate()
    # A controlled electric test and a controlled magnetic test, both
    # strength1/64. These are additional complete comparison settings.
    # Exact algebra gives electric variance <1/2 and magnetic variance <2.
    s=F(1,64)*(bounds.root(F(1,2))+bounds.root(2))
    history=bounds.cutoff_bound(3,apparatus_radius=s,histories=True,order=24)
    tables=[bounds.cutoff_bound(k,order=max(24,k+20)) for k in (1,2,3,4,5,6,8,10,16,32)]
    return {'schema':'oph.fixed_native_wilson_photon_cutoff.v1',
        'parent_finite_source_sha256':PARENT_SHA,
        'construction':{'changed_parameter':'total photon cap K only',
            'native_constituents':32,'complete_fermion_code_dimension':16**8,
            'fixed_five_photon_frequencies_over_m':['2','2','2','sqrt(6)','sqrt(6)'],
            'finite_Hph':'m [dGamma(nu)+Tr(nu)/2], with positive boundary completion',
            'literal_finite_link':'exp(i g P_K Q_e P_K)',
            'compression_of_unbounded_exponential_assumed':False,
            'same_native_trial_state_reference_density_Gauss_and_clock':True,
            'unbounded_oscillator_is_mathematical_comparison_only':True,
            'source_mass_or_empirical_coupling_changed':False},
        'cap3_complete_state_and_record_bound':cut,
        'cap3_pair_photon_transfer':record_transfer(saved,cut),
        'cutoff_table':tables,
        'two_bounded_field_test_history':{'electric_strength':'1/64','magnetic_strength':'1/64',
            'all_control_outcomes_retained':True,'comparison':history,
            'finite_Hadamard_pointer_labels':4,
            'additional_common_free_evolution_after_readout':'must be included in total autonomous time if terminal source is reused',
            'is_original_parent_count_experiment':False},
        'nonincident_observer_order_comparison':{
            'electric_edge_index':0,'magnetic_face_index':3,
            'same_source_wait_before_tests':'1/(64m)',
            'each_Weyl_test_strength':'1/64',
            'complete_joint_record_TV_between_two_orders_upper':str(bounds.rounded(2*F(history['complete_record_TV_upper']))),
            'all_Hadamard_records_and_source_disturbance_retained':True,
            'no_autonomous_wait_between_ordered_test_factors':True,
            'separated_time_propagation_cone_claimed':False},
        'exact_local_comparison_algebra':algebra_data,
        'capacity':{'photon_dimensions':{str(k):comb(k+5,5) for k in (1,2,3,4,8,16,32)},
            'finite_count_pointer_dimension':'290(K+1)+1',
            'comparison_with_overflow_pointer_dimension':'290(K+2)+1',
            'full_capacity_formula_nats':'32*(256ln2+8ln5+24ln119+ln(2^31+1))+7ln9+ln65+ln binom(K+5,5)+ln[290(K+1)+1]',
            'free_photon_raw_norm_upper':'m*(3+sqrt(6)+sqrt(6)*K)',
            'fixed_photon_free_full_source_bound':'parent H_M, independent of K',
            'bounded_oriented_Peierls_sum_norm_upper':'48m',
            'uniform_infinite_controller_cost_claimed':False},
        'boundaries':{'arbitrary_uncontrolled_photon_injecting_instruments':False,
            'uniform_spatial_propagation_cone':False,'growing_volume_or_UV_limit':False,
            'full_native_complement_local_field_identification':False,
            'interacting_Lorentz_QED_or_SM_identification':False}}


def certificate():
    result=calculate();saved=source()
    inputs=dict(saved['input_sha256'])
    inputs.update({str(Path(PARENT).parent/k):v for k,v in saved['owned_source_sha256'].items()})
    inputs[PARENT]=PARENT_SHA
    # Directly used exact cube providers are also present in parent pins.
    for relative in ('model/native_wilson_maxwell/field.py','model/native_wilson_maxwell/gauss.py'):
        assert relative in inputs
    seeds=set(inputs)|{str(RESEARCH.relative_to(RER)),'model/particle_common_contract/provenance.py'}
    seeds.update(str(path.relative_to(RER)) for path in HERE.iterdir() if path.suffix in ('.py','.md'))
    result['input_sha256']={name:digest(RER/name) for name in sorted(import_closure(seeds))}
    result['owned_source_sha256']={p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('immutable research comparison differs')
    ignored={'input_sha256','owned_source_sha256'}
    science=lambda value:{k:v for k,v in value.items() if k not in ignored}
    if science(result)!=science(json.loads(RESEARCH.read_text())):
        raise ValueError('public cutoff science differs from frozen research')
    result['independent_research_comparison']={'reference_sha256':RESEARCH_SHA,
        'all_scientific_fields_equal':True,'ignored_only_provenance_fields':sorted(ignored)}
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args();value=certificate();path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    if args.check:
        if json.loads(path.read_text())!=value:raise SystemExit('receipt differs')
        print('PASS: fixed-source photon cutoff and complete-record receipt')
    if not(args.write or args.check):print(json.dumps(value,indent=2,sort_keys=True))


if __name__=='__main__':main()
