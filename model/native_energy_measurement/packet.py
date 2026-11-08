"""Certify every emitted native word of the fair transfer-energy instrument."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from model.native_energy_identification import packet as energy_packet
from model.native_transfer_control import packet as control_packet
from model.native_effect_synthesis import packet as synthesis_packet
from model.native_noncentral_reference import packet as reference_packet
from model.native_phase_experiment import provenance as phase_provenance
from . import measurement as m, spectral, reference

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_energy_measurement/'+name for name in (
    '__init__.py','spectral.py','measurement.py','reference.py','packet.py','proposal.json',
    'test_measurement.py','PROOF.md','README.md'))+('model/native_phase_experiment/results.json',)


def canonical(value):return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=set(energy_packet.source_bindings())|set(control_packet.source_bindings())|set(reference_packet.source_bindings())|set(phase_provenance.source_hashes())|set(SOURCE_FILES)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def produce():
    before=source_bindings();certificate=spectral.certify()
    m.syn.certify_primitive.cache_clear()
    initial=np.zeros((66,2),complex);initial[0,0]=initial[1,1]=1/np.sqrt(2)
    expected,residual=m.ideal_probabilities(initial)
    current=initial.copy();words=[];observed=[]
    for index in range(m.STAGES):
        word=m.word(index)
        checked=m.syn.certify_word(word)
        admitted=m.syn.interpreter.source_word(word,9)
        stop=m.syn.interpreter.execute_word(admitted,np.zeros((1,1)))
        if stop['status']!='RESOURCE_STOP' or stop['pulses_executed']!=0:
            raise ArithmeticError('the generic coherent evaluator limit changed')
        current,abort=m.menus.successful_word(word,current)
        observed.append(abort)
        words.append({'stage':index,'word_sha256':synthesis_packet.word_digest(word),
                      'distinct_primitive_checks':checked,'resources':admitted['declared_resources'],
                      'observed_padding_energy_flags':1,'internal_padding_rows':446,
                      'generic_host_status':'RESOURCE_STOP','generic_host_pulses_executed':0})
    if sum(np.linalg.norm(a-b) for a,b in zip(observed,expected,strict=True))+np.linalg.norm(current-residual)>1e-7:
        raise ArithmeticError('the direct stored-product energy diagnostic failed')
    # Every possible conditional basis reset uses these same eighteen pulses.
    for bit in range(9):m.syn.certify_word(m.syn.tomo.route_word((('X',bit),)))
    routes=[]
    for terminal in range(66):
        word=m.reset_word(terminal)
        operations=tuple(('X',bit) for bit in range(9)
                         if (m.syn.tomo.codes()[terminal]^m.syn.tomo.codes()[0])>>bit&1)
        if m.syn.tomo.permutation(operations)[m.syn.tomo.codes()[terminal]]!=m.syn.tomo.codes()[0]:
            raise ArithmeticError('a retained-record conditional reset route failed')
        routes.append(len(word))
    for initial_label in range(66):
        run=m.bounds.detector.run(initial_label,tuple(j%2 for j in range(11)))
        weight,mapping,_=m.preparation_branch(run['slots'])
        if weight!=F(1,2**11) or mapping!={initial_label:run['final']}:
            raise ArithmeticError('the original fair preparation instrument lost its source map')
    value={'schema':'oph.native_energy_measurement.v1','source_bindings':before,
           'spectral_certificate':certificate,
           'complete_instrument_budget':m.budget(),
           'optional_preparation':m.preparation_budget(),
           'counted_energy_label_intervals':[[str(x) for x in pair] for pair in spectral.energy_labels()],
           'all_native_words':words,'all_conditional_reset_pulse_counts':routes,
           'finite_resolution':m.resolution(F(1,1000)),
           'native_prepared_complete_energy_law':m.prepared_reference_histogram(),
           'reused_count_derived_reference_and_modular_clock':reference.certificate(),
           'coarse_histogram_sampling':m.sampling_budget(preparation_rejection=True),
           'direct_stored_product_diagnostic':{'reference_dimension':2,
                'complete_padding_reference_blocks_checked':66,
                'all_success_prefixes_executed':True,'aggregate_diagnostic_gate':'1/10000000',
                'no_normalization_or_postselection':True,
                'floating_propagation_is_exact_arithmetic_certificate':False},
           'scope':{'complete_compiled_finite_energy_instrument':True,
                'every_stored_success_prefix_evaluated_on_reference_control':True,
                'full_expanded_generic_instrument_executed':False,
                'all_emitted_stored_primitives_certified':True,
                'arbitrary_source_reference_in_instrument_bound':True,
                'fair_repair_and_coherent_service_source_charts_equal':True,
                'all66_spectral_labels_and_final_residual_retained':True,
                'energy_detection_preserves_the_source_after_terminal_flag':False,
                'rounded_instrument_is_exactly_trace_preserving':False,
                'single_monolithic_program_admitted':False,
                'native_repair_automatically_issues_coherent_words':False,
                'growing_all_number_coherent_embedding_derived':False,
                'physical_pulse_attempt_clock_identified':False,
                'physical_particle_mass_or_relativistic_dynamics_identified':False}}
    if before!=source_bindings():raise ArithmeticError('source bytes changed during complete energy execution')
    return value


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('complete energy measurement receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write',action='store_true');action.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Complete native energy measurement receipt written.')
    else:
        verify();print('Complete native energy measurement receipt verified.')
