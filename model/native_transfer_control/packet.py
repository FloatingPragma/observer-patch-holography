"""Replay source-derived controlled transfer-ground experiments."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from tools import strict_json
from model.native_coexistent_effects import packet as phase_packet
from model.native_effect_synthesis import packet as synthesis_packet
from . import control as c, spectral

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_transfer_control/'+n for n in (
    '__init__.py','spectral.py','control.py','packet.py','proposal.json','test_control.py','PROOF.md','README.md'))+(
    'model/native_quantum_history/README.md','model/native_quantum_history/construction.py',
    'model/native_fixed_algebra/README.md','tools/strict_json.py')

# Explicit transitive import custody; no session-dependent binding discovery.
IMPORT_FILES=(
    'model/calibrated_measurement/__init__.py',
    'model/calibrated_measurement/bounds.py',
    'model/clock_defect_energy/__init__.py',
    'model/clock_defect_energy/contract.py',
    'model/clock_defect_energy/defect_history.py',
    'model/clock_defect_energy/joins_contract.py',
    'model/common_source/chiral_candidate/__init__.py',
    'model/common_source/chiral_candidate/action.py',
    'model/common_source/chiral_candidate/gaussian_source.py',
    'model/common_source/halfline.py',
    'model/common_source/mass_response/__init__.py',
    'model/common_source/mass_response/history.py',
    'model/common_source/mass_response/model.py',
    'model/common_source/native_code_transport.py',
    'model/common_source/native_gauge_object.py',
    'model/common_source/native_matter_attachment.py',
    'model/common_source/native_sm_current.py',
    'model/common_source/parameter_selection.py',
    'model/common_source/particle_identification/__init__.py',
    'model/common_source/particle_identification/family_attachment.py',
    'model/common_source/record_code.py',
    'model/common_source/repair_field_attachment.py',
    'model/defect_protection/__init__.py',
    'model/defect_protection/protection.py',
    'model/exchange_sector/__init__.py',
    'model/gauge_defect_transport/__init__.py',
    'model/gauge_defect_transport/transport.py',
    'model/interacting_defect_response/__init__.py',
    'model/interacting_defect_response/algebra.py',
    'model/native_equilibration/common_variation.py',
    'model/native_mass_channels/__init__.py',
    'model/native_mass_channels/algebra.py',
    'model/native_modular_reference/history_transport.py',
    'model/native_modular_reference/learned_clock.py',
    'model/native_modular_reference/phase_inference.py',
    'model/native_modular_reference/phase_readout.py',
    'model/native_modular_reference/reference_tomography.py',
    'model/native_quantum_gauss/__init__.py',
    'model/native_quantum_gauss/algebra.py',
    'model/observer/adaptive_record_recovery.py',
    'model/observer/native_loop_readout.py',
    'model/particle_gauge_band/__init__.py',
    'model/particle_gauge_band/quantum_link_completion.py',
    'model/radial_defect_gap/__init__.py',
    'model/radial_defect_gap/bounds.py',
    'model/scalar_observations/experimental_readout/stats.py',
)


def source_bindings():
    paths=set(SOURCE_FILES)|set(IMPORT_FILES)|set(phase_packet.source_bindings())|set(c.response.source_bindings())
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def word_record(word):
    checked=c.syn.certify_word(word); admitted=c.syn.interpreter.source_word(word,9)
    stopped=c.syn.interpreter.execute_word(admitted,np.zeros((1,1)))
    if stopped['status']!='RESOURCE_STOP' or stopped['pulses_executed']!=0:
        raise ArithmeticError('the unchanged generic word host boundary failed')
    return {'word_sha256':synthesis_packet.word_digest(word),
            'native_pulses':len(word),'distinct_primitive_checks':checked,
            'generic_host_status':'RESOURCE_STOP','generic_host_pulses_executed':0}


def produce():
    before=source_bindings(); cert=spectral.certify(); data=c.response.data()
    _,pi=spectral.source_intervals(data); omega=np.sqrt(np.array(pi,dtype=float))
    run=c.detector.run(17,tuple(i%2 for i in range(11)))
    weight,mapping,records=c.detector.branch(run['slots'])
    if mapping!={17:run['final']} or weight!=run['weight']:
        raise ArithmeticError('the complete preparation record decoder changed')
    c.syn.certify_primitive.cache_clear()
    ground=c.ground_word(run['final']); ground_check=word_record(ground)
    # All possible conditional X routes use these same nine checked native gates.
    for bit in range(9):c.syn.certify_word(c.syn.tomo.route_word((('X',bit),)))
    routes=[]
    for terminal,code in enumerate(c.syn.tomo.codes()):
        operations=tuple(('X',bit) for bit in range(9) if code>>bit&1)
        remaining=code
        for _,bit in operations:remaining^=1<<bit
        if remaining!=0 or len(operations)>9:raise ArithmeticError('a conditional source reset route failed')
        routes.append({'terminal':terminal,'routing_pulses':2*len(operations)})
    initial=np.zeros((66,1),complex);initial[run['final'],0]=1
    prepared,padding=c.menus.successful_word(ground,initial)
    if np.linalg.norm(prepared[:,0]-omega)>=1e-8:
        raise ArithmeticError('the executed source-ground preparation exceeds its diagnostic gate')
    current=[]; kick_checks=[]
    for enabled in (False,True):
        word=c.intervention_word(enabled); kick_checks.append(word_record(word))
        output,abort=c.menus.successful_word(word,prepared);current.append(output[:,0])
        if float(np.trace(abort).real)>=1e-8:raise ArithmeticError('the intervention padding diagnostic failed')
    current=np.column_stack(current); phase_checks=[]
    for index,word in enumerate(c.evolution_words()):
        phase_checks.append({'index':index,**word_record(word)})
        current,abort=c.menus.successful_word(word,current)
        if float(np.trace(abort).real)>=1e-8:raise ArithmeticError('a transfer phase padding diagnostic failed')
    P=np.array(data['P'],float); ps=np.array(pi,float)
    S=P*np.sqrt(ps)[None,:]/np.sqrt(ps)[:,None]
    vals,basis=np.linalg.eigh(S);unitary=(basis*np.exp(1j*float(c.TIME)*np.log(vals)))@basis.T
    ideal=np.column_stack((omega,omega.astype(complex)));ideal[c.TARGET,1]*=1j
    ideal=unitary@ideal
    if np.linalg.norm(current-ideal)>=1e-7:
        raise ArithmeticError('the executed complete successful word differs from its derived transfer target')
    increment=abs(current[c.TARGET,1])**2-abs(current[c.TARGET,0])**2
    if increment<=1/2000:raise ArithmeticError('the executed source score diagnostic lost its nonzero margin')
    # The final E instrument is complete and its decoded premeasurement label
    # supplies the configuration projector, rather than its disturbed terminal label.
    for initial_label in range(66):
        final_run=c.detector.run(initial_label,(1,)*11)
        if c.readout_score(final_run['slots'])!=int(initial_label==c.TARGET):
            raise ArithmeticError('the final measured score is not the intended source effect')
    value={'schema':'oph.native_transfer_control.v1','source':c.response.source(),
           'spectral_residual_certificate':cert,'complete_protocol_budget':c.protocol_budget(),
           'ground_response':c.witness(),'complete_outcomes':c.outcomes(),
           'actual_compilation':{'ground':ground_check,'all_terminal_routes':routes,'interventions':kick_checks,
                                 'transfer_phases':phase_checks,
                                 'successful_product_diagnostic_gate':'1/10000000',
                                 'ground_diagnostic_gate':'1/100000000',
                                 'one_live_preparation_path_diagnostic':list(run['slots']),
                                 'full_preparation_instrument_bound_is_analytic':True,
                                 'normalized_success_postselection_used_in_claim':False,
                                 'floating_residual_digits_serialized':False},
           'scope':{'controlled_source_ground_preparation':True,
                    'record_derived_transfer_unitary_compiled':True,
                    'complete_preparation_intervention_evolution_readout_protocol':True,
                    'arbitrary_untouched_reference_in_complete_bound':True,
                    'preparation_archive_retained_but_can_correlate_with_reference':True,
                    'all_emitted_stored_primitive_checks_pass':True,
                    'rounded_maps_asserted_exactly_trace_preserving':False,
                    'floating_propagation_certified_or_hardware_error_included':False,
                    'generic_evaluator_limits_changed':False,
                    'bare_repair_is_autonomous_unitary':False,
                    'working_mixed_reference_equals_transfer_ground':False,
                    'native_modular_and_transfer_clocks_identical':False,
                    'growing_coherent_service_embedding':False,
                    'physical_vacuum_relativistic_clock_or_particle_masses':False},
           'source_bindings':before}
    if source_bindings()!=before:raise ArithmeticError('source inputs changed during transfer replay')
    return value


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native transfer-control receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native transfer-control receipt written.')
    else:
        verify();print('Native transfer-control receipt verified.')
