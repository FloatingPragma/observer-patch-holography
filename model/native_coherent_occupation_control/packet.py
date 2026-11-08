"""Replay exact source controls, native Bell handoff and literal-word obstructions."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from model.native_fair_local_preparation import packet as inherited
from . import bell, derived_control as beam, diagonal_control as phase
from . import literal_angles, literal_controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'


def encode(value):
    if type(value)is F:return str(value)
    if type(value)in (tuple,list):return [encode(x)for x in value]
    if type(value)is dict and all(type(k)is str for k in value):
        return {k:encode(v)for k,v in value.items()}
    if value is None or type(value)in (str,int,bool):return value
    raise TypeError('exact typed receipt values required')


def canonical(value):
    return json.dumps(encode(value),sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths=set(inherited.source_bindings())
    paths.update(p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in ('.py','.md'))
    pending=list(paths);visited=set()
    def add(module):
        if not module.startswith(('model.','tools.')):return
        parts=module.split('.');stem='/'.join(parts)
        candidates=[stem+'.py',stem+'/__init__.py']+['/'.join(parts[:i])+'/__init__.py'for i in range(1,len(parts))]
        for name in candidates:
            if (ROOT/name).is_file()and name not in paths:paths.add(name);pending.append(name)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name);module=name[:-3].replace('/','.')
        package=module.removesuffix('.__init__')if module.endswith('.__init__')else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(target,package)if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in sorted(paths)}


def stable_bell_certificate(level):
    result=beam.coupled_bell(level)
    # Retain exact analytic bounds and the pass/fail checks on actual stored
    # arrays. Their platform-dependent last bits are diagnostics, not labels
    # in the cross-platform mathematical receipt.
    omitted=('stored_success_mass','stored_complete_output_mass','stored_conditional_Bell_exact')
    value={k:v for k,v in result.items()if k not in omitted}
    value.update(stored_success_mass_certified_above_4999_over_10000=True,
                 stored_complete_mass_checked_without_normalization=True,
                 stored_conditional_Bell_equality_claimed=False,
                 stored_trace_preservation_claimed=False)
    return value


def stable_phase_certificate(kind):
    result=phase.execute(phase.declare(kind,F(1,2)))
    fields=('declaration','operator_error_upper','complete_halfdiamond_upper',
            'complete_map_gain_upper','exact_stored_number_preservation')
    return {**{k:result[k]for k in fields},'complete_complement_rows':2,
            'full_six_level_response_passes_directed_source_enclosure':True,
            'arbitrary_precision_logical_synthesis_claimed':False}


def produce():
    before=source_bindings()
    beam.source.specification.cache_clear()
    beam.compiler._source.cache_clear()
    beam.compiler._native.cache_clear()
    literal_controls.primitive.cache_clear()
    exact_native=[]
    for level in (0,1):
        for reference in (False,True):
            value,_=bell.branch_run(level,reference);exact_native.append(value)
        if not bell.charge_block_check(level):raise ArithmeticError('source charge-block identity failed')
    value={'schema':'oph.native_coherent_occupation_control.v1','source_bindings':before,
           'exact_six_level_exchange_generator':beam.exact_identity(),
           'exact_six_level_phase_generators':{kind:phase.exact_identity(kind)for kind in phase.KINDS},
           'complete_native_even_parity':exact_native,
           'executed_addressed_Bell_services':[stable_bell_certificate(level)for level in (0,1)],
           'executed_diagonal_services':[stable_phase_certificate(kind)for kind in phase.KINDS],
           'literal_word_angle_obstruction':literal_angles.certificate(),
           'literal_full_source_control_obstruction':literal_controls.certificate(),
           'complete_22_pulse_service_bounds':bell.beam_execution_certificate(),
           'ideal_balanced_retry_budget':bell.retry_budget(8,16),
           'cell_charge_Dicke_fidelity_ceiling_V240_q2':'6/1195',
           'cell_charge_Dicke_trace_distance_lower_V240_q2':'1189/1195',
           'scope':{'explicit_versioned_M1_coherent_service':True,
                    'frozen_six_kind_native_word_schema_changed':False,
                    'same_exact_twelve_source_generators':True,
                    'occupation_wire_coherent_loading_premise_retained':True,
                    'known_A0_B1_preparation_interface_retained':True,
                    'fixed_q_handoff_without_projection_or_renormalization':True,
                    'complete_failure_reference_and_record_outputs_retained':True,
                    'pulse_and_native_attempt_counters_distinct':True,
                    'bare_cross_cell_charge_coherence_prepared':False,
                    'growing_transfer_reference_or_ground_prepared':False,
                    'universal_logical_network_executed':False,
                    'relativistic_field_or_particle_mass_identified':False}}
    if source_bindings()!=before:raise ArithmeticError('source bytes changed during coherent-control replay')
    return encode(value)


def check_payload(retained,fresh):
    if canonical(retained)!=canonical(fresh):raise ArithmeticError('coherent-control receipt differs from fresh exact replay')
    return True


def verify(path=OUTPUT):
    return check_payload(strict_json.load(Path(path)),produce())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        value=produce();OUTPUT.write_text(canonical(value),encoding='utf-8')
        check_payload(strict_json.load(OUTPUT),value)
        print('Native coherent occupation-control receipt written.')
    else:
        verify();print('Native coherent occupation-control receipt verified.')


if __name__=='__main__':main()
