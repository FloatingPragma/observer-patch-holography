"""Reproduce complete finite native clock observables and their source bounds."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from tools import strict_json
from . import observables as o

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
SOURCE_FILES=tuple('model/native_clock_observables/'+name for name in
    ('__init__.py','observables.py','packet.py','PROOF.md','README.md','test_observables.py'))+(
    'model/native_program/interpreter.py',
    'model/local_repair_family/mixing.py',
    'model/native_modular_reference/phase_readout.py',
    'model/native_modular_reference/test_phase_readout.py',
    'Lean/EventAlgebra/Robertson.lean',
    'Lean/EventAlgebra/OperationalBornBridge.lean',
    'tools/strict_json.py',
    'model/requirements.txt')


def source_bindings():
    result={path:digest for path,digest in o.learned_flow.source_sha256().items()
            if not path.startswith('docs/')}
    result.update({path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in SOURCE_FILES})
    return dict(sorted(result.items()))


def execution_check(d, arrays, protocols):
    rng=np.random.default_rng(291)
    v=rng.normal(size=132)+1j*rng.normal(size=132);v/=np.linalg.norm(v)
    rho=np.outer(v,v.conj())
    for selected in o.AXES:
        for corrected in (False,True):
            actual=o.compiled_output(d,rho,selected,arrays,corrected=corrected)
            ideal=o.ideal_output_with_same_exports(d,rho,selected,corrected=corrected)
            if actual.keys()!=ideal.keys():
                raise ArithmeticError('complete execution lost an output address')
            error=sum(np.linalg.norm(actual[k]-ideal[k],'nuc') for k in actual)/2
            label='corrected_worst_branch' if corrected else 'raw'
            allowed=float(o.F(protocols[selected][label]['complete_halfdiamond_upper']))
            mass=sum(np.trace(block).real for block in actual.values())
            if error>allowed+1e-13 or abs(mass-1)>2*allowed+1e-13:
                raise ArithmeticError('literal execution exceeds its complete channel allowance')
    return {'untouched_reference_dimension':2,'axis_and_correction_cases':6,
            'all_raw_slots_and_address_export_aborts_retained':True,
            'all_executed_cases_within_analytic_bound_plus_rounding_slack':True,
            'diagnostic_float_rounding_slack':'1/10000000000000',
            'arbitrary_reference_proof':'full-carrier exact branch identities and complete operator bounds'}


def produce():
    bindings=source_bindings();d=o.data();arrays=o.source_arrays()
    compiler=o.compile_certificate(arrays);protocols=o.protocol_certificate(compiler)
    result={'schema':'oph.native_clock_observables.v1',
            'source':o.read_program(),
            'instruments':o.instrument_certificate(d),
            'uncertainty':o.uncertainty_certificate(d),
            'modular_clock':o.modular_clock_certificate(d),
            'preparation':o.preparation_certificate(),
            'compiler':compiler,'protocols':protocols,
            'executed_complete_reference_test':execution_check(d,arrays,protocols),
            'scope':{'coherent_services_and_code_resources_declared':True,
                     'feedback_protocol_composition_declared':True,
                     'bare_repair_grammar_implements_coherent_protocol':False,
                     'passive_repair_derives_coherent_primitives':False,
                     'growing_capacity_service_embedding_proved':False,
                     'locality_is_completed_two_port_algebra':True,
                     'spatial_locality_of_intermediate_compiler_pulses_proved':False,
                     'universal_operational_Born_premises_supplied':False,
                     'physical_CCR_or_particle_identification':False},
            'source_bindings':bindings}
    if source_bindings()!=bindings:
        raise ArithmeticError('native clock observable sources changed during replay')
    return result


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if strict_json.load(Path(path))!=produce():
        raise ValueError('native clock observables receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native clock observables receipt written.')
    else:
        verify();print('Native clock observables receipt verified.')
