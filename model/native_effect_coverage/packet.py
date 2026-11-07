"""Reproduce complete weighted-native effect and preparation coverage."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from tools import strict_json
from . import detector as d, coverage as c
from model.native_causal_dataflow import dataflow

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_effect_coverage/'+name for name in (
    '__init__.py','detector.py','coverage.py','packet.py','test_coverage.py','PROOF.md','README.md'))+(
    'model/native_program/interpreter.py','model/local_repair_family/local_bias.py',
    'model/local_repair_family/modular_clock.py','model/native_causal_dataflow/dataflow.py',
    'model/observer/native.py','Lean/EventAlgebra/OperationalBornBridge.lean',
    'Lean/EventAlgebra/OperationalAdditivityBoundary.lean',
    'Lean/EventAlgebra/FiniteBuschGleason.lean','tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|{p for p in c.t.SOURCE_PATHS if not p.startswith('docs/')}
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def execution_certificate():
    rows=[]
    for initial,coins in ((0,(0,)*11),(65,(1,)*11),(29,tuple(i%2 for i in range(11)))):
        record=d.run(initial,coins)
        density={(initial,0,initial,0):F(1)}
        output=d.continue_raw(record['slots'],density)
        if output['density']!=d.apply_path(record['slots'],density):
            raise ArithmeticError('continuation kernel and exact fine instrument disagree')
        authenticated=dataflow.replay_trace(d.program(),d.local.basis()[initial],record['records'])
        if len(authenticated)!=11 or output['archive']!=record['archive'] or output['memory']!=record['FIFO']:
            raise ArithmeticError('an actual record/FIFO continuation was omitted')
        rows.append({'initial':initial,'terminal':record['final'],'slots':list(record['slots']),
                     'weight':str(record['weight']),'authenticated_events':len(authenticated)})
    stopped=d.interpreter.execute(d.program(),{(0,0,0,0):F(1)})
    if stopped.status!='RESOURCE_STOP' or stopped.attempted_steps!=3 or len(stopped.branches)!=25**3:
        raise ArithmeticError('the stock expanded evaluator stop boundary changed')
    return {'complete_streamed_words':rows,'expanded_host_status':stopped.status,
            'expanded_attempts_executed':stopped.attempted_steps,'expanded_retained_blocks':len(stopped.branches),
            'host_stop_relabelled_as_success':False,
            'lazy_branch_query_is_normalized_outcome_sampler':False}


def compiled_execution_certificate(coherent):
    rng=np.random.default_rng(871)
    factor=rng.normal(size=(66,2))+1j*rng.normal(size=(66,2));factor/=np.linalg.norm(factor)
    work=np.zeros((512,2),complex);work[list(c.t.codes())]=factor
    pairs=((0,1),(0,65),(32,65));checked=0;structural_abort=False
    for a,b in pairs:
        for kind in ('mixer','y_analyzer'):
            actual=c.compiler.execute(c.t.word(a,b,kind),9,work)
            ideal=c.t.ideal(a,b,kind)@work
            # Complete address readout, with padding coarse-grained to abort.
            def records(array):
                result=[np.outer(array[code],array[code].conj()) for code in c.t.codes()]
                padding=sorted(set(range(512))-set(c.t.codes()))
                result.append(array[padding].T@array[padding].conj())
                return result
            raw,target=records(actual),records(ideal)
            error=sum(np.linalg.norm(x-y,'nuc') for x,y in zip(raw,target,strict=True))/2
            allowed=float(F(coherent['uniform_words']['analyzer']['complete_halfdiamond_upper']))
            if error>allowed+1e-13:
                raise ArithmeticError('literal source execution exceeds the complete source allowance')
            structural_abort |= np.trace(target[-1]).real>1e-5
            source=np.zeros(512,complex);source[c.t.codes()[a]]=1
            prepared=c.compiler.execute(c.preparation_word(a,b,kind),9,source)
            restricted=prepared[list(c.t.codes())]
            desired=c.ideal_preparation(a,b,kind)
            if np.linalg.norm(np.outer(restricted,restricted.conj())-desired,'nuc')/2>allowed+1e-13:
                raise ArithmeticError('literal inverse analyzer has the wrong prepared state')
            checked+=1
    if not structural_abort:
        raise ArithmeticError('execution diagnostic failed to retain structural padding output')
    return {'analyzer_and_inverse_preparation_pairs_checked':checked,
            'reference_dimension':2,'full_work_carrier_dimension':512,
            'structural_ideal_padding_abort_observed':True,
            'all_cases_within_theorem_bound_plus_diagnostic_rounding_slack':True,
            'diagnostic_rounding_slack':'1/10000000000000',
            'full_primitive_failure_mass_bounded_separately':True}


def produce():
    bindings=source_bindings();coherent=c.coherent_certificate()
    value={'schema':'oph.native_effect_coverage.v1','native_detector':d.certificate(),
           'continuation_and_host_boundary':execution_certificate(),
           'coherent_services':coherent,'preparation':c.preparation_certificate(),
           'effect_coverage':c.coverage_certificate(),'stable_actual_span':c.stable_span_certificate(coherent),
           'native_diagonal_coexistence':c.dyadic_diagonal_refinement(
               tuple(F(i%5,8) for i in range(66)),tuple(F((i+1)%3,8) for i in range(66)),3),
           'literal_source_execution':compiled_execution_certificate(coherent),
           'scope':{'source_law_changed':False,'same_fixed_weighted_H66_source':True,
                    'declared_coherent_services_and_work_resources_used':True,
                    'bounded_controller_composition_is_bare_repair_grammar':False,
                    'expanded_interpreter_completes_full_protocol':False,
                    'universal_Born_or_all_effect_coexistence_supplied':False,
                    'capacity_family_coherent_services_embedded':False,
                    'physical_QFT_or_particle_identification':False},
           'source_bindings':bindings}
    if source_bindings()!=bindings:
        raise ArithmeticError('effect coverage sources changed during replay')
    return value


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if strict_json.load(Path(path))!=produce():
        raise ValueError('native effect coverage receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native effect coverage receipt written.')
    else:
        verify();print('Native effect coverage receipt verified.')
