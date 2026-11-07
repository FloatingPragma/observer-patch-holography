"""Replay finite native local-response and dependency-cone certificates."""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import cone as c,bounds as b
from model.native_causal_dataflow import dataflow
from model.native_record_causality import causality
from model.native_program import receipt as program_receipt

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple('model/native_local_response/'+name for name in (
    '__init__.py','cone.py','bounds.py','packet.py','test_local_response.py','PROOF.md','README.md'))+(
    'model/native_program/interpreter.py','model/weighted_source_family/family.py',
    'model/weighted_source_family/PROOF.md','model/native_record_causality/causality.py',
    'model/native_record_causality/PROOF.md','model/native_causal_dataflow/dataflow.py',
    'model/local_repair_family/local_bias.py','model/observer/native.py',
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
    'tools/strict_json.py','model/requirements.txt')


def source_bindings():
    paths=set(SOURCE_FILES)|set(c.native.geometry.SOURCE_PATHS)
    paths.update(program_receipt.source_bindings())
    paths.update(('model/observer/source_discrimination/circuit.py',
                  'model/observer/source_scalar/compiler.py'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in sorted(paths) if not p.startswith('docs/')}


def fixed_program(edges,*,federation=False,overflow='export_oldest'):
    commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(e)},'law':'mixture','repeat':1} for e in edges]
    options={'commands':commands,'capacity':1,'overflow':overflow}
    if federation:options.update(federation={'level':0,'total':2},biased=False)
    return c.source.program(**options)


def source_attachment():
    marker={(i,i,i,i):F(1) for i in range(66)};checked=0
    for edge in c.source.local.SEAMS:
        value=fixed_program((edge,));family=c.source_family(value)
        actual=c.source._outcomes(value,value['commands'][0],marker)
        lifted=c.fixed_instrument(family,c.source_masks(value,marker),edge)
        if len(actual)!=25 or len(lifted)!=25:
            raise ArithmeticError('a complete native address was omitted')
        for original,(record,density) in zip(actual,lifted,strict=True):
            if original.record!=record or c.source_masks(value,original.density)!=density:
                raise ArithmeticError('literal H66 endpoint attachment differs')
            checked+=1
    graph=c.native.geometry.federation(0);edge=graph.inter_seams[0]
    outside=tuple(v for v in range(graph.ports) if v not in edge)[:2]
    pairs=(outside,(edge[0],outside[0]),(edge[1],outside[0]),edge)
    masks=tuple(sum(1<<v for v in pair) for pair in pairs)
    marker={(mask,i,mask,i):F(1) for i,mask in enumerate(masks)}
    value=fixed_program((edge,),federation=True);family=c.source_family(value)
    actual=c.source._outcomes(value,value['commands'][0],marker)
    lifted=c.fixed_instrument(family,marker,edge)
    if any(x.record!=record or x.density!=rho for x,(record,rho) in zip(actual,lifted,strict=True)):
        raise ArithmeticError('the actual fair federation endpoint lift differs')
    return {'weighted_H66_fixed_raw_maps_checked':checked,'source_basis_columns_per_map':66,
            'fair_federation_endpoint_patterns':4,'fair_federation_complete_raw_maps':25,
            'arbitrary_complex_matrix_units_by_partial_injection_identity':True,
            'weighted_federation_admission_claimed':False}


def actual_comparison(federation=False):
    edges=((0,5),(12,13),(0,1)) if federation else ((0,5),(2,3),(0,1))
    value=fixed_program(edges,federation=federation);family=c.source_family(value)
    states=((0,6),(1,6))
    masks=tuple(sum(1<<v for v in row) for row in states)
    if federation:indices=masks
    else:
        basis=c.source.local.basis()
        indices=tuple(basis.index(tuple(int(v in row) for v in range(12))) for row in states)
    rho={(indices[i],i,indices[j],j):F(1,2) for i in (0,1) for j in (0,1)}
    observed=(0,2);final=(0,1)
    result=c.source.execute(value,rho)
    actual=c.observe_history(value,result,final,observed)
    source_rho=c.source_masks(value,rho)
    full=c.instrument(family,source_rho,edges,final,observed)
    reduced=c.instrument(family,source_rho,edges,final,observed,reduce=True)
    if actual!=full or full!=reduced or len(actual)!=625:
        raise ArithmeticError('complete native local/reference history differs from its cone reduction')
    return {'carrier':value['carrier'],'program_digest':c.source.digest(value),
            'source_raw_histories':len(result.branches),'observed_raw_histories':len(actual),
            'observed_zero_blocks':sum(not block for block in actual.values()),
            'source_FIFO_capacity':1,'actual_exported_records_per_path':2,
            'unselected_record_indices_explicitly_traced':[1],
            'observed_indices':list(observed),'final_vertices':list(final),
            'exact_initial_cone':list(c.backward(edges,final,observed)['initial_support']),
            'actual_expanded_status':result.status,'arbitrary_reference_identity_supported':True}


def halo_certificate():
    census=[]
    for level in (0,1,2):
        graph=c.native.geometry.federation(level)
        rows=[b.halo_geometry(level,cell) for cell in range(graph.cells)]
        if len(graph.seams)!=630*4**level or graph.max_degree!=6:
            raise ArithmeticError('literal federation rate/degree changed')
        census.append({'level':level,'cells_checked':len(rows),'ports':graph.ports,
                       'seams':len(graph.seams),'maximum_degree':graph.max_degree,
                       'halo_ports_per_cell':3,'outward_edges_per_halo_port':5,
                       'all_halos_have_distinct_cell_ownership_and_no_internal_halo_edge':True})
    family=c.native.federation_family(0,(1,)*12);halo=b.halo_geometry(0,0)
    controls=[]
    for attempts in (0,1,2,3,6,12):
        exact=b.escape(family,halo['root'],halo['window'],attempts)
        reduced=b.halo_global_error(630,attempts)
        if exact['escape_probability']!=reduced['cone_escape_probability']:
            raise ArithmeticError('full source cone and four-state halo chain differ')
        controls.append({'attempts':attempts,'exact_escape':str(exact['escape_probability']),
                         'subset_states':exact['maximum_subset_states']})
    return {'literal_geometry':census,'full_subset_vs_four_state_controls':controls,
            'global_attempt_examples':[
                {'level':level,'attempts':30,**{k:str(v) if type(v) is F else v
                  for k,v in b.halo_global_error(630*4**level,30).items()}}
                for level in (0,1,2)],
            'local_clock_examples':[{k:str(v) if type(v) is F else v for k,v in
                                     b.halo_internal_error(n).items()} for n in (0,1,2,4,8)],
            'global_union_bound':'min(1,15*binomial(T,2)/M^2)',
            'local_clock_union_bound':'min(1,k*(k+1)/120)',
            'local_clock_kernel':[[str(x) for x in row] for row in b.internal_block()],
            'bounded_Abel_record_response':[
                {k:[str(x) for x in v] if type(v) is tuple else str(v) if type(v) is F else v
                 for k,v in b.abel_response_bound(q).items()} for q in (F(0),F(1,4),F(1,2),F(3,4))],
            'exact_clock_elapsed_attempts_are_retained_not_conditioned_away':True,
            'no_vanishing_error_claim_at_fixed_local_repair_count':True}


def boundaries():
    adaptive=dataflow.certificate()
    try:c.audit_program(dataflow.feedback_program())
    except ValueError:pass
    else:raise ArithmeticError('adaptive FIFO influence was silently removed')
    halted=fixed_program(((0,5),(0,1)),overflow='halt')
    try:c.audit_program(halted)
    except ValueError:pass
    else:raise ArithmeticError('intrinsic halt was silently removed')
    stopped=c.source.execute(fixed_program(((0,5),(2,3),(0,1),(1,7))),{(0,0,0,0):F(1)})
    if stopped.status!='RESOURCE_STOP' or stopped.attempted_steps!=3:
        raise ArithmeticError('actual expanded-tree host stop changed')
    return {'actual_adaptive_receiver_probabilities':adaptive['receiver_probabilities'],
            'controller_reads_are_included_in_actual_causal_dataflow':True,
            'adaptive_selector_rejected_by_spatial_bound':True,
            'intrinsic_FIFO_halt_rejected':True,
            'four_step_expanded_status':stopped.status,
            'four_step_expanded_attempts':stopped.attempted_steps,
            'host_stop_is_not_relabelled_as_completed_kernel_family':True,
            'admission_attempt_limit':c.source.MAX_ATTEMPTS,
            'expanded_attempt_limit':c.source.MAX_EVALUATED_STEPS,
            'expanded_branch_limit':c.source.MAX_EVALUATED_BRANCHES,
            'exact_probability_evaluator_attempt_limit':256,
            'full_remote_archive_observation_covered':False,
            'reading_source_or_scheduler_specifications_forbidden':False}


def produce():
    bindings=source_bindings()
    result={'schema':'oph.native_local_response.v1','source_attachment':source_attachment(),
            'complete_local_history_comparisons':[actual_comparison(),actual_comparison(True)],
            'literal_cell_halo':halo_certificate(),'controller_and_evaluator_boundaries':boundaries(),
            'scope':{'fixed_input_independent_selector_required':True,
                     'arbitrary_correlated_references_in_exact_cone_reduction':True,
                     'initial_state_comparison_requires_joint_window_reference_bound':True,
                     'unreferenced_product_limit_promoted_to_correlated_reference_limit':False,
                     'finite_volume_comparison_keeps_original_selector_rates':True,
                     'weighted_federation_is_admitted':False,
                     'physical_spatial_dimension_or_Lorentz_clock_derived':False,
                     'particle_mass_or_interacting_QFT_claimed':False},
            'source_bindings':bindings}
    if source_bindings()!=bindings:raise ArithmeticError('local-response sources changed during replay')
    return result


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native local-response receipt differs from source replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native local-response receipt written.')
    else:
        verify();print('Native local-response receipt verified.')
