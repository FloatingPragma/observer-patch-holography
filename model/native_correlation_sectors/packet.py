"""Replay the exact native correlation-sector and reflected-record evidence."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..native_program import receipt as native_receipt
from ..federated_observer import scaling
from . import sectors as s, records as r

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')
SOURCE_FILES=tuple(sorted(set((*scaling.SOURCE_PATHS,
    'model/federated_observer/SCALING.md','model/native_equilibration/README.md',
    'model/native_equilibration/construction.py','model/observer/QUANTUM_MIXING.md',
    'model/observer_time_fields/dual.py','model/observer_time_fields/REFINEMENT.md',
    'model/native_effective_fields/construction.py','model/native_effective_fields/README.md',
    'model/weighted_source_family/__init__.py','model/weighted_source_family/family.py','model/weighted_source_family/PROOF.md',
    'model/observer/source_discrimination/circuit.py','model/observer/source_scalar/compiler.py',
    'tools/strict_json.py','model/requirements.txt',
    *(f'model/native_correlation_sectors/{name}' for name in (
        '__init__.py','sectors.py','records.py','packet.py','test_sectors.py','test_packet.py',
        'README.md','PROOF.md'))))))


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    names=set(SOURCE_FILES)|set(native_receipt.source_bindings())
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names) if not name.startswith('docs/')}


def attachment():
    interpreter=s.source.interpreter
    states=s.source.local.basis()
    masks=tuple(sum(value << i for i,value in enumerate(x)) for x in states)
    by_mask={mask:i for i,mask in enumerate(masks)}
    edge=s.source.local.SEAMS[0]
    commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(edge)},'law':'E','repeat':1}]
    program=interpreter.program(biased=False,commands=commands)
    marker={(i,i,i,i):F(1) for i in range(66)}
    actual=interpreter._outcomes(program,commands[0],marker)
    predicted=r.raw_maps(edge,12,masks)
    if len(actual)!=len(predicted):
        raise ArithmeticError('raw E address alphabet changed')
    for original,(flag,weight,mapping,_,_) in zip(actual,predicted,strict=True):
        expected={(by_mask[j],by_mask[i],by_mask[j],by_mask[i]):weight for i,j in mapping.items() if weight}
        actual_flag=('hold',) if original.record.kind=='hold' else (
            original.record.law,original.record.coin,original.record.raw)
        if actual_flag!=flag or original.density!=expected:
            raise ArithmeticError('full H66 E map differs from the native interpreter')
    graph,first,second,_=s.probe_data(0)
    commands=[{'op':'repair','selector':{'kind':'fixed','seam':list(e)},'law':'E','repeat':1}
              for e in (first,second)]
    federation=interpreter.program(federation={'level':0,'total':2},commands=commands)
    # All four local input patterns are embedded in the same actual total-two carrier.
    other=tuple(i for i in range(len(graph.weights)) if i not in first)[:2]
    cases=(other,(first[0],other[0]),(first[1],other[0]),first)
    masks=tuple(sum(1 << i for i in x) for x in cases)
    marker={(mask,i,mask,i):F(1) for i,mask in enumerate(masks)}
    actual=interpreter._outcomes(federation,commands[0],marker)
    for original,(flag,weight,mapping,_,_) in zip(actual,r.raw_maps(first,len(graph.weights),masks),strict=True):
        expected={(j,masks.index(i),j,masks.index(i)):weight for i,j in mapping.items() if weight}
        if original.density!=expected:
            raise ArithmeticError('full federation E map differs from the native interpreter')
    return {'H66_fair_raw_maps':25,'H66_source_columns':66,
            'federation_raw_maps':25,'federation_endpoint_patterns':4,
            'program_digest':interpreter.digest(r.probe_program(0,2,6)),
            'complete_probe_attempts':10,
            'kernel_attachment_is_full_expanded_history_execution':False,
            'expanded_history_and_FIFO_resource_stops_preserved':True}


def controls():
    rows=[]
    for ports,seams in ((4,((0,1),(1,2),(2,3))),
                        (6,((0,1),(0,5),(1,2),(2,3),(3,4),(4,5)))):
        family=s.source.Family((1,)*ports,seams,(1,)*len(seams))
        h=s.contrast((0,1),(2,3),ports)
        for occupied in range(ports+1):
            values=h;denominator=1
            for t in range(3):
                expected=(s.norm_factor(ports,2,occupied)*s.dot(h,values)/denominator
                          if 2<=occupied<=ports-2 else F(0))
                actual=r.direct_probe(family,occupied,t,(0,1),(2,3))
                if actual!=expected:
                    raise ArithmeticError('complete native records violate primitive response identity')
                rows.append({'ports':ports,'occupied':occupied,'attempts':t,'reflected_score':str(actual)})
                values=s.counts_step(family,values,2);denominator*=2*len(seams)
    path=s.source.Family((1,)*4,((0,1),(1,2),(2,3)),(1,)*3)
    naive=r.direct_probe(path,2,0,(0,1),(2,3),reflected=False)
    if naive!=0:
        raise ArithmeticError('incoming-past backaction control changed')
    return {'exact_control_rows':rows,'naive_incoming_past_t0':'0','reflected_past_t0':'2/3'}


def produce():
    before=source_bindings();levels=[]
    for level in (0,1,2):
        row=s.autocorrelation(level)
        floor=scaling.transfer_floor(level)
        if floor!=F(10,21):
            raise ArithmeticError('owning federation transfer positivity changed')
        ports=row['ports']
        row['variance_by_declared_occupation']={str(q):str(4*s.norm_factor(ports,2,q))
                                                for q in (2,4,ports//3,ports//2)}
        levels.append(row)
    value={'schema':'oph.native_correlation_sectors.v1',
           'source_bindings':before,'native_attachment':attachment(),
           'complete_reflected_record_identity':r.pair_certificate(),
           'finite_record_controls':controls(),'evaluated_federations':levels,
           'all_domain_identities':{
               'primitive_dimension':'binomial(V,r)-binomial(V,r-1)',
               'inclusion_squared_norm':'binomial(V-2r,q-r)',
               'uniform_variance':'4*binomial(V-4,q-2)/binomial(V,q)',
               'normalized_probe':'h^T P2^t h/4',
               'normalized_preparation_error':'2*epsilon/variance',
               'degree_range':'0<=r<=min(q,V-q)',
               'positive_transfer_floor_reused':'10/21'},
           'scope':{
               'source':'fair fixed-selector native binary population/history transfer',
               'orthogonal_blocks_exhaust_fixed_occupation':True,
               'degree_two_eigenvalues_disjoint_from_degree_one':False,
               'raw_quantum_branches_intertwine_all_occupations':False,
               'normalized_probe_independent_of_occupation_on_fixed_graph':True,
               'normalized_probe_independent_of_capacity':False,
               'initial_stationary_preparation_assumed_or_separately_bounded':True,
               'arbitrary_adaptive_selector_covered':False,
               'density_or_physical_clock_selected':False,
               'physical_particle_or_gauge_identification':False,
               'new_canonical_premise_discharged':False}}
    if before!=source_bindings():
        raise ArithmeticError('source bytes changed during replay')
    return value


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('native correlation-sector receipt differs from exact replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write',action='store_true');actions.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()));print('Native correlation-sector receipt written.')
    else:
        verify();print('Native correlation-sector receipt verified.')
