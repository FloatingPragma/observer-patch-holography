"""Exact source energy and local-clock receipt with authenticated inputs."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.completion import native_measurement_parent_adapter as parent_custody
from . import energy
from ..native_exchange_clock import packet as exchange_owner
from ..native_energy_measurement import packet as measurement_owner
from ..native_energy_measurement import measurement, spectral
from model.observer.quantum_mixing import seams
from model.federated_observer import source

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_local_energy_clock/'+name for name in (
    '__init__.py','energy.py','packet.py','test_energy.py','test_packet.py',
    'PROOF.md','README.md'))
# Full 66-word replay established at 97633a11f; audited metadata-only renewals
# preserve its scientific payload. A live input check remains mandatory.
AUDITED_MEASUREMENT_SHA256 = '6c4f66f29e34eac29bac5d7ebe1c857606888eede03113761556c41eb1b0ffeb'


def canonical(value):
    return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def source_bindings():
    paths=set(exchange_owner.source_bindings())|set(measurement_owner.source_bindings())|set(SOURCE_FILES)|set(parent_custody.SOURCE_FILES)
    paths.add(measurement_owner.OUTPUT.relative_to(measurement_owner.ROOT).as_posix())
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def owning_measurement():
    """Authenticate the prior full replay and recompute its consumed bounds."""
    try:
        digest=parent_custody.authenticate(measurement_owner.OUTPUT,
                                          AUDITED_MEASUREMENT_SHA256,
                                          measurement_owner.source_bindings())
    except ValueError as error:
        raise ValueError('audited complete measurement receipt bytes differ: '+str(error)) from error
    value=strict_json.load(measurement_owner.OUTPUT)
    expected={'source_bindings':measurement_owner.source_bindings(),
              'spectral_certificate':spectral.certify(),
              'complete_instrument_budget':measurement.budget(),
              'counted_energy_label_intervals':[[str(x) for x in pair] for pair in spectral.energy_labels()]}
    if any(canonical(value.get(key))!=canonical(wanted) for key,wanted in expected.items()):
        raise ValueError('owning complete measurement or source bounds differ')
    words=value.get('all_native_words',[])
    if ([row.get('stage') for row in words]!=list(range(66))
            or any(type(row.get('stage')) is not int for row in words)
            or not value['scope']['all_emitted_stored_primitives_certified']
            or value['scope']['rounded_instrument_is_exactly_trace_preserving']):
        raise ValueError('complete measured-word evidence differs')
    return {'owning_receipt_sha256':digest,
            'source_bindings_and_consumed_exact_bounds_recomputed':True,
            'owning_full_66_word_execution_repeated_by_this_packet':False}


def source_controls():
    edges=tuple(sorted(seams()));edge_set=set(edges)
    from itertools import combinations
    triangles=tuple(t for t in combinations(range(12),3)
                    if all(tuple(sorted(e)) in edge_set for e in combinations(t,2)))
    if len(triangles)!=20 or any(sum(set(edge)<=set(t) for t in triangles)!=2 for edge in edges):
        raise ArithmeticError('native internal triangle cover differs')
    rows=[]
    for level in (0,1):
        family=source.federation(level)
        counts=energy.graph_data(family.ports,tuple(sorted(family.seams)))
        if counts!=(630*4**level,6,5400*4**level):
            raise ArithmeticError('literal source graph counts differ')
        rows.append({'level':level,'ports':family.ports,'seams':counts[0],
                     'maximum_degree':counts[1],'ordered_adjacent_seam_pairs':counts[2]})
    return {'native_internal_seams':len(edges),'triangles':len(triangles),
            'each_internal_seam_in_two_triangles':True,'literal_federations':rows}


def serial(value):
    if type(value) is F:
        return str(value)
    if type(value) is dict:
        return {k:serial(v) for k,v in value.items()}
    if type(value) in (list,tuple):
        return [serial(v) for v in value]
    if value is None or type(value) in (str,int,bool):
        return value
    raise TypeError('exact scalar or container required')


def compute():
    rows=[]
    for level in (0,6,12,24,100):
        size=4**level
        row=energy.local_clock_bound(240*size,630*size,5400*size,80*size,6,1,F(1))
        rows.append({'level':level,'ports':240*size,'seams':630*size,
                     'total':80*size,'time':'1',**row})
    histories=[]
    for level in (12,24,100):
        n=4**level
        histories.append({'level':level,**energy.finite_history_bound(
            240*n,630*n,5400*n,80*n,6,(1,2,1),
            (('first','second','reserved_zero'),)*3,(F(1),)*3)})
    return {'schema':'oph.native_local_energy_clock.v1',
        'source_reference':'rho=P/Tr(P), P=I-G/(2M)',
        'source_transfer_clock':'t=N/(2M), H=-2M log(P)',
        'density_limit_controlled_clock':'G/[1-p(1-p)]',
        'exact_H66':energy.exact_energy_certificate(12,tuple(sorted(seams())),2),
        'complete_energy_detector_heating':energy.h66_detector_heating_bound(),
        'rational_local_clock_bounds':rows,
        'complete_fixed_finite_controlled_history_bounds':histories,
        'scope':{'bound_norm':'sqrt(Tr(rho X*X))',
                 'fixed_local_support_bound':True,
                 'complete_fixed_number_sector':True,
                 'native_unread_repair_is_controlled_unitary':False,
                 'whole_growing_reference_preparation_supplied':False,
                 'operator_norm_or_diamond_convergence':False,
                 'relativistic_or_vacuum_limit_supplied':False,
                 'symmetric_physical_stress_supplied':False,
                 'large_level_rows_are_arithmetic_bounds_not_expanded_executions':True}}


def produce():
    before=source_bindings()
    value=serial(compute())
    value['source_controls']=source_controls()
    value['owning_measurement_attachment']=owning_measurement()
    value['source_bindings']=before
    if source_bindings()!=before:
        raise ArithmeticError('source bytes changed during energy-clock replay')
    return value


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path)))!=canonical(produce()):
        raise ValueError('source energy-clock receipt differs from replay')
    return True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write',action='store_true')
    actions.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Exact local energy and controlled clock receipt written.')
    else:
        verify()
        print('Exact local energy and controlled clock receipt verified.')
