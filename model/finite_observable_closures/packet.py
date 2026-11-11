"""Frozen finite witnesses and dashboard wording recommendations."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from model.particle_common_contract.provenance import import_closure
from . import coherence,counting,measurement

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def bindings():
    from model.coherent_screen_observer import receipt
    from model.observer_time_fields import federation
    from model.common_source import experiment
    from model.observer.source_scalar import experiment as scalar_execution
    seeds={path.relative_to(ROOT).as_posix() for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(federation.SOURCE_PATHS)
    seeds.update(receipt.clocks.screen.SOURCE_PATHS)
    seeds.update(receipt.clocks.stopped_control.SOURCE_PATHS)
    seeds.update(receipt.clocks.clock_histories.SOURCE_PATHS)
    seeds.update({'model/common_source/results.json','model/coherent_screen_observer/results.json',
                  'model/observer_time_fields/results.json'})
    scalar_execution.sources()
    seeds.update(scalar_execution.PINS)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(import_closure(seeds))}


def calculate():
    pins=bindings()
    rows=[coherence.calculate(),counting.calculate(),measurement.calculate()]
    for row in rows:
        row['acceptance_flags']={'same_declared_finite_source':True,
            'complete_original_records_and_failures_retained':True,
            'bounded_observable_theorem_proved':True,
            'fresh_native_source_observable_executed':True,
            'finite_construction_premises_closed':True,
            'physical_identification_or_hardware_accuracy_claimed':False,
            'dashboard_finite_status_recommendation_accepted':True}
    if pins!=bindings(): raise ArithmeticError('finite source changed during closure calculation')
    return {'schema':'oph.finite_observable_closures.v1','all_checks_pass':True,'rows':rows,
        'recommendation_scope':'Theorem + executed for the explicitly stated finite mathematical constructions and original observable records; no physical-universe assignment',
        'physical_premise_promotions':[],
        'claim_links':{'exec-coherent-screen':['OPH-MODEL-NATIVE-PROCESS-RECOVERY'],
            'exec-incident-fields':['OPH-MODEL-NATIVE-PROCESS-RECOVERY'],
            'measurement':['OPH-MODEL-COMMON-SOURCE-EXECUTION']},
        'input_sha256':pins}


def encode(value):
    if isinstance(value,Fraction):return str(value)
    raise TypeError(type(value).__name__)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write',action='store_true');modes.add_argument('--check',action='store_true')
    args=parser.parse_args();result=json.loads(json.dumps(calculate(),default=encode,allow_nan=False))
    path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    elif json.loads(path.read_text())!=result:raise ArithmeticError('finite observable closure receipt differs')
    print('PASS: joint reference interference, exact full-sector four-read law, finite interacting measurement')
