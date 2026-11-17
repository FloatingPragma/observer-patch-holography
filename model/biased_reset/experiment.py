"""Source-bound continuation and reset-optimality receipt."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from ..moving_defect.experiment import SOURCES,compare
from . import reset
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUTPUT=HERE/'results.json'


def source_bindings():
    paths=set(SOURCES)|{p.relative_to(ROOT).as_posix() for p in HERE.glob('*.py')}|{
        'model/biased_reset/README.md','model/completion/SAME_PROCESS_INTERFACES.md'}
    return bound_pins({name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)})


def run():
    fair=reset.renewal_certificate(True)
    if fair['maximum_mean_excursion_attempts']!='150' or fair['stationary_resets_per_native_attempt']!='1/150':
        raise ArithmeticError('fair-source control no longer reproduces its own law')
    return {'schema':'oph.biased_reset.v1','owner_issues':[1139,1140,1137],
            'scope':'ACTIVE_CONTINUATION_OF_THE_FIXED_ENDPOINT_BIASED_CODE_FAMILY',
            'source_bindings':bound_pins(source_bindings()),
            'renewal':reset.renewal_certificate(),
            'finite_horizon':reset.bellman_certificate(),
            'complete_service':reset.exact_service_certificate(),
            'continuing_resources':reset.resource_summary(),
            'continuing_clock_boundary':reset.clock_boundary_certificate(),
            'source_execution':reset.admitted_source_report(),
            'fair_control':{'maximum_mean_excursion_attempts':'150','stationary_resets_per_native_attempt':'1/150',
                            'mean_boundary_time_on_distance_two_cores':'120',
                            'fair_and_biased_fixed_sources_identified':False},
            'retained_obstructions':{
                'all_moving_codes_clock_invariant':False,
                'fixed_reference_clock_invariant_cores':reset.t.clock_compatibility()['invariant_core_count'],
                'reference_like_full_exterior_obtained':False,
                'minimum_full_exterior_distance':min((r['sharp_q1_full_exterior_distance_minimum'] for r in reset.t.exterior_certificate()['rows']),key=reset.F),
                'topological_particle_or_physical_mass_assigned':False},
            'integrated_gate_accepted':False,'physical_comparison_count_added':0}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    modes=p.add_mutually_exclusive_group(required=True)
    modes.add_argument('--write',nargs='?',const=OUTPUT,type=Path)
    modes.add_argument('--check',nargs='?',const=OUTPUT,type=Path)
    args=p.parse_args();result=run()
    if args.write is not None:
        args.write.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n',encoding='utf8')
    if args.check is not None:
        compare(strict_json.load(args.check),result)
    print(json.dumps({'scope':result['scope'],'reset_target':result['renewal']['target_core'],
                      'asymptotic_reset_rate':result['renewal']['stationary_resets_per_native_attempt'],
                      'integrated_gate_accepted':False},indent=2))


if __name__=='__main__':
    main()
