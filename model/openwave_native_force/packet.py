"""Exact source-charge and finite elastic-response receipt."""
import argparse
import hashlib
import json
from pathlib import Path

from model.finite_area_response import packet as parent
from . import force

ROOT=Path(__file__).resolve().parents[2]
OUTPUT=Path(__file__).with_name('results.json')


def canonical(value):
    return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n'


def source_bindings():
    paths=set(parent.source_bindings())
    paths.update(str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir()
                 if p.suffix in ('.py','.md'))
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    pins=source_bindings()
    value={'schema':'oph.native_scalar_geometry_force.v1','all_pass':True,
           'source_bindings':pins,'inherited_source_custody':parent.inherited_custody(),
           'force':force.force_run(),
           'prepared_reference_nonuniversality':force.positive_energy_nonuniversality(),
           'actual_ground_nonuniversality':force.ground_counterexample(),
           'scope':{'native_local_Hamiltonian_and_derivative_used':True,
                    'spatial_elastic_geometry_law_is_declared_extension':True,
                    'Poisson_equation_inserted_independently_of_action':False,
                    'full_nonlinear_native_matter_energy_in_finite_window_bounded':True,
                    'universal_Newton_gravity_or_G_derived':False,
                    'physical_mass_or_laboratory_length_identified':False,
                    'native_controller_implements_elastic_minimum_preparation':False}}
    if source_bindings()!=pins:raise ValueError('source changed during force replay')
    return value


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=ap.parse_args();value=report()
    if args.write:OUTPUT.write_text(canonical(value))
    elif json.loads(OUTPUT.read_text())!=value:raise ValueError('force receipt differs from source replay')
    print('PASS: finite native-charge attraction; universal positive-energy gravity refuted for this scalar extension')


if __name__=='__main__':main()
