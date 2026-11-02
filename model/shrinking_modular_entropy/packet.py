"""Reproduce the natural-shell entropy theorem's finite certificates."""
import argparse
import hashlib
import json
from pathlib import Path
from tools import strict_json
from . import algebra, bounds, controls

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUTPUT=HERE/'results.json'


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def source_bindings():
    paths={p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({
        'model/__init__.py', 'tools/strict_json.py',
        'model/shrinking_modular_entropy/README.md',
        'paper/tex_fragments/SHRINKING_MODULAR_ENTROPY.tex',
        'model/lattice_modular_limit/README.md',
        'model/lattice_weyl_response/README.md',
        'model/lattice_weyl_response/algebra.py',
        'model/weyl_gravity_response/README.md',
        'model/weyl_gravity_response/algebra.py',
        'model/weyl_gravity_response/bounds.py',
        'model/weyl_geometric_modular/README.md',
        'model/continuum_field_net/cells.py',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
    })
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def theorem():
    return {
        'source_U':0, 'source_valleys_retained':8,
        'natural_exterior_shell_logarithm_control':True,
        'actual_cell_norm_used':True,
        'uniform_covariance_gap_assumed':False,
        'explicit_positive_symbolic_refinement_schedule':True,
        'ordinary_finite_lattice_entropy_differences':True,
        'normalized_shrinking_matter_error_tends_to_zero':True,
        'observer_frame_preparation_matrix_invertible':True,
        'finite_lattice_regions':'ordinary rest balls for every observer preparation',
        'finite_observer_rows':'distinct Lorentz-transformed exterior shell preparations',
        'finite_lattice_boost_symmetry_constructed':False,
        'single_lattice_family_of_boosted_sharp_algebras_constructed':False,
        'absolute_continuum_entropy_computed':False,
        'nonzero_U_regional_entropy_limit_proved':False,
        'area_Newton_or_stationarity_input_supplied':False,
        'Einstein_endpoint_closed':False,
        'native_entropy_circuit_executed':False,
        'astronomical_refinement_family_numerically_executed':False,
    }


def custody():
    return {
        'current_narrow_mathematical_sources_checked':True,
        'predecessor_full_repository_receipts_replayed':False,
        'historical_receipts_reinterpreted':False,
        'new_proof_requires_no_runtime_Git_history':True,
        'scope':'new analytic proof, exact finite controls and recorded numerical diagnostics only',
    }


def validate(value):
    if type(value) is not dict or set(value)!={'schema','finite_attachment','new_checks_pass',
            'full_repository_replay_performed','theorem','bounds','algebra','controls','custody',
            'source_bindings'}:
        raise ArithmeticError('unexpected entropy certificate root fields')
    if value.get('schema')!='oph.shrinking_modular_entropy.v1':
        raise ArithmeticError('wrong entropy certificate schema')
    if value.get('new_checks_pass') is not True:
        raise ArithmeticError('new direct checks did not pass')
    if value.get('full_repository_replay_performed') is not False:
        raise ArithmeticError('unperformed repository replay cannot be promoted')
    expected={'finite_attachment':'FA101','theorem':theorem(),'bounds':bounds.certificate(),
              'algebra':algebra.certificate(),'controls':controls.certificate(),
              'custody':custody(),'source_bindings':source_bindings()}
    for key, item in expected.items():
        if canonical(value.get(key))!=canonical(item):
            raise ArithmeticError('entropy certificate differs: '+key)


def build():
    pins=source_bindings()
    value={'schema':'oph.shrinking_modular_entropy.v1','finite_attachment':'FA101','new_checks_pass':True,
           'full_repository_replay_performed':False,'theorem':theorem(),
           'bounds':bounds.certificate(),'algebra':algebra.certificate(),
           'controls':controls.certificate(),'custody':custody(),'source_bindings':pins}
    if pins!=source_bindings():
        raise ArithmeticError('source changed during entropy replay')
    validate(value)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args();value=build()
    if args.write:
        args.write.write_text(canonical(value))
    else:
        retained=strict_json.load(args.check);validate(retained)
        if canonical(retained)!=canonical(value):
            raise ArithmeticError('cold entropy replay differs')
    print(json.dumps({'new_checks_pass':True,'full_repository_replay_performed':False}))


if __name__=='__main__':main()
