"""Current source bindings and exact finite checks for energy polarization."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import algebra, bounds, samples

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent


def serialize(value):
    if type(value)is F:return str(value)
    if isinstance(value,dict):return {k:serialize(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v)for v in value]
    return value


def bindings():
    paths=set(HERE.glob('*.py'));queue=list(paths)
    def include(path):
        if path.is_file()and path not in paths:paths.add(path);queue.append(path)
    while queue:
        path=queue.pop();package='.'.join(path.relative_to(ROOT).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names=[]
            if isinstance(node,ast.Import):names=[a.name for a in node.names]
            elif isinstance(node,ast.ImportFrom):
                name='.'*node.level+(node.module or '')
                base=importlib.util.resolve_name(name,package)if node.level else name
                names=[base]+[base+'.'+a.name for a in node.names]
            for name in names:
                parts=name.split('.');include(ROOT.joinpath(*parts).with_suffix('.py'))
                for j in range(1,len(parts)+1):include(ROOT.joinpath(*parts[:j])/'__init__.py')
    paths.update(ROOT/p for p in (
        'model/spin_canonical_field/README.md','model/spin_canonical_energy/README.md',
        'model/spin_current_spectrum/README.md','model/spin_boost_stress/PROOF.md',
        'model/spin_charge_dispersion/PROOF.md','model/protected_defect_catalog/PROOF.md',
        'model/source_energy_effort/PROOF.md','model/native_spin_clock_join/energy.py',
        'model/native_ordered_ground/source.py','model/protected_spin_limit/README.md',
        'Lean/QFT/CoarseEnergyWard.lean','Lean/Dynamics/WardLimitManifest.lean',
        'paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex'))
    paths.update((HERE/'README.md',HERE/'PROOF.md'))
    return {p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(paths)}


def claims():
    return {'same_ordered_source_state_clock_retained':True,
            'compact_preparation_increment_cutoff_independence':True,
            'source_weighted_twist_limit_with_explicit_errors':True,
            'actual_translation_energy_affine_law':True,
            'exact_directional_current_scalar_field_remainder':True,
            'Schwartz_weighted_twist_tail_control':True,
            'complete_charged_preparation_truncation_bound':True,
            'annular_test_assumed_compact_spatial_support':False,
            'coherent_charged_energy_replaced_by_mixture':False,
            'source_kappa_numerical_value_selected':False,
            'weighted_source_energy_identified_as_continuum_boost':False,
            'current_remainder_proved_to_vanish':False,
            'vacuum_and_full_history_boost_covariance_proved':False,
            'first_energy_escape_removed':False,
            'Lorentz_stress_mass_shell_or_Einstein_dynamics_proved':False,
            'growing_native_tape_executed':False}


def build():
    before=bindings()
    value={'schema':'oph.spin_energy_polarization.v1','all_pass':True,
        'source_bindings':before,'exact_source_algebra':algebra.exact_checks(),
        'translated_rational_source_meshes':serialize([samples.control(ch,n)for ch in ('Q','P')for n in (2,4,8,16)]),
        'general_mesh_bounds_per_kappa':serialize([bounds.mesh(ch,F(1,2**n),F(2),F(1),F(8),F(64),F(512))for ch in ('Q','P')for n in (2,4,8)]),
        'symmetric_Weyl_polarization_controls':serialize([bounds.mixed_kick(F(1,2**n),F(2),F(1),F(8),F(64),F(1),F(8),F(64),F(512))for n in (2,4,8)]),
        'weighted_Schwartz_tail_bounds_per_kappa':serialize([bounds.schwartz_tail(ch,F(l),6,F(1))for ch in ('Q','P')for l in (2,4,16)]),
        'complete_charge_truncation_examples':serialize([bounds.charged_truncation(F(l),6,F(1))for l in (2,4,16)]),
        'scalar_action_residual_bound':serialize(bounds.scalar_field_remainder(F(1,16),F(2),F(8))),
        'finite_complete_score_controls':serialize([bounds.finite_readout(a,F(-3),F(1,1000),F(1,1000))for a in ((),(F(1,4),F(-1,8)))]),
        'claims':claims(),
        'inherited_producer_scope':{'source_twist_algebra_recomputed':True,
            'analytic_ground_and_field_premises_reused_explicitly':True,
            'all_prior_native_experiments_replayed':False}}
    if before!=bindings():raise ValueError('source changed during polarization checks')
    return value


def validate(value):
    expected=build()
    if json.dumps(value,sort_keys=True,allow_nan=False)!=json.dumps(expected,sort_keys=True,allow_nan=False):
        raise ValueError('source, exact control or physical scope differs')
    return expected


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',nargs='?',type=Path,const=HERE/'results.json')
    g.add_argument('--check',nargs='?',type=Path,const=HERE/'results.json')
    a=p.parse_args()
    if a.write:a.write.write_text(json.dumps(build(),indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:validate(strict_json.load(a.check))
    print('PASS: source energy polarization, exact current remainder and complete cutoff bounds')


if __name__=='__main__':main()
