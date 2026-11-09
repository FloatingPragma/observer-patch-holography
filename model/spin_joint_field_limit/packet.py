"""Replay the prepared positive-field family and its source concentration bounds."""
import argparse
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

from tools import strict_json
from . import algebra, bounds, clock_comparison, construction, force, joint_symmetry, ward, uniform_mode, uniform_spectral
from . import arctan_bands, directed_bands, quantum_force_slack, speed_concentration

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=HERE/'results.json'
PARENTS=('model/spin_twist_closure/FINITE_FIELD_GAP.md',
         'model/spin_twist_closure/PROOF.md',
         'model/spin_twist_closure/results.json',
         'model/spin_fixed_field_locality/PROOF.md',
         'model/spin_fixed_field_locality/results.json',
         'model/spin_canonical_field/README.md',
         'model/spin_canonical_field/results.json',
         'model/spin_canonical_energy/README.md',
         'model/spin_relaxed_response/README.md',
         'model/spin_shell_diagnostic/README.md',
         'model/native_joint_spin_spectroscopy/README.md',
         'model/spin_current_spectrum/README.md',
         'model/spin_dynamical_memory/RESOLVENT_CLOSURE.md',
         'model/spin_goldstone_sector/README.md',
         'model/spin_goldstone_sector/algebra.py',
         'model/spin_charge_dispersion/PROOF.md',
         'model/spin_charge_dispersion/DIMER_BOUND.md',
         'model/spin_charge_dispersion/order.py',
         'model/spin_charge_dispersion/results.json',
         'Lean/SourceLaw/SpinMarginalForce.lean',
         'model/native_spin_observer/README.md',
         'model/requirements.txt')


def serialize(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {key:serialize(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(item) for item in value]
    return value


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def bindings():
    legacy='evidence/observer_dynamics_20260925/oph-physics-sim/'
    names=set(PARENTS)|{legacy+'oph_exact/carrier.py',legacy+'oph_exact/__init__.py'}
    names.update((uniform_spectral.WITNESS,uniform_spectral.zero_field.WITNESS))
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py','.md'))
    pending=list(names);visited=set()
    def add(module):
        parts=module.split('.')
        prefix=legacy if parts[0] in ('oph_exact','oph_fpe') else ''
        candidates=[prefix+'/'.join(parts)+'.py',prefix+'/'.join(parts)+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name);pending.append(name)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name)
        module=name.removeprefix(legacy)[:-3].replace('/','.')
        package=module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(target,package) if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before=bindings()
    value={'schema':'oph.spin_joint_field_limit.v1','source_bindings':before,
        'actual_local_source_algebra':algebra.local_certificate(),
        'finite_selected_ground_controls':algebra.finite_ground_controls(),
        'joint_family':[bounds.family(n) for n in (2,8,32,128)],
        'order_concentration':[bounds.order_concentration(n,1,1) for n in (2,8,32,128)],
        'bond_concentration':[bounds.bond_concentration(n,1,1) for n in (2,8,32,128)],
        'bond_history_transport':[bounds.bond_transport(n,1,1,1,1) for n in (32,128,512)],
        'transverse_energy':[bounds.field_energy(n,kind,1,1)
                              for n in (2,8,32) for kind in ('Q','P')],
        'complete_local_preparation':[bounds.local_preparation(n,1,1,F(1,10)) for n in (2,3)],
        'selection_clock_comparison_and_native_clock':clock_comparison.build(),
        'complete_prepared_local_history':construction.probe_history(2,
            (('X',(0,0,0)),('Z',(1,0,0))),(F(1,32),F(1,16),F()),F(1,10)),
        'second_canonical_pair_and_full_internal_circle':joint_symmetry.certificate(),
        'same_clock_filtered_current_and_complete_force_readout':force.certificate(),
        'source_Ward_susceptibilities_and_force_equations':ward.certificate(),
        'source_quantum_correction_to_mean_spin_force':quantum_force_slack.certificate(),
        'uncut_speed_concentration_from_source_susceptibility':speed_concentration.certificate(),
        'uniform_mode_original_weight_and_concentration':uniform_mode.certificate(),
        'complete_finite_uniform_spectral_readout':uniform_spectral.certificate(),
        'bounded_joint_charged_energy_momentum_bands':arctan_bands.certificate(),
        'directed_joint_charged_bands':directed_bands.certificate(),
        'scope':{'source_field_sequence':'p=n^-2, B=J*n^-5, n>=2',
            'sequence_is_a_chosen_source_preparation_prescription':True,
            'sequence_forced_by_unread_repair':False,
            'state':'the unique positive-B infinite-volume source ground at each n',
            'clock':'K_B/(Jp) at each n',
            'finite_positive_field_preparation_and_local_accuracy_computable':True,
            'large_preparation_circuits_executed':False,
            'complete_native_failure_records_retained':True,
            'positive_field_order_and_bond_variance_bounds_quantitative':True,
            'direct_single_pair_regular_positive_energy_history_subsequence':True,
            'second_canonical_pair_on_same_joint_history_subsequence':True,
            'full_internal_U1_on_both_canonical_pairs':True,
            'fixed_scale_complete_prepared_H_B_probe_history_constructed':True,
            'same_clock_ground_vector_filtered_current_O_p_squared_bound':True,
            'complete_finite_original_weight_force_defect_readout_constructed':True,
            'source_Ward_reduction_of_two_force_form_defects':True,
            'strict_source_bond_minus_squared_order_margin_proved':True,
            'unfiltered_microscopic_mean_spin_force_factorization_excluded':True,
            'quantum_force_margin_excludes_corrected_infrared_shell':False,
            'uncut_original_weight_speed_tail_from_one_source_susceptibility':True,
            'complete_finite_field_inverse_score_error_bound':True,
            'actual_source_speed_concentration_defect_vanishes':False,
            'original_weight_energy_tube_bound_from_filtered_Ward_data':True,
            'uniform_mode_weighted_centre_and_variance_concentrate_at_zero':True,
            'finite_source_uniform_order_floor_computable':True,
            'finite_spectral_table_all_original_records_certified':True,
            'bounded_charged_joint_band_original_weight_uniform_in_momentum_scale':True,
            'finite_source_monotone_creator_commutator_error_controlled':True,
            'uniform_mode_nonvanishing_original_continuum_residue_proved':False,
            'fixed_cube_spectrum_identified_with_growing_volume_field':False,
            'positive_full_pair_slack_excludes_sharp_subbranches':False,
            'filtered_linear_force_residual_proved_o_p_squared':False,
            'mean_magnetization_or_bond_convergence_rate_proved':False,
            'joint_continuum_subsequence_uniqueness_proved':False,
            'equality_to_ordered_zero_field_state_clock_limit_proved':False,
            'closed_linear_field_equation_proved':False,
            'sharp_particle_mass_shell_proved':False,
            'Lorentz_covariance_proved':False,'physical_particle_assignment':None}}
    if bindings()!=before:raise RuntimeError('source inputs changed during replay')
    return serialize(value)


def validate(value):
    if canonical(value)!=canonical(build()):
        raise ValueError('joint family source, bounds, records or scope differs')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args()
    if args.write:args.write.write_text(canonical(build()))
    else:validate(strict_json.load(args.check))
    print('PASS: prepared positive-field source, concentration and scoped histories')
