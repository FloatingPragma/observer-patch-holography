"""Replay exact source twists, complete coefficient scores and closure controls."""
import argparse
import ast
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

from tools import strict_json
from . import algebra, bounds, polarization, readout, state_bounds

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PARENTS = ('model/protected_spin_limit/README.md',
           'model/protected_spin_limit/algebra.py',
           'model/protected_spin_limit/bounds.py',
           'model/protected_spin_limit/experiment.py',
           'model/protected_spin_limit/results.json',
           'model/spin_goldstone_sector/README.md',
           'model/spin_canonical_field/README.md',
           'model/spin_canonical_field/results.json',
           'model/spin_canonical_energy/README.md',
           'model/spin_canonical_energy/results.json',
           'model/spin_shell_diagnostic/README.md',
           'model/spin_current_spectrum/README.md',
           'model/spin_energy_polarization/PROOF.md',
           'model/spin_charge_dispersion/DIMER_BOUND.md',
           'model/native_spin_observer/README.md',
           'model/spin_operational_position/PROOF.md',
           'model/spin_operational_position/generic.py',
           'model/native_ordered_ground/source.py',
           'model/requirements.txt')


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def bindings():
    legacy = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
    names = set(PARENTS)|{legacy+'oph_exact/carrier.py',legacy+'oph_exact/__init__.py'}
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py','.md'))
    pending = list(names)
    visited = set()
    def add(module):
        parts = module.split('.')
        prefix = legacy if parts[0] in ('oph_exact','oph_fpe') else ''
        candidates = [prefix+'/'.join(parts)+'.py',prefix+'/'.join(parts)+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removeprefix(legacy)[:-3].replace('/','.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target,package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before = bindings()
    value = {'schema':'oph.spin_twist_closure.v1','source_bindings':before,
        'actual_source_frame':algebra.source_frame_certificate(),
        'local_rotation_algebra':algebra.local_certificate(),
        'actual_source_force':algebra.source_force_certificate(),
        'actual_source_reset_identity':algebra.reset_energy_certificate(),
        'scalar_total_energy_controls':[bounds.ground_energy_cube(n) for n in (4,16,64,256)],
        'selected_ground_region_energy_controls':[state_bounds.regional_second_moment(n)
                                                   for n in (4,16,64,256)],
        'actual_source_energy_only_phase_control':state_bounds.source_phase_control(),
        'actual_finite_gap_boundary_controls':[state_bounds.source_gap_control(n) for n in (2,3,4)],
        'vanishing_scalar_polarization_residual':polarization.exact_checks(),
        'actual_product_state_residual_countercontrol':polarization.product_source_control(),
        'scalar_residual_history_bounds':[polarization.compact_bounds(F(1,n),1,1,1)
                                          for n in (16,256,4096)],
        'hidden_force_control':algebra.hidden_force_control(),
        'compact_test_controls':[bounds.compact_twist(F(1,n),1,1,axis)
                                 for n in (4,16,64) for axis in ('P','Q')],
        'bond_phase_controls':[bounds.phase_transport(F(1,n),1,1,1,1)
                               for n in (16,256,4096)],
        'cubic_bond_overlap':bounds.bond_kick(6,F(1,64),1),
        'complete_coefficient_readout':readout.plan((F(1,3),F(-1,5)),2),
        'scope':{
            'same_ordered_spin_source_state_and_clock':True,
            'Q_and_P_single_direction_BCH_remainder_vanishes_in_operator_norm':True,
            'actual_bond_coefficients_have_complete_finite_score_instruments':True,
            'coefficient_phase_transport_through_finite_time_histories':True,
            'coefficient_phases_can_be_adjoined_as_central_source_variables':True,
            'aggregate_bond_energy_limit_is_scalar_without_ergodicity':True,
            'Q_energy_polarization_residual_vanishes_on_actual_finite_histories':True,
            'separate_unbounded_boost_generator_or_commutator_limits_proved':False,
            'selected_state_all_bond_limits_proved_scalar':False,
            'original_minimal_field_space_proved_to_reduce_added_coefficients':False,
            'original_equal_time_canonical_algebra_proved_dynamically_closed':False,
            'first_energy_moment_limit_interchange_proved':False,
            'actual_filtered_force_proved_to_vanish':False,
            'sharp_dispersion_or_Lorentz_covariance_proved':False,
            'physical_particle_assignment':None}}
    if bindings() != before:
        raise RuntimeError('source twist inputs changed during replay')
    return bounds.serialize(value)


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('source twist certificate, complete records or scope differs')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args = parser.parse_args()
    if args.write:
        args.write.write_text(canonical(build()))
    else:
        validate(strict_json.load(args.check))
    print('PASS: exact source twists, coefficient readout and scoped closure controls')
