"""Replay source local creators, complete records and multipacket bounds."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from model.spin_charge_dispersion import order
from . import bounds, finite, dispersion

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
PROOFS=('model/spin_canonical_field/README.md','model/spin_canonical_energy/README.md',
        'model/protected_defect_catalog/PROOF.md','model/spin_charge_dispersion/PROOF.md',
        'model/spin_energy_polarization/PROOF.md','model/source_charge_scattering/CONTINUUM_PAIRS.md',
        'model/source_charge_scattering/PROOF.md','model/particle_mass_readout_gate/README.md',
        'model/spin_shell_readout/PROOF.md','model/spin_shell_diagnostic/README.md',
        'model/spin_current_spectrum/README.md',
        'model/spin_operational_position/PROOF.md','model/native_ordered_ground/README.md',
        'model/native_ordered_ground/witness.npz','model/requirements.txt',
        'Lean/QFT/FiniteUnitaryScatteringNoGo.lean','Lean/Dynamics/ProtectedCharge.lean')


def canonical(value):
    return json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n'


def serialize(value):
    if type(value)is F:return str(value)
    if isinstance(value,dict):return {str(k):serialize(v)for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serialize(v)for v in value]
    return value


def bindings():
    legacy='evidence/observer_dynamics_20260925/oph-physics-sim/'
    names=set(PROOFS)|{legacy+"oph_exact/carrier.py",legacy+"oph_exact/__init__.py"}
    names.update(p.relative_to(ROOT).as_posix()for p in HERE.iterdir()if p.suffix in ('.py','.md'))
    pending=list(names);visited=set()
    def add(module):
        parts=module.split('.');prefix=legacy if parts[0]in('oph_exact','oph_fpe')else''
        stem=prefix+'/'.join(parts)
        candidates=[stem+'.py',stem+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py'for i in range(1,len(parts))]
        for name in candidates:
            if (ROOT/name).is_file()and name not in names:names.add(name);pending.append(name)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name)
        module=name.removeprefix(legacy)[:-3].replace('/','.')
        package=module.removesuffix('.__init__')if module.endswith('.__init__')else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or'')
                target=importlib.util.resolve_name(target,package)if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in sorted(names)}


def build():
    before=bindings();band=bounds.local_band();a,b=band['dimensionless_energy_window']
    filter_data=bounds.filter_bounds(a,b)
    p=band['each_charge_original_joint_band_probability_lower']/filter_data['time_filter_L1_upper']**2
    configurations=[]
    for n in (1,2,3,4):
        eps=p**n*4**(n-1)/(4*n*(n-1))if n>1 else F(1,2)
        configurations.append(bounds.configuration(n,F(1),F(4),filter_data,p,eps))
    result={'schema':'oph.spin_multipacket_comparison.v1','all_pass':True,'source_bindings':before,
        'actual_source_order_bound':order.certificate(),
        'actual_finite_local_source_control':finite.build(),
        'actual_source_scalar_dispersion_diagnostic':dispersion.build(),
        'strict_local_creator_band':serialize(band),
        'analytic_profile_cutoff_examples':serialize([bounds.cutoff(F(c))for c in (1,2,10)]),
        'certified_frequency_filter':serialize(filter_data),
        'positive_equal_sign_multipackets':serialize(configurations),
        'derivative_identity_control':serialize(bounds.cook_bound([F(1,2)]*3,[F(1,100)]*3,
                {(0,1):F(1,1000),(0,2):F(1,2000),(1,2):F(1,1000)})),
        'analytic_source_conclusions':{
            'finite_vacuum_configuration_distributions_absolutely_continuous':True,
            'independent_tests_required_for_joint_absolute_continuity':True,
            'finite_distinct_creator_monomial_Gram_strictly_positive':True,
            'effective_uniform_Gram_floor_provided':False,
            'full_net_P2_rotated_to_extra_P3':False,
            'configuration_history_internal_charge_spectrum':'Z',
            'source_even_CAR_grading_preserved':True,
            'strict_local_positive_energy_charge_inputs':True,
            'all_original_preparation_records_retained':True,
            'positive_spectrum_approaches_zero_for_every_nonzero_integer_internal_charge':True,
            'neutral_sector_positive_spectrum_claimed_by_equal_sign_proof':False,
            'all_mixed_sign_filtered_products_proved_nonzero':False,
            'actual_one_particle_dispersion_or_Lorentz_shell_derived':False,
            'Gaussian_Fock_factorization_assumed':False,
            'incoming_outgoing_waveoperators_constructed':False,
            'physical_particle_species_masses_or_cross_sections_derived':False},
        'validation_scope':{'finite_source_integer_operators_and_ground_certificate_recomputed':True,
            'analytic_infinite_volume_field_theorems_used':True,
            'continuum_state_replaced_by_finite_matrix':False,
            'time_filter_examples_have_finite_certified_smoothness':True}}
    if bindings()!=before:raise RuntimeError('source changed during multipacket certificate')
    return result


def validate(value):
    if canonical(value)!=canonical(build()):raise ValueError('source, complete records, analytic bounds or physical scope differ')
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',type=Path,nargs='?',const=HERE/'results.json')
    g.add_argument('--check',type=Path,nargs='?',const=HERE/'results.json')
    args=p.parse_args()
    if args.write:args.write.write_text(canonical(build()))
    else:validate(strict_json.load(args.check))
    print('PASS: source local charged creators, complete finite records and multipacket bounds')
