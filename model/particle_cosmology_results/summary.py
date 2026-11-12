"""Keep calibrated readouts, native calculations and changed laws distinct."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TASKS = (
    ('T1', 'particle_joint_mass_identity', 'Joint particle identity and GeV coefficients', 'calibrated_quadratic'),
    ('T2', 'lepton_compiler_compression', 'Exact compiler compression and complete-output error', 'neutral_lepton_compiler'),
    ('T3', 'lepton_service_compiler', 'Native lepton service words and complete archives', 'native_coherent_lepton_cycle'),
    ('T4', 'native_thermal_admission', 'Stationary classes and thermal obstruction', 'native_closed_cell'),
    ('T5', 'native_freezeout_response', 'Absorption distributions and lifetimes', 'native_closed_cell'),
    ('T6', 'particle_thermal_history', 'Finite Gibbs states and conditional free-gas history', 'calibrated_quadratic_and_declared_continuum'),
    ('T7', 'correlated_hadron_spectrum', 'Correlated hadron fits and paired GeV errors', 'archived_quenched_Wilson'),
    ('T8', 'source_flavour_identifiability', 'Flavour parameter quotient and CP readout', 'native_family_action'),
    ('T9', 'frequency_yukawa_vertex', 'Connected finite-frequency Yukawa vertex', 'local_320_state_unit_source'),
    ('T10', 'native_thermal_completion', 'Explicit reversible bath and cooling laws', 'declared_bath_extension'),
)
JOIN = 'thermal_particle_records'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(name: str, root: Path = ROOT) -> tuple[dict, dict]:
    path = root / 'model' / name / 'results.json'
    result = json.loads(path.read_text())
    if result.get('schema') != f'oph.{name}.v1':
        raise ValueError(f'{name}: schema mismatch')
    structural_schemas = {'frequency_yukawa_vertex', 'native_thermal_completion'}
    if name not in structural_schemas and result.get('all_checks_pass') is not True:
        raise ValueError(f'{name}: failed scientific checks')
    if name == 'frequency_yukawa_vertex':
        pointer = result.get('complete_invariant_pointer', {})
        if (result.get('all_one_point_functions_zero_exact') is not True
                or pointer.get('four_positive_effects_sum_identity_exact') is not True
                or pointer.get('record_count') != 17
                or len(result.get('complete_first_frequency_tensor', [])) != 72
                or not pointer.get('total_current_commutators_exact')
                or not all(x is True for x in pointer['total_current_commutators_exact'])):
            raise ValueError(f'{name}: incomplete final vertex certificate')
    if name == 'native_thermal_completion':
        rows = result.get('temperature_rows', [])
        if len(rows) != 4 or any(
                row.get('detailed_balance_exact') is not True
                or row.get('stationarity_exact') is not True
                or row.get('conjugation_exact') is not True
                or row.get('complete_proposal_records') != 17352 for row in rows):
            raise ValueError(f'{name}: incomplete exact bath certificate')
    bindings = result.get('source_bindings', result.get('input_sha256'))
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError(f'{name}: missing source ancestry')
    for relative, expected in bindings.items():
        target = (root / relative).resolve()
        if not target.is_relative_to(root.resolve()) or digest(target) != expected:
            raise ValueError(f'{name}: changed input {relative}')
    return result, {**bindings, str(path.relative_to(root)): digest(path)}


def require_false(value: object, label: str) -> None:
    if value is not False:
        raise ValueError(f'Unsupported physical promotion: {label}')


def certificate(root: Path = ROOT) -> dict:
    receipts, bindings, rows = {}, {}, []
    for task, name, title, domain in (*TASKS, ('join', JOIN, 'Thermal particle-number records', 'calibrated_quadratic')):
        result, pins = authenticate(name, root)
        for key, value in pins.items():
            if key in bindings and bindings[key] != value:
                raise ValueError(f'Inconsistent shared source: {key}')
            bindings[key] = value
        receipts[name] = result
        rows.append({'task': task, 'package': f'model/{name}', 'title': title,
                     'source_domain': domain, 'receipt_sha256': pins[f'model/{name}/results.json']})
    joint = receipts['particle_joint_mass_identity']
    require_false(joint['scope']['physical_Standard_Model_particle_dictionary_derived'], 'particle dictionary')
    require_false(joint['scope']['calibrated_masses_are_independent_predictions'], 'calibrated masses')
    require_false(receipts['source_flavour_identifiability']['scope']['source_selects_all_flavour_parameters'], 'flavour selection')
    if receipts['source_flavour_identifiability']['physical_premise_promotions']:
        raise ValueError('Unsupported physical promotion in flavour receipt')
    require_false(receipts['frequency_yukawa_vertex']['physical_Yukawa_or_on_shell_mass_identification'], 'finite vertex')
    require_false(receipts['native_thermal_completion']['scope']['unchanged_native_repair_claimed_thermal'], 'changed bath law')
    require_false(receipts['particle_thermal_history']['physical_admission']['physical_abundances_or_freezeout_predicted'], 'cosmological abundances')
    targets = joint['conditional_GeV_targets']
    if len(targets) != 17 or any(row['native_prediction'] for row in targets):
        raise ValueError('Incomplete or misclassified calibrated mass menu')
    # The actual common source, state and detector transfer is proved by this
    # consumer. A common file hash alone is not a physical transfer theorem.
    join_pins = receipts[JOIN].get('source_bindings', receipts[JOIN].get('input_sha256', {}))
    for parent in ('particle_joint_mass_identity', 'particle_thermal_history'):
        path = f'model/{parent}/results.json'
        if join_pins.get(path) != bindings[path]:
            raise ValueError(f'Thermal detector omits its actual parent receipt: {path}')
    frequency = receipts['frequency_yukawa_vertex']
    if frequency['complete_invariant_pointer']['record_count'] != 17:
        raise ValueError('Finite vertex preparation/pointer records incomplete')
    hadron = receipts['correlated_hadron_spectrum']
    for path in HERE.iterdir():
        if path.suffix in {'.py', '.md'}:
            bindings[str(path.relative_to(ROOT))] = digest(path)
    return {
        'schema': 'oph.particle_cosmology_results.v1',
        'source_bindings': dict(sorted(bindings.items())),
        'all_checks_pass': True, 'tasks': rows[:10], 'same_source_join': rows[10],
        'conditional_GeV_targets': targets,
        'compiled_lepton_GeV_readout': joint['compiled_lepton_joint']['families'],
        'compiler_resource_examples': receipts['lepton_compiler_compression']['precision_resource_examples'],
        'hadron_GeV_central': hadron['gev']['primary']['central'],
        'hadron_GeV_statistical': hadron['gev']['primary']['statistical'],
        'hadron_GeV_fit_sensitivity': hadron['gev_fit_sensitivity'],
        'finite_vertex_frequency_rows': frequency['frequency_rows'],
        'source_relations': [
            'Joint identity and finite Gibbs number records share the calibrated quadratic action; the detector consumer supplies the state and instrument transfer.',
            'Native absorption and thermal obstruction use the unchanged closed-cell kernel. The thermal bath changes its law on the same carrier.',
            'The finite-cap Gibbs state and continuum free gas have different state spaces and population laws.',
            'Compiler compression preserves the specified full Fock source and archives. Native service words implement a separately specified coherent cycle.',
            'The local finite-frequency vertex, family-action quotient and quenched Wilson archive retain their separate preparations, clocks and parameter ancestry.',
        ],
        'physical_premise_promotions': [],
        'scope': 'Authenticated finite results and explicit calibrated readouts; no theorem identifying the full observed particle spectrum, physical Yukawas, or early-Universe abundances.',
    }


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    path = HERE / 'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    elif json.loads(path.read_text()) != result:
        raise SystemExit('Particle/cosmology overview differs from its authenticated inputs')
    print('PASS: ten finite results, calibrated GeV ancestry, and thermal detector join')
