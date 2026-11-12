"""Keep calibrated readouts, native calculations and changed laws distinct."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
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
    ('T11', 'lepton_clock_energy', 'Complete native-event timeline and energy-unit freedom', 'native_clock_with_coherent_lepton_services'),
    ('T12', 'lepton_archive_compression', 'Minimal lossless quantum archives and native compressors', 'native_coherent_lepton_cycle'),
    ('T13', 'quantum_sector_measure', 'Quantum electric sectors and complete stationary family', 'native_quantum_Gauss_link'),
    ('T14', 'finite_particle_spectroscopy', 'Finite-time energy detector and conditional GeV intervals', 'calibrated_quadratic'),
    ('T15', 'neutrino_mass_observables', 'Conditional mass norms and Majorana-phase readout', 'calibrated_Weinberg_kernel'),
    ('T16', 'native_flavour_transport', 'Coherent family transport and complete CP contrast', 'declared_native_family_transport'),
    ('T17', 'quantum_thermal_relaxation', 'Complete quantum bath and source relaxation', 'calibrated_quadratic_with_declared_bath'),
    ('T18', 'expanding_quadratic_particles', 'Particle production under a declared expansion', 'calibrated_quadratic_with_declared_expansion'),
    ('T19', 'native_bath_mixing', 'Certified component gaps and cooling lag', 'declared_bath_extension'),
    ('T20', 'paired_hadron_quark_scale', 'Paired potential scales and bare Wilson quark coefficients', 'archived_quenched_Wilson'),
)
JOIN = 'thermal_particle_records'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(name: str, root: Path = ROOT) -> tuple[dict, dict]:
    path = root / 'model' / name / 'results.json'
    result = json.loads(path.read_text())
    if result.get('schema') != f'oph.{name}.v1':
        raise ValueError(f'{name}: schema mismatch')
    structural_schemas = {'frequency_yukawa_vertex', 'native_thermal_completion', 'quantum_sector_measure'}
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
    if name == 'quantum_sector_measure':
        history = result.get('complete_histories', [])
        fixed = result.get('stationarity', {})
        if (result.get('physical_promotions') != []
                or result.get('closure', {}).get('scalar_sector_adjoint_closure') is not False
                or len(result.get('closure', {}).get('family_kernels', [])) != 5
                or fixed.get('complex_fixed_operator_dimension') != 7
                or fixed.get('family_weights_selected_by_channel') is not False
                or len(fixed.get('pair_fixed_equations', [])) != 25
                or len(history) != 5 or any(len(row.get('complete_records', {})) != 8
                    or sum(map(Fraction, row['complete_records'].values())) != 1 for row in history)):
            raise ValueError(f'{name}: incomplete quantum sector certificate')
    if name == 'finite_particle_spectroscopy':
        confidence = result.get('uniform_confidence', {})
        if (confidence.get('coverage_domain') != 'normalized one-excitation inputs at detector entrance'
                or confidence.get('all_source_unconditional_65_over72_coverage_claimed') is not False
                or confidence.get('vacuum_many_excitation_outcomes_supply_mass_interval') is not False
                or confidence.get('all_input_bad_phase_effect_bound') != 'F_bad <= (7/72) P_one'
                or not confidence.get('complete_original_probe_bound')):
            raise ValueError(f'{name}: confidence domain or complete-probe pullback lost')
    if name == 'expanding_quadratic_particles':
        protocols = result.get('protocols', [])
        if len(protocols) != 7:
            raise ValueError(f'{name}: incomplete expansion protocols')
        for protocol in protocols:
            preparations = protocol.get('preparations', [])
            if len(preparations) != 2:
                raise ValueError(f'{name}: incomplete source preparations')
            for preparation in preparations:
                law = preparation['complete_number_readout']
                probabilities = law['all181_total_number_probabilities']
                mean = law['expected_total_number']
                species = law['complete12_species_number_laws']
                for row in species:
                    values = row.get('probabilities', [])
                    species_mean = row['mean_number']
                    if (not 2 <= len(values) <= 181
                            or any(not math.isfinite(p) or not 0 <= p <= 1 for p in values)
                            or abs(math.fsum(values)-1) > 1e-10
                            or not math.isfinite(species_mean) or not 0 <= species_mean <= len(values)-1
                            or abs(species_mean-math.fsum(n*p for n, p in enumerate(values))) > 1e-9):
                        raise ValueError(f'{name}: invalid normalized species count law')
                if (len(probabilities) != 181 or len(species) != 12
                        or any(not math.isfinite(p) or not 0 <= p <= 1 for p in probabilities)
                        or abs(math.fsum(probabilities)-1) > 1e-10
                        or not math.isfinite(mean) or not 0 <= mean <= 180
                        or abs(mean-math.fsum(n*p for n, p in enumerate(probabilities))) > 1e-9
                        or abs(mean-math.fsum(row['mean_number'] for row in species)) > 1e-9):
                    raise ValueError(f'{name}: invalid normalized source count law')
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
    require_false(receipts['lepton_clock_energy']['clock_energy_result']['absolute_GeV_or_SI_scale_derived'], 'absolute clock scale')
    require_false(receipts['lepton_archive_compression']['scope']['original_record_histories_erased_or_coarsened'], 'archive erasure')
    for name in ('neutrino_mass_observables', 'native_flavour_transport'):
        if receipts[name].get('physical_promotions') != []:
            raise ValueError(f'Unsupported physical promotion in {name} receipt')
    if receipts['neutrino_mass_observables']['mass_observables']['all_masses_are_conditional_calibration'] is not True:
        raise ValueError('Neutrino calibration ancestry lost')
    require_false(receipts['quantum_thermal_relaxation']['scope']['physical_clock_or_Universe_abundances_derived'], 'quantum bath cosmology')
    require_false(receipts['quantum_thermal_relaxation']['scope']['complete_archive_claimed_to_forget_initial_information'], 'quantum archive mixing')
    require_false(receipts['finite_particle_spectroscopy']['scope']['planned_large_registers_executed'], 'planned precision execution')
    require_false(receipts['finite_particle_spectroscopy']['scope']['neutrino_species_resolved_by_executed_five_qubit_pointer'], 'unresolved neutrino detector')
    adaptive = receipts['finite_particle_spectroscopy']['adaptive_pointer']
    if (adaptive['exact_all_input_instrument_equivalence'] is not True
            or adaptive['reusable_quantum_qubits'] != 1
            or adaptive['retained_original_classical_bits'] != 5
            or adaptive['complete_source_reference_output_preserved'] is not True):
        raise ValueError('Adaptive phase-pointer equivalence or record retention lost')
    require_false(adaptive['classical_history_erased'], 'adaptive record erasure')
    require_false(adaptive['large_n_adaptive_precision_plans_executed'], 'adaptive precision planning')
    require_false(receipts['paired_hadron_quark_scale']['scope']['identified_up_down_or_strange_mass'], 'bare quark identification')
    require_false(receipts['paired_hadron_quark_scale']['scope']['renormalized_MSbar_mass'], 'bare quark scheme')
    require_false(receipts['native_bath_mixing']['scope']['global_unique_equilibrium_claimed'], 'component weights')
    require_false(receipts['expanding_quadratic_particles']['scope']['Einstein_backreaction_or_physical_cosmic_time_derived'], 'supplied expansion')
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
    paired = receipts['paired_hadron_quark_scale']
    bare = []
    for scale in ('rho', 'string', 'r0_benchmark'):
        key = f'bare_Wilson_GeV/{scale}/ratio_tuned'
        primary = paired['conditional_coordinates']['primary']
        bare.append({'scale': scale, 'bare_Wilson_GeV': primary['central'][key],
                     'conditional_statistical_error_GeV': primary['statistical'][key],
                     'separate_fit_allowance_GeV': paired['crossed_fit_sensitivity'][key],
                     'identified_physical_quark_mass': False})
    for path in HERE.iterdir():
        if path.suffix in {'.py', '.md'}:
            bindings[str(path.relative_to(ROOT))] = digest(path)
    return {
        'schema': 'oph.particle_cosmology_results.v1',
        'source_bindings': dict(sorted(bindings.items())),
        'all_checks_pass': True, 'tasks': rows[:len(TASKS)], 'same_source_join': rows[len(TASKS)],
        'conditional_GeV_targets': targets,
        'compiled_lepton_GeV_readout': joint['compiled_lepton_joint']['families'],
        'compiler_resource_examples': receipts['lepton_compiler_compression']['precision_resource_examples'],
        'hadron_GeV_central': hadron['gev']['primary']['central'],
        'hadron_GeV_statistical': hadron['gev']['primary']['statistical'],
        'hadron_GeV_fit_sensitivity': hadron['gev_fit_sensitivity'],
        'finite_vertex_frequency_rows': frequency['frequency_rows'],
        'conditional_neutrino_mass_observables_GeV': receipts['neutrino_mass_observables']['mass_observables'],
        'lepton_archive_resources': receipts['lepton_archive_compression']['resources'],
        'conditional_bare_Wilson_coordinates_GeV': bare,
        'finite_spectroscopy': {'parallel_pointer_qubits': receipts['finite_particle_spectroscopy']['finite_pointer_qubits'],
            'equivalent_adaptive_pointer': adaptive,
            'complete_records': receipts['finite_particle_spectroscopy']['record_count'],
            'ideal_coverage_lower': receipts['finite_particle_spectroscopy']['uniform_confidence']['coverage_lower_exact'],
            'coverage_domain': receipts['finite_particle_spectroscopy']['uniform_confidence']['coverage_domain'],
            'all_input_bad_phase_effect_bound': receipts['finite_particle_spectroscopy']['uniform_confidence']['all_input_bad_phase_effect_bound'],
            'complete_original_probe_bound': receipts['finite_particle_spectroscopy']['uniform_confidence']['complete_original_probe_bound'],
            'vacuum_many_excitation_outcomes_supply_mass_interval': False,
            'all_source_unconditional_65_over72_coverage_claimed': False,
            'conditional_targets_and_separate_precision_plans': receipts['finite_particle_spectroscopy']['conditional_GeV_menu']},
        'source_relations': [
            'Joint identity and finite Gibbs number records share the calibrated quadratic action; the detector consumer supplies the state and instrument transfer.',
            'Native absorption and thermal obstruction use the unchanged closed-cell kernel. The thermal bath changes its law on the same carrier.',
            'The finite-cap Gibbs state and continuum free gas have different state spaces and population laws.',
            'Compiler compression preserves the specified full Fock source and archives. Native service words implement a separately specified coherent cycle.',
            'The local finite-frequency vertex, family-action quotient and quenched Wilson archive retain their separate preparations, clocks and parameter ancestry.',
            'Finite spectroscopy, quantum thermal collisions and prescribed expansion consume the calibrated quadratic source with different declared apparatus and interventions.',
            'Neutrino mass norms consume the actual calibrated Weinberg kernel and weak current; the Majorana-phase envelope retains the unknown relative phase.',
            'Coherent native-family transport changes the down vertex architecture without changing its singular values; it does not select empirical CKM coordinates.',
            'Quantum electric sectors belong to the configured Gauss-link representation; their stationary family weights are not a physical running coupling.',
            'Lepton archive compression preserves every original return and reference. The native-event timeline retains an independent energy-unit calibration freedom.',
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
    print('PASS: twenty finite results, calibrated GeV ancestry, and thermal detector join')
