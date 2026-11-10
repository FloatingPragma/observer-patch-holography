"""A compact consumer of independently calculated candidate identity tests."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from model.particle_common_contract import contract
from . import catalogue, paths, steps

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
INPUTS = {
    'dictionary': 'model/standard_model_source_map/candidates/correspondence.json',
    'charged': 'model/particle_candidate_charged/results.json',
    'composites': 'model/particle_candidate_composites/results.json',
    'bosons': 'model/particle_candidate_bosons/results.json',
    'charged_propagation': 'model/particle_candidate_charged_propagation/results.json',
    'composite_spectrum': 'model/particle_candidate_composite_spectrum/results.json',
    'scalar_response': 'model/particle_candidate_scalar_gauge_response/results.json',
    'endpoint_threshold': 'model/particle_candidate_endpoint_threshold/results.json',
    'neutral_dispersion': 'model/particle_candidate_neutral_dispersion/results.json',
    'scalar_spectrum': 'model/particle_candidate_scalar_spectrum/results.json',
}
SCHEMAS = {
    'dictionary': 'oph.particle_candidate_correspondence.v1',
    'charged': 'oph.particle_candidate_charged_certificate.v1',
    'composites': 'oph.particle_candidate_composites_certificate.v1',
    'bosons': 'oph.particle_candidate_bosons_certificate.v1',
    'charged_propagation': 'oph.particle_candidate_charged_propagation.v1',
    'composite_spectrum': 'oph.particle_candidate_composite_spectrum.v1',
    'scalar_response': 'oph.particle_candidate_scalar_gauge_response.v1',
    'endpoint_threshold': 'oph.particle_candidate_endpoint_threshold.v1',
    'neutral_dispersion': 'oph.particle_candidate_neutral_dispersion_certificate.v1',
    'scalar_spectrum': 'oph.particle_candidate_scalar_spectrum.v1',
}
SOURCE_INPUTS = tuple('model/particle_candidate_confirmation/' + name for name in
                     ('__init__.py', '__main__.py', 'model.py', 'test_model.py',
                      'PROOF.md', 'README.md', 'steps.py', 'catalogue.py', 'paths.py')) + tuple(
                          row['path'] for row in paths.REUSE) + (
                              'model/common_spin_protection/README.md',
                              'model/cube_eta_protection/README.md')
SOURCE_TESTS = {
    'charged_leptons': 'A Higgs-conjugate composite gives a charged weak/color-singlet dressing with positive original dressed records at three endpoint separations and an admitted canonical energy/momentum kernel; its populated particle branch and flavor are unproved.',
    'up_type_quarks': 'Exact native quark representations have no onsite generation multiplicity; a common bare Yukawa attachment is constrained by native charges.',
    'down_type_quarks': 'Exact native quark representations have no onsite generation multiplicity; a common bare Yukawa attachment is constrained by native charges.',
    'active_neutrinos': 'The native weak doublet differs from the tested neutral singlet trion; a direct current-preserving identification fails.',
    'photon': 'Finite polar and axial loop channels couple; the isolated A1u line is dark to vector probes. The physical field dictionary and masslessness are unproved.',
    'gluon': 'Eight native color directions exist; an undressed local adjoint current has zero physical Gauss compression.',
    'W_plus': 'Native charged weak directions exist; neutral loop multiplets do not identify a charged W field.',
    'W_minus': 'Native charged weak directions exist; neutral loop multiplets do not identify a charged W field.',
    'Z': 'A neutral source-vector channel couples; its weak/hypercharge mixture and physical Z attachment are unproved.',
    'Higgs': 'A neutral radial scalar made from the Higgs-conjugate composite has uniform elastic detection and a separate resolved positive original energy-transfer window; a surviving inelastic Higgs branch and common physical phase are unproved.',
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_result_digests(value):
    if isinstance(value, dict):
        if 'result_sha256' in value and isinstance(value.get('result'), dict):
            if contract.digest(value['result']) != value['result_sha256']:
                raise ValueError('candidate consequence digest differs')
        for child in value.values():
            validate_result_digests(child)
    elif isinstance(value, list):
        for child in value:
            validate_result_digests(child)


def read_inputs():
    values, pins, cache = {}, {}, {}
    for role, name in INPUTS.items():
        path = ROOT / name
        value = json.loads(path.read_text())
        if value.get('schema') != SCHEMAS[role]:
            raise ValueError('candidate producer schema differs: ' + role)
        ancestry = value.get('source_sha256')
        if not ancestry:
            raise ValueError('missing candidate producer custody: ' + name)
        for relative, expected in ancestry.items():
            source = (ROOT / relative).resolve()
            if not source.is_relative_to(ROOT) or not source.is_file():
                raise ValueError('invalid candidate producer path: ' + relative)
            if relative not in cache:
                cache[relative] = sha(source)
            actual = cache[relative]
            if actual != expected:
                raise ValueError('stale candidate producer: ' + relative)
            if relative in pins and pins[relative] != expected:
                raise ValueError('conflicting candidate producer: ' + relative)
            pins[relative] = expected
        validate_result_digests(value)
        if role != 'dictionary':
            validate_primary_binding(value)
        else:
            from model.standard_model_source_map.candidates.model import validate as validate_dictionary
            validate_dictionary(value)
        values[role] = value
        pins[name] = sha(path)
    return values, pins


@lru_cache(maxsize=1)
def expected_bindings():
    members = [contract.member(n) for n in (0, 1, 2, 8)]
    return {m['member_id']: {'member_id': m['member_id'], **m['ids']} for m in members}


def validate_primary_binding(value):
    expected_by_id = expected_bindings()
    bounds = value.get('members', [value.get('selected_member_bound')])
    observed = []
    for bound in bounds:
        result = bound['result']
        expected = expected_by_id.get(result.get('member_id'))
        if expected is None or any(result.get(key) != item for key, item in expected.items()):
            raise ValueError('candidate calculation changed source, state, reference or clock')
        if any(result.get(key) for key in ('physical_particle_map', 'physical_mass_table_GeV',
                                           'physical_species_confirmed', 'physical_species_map', 'mass_predictions_GeV')):
            raise ValueError('candidate calculation requires physical interpretation review')
        observed.append(result['member_id'])
    if next(iter(expected_by_id)) not in observed or len(observed) != len(set(observed)):
        raise ValueError('candidate calculation lacks a distinct member-zero binding')


def target_family(name):
    if name in ('electron', 'muon', 'tau'):
        return 'charged_leptons'
    if name in ('up', 'charm', 'top'):
        return 'up_type_quarks'
    if name in ('down', 'strange', 'bottom'):
        return 'down_type_quarks'
    if name.endswith('_neutrino'):
        return 'active_neutrinos'
    return name


def supporting_receipts(name):
    if target_family(name) == 'active_neutrinos':
        return [INPUTS['composites'], INPUTS['composite_spectrum'], INPUTS['neutral_dispersion']]
    if target_family(name) == 'charged_leptons':
        return [INPUTS['charged_propagation'], INPUTS['endpoint_threshold']]
    if name == 'Higgs':
        return [INPUTS['charged'], INPUTS['scalar_response'], INPUTS['scalar_spectrum']]
    return []


def calculated_summary(values):
    charged = next(b['result'] for b in values['charged']['members']
                   if b['result']['member_id'] == next(iter(expected_bindings())))
    composites = values['composites']['selected_member_bound']['result']
    bosons = values['bosons']['selected_member_bound']['result']
    selected = lambda role: next(b['result'] for b in values[role]['members']
                                if b['result']['member_id'] == next(iter(expected_bindings())))
    propagation = selected('charged_propagation')
    spectrum = values['composite_spectrum']['selected_member_bound']['result']
    scalar = selected('scalar_response')
    endpoints = values['endpoint_threshold']['selected_member_bound']['result']
    neutral = selected('neutral_dispersion')
    inelastic = selected('scalar_spectrum')
    return {
        'charged': {
            'local_coefficient_commutant_dimension': charged['native_coefficient_control']['gauge_and_each_Wilson_coefficient_commutant_dimension'],
            'one_sharp_H_dagger_for_all_bare_Yukawas': charged['native_number_Ward_control']['one_sharp_native_number_H_dagger_for_all_three'],
            'scalar_and_dressed_operator': charged['scalar_and_dressed_odd_control'],
            'original_coherent_bare_dressed_join': charged['original_coherent_bare_dressed_join'],
            'original_compensated_dressed_energy_window': charged['original_compensated_dressed_energy_window'],
        },
        'composites': {
            'native_counter_coefficients': composites['complete_detector']['S_native_coefficients'],
            'cyclic_radix': composites['complete_detector']['cyclic_radix'],
            'signed_family_increments': {r['family'] + ('+' if r['sign'] > 0 else '-'): r['S_increment']
                                         for r in composites['candidate_rows']},
            'multiple_excitation_alias': composites['actual_multiple_excitation_alias'],
        },
        'bosons': {
            'multiplet_dimensions': bosons['actual_loop_symmetry']['multiplet_dimensions'],
            'complete_menu': bosons['actual_complete_loop_menu'],
            'full_source_spectral_bounds': bosons['full_source_interacting_multiplet_bounds'],
            'original_row_weight': bosons['each_row_original_reported_weight_strict_lower_factored'],
            'forbidden_record_categories': bosons['forbidden_original_record_categories'],
        },
        'charged_propagation': {
            'filtered_original_window': propagation['filtered_original_window'],
            'charged_joint_momentum_plan': propagation['charged_joint_momentum_plan'],
        },
        'composite_spectrum': {
            'component_count': len(spectrum['component_rows']),
            'common_original_spatial_energy_contrast': spectrum['common_original_spatial_energy_contrast'],
            'complete_instrument': spectrum['complete_instrument'],
            'clock_histogram_executed': spectrum['complete_final_energy_clock']['clock_histogram_executed'],
            'original_odd_neighboring_two_point_m1_value_or_sign_computed':
                spectrum['original_odd_neighboring_two_point_m1_value_or_sign_computed'],
        },
        'endpoint_threshold': {
            'separations': [{key: row[key] for key in (
                'geometry', 'onsite_electric_commutator_center', 'centered_second_moment_upper',
                'reported_transfer_interval', 'original_filtered_recorded_simple_strict_lower',
                'complete_probe_outcomes', 'finite_pair_interval_is_single_endpoint_mass')}
                for row in endpoints['separations']],
            'pair_only_mass_ambiguity': endpoints['pair_only_mass_ambiguity'],
            'charged_asymptotic_threshold_proved': endpoints['charged_asymptotic_threshold_proved'],
        },
        'neutral_dispersion': {
            'native_increment': neutral['native_increment'],
            'channels': [{key: row[key] for key in (
                'probe_irrep', 'dimension', 'centroid_strict_interval',
                'simple_centroid_span_strict_upper', 'within_multiplet_centroid_ratio_strict_upper',
                'neighbor_first_moment_value_or_sign_computed', 'centroids_are_isolated_particle_energies')}
                for row in neutral['channels']],
            'finite_Walsh_labels_are_conserved_operational_momentum':
                neutral['finite_Walsh_labels_are_conserved_operational_momentum'],
        },
        'scalar_spectrum': {key: inelastic[key] for key in (
            'source_spacing', 'source_width', 'source_density_floor',
            'original_ideal_positive_transfer_interval', 'original_ideal_positive_transfer_weight_lower',
            'reported_positive_transfer_interval', 'original_reported_positive_transfer_weight_lower',
            'uniform_elastic_weight_lower', 'member_zero_simple_elastic_strict_lower',
            'positive_cutoff_is_capacity_independent_in_source_clock',
            'uniform_compact_positive_energy_weight_proved', 'isolated_scalar_pole_or_mass_gap_proved')},
        'scalar_response': {
            'local_radial_spectrum': scalar['local_radial_control']['spectrum'],
            'uniform_original_total_success_lower': scalar['complete_radial_records']['uniform_original_total_success_lower'],
            'member_zero_simple_strict_lower': scalar['complete_radial_records']['member_zero_simple_strict_lower'],
            'uniform_positive_energy_success_lower_proved': scalar['complete_radial_records']['uniform_positive_energy_success_lower_proved'],
            'original_radial_dynamics': scalar['original_radial_dynamics'],
        },
    }


def build():
    values, pins = read_inputs()
    dictionary = values['dictionary']
    if dictionary['physical_particle_map'] or dictionary['mass_predictions_GeV']:
        raise ValueError('candidate input requires a physical interpretation review')
    rows = []
    for target in dictionary['target_rows']:
        name = target['target']
        role = 'charged' if target['source_component_ids'][0] in {
            'charged_lepton_left', 'neutrino_left', 'up_left', 'down_left'} else 'bosons'
        rows.append({
            'target': name,
            'family': target_family(name),
            'source_component_ids': target['source_component_ids'],
            'actual_primary_excitation_ids': target['actual_primary_excitation_ids'],
            'test_receipt': INPUTS[role],
            'supporting_receipts': supporting_receipts(name),
            'source_test_summary': SOURCE_TESTS[target_family(name)],
            'confirmation_steps': steps.candidate_steps(name),
            'scope': 'source representation, conserved-charge and detector selection tests',
            'physical_identity_status': 'unresolved',
            'source_flavor_index': None,
            'source_derived_mass_GeV': None,
        })
    result = {
        'schema': 'oph.particle_candidate_confirmation.v1',
        'target_rows': rows,
        'source_excitation_ids': [r['id'] for r in dictionary['source_excitation_rows']],
        'protected_catalogue_ids': [r['id'] for r in dictionary['protected_catalogue_rows']],
        'evidence': {role: {'path': path, 'sha256': sha(ROOT / path),
                            'schema': values[role]['schema']}
                     for role, path in INPUTS.items()},
        'calculated_source_results': calculated_summary(values),
        'confirmation_plan': steps.plan(),
        'labelled_catalogue': catalogue.build(dictionary),
        'protected_particle_paths': paths.build(dictionary),
        'physical_particle_map': [],
        'mass_predictions_GeV': [],
        'fit_parameters_used': [],
        'physical_gate_promotions': [],
        'source_sha256': {**pins, **{p: sha(ROOT / p) for p in SOURCE_INPUTS}},
    }
    result['content_sha256'] = content_digest(result)
    validate(result, authenticated=(values, pins))
    return result


def content_digest(value):
    return contract.digest({k: v for k, v in value.items() if k not in ('source_sha256', 'content_sha256')})


def validate(value, *, authenticated=None):
    targets = value['target_rows']
    names = [row['target'] for row in targets]
    expected = {'electron', 'muon', 'tau', 'up', 'charm', 'top', 'down', 'strange',
                'bottom', 'electron_neutrino', 'muon_neutrino', 'tau_neutrino',
                'photon', 'gluon', 'W_plus', 'W_minus', 'Z', 'Higgs'}
    if len(names) != 18 or set(names) != expected:
        raise ValueError('candidate identity coverage differs')
    if any(value[key] for key in ('physical_particle_map', 'mass_predictions_GeV',
                                  'fit_parameters_used', 'physical_gate_promotions')):
        raise ValueError('source identity tests do not certify physical particles')
    original = json.loads((ROOT / INPUTS['dictionary']).read_text())
    originals = {r['target']: r for r in original['target_rows']}
    for row in targets:
        if (row['physical_identity_status'] != 'unresolved'
                or row['source_flavor_index'] is not None
                or row['source_derived_mass_GeV'] is not None):
            raise ValueError('candidate identity promoted beyond the calculated tests')
        expected_family = target_family(row['target'])
        expected_role = 'charged' if expected_family in (
            'charged_leptons', 'up_type_quarks', 'down_type_quarks', 'active_neutrinos') else 'bosons'
        if row['family'] != expected_family or row['test_receipt'] != INPUTS[expected_role]:
            raise ValueError('candidate identity test has the wrong source owner')
        source = originals[row['target']]
        if any(row[key] != source[key] for key in ('source_component_ids', 'actual_primary_excitation_ids')):
            raise ValueError('candidate source assignments differ from the original dictionary')
        if row.get('confirmation_steps') != steps.candidate_steps(row['target']):
            raise ValueError('candidate confirmation requirements differ')
        if row.get('supporting_receipts') != supporting_receipts(row['target']):
            raise ValueError('candidate supporting evidence differs')
        if row['source_test_summary'] != SOURCE_TESTS[row['family']]:
            raise ValueError('candidate identity summary differs from its tested domain')
    if value.get('confirmation_plan') != steps.plan():
        raise ValueError('candidate confirmation plan differs')
    if value.get('labelled_catalogue') != catalogue.build(original):
        raise ValueError('candidate labelled catalogue differs')
    if value.get('protected_particle_paths') != paths.build(original):
        raise ValueError('protected-object particle paths differ')
    if len(set(value['source_excitation_ids'])) != 16 or len(set(value['protected_catalogue_ids'])) != 11:
        raise ValueError('candidate source IDs were lost or conflated')
    if value.get('content_sha256') != content_digest(value):
        raise ValueError('candidate identity result content differs')
    values, pins = read_inputs() if authenticated is None else authenticated
    dictionary = values['dictionary']
    expected_evidence = {role: {'path': path, 'sha256': sha(ROOT / path),
                                 'schema': values[role]['schema']}
                         for role, path in INPUTS.items()}
    expected_pins = {**pins, **{p: sha(ROOT / p) for p in SOURCE_INPUTS}}
    if (value['calculated_source_results'] != calculated_summary(values)
            or value['evidence'] != expected_evidence
            or value['source_sha256'] != expected_pins
            or value['source_excitation_ids'] != [r['id'] for r in dictionary['source_excitation_rows']]
            or value['protected_catalogue_ids'] != [r['id'] for r in dictionary['protected_catalogue_rows']]):
        raise ValueError('candidate identity result differs from authenticated producers')


def render(value):
    facts = value['calculated_source_results']
    join = facts['charged']['original_coherent_bare_dressed_join']
    window = facts['charged']['original_compensated_dressed_energy_window']
    filtered = facts['charged_propagation']['filtered_original_window']
    scalar = facts['scalar_response']
    lines = ['# Candidate identity tests', '',
             'These tests distinguish source excitations and constrain their physical',
             'interpretations. No named particle or GeV mass is certified by this table.', '',
             'The source calculations give:', '',
             '- A complete ' + str(join['complete_probe_outcomes']) + '-record compensated '
             'bare/dressed interference experiment, with positive original contrast '
             '`' + join['positive_original_contrast_multiplier_of_density_floor'] + ' * c`. '
             'The transition frequency is unresolved and can include elastic response.',
             '- A separate complete dressed menu has original reported weight greater than `'
             + window['original_dressed_recorded_simple_strict_lower'] + '` in the source pair-transfer '
             'interval `[' + ', '.join(window['reported_transfer_interval']) + ']`. '
             'The compensating endpoint and flux are retained; this is not a single-particle rest mass.',
             '- A local coefficient commutant of dimension ' + str(facts['charged']['local_coefficient_commutant_dimension'])
             + ', with no extra onsite three-generation factor.',
             '- A coherent ' + str(facts['composites']['cyclic_radix']) + '-state increment counter '
             'for the six signed composite interventions plus failure, with an explicit multiple-excitation alias.',
             '- Six interacting gauge multiplets: ' + ', '.join('`' + label + '`' for label in facts['bosons']['multiplet_dimensions'])
             + '. Their complete source menu retains the original rare-sector probability.',
             '- A time-filtered dressed pair has original reported probability above `'
             + filtered['original_filtered_recorded_simple_strict_lower'] + '` in `['
             + ','.join(filtered['reported_transfer_interval']) + ']` at member zero. '
             'Its seed admits a canonical energy/momentum kernel. '
             'The distinct bulk aperture is inactive at member zero, and no limiting residue is transferred.',
             '- All ' + str(facts['composite_spectrum']['component_count'])
             + ' neutral composite components have nonzero original spatial/final-energy '
             'records, with explicit rare-state probability floors. Final energy is not transfer energy.',
             '- The neutral radial scalar detector has total success at least `'
             + scalar['uniform_original_total_success_lower'] + '` throughout the family and above `'
             + scalar['member_zero_simple_strict_lower'] + '` at member zero. Its separate finite dynamical '
             'response is nonzero; uniform total detection includes elastic records.',
             '- At endpoint path lengths 1, 2 and 3, complete filtered dressed records have '
             'probabilities greater than `4/25`, `2/25` and `7/250`. Their broad positive '
             'transfer bands overlap. The pair data do not fix an endpoint mass.',
             '- For the weak-pair/singlet neutral trion, the four Walsh spectral centroids '
             'within either G1 or H lie in `(16/3,38/3)`, span less than `5/9`, and have '
             'ratio less than `53/48`. These are mean energies, not rest masses; the signed '
             'neighbor coefficient has no calculated value or sign.',
             '- The radial scalar has positive original transfer weight in `[1/40,D]` '
             'at least `c*9/(128*a^2*D^2)`, where `rho >= c I` and D is the full energy '
             'enclosure width. Complete finite clocks retain at least three quarters of '
             'this bound in the expanded band. Its uniform elastic floor is `3481/160000`; '
             'the inelastic floor can vanish with capacity and the enormous clock histogram '
             'is not executed.', '',
             '| Target | Shared source test | Calculated discriminator | Physical identity |',
             '| --- | --- | --- | --- |']
    for row in value['target_rows']:
        owner = row['test_receipt'].split('/')[1]
        support = ''.join(' + [support](../' + p.split('/')[1] + '/PROOF.md)'
                          for p in row['supporting_receipts'])
        lines.append('| ' + row['target'] + ' | [' + row['family'] + '](../' + owner
                     + '/PROOF.md)' + support + ' | ' + row['source_test_summary']
                     + ' | Unresolved |')
    lines += ['', 'The neutral composite families have a separate '
              '[conserved-charge identity test](../particle_candidate_composites/PROOF.md). '
              'They are not inserted into the elementary-particle target rows.', '',
              'A physical confirmation needs the same source current and spin dictionary, '
              'a populated propagating sector with controlled dispersion, and matching '
              'interaction records. Generation names require a source flavor discriminator. '
              'A unit conversion supplies none of these tests.', '',
              '[Candidate confirmation steps](CONFIRMATION_STEPS.md); '
              '[per-object mass paths](MASS_PATHS.md); '
              '[labelled defect and excitation families](CATALOGUE.md); '
              '[Calculated source tests and their domains](README.md); '
              '[candidate map](../standard_model_source_map/candidates/CANDIDATES.md).', '']
    return '\n'.join(lines)
