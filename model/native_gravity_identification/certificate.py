"""Replay the finite entropy attachment and its geometric identification tests."""
import argparse
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
from math import comb

from tools import strict_json
from . import entropy as e

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
SOURCE_FILES = tuple('model/native_gravity_identification/' + name for name in (
    '__init__.py', 'entropy.py', 'certificate.py', 'test_entropy.py', 'PROOF.md', 'README.md')) + (
    'model/native_capacity_restriction/restriction.py',
    'model/native_capacity_restriction/PROOF.md',
    'model/observer_time_fields/canonical_sector.py',
    'model/native_fair_local_preparation/preparation.py',
    'model/native_fair_local_preparation/PROOF.md',
    'model/native_action/screen_compatibility.py',
    'model/native_equilibration/construction.py',
    'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
    'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
    'tools/strict_json.py', 'model/requirements.txt')


def source_bindings():
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in SOURCE_FILES}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n'


def expr(terms):
    return [{'argument': str(argument), 'coefficient': str(coefficient)} for argument, coefficient in terms]


def row(value):
    return {key: str(item) if type(item) is F else item for key, item in value.items()}


def small_control():
    """Independent binomial sums and every regional number block at small sizes."""
    cases = matrix_units = 0
    for population in range(1, 13):
        for occupied in range(population + 1):
            for retained in range(population + 1):
                rows, bulk, edge, full = e.finite_entropy_data(population, occupied, retained)
                if sum(item['probability'] for item in rows) != 1:
                    raise ArithmeticError('regional outcome probabilities do not sum to one')
                if e.combine((1, bulk), (1, edge)) != full:
                    raise ArithmeticError('finite central entropy premise failed')
                # The full regional algebra has every matrix unit inside each
                # supported number block. Each commutator has coefficient w-w.
                for item in e.supported(rows):
                    d = item['dimension']
                    if item['eigenvalue'] != F(item['exterior_dimension'], comb(population, occupied)):
                        raise ArithmeticError('reference is not scalar on its entire number block')
                    matrix_units += d*d
                cases += 1
    return {'all_population_occupation_region_cases': cases,
            'number_preserving_regional_matrix_units': matrix_units,
            'impossible_number_outcomes_retained': True,
            'bulk_edge_log_identity_exact': True,
            'modular_commutators_exactly_zero': True}


def report():
    bindings = source_bindings()
    large = []
    for level in (0, 6, 200):
        population = e.restriction.federation_population(level)
        for occupied in (2, population // 3):
            rows = e.blocks(population, occupied, 12)
            floor = min(item['eigenvalue'] for item in e.supported(rows))
            error = F(1, 2**64)
            stable = error < floor
            large.append({'level': level, 'population': str(population), 'occupied': str(occupied),
                          'all_number_outcomes': [row(item) for item in rows],
                          'positive_support_minimum': str(floor),
                          'modular_action_on_supported_number_algebra': 'identity',
                          'assumed_preparation_half_trace_error': str(error),
                          'spectral_stability_condition_passes': stable,
                          'flow_error_per_unit_observable_at_modular_time_one':
                              str(e.modular_flow_allowance(floor, error, F(1))) if stable else None,
                          'global_canonical_preparation_inferred': False})
    rows, bulk, edge, full = e.finite_entropy_data(240, 80, 12)
    variation = e.sector_redistribution(240, 80, 12)
    if variation['edge_ratio'] == 1 or e.combine((1, variation['edge_derivative']), (1, variation['bulk_derivative'])):
        raise ArithmeticError('finite entropy or static-area countercontrol failed')
    tangent = tuple(F(0) for _ in rows)
    if e.first_law_expression(rows, tangent):
        raise ArithmeticError('within-number traceless directions have a modular response')
    # In a two-dimensional subspace of the k=1 block, X swaps two basis
    # vectors. ||X||=1 and ||X||_HS^2=2; choose t=m/4.
    floor = min(item['eigenvalue'] for item in e.supported(rows))
    t = floor / 4
    remainder = e.first_law_remainder(floor, F(2), t, F(1))
    value = {'schema': 'oph.native_gravity_identification.v1',
             'source_sha256': bindings, 'finite_exact_controls': small_control(),
             'growing_source_references': large,
             'same_state_entropy_split': {
                 'population': 240, 'occupied': 80, 'region_ports': 12,
                 'all_number_outcomes': [row(item) for item in rows],
                 'bulk_entropy': expr(bulk), 'edge_entropy': expr(edge), 'full_entropy': expr(full),
                 'full_reference_unchanged_by_sector_product_completion': True,
                 'central_Shannon_entropy_in_bulk': True,
                 'full_entropy_deficit': 'D(rho || tau_(V,q))',
                 'sector_reweighted_entropy_deficit': 'D(p || hypergeometric_number_law)'},
             'same_support_first_law': {
                 'nonconstant_sector_traceless_tangent_slope': '0',
                 'offdiagonal_tangent_hs_squared': '2', 'offdiagonal_tangent_norm': '1',
                 'positive_parameter': str(t), 'relative_entropy_upper': str(remainder),
                 'actual_fair_preparation_error_consumable_with_support_floor': True},
             'static_cut_normalization_control': {
                 'sector_probability_tangent': [str(x) for x in variation['tangent']],
                 'positive_parameter_radius': str(variation['positive_parameter_radius']),
                 'edge_derivative': expr(variation['edge_derivative']),
                 'bulk_derivative': expr(variation['bulk_derivative']),
                 'edge_ratio_exact': str(variation['edge_ratio']),
                 'fixed_graph_cut_derivative': '0',
                 'constant_cut_times_any_coefficient_matches_edge_variation': False},
             'scope': {
                 'canonical_reference_parameters_read_from_source_specification': True,
                 'finite_entropy_split_and_first_law': True,
                 'sector_scalar_reference_modular_action_derived': True,
                 'nontrivial_native_attempt_dynamics_identified_with_modular_flow': False,
                 'regional_canonical_marginal_prepared_from_arbitrary_input': False,
                 'all_reference_correlations_factored': False,
                 'joint_archive_modular_flow_identified': False,
                 'responsive_geometry_or_proper_volume_derived': False,
                 'same_stress_G2_and_G3_derived': False,
                 'geometric_G4_area_and_stationarity_derived': False,
                 'controlled_G5_and_universal_positive_G6_derived': False,
                 'Einstein_equation_or_OPH_program_no_go': False},
             'custody': {'direct_source_and_proof_hashes': True,
                         'full_import_closure_or_inherited_producer_replay_claimed': False}}
    if bindings != source_bindings():
        raise ArithmeticError('source bindings changed during the finite replay')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    if args.write:
        OUTPUT.write_text(canonical(value))
    if args.check and canonical(strict_json.load(OUTPUT)) != canonical(value):
        raise ArithmeticError('finite gravity identification receipt is stale')
    print(json.dumps({'finite_attachment_passes': True, 'geometric_gravity_premises_discharged': False}))


if __name__ == '__main__':
    main()
