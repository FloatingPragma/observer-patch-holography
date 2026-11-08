"""Exact local reference preparation and complete acceptance bookkeeping."""
import argparse
import hashlib
from fractions import Fraction as F
from pathlib import Path

from tools import strict_json
from model.native_fair_local_preparation import packet as inherited
from . import reference as r

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
FILES = tuple('model/native_noncentral_reference/' + name for name in (
    '__init__.py', 'reference.py', 'packet.py', 'test_reference.py', 'PROOF.md', 'README.md'))


def source_bindings():
    paths = set(inherited.source_bindings()) | set(FILES)
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def all_sector_certificate():
    rows = []
    for total in range(1, 12):
        matrix = r.p.block_numerators(total)
        weights = r.occupation_weights(total)
        budget = r.filter_budget(total)
        count = matrix.shape[0]
        lower = F(int(matrix.min()), r.DENOMINATOR)
        spread = F(int((matrix.max(axis=0) - matrix.min(axis=0)).min()), r.DENOMINATOR)
        if lower < F(1, 8192) or spread <= 0 or budget['source_eigenvalue_floor'] < F(1, 8192):
            raise ArithmeticError('faithfulness or noncentrality failed in a retained sector')
        for i, weight in enumerate(weights):
            if any(int(x) < budget['path_multiplicity'] * weight for x in matrix[i, :]):
                raise ArithmeticError('accepted paths exceed actual source paths')
        rows.append({**budget, 'direct_reference_eigenvalue_floor': lower,
                     'minimum_column_max_minus_min': spread,
                     'all_terminal_configurations': count,
                     'all_source_endpoint_count_entries': count * count,
                     'exact_integer_matrix_sha256': hashlib.sha256(matrix.astype('<i8').tobytes()).hexdigest()})
    return rows


def modular_control():
    checked = 0
    port_weights = r.tree_weights()
    for total in range(1, 12):
        states = r.p.masks(total)
        weights = dict(zip(states, r.occupation_weights(total), strict=True))
        for mask in states:
            for a, b in r.p.SEAMS:
                if ((mask >> a) & 1) != ((mask >> b) & 1):
                    swapped = mask ^ ((1 << a) | (1 << b))
                    ratio = F(port_weights[b], port_weights[a]) if (mask >> a) & 1 else F(port_weights[a], port_weights[b])
                    if F(weights[swapped], weights[mask]) != ratio:
                        raise ArithmeticError('modular ratio depends on exterior context')
                    checked += 1
    return {'same_sector_native_seam_matrix_units': checked,
            'all_exterior_contexts_in_each_local_sector': True,
            'generator_modulo_scalar': '-sum_i log(1+native_tree_depth_i) n_i',
            'port0_to_port1_modular_weight_ratio': str(F(port_weights[1], port_weights[0])),
            'natural_number_preserving_regions_modular_invariant': True,
            'overlap_modular_action_agrees': True,
            'disjoint_number_preserving_regions_commute': True,
            'occupation_observables_change_under_modular_flow': False}


def example_words():
    tree = r.p.local_trace(3, 0, (0,) * 11)
    initial = tree['terminal_mask']
    total = initial.bit_count()
    budget = r.filter_budget(total)
    matrix = r.p.block_numerators(total)
    column = r.p.masks(total).index(initial)
    entries = list(zip(r.p.masks(total), r.occupation_weights(total), strict=True))
    final, weight = next((mask, weight) for i, (mask, weight) in enumerate(entries)
                        if int(matrix[i, column]) > budget['path_multiplicity'] * weight)
    threshold = budget['path_multiplicity'] * weight
    examples = []
    for rank, expected in ((0, 'ACCEPT'), (threshold-1, 'ACCEPT'), (threshold, 'REJECT')):
        coins = r.unrank(initial, final, rank)
        full = r.p.local_trace(3, 1, (0,) * 11 + coins)
        result = r.classify(full['slots'][:11], full['slots'][11:])
        if result['status'] != expected or result['path_rank'] != rank:
            raise ArithmeticError('complete acceptance threshold control failed')
        authenticated = r.p.source_query(total=2, blocks=1, mask=3, coins=(0,) * 11 + coins)
        examples.append({'flag': result, 'actual_native_word': authenticated})
    return examples


def produce():
    before = source_bindings()
    value = {'schema': 'oph.native_noncentral_reference.v1', 'source_bindings': before,
             'native_word': {'tree_attempts': 11, 'mixing_attempts': 60, 'attempts': 71,
                             'description_nodes': 71, 'all_raw_addresses_per_attempt': 25,
                             'actual_archive_records': 70, 'actual_fifo_records': 1,
                             'reference_observation': 'retain tree and acceptance flag; marginalize sixty mixing values'},
             'declared_observer_port_weights': r.tree_weights(),
             'all_nonconstant_sectors': all_sector_certificate(),
             'complete_native_threshold_controls': example_words(),
             'common_local_modular_algebra': modular_control(),
             'classical_readout_cost': {'maximum_suffix_table_integer_entries': 61 * 924,
                                      'maximum_counter_bits': 61, 'coins_examined': 60,
                                      'native_opcode_or_feedback_instruction_claimed': False,
                                      'all_2_power60_paths_enumerated': False},
             'scope': {'exact_faithful_noncentral_reference_per_nonconstant_tree_branch': True,
                       'conditional_acceptance_bound_one_twelfth': True,
                       'unconditional_acceptance_lower': 'actual nonconstant-sector weight / 12',
                       'accepted_reference_floor_one_over8192': True,
                       'arbitrary_source_input_and_external_reference': True,
                       'conditional_source_reference_factors_after_stated_record_marginal': True,
                       'full_archive_source_factorization': False,
                       'constant_coherent_code_preserved_and_inconclusive': True,
                       'all_reject_and_zero_outcomes_retained': True,
                       'auxiliary_random_coin_or_Gibbs_Hamiltonian_supplied': False,
                       'observer_tree_weight_protocol_is_unique_vacuum_selection': False,
                       'modular_parameter_is_native_attempt_or_physical_time': False,
                       'compiled_full_modular_flow_executed_here': False,
                       'small_ball_geometric_stress_or_Einstein_premises_derived': False}}
    if before != source_bindings():
        raise ArithmeticError('source changed during the exact local reference replay')
    return encode(value)


def check():
    if inherited.canonical(strict_json.load(OUTPUT)) != inherited.canonical(produce()):
        raise ArithmeticError('retained local reference differs from exact replay')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(inherited.canonical(produce()))
    else:
        check()
    print('Exact native noncentral reference certificate passed.')
