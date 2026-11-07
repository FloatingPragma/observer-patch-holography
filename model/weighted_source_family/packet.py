"""Reproduce the declared weighted native-family certificate."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

import sympy as s

from tools import strict_json
from ..native_program import receipt as native_receipt
from . import family as f


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
SOURCE_FILES = tuple(sorted(set((*f.geometry.SOURCE_PATHS,
    'Lean/Thermodynamics/FiniteConditionalRepair.lean',
    'Lean/ObserverPatchHolography/Locality/NoSignalling.lean',
    'tools/strict_json.py',
    *(f'model/weighted_source_family/{name}' for name in (
        '__init__.py', 'family.py', 'packet.py', 'test_family.py', 'test_packet.py',
        'README.md', 'PROOF.md'))))))


def source_bindings():
    result = {name: digest for name, digest in native_receipt.source_bindings().items()
              if not name.startswith('docs/')}
    result.update({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                   for name in SOURCE_FILES})
    return dict(sorted(result.items()))


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def native_attachment():
    family = f.local_family()
    masks = [sum(value << i for i, value in enumerate(x)) for x in f.local.basis()]
    marker = {(i, i, i, i): F(1, 66) for i in range(66)}
    initial = {(masks[i], i, masks[i], i): value for (i, _, _, _), value in marker.items()}
    program = f.interpreter.program()
    expected = f.interpreter._outcomes(program, program['commands'][0], marker)
    actual = f.attempt(family, initial)
    encoded = []
    if len(expected) != len(actual) or len(actual) != 721:
        raise ArithmeticError('complete native raw-outcome alphabet changed')
    for native, (record, density) in zip(expected, actual, strict=True):
        translated = {(masks[i], a, masks[j], b): value
                      for (i, a, j, b), value in native.density.items()}
        if record != native.record or density != translated:
            raise ArithmeticError('weighted local source attachment changed')
        encoded.append([record.__dict__, [[list(k), str(v)] for k, v in sorted(density.items())]])
    proposed = f.interpreter.program(federation={'level':0, 'total':2})
    proposed['parameters']['port_weights'] = [25, 25]+[8]*10
    proposed['declared_resources'] = f.interpreter._native_cost(proposed)
    try:
        f.interpreter.admit(proposed)
    except ValueError as error:
        rejection = str(error)
        if 'not an invented bias lift' not in rejection:
            raise ArithmeticError('weighted federation rejection changed') from error
    else:
        raise ArithmeticError('weighted federation was silently admitted')
    return {'H66_complete_raw_outcomes': len(actual),
            'zero_raw_outcomes': sum(not density for _, density in actual),
            'complete_partial_injection_attachment_sha256': hashlib.sha256(canonical(encoded).encode()).hexdigest(),
            'all_complex_matrix_units_and_untouched_references': True,
            'connected_weighted_extension_admitted_by_native_interpreter': False,
            'actual_admission_rejection': rejection}


def finite_certificate():
    balance = restrictions = modular = 0
    for weights in ((1, 1), (2, 3, 5), (25, 8, 25, 8)):
        seams = tuple((i, i+1) for i in range(len(weights)-1))
        family = f.Family(weights, seams, tuple(range(1, len(seams)+1)))
        for m in range(len(weights)+1):
            pi = f.canonical(weights, m)
            columns = {x: f.population_column(family, x) for x in pi}
            for x, edge, coin in product(pi, family.seams, (0, 1)):
                y, current, _ = f.image(x, edge, coin, len(weights))
                if pi[y]/pi[x] != f.modular_ratio(family, edge, current):
                    raise ArithmeticError('native modular current gap changed')
                modular += 1
            for x, y in product(pi, repeat=2):
                if pi[x]*columns[x].get(y, 0) != pi[y]*columns[y].get(x, 0):
                    raise ArithmeticError('native weighted detailed balance failed')
                balance += 1
            for bits in range(1 << len(weights)):
                window = tuple(i for i in range(len(weights)) if bits >> i & 1)
                trace = defaultdict(F)
                for mask, p in pi.items():
                    trace[sum(((mask >> i) & 1) << j for j, i in enumerate(window))] += p
                if dict(trace) != f.restricted_state(weights, m, window):
                    raise ArithmeticError('weighted canonical partial trace failed')
                restrictions += 1
    generators = []
    for a, b in ((1, 1), (25, 8), (3, 7)):
        H = f.local_generator(a, b)
        omega = s.Matrix([1, s.sqrt(a), s.sqrt(b), s.sqrt(a*b)])
        if H != H.T or s.simplify(H*H-H) != s.zeros(4) or s.simplify(H*omega) != s.zeros(4, 1):
            raise ArithmeticError('positive native population generator changed')
        generators.append({'endpoint_weights':[a, b], 'singleton_block':[
            [str(H[i, j]) for j in (1, 2)] for i in (1, 2)],
            'orthogonal_projector':True, 'square_root_stationary_vector_annihilated':True})
    return {'exact_detailed_balance_pairs':balance, 'enumerated_full_window_partial_traces':restrictions,
            'exact_modular_current_ratios':modular,
            'modular_gap_formula':'current*log(w_a/w_b)',
            'modular_gap_independent_of_global_size_and_number':True,
            'fixed_product_reference_modular_partial_trace_covariance':True,
            'modular_covariance_is_executed_coherent_service':False,
            'weighted_canonical_number_mixture':True,
            'restriction_nesting_by_elementary_symmetric_convolution':True,
            'fixed_fugacity_product_restriction_compatible':True,
            'local_population_generators':generators,
            'generator_is_native_coherent_or_physical_Hamiltonian':False}


def modular_covariance():
    weights, z = (25, 25, 8, 8), F(1, 8)
    mixture = f.grand_canonical(weights, z)
    pi = {mask: p*mixture[m] for m in range(5) for mask, p in f.canonical(weights, m).items()}
    nonzero = zero = 0
    def retain(mask, sites):
        return sum(((mask >> i) & 1) << j for j, i in enumerate(sites))
    for selection in range(16):
        inside = tuple(i for i in range(4) if selection >> i & 1)
        outside = tuple(i for i in range(4) if not selection >> i & 1)
        local = defaultdict(F)
        for mask, probability in pi.items():
            local[retain(mask, inside)] += probability
        for left, right in product(range(16), repeat=2):
            if retain(left, outside) != retain(right, outside):
                zero += 1
                continue
            if pi[left]/pi[right] != local[retain(left, inside)]/local[retain(right, inside)]:
                raise ArithmeticError('product modular partial-trace covariance failed')
            nonzero += 1
    return {'declared_weights':list(weights), 'declared_fugacity':str(z),
            'nonzero_partial_trace_matrix_units':nonzero,
            'zero_partial_trace_matrix_units':zero,
            'exact_probability_ratios_determine_all_real_modular_times':True,
            'arbitrary_correlated_input_and_untouched_reference':True,
            'executed_coherent_clock_service':False}


def source_federations():
    records = []
    for level in (0, 1, 2):
        family = f.federation_family(level)
        graph = f.geometry.federation(level)
        incident_cells = Counter(i//12 for edge in graph.inter_seams for i in edge)
        if incident_cells != {cell:3 for cell in range(graph.cells)}:
            raise ArithmeticError('a source cell does not have three crossing seams')
        window = tuple(range(12))
        inside = tuple(e for e in family.seams if set(e) <= set(window))
        crossing = f.boundary(family, window)
        if len(inside) != 30 or len(crossing) != 3:
            raise ArithmeticError('source cell boundary census changed')
        sectors = f.window_law(family.weights, 2, window)
        bounds = f.cut_bounds(family, window, attempts=4, internal_repairs=3)
        if bounds['no_boundary_through_kth_internal'] != F(10, 11)**3:
            raise ArithmeticError('local repair clock comparison changed')
        records.append({'level':level, 'ports':len(family.weights), 'seams':len(family.seams),
            'repeated_declared_port_word':list(family.weights[:12]),
            'source_seam_sha256':hashlib.sha256(canonical(family.seams).encode()).hexdigest(),
            'window_ports':list(window), 'internal_seams':len(inside), 'boundary_seams':len(crossing),
            'all_cells_boundary_census_verified':True,
            'declared_global_occupation':2, 'full_local_number_law':[str(p) for p in sectors],
            'conditional_two_load_reference_equals_H66':True,
            'clock_comparison':{key:str(value) for key, value in bounds.items()},
            'comparison_global_attempts':4, 'comparison_internal_repairs':3})
    # An explicit positive-density control, rather than a selected occupation law.
    family = f.federation_family(0)
    density = f.window_law(family.weights, 40, tuple(range(12)))
    if sum(density) != 1 or not all(density):
        raise ArithmeticError('finite-density weighted number mixture changed')
    return {'source_levels':records,
            'declared_finite_density_control':{'ports':240, 'global_occupation':40,
                'full_cell_number_law':[str(p) for p in density]},
            'all_local_number_sectors_retained':True,
            'conditioning_is_not_deterministic_preparation':True}


def produce():
    bindings = source_bindings()
    result = {'schema':'oph.weighted_source_family.v1',
              'native_attachment':native_attachment(),
              'finite_construction':finite_certificate(),
              'algebraic_product_modular_restriction':modular_covariance(),
              'connected_source_graph_examples':source_federations(),
              'scope':{'extension_graph_and_positive_weights_are_explicit_inputs':True,
                       'same_literal_endpoint_translation_and_raw_C_E_codes':True,
                       'fixed_independent_selector_for_history_bounds':True,
                       'arbitrary_correlated_source_exterior_reference_allowed_in_cut_bound':True,
                       'source_selected_weighted_capacity_trajectory':False,
                       'unique_stationary_population_not_unique_coherent_state':True,
                       'unique_coherent_stationary_state_proved':False,
                       'uniform_preparation_or_mixing_bound':False,
                       'local_fixed_two_sector_autonomous_in_connected_source':False,
                       'boundary_event_bound_is_exact_TV_discrepancy':False,
                       'physical_clock_or_capacity_selected':False,
                       'coherent_service_embedding_derived':False,
                       'interacting_relativistic_QFT_or_measured_particle_identification':False},
              'source_bindings':bindings}
    if source_bindings() != bindings:
        raise ArithmeticError('weighted family sources changed during replay')
    return result


def verify(path=OUTPUT):
    if strict_json.load(Path(path)) != produce():
        raise ValueError('weighted native-family receipt differs from exact replay')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', action='store_true')
    actions.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print('Weighted native-family receipt written.')
    else:
        verify()
        print('Weighted native-family receipt verified.')
