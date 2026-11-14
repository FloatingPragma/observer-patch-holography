"""A faithful reference from a bounded endpoint-only C/E law adjustment.

The supplied positive port weights select a product-pair reference before the
update is defined. Endpoint heat-bath coins preserve it. This finite control
does not select the weights, implement modular flow, or freeze the complete
physical hypothesis. The C/E record partitions and endpoint map are retained.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import hashlib
import json
from math import gcd, prod
from pathlib import Path

from tools import strict_json
from ..observer import native, record_recovery as records
from ..native_modular_reference import common_caps

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CONTRACT = HERE/'campaign.json'
OUTPUT = HERE/'results.json'
CAMPAIGN = 'OPH-LOCAL-PLACEMENT-FAMILY-2026-10-28-v1'
SEAMS = records.seams()
OLD_WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))


def _probability(value, floor=F()):
    if type(value) is not F or not floor <= value <= 1 or value.denominator > 256:
        raise ValueError('exact bounded rational probability required')


@dataclass(frozen=True)
class Parameters:
    port_weights: tuple[int, ...] = (1,)*12
    seam_weights: tuple[int, ...] = (1,)*30
    c_weight: F = F(1, 2)
    repair_probability: F = F(1)

    def __post_init__(self):
        if (type(self.port_weights) is not tuple or len(self.port_weights) != 12
                or any(type(x) is not int or not 1 <= x <= 32 for x in self.port_weights)
                or gcd(*self.port_weights) != 1
                or max(self.port_weights) > 4*min(self.port_weights)):
            raise ValueError('twelve primitive integer port weights in1..32 with ratio at most4 required')
        if (type(self.seam_weights) is not tuple or len(self.seam_weights) != 30
                or any(type(x) is not int or not 1 <= x <= 16 for x in self.seam_weights)):
            raise ValueError('thirty fixed positive integer seam weights in1..16 required')
        _probability(self.c_weight)
        _probability(self.repair_probability, F(1, 16))

    def description(self):
        return {'port_weights': list(self.port_weights), 'seam_weights': list(self.seam_weights),
                'c_weight': str(self.c_weight), 'repair_probability': str(self.repair_probability)}


@lru_cache(None)
def basis(binary=True):
    if type(binary) is not bool:
        raise ValueError('Boolean shell selector required')
    return tuple(native.states(binary=binary))


def reference(parameters):
    """Compatibility target: the fixed two-load conditional product measure."""
    unnormalized = tuple(prod(w**n for w, n in zip(parameters.port_weights, x, strict=True))
                         for x in basis())
    z = sum(unnormalized)
    return tuple(F(w, z) for w in unnormalized)


def coin_probability(parameters, edge, coin):
    edge = records.seam(edge)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError('binary integer coin required')
    a, b = edge
    # With one endpoint load the native map places it at a for coin1, b for coin0.
    return F(parameters.port_weights[a if coin else b],
             parameters.port_weights[a]+parameters.port_weights[b])


@lru_cache(None)
def raw_maps(law, edge, binary=True):
    """Partial injections with all raw C/E slots, including zero maps.

    The actual Kraus map is sqrt(p_coin) times the returned partial injection.
    No source load is measured to choose the coin: its bias depends only on
    the requested seam's two fixed weights.
    """
    law, edge = records.candidate(law), records.seam(edge)
    states = basis(binary)
    index = {x: i for i, x in enumerate(states)}
    slots = {(c, a): [] for c in (0, 1) for a in (range(-1, 2) if law == 'C' else range(-2, 3))}
    for i, x in enumerate(states):
        for coin in (0, 1):
            y, rec = records.recorded_step(law, x, edge, coin)
            slots[coin, rec[-1]].append((i, index[y]))
    return tuple((coin, value, tuple(pairs)) for (coin, value), pairs in slots.items())


def seam_kernel(parameters, edge):
    states = basis(); index = {x: i for i, x in enumerate(states)}
    rows = [[F() for _ in states] for _ in states]
    for i, x in enumerate(states):
        for coin in (0, 1):
            y, _ = records.recorded_step('C', x, edge, coin)
            rows[i][index[y]] += coin_probability(parameters, edge, coin)
    return tuple(map(tuple, rows))


def kernel(parameters):
    rows = [[F() for _ in basis()] for _ in basis()]
    total = sum(parameters.seam_weights)
    for edge, weight in zip(SEAMS, parameters.seam_weights, strict=True):
        factor = parameters.repair_probability*F(weight, total)
        for i, row in enumerate(seam_kernel(parameters, edge)):
            for j, value in enumerate(row):
                rows[i][j] += factor*value
    for i, row in enumerate(rows):
        row[i] += 1-parameters.repair_probability
    return tuple(map(tuple, rows))


def exact_certificate(parameters):
    pi = reference(parameters)
    quantum = {}
    for law in ('C', 'E'):
        raw_slots = nonzero = checked_columns = 0
        modular_ratios = set()
        for edge in SEAMS:
            output = [F()]*66
            for binary in (True, False):
                effects = [F()]*len(basis(binary))
                for coin, value, pairs in raw_maps(law, edge, binary):
                    p = coin_probability(parameters, edge, coin)
                    assert F(1, 5) <= p <= F(4, 5)
                    assert len({i for i, _ in pairs}) == len(pairs)
                    assert len({j for _, j in pairs}) == len(pairs)
                    for i, j in pairs:
                        effects[i] += p
                        if binary:
                            output[j] += p*pi[i]
                    if binary:
                        raw_slots += 1; nonzero += bool(pairs); checked_columns += len(pairs)
                        ratios = {pi[j]/pi[i] for i, j in pairs}
                        assert len(ratios) <= 1
                        current = records.record_current(law, (*edge, coin, value))
                        expected = F(parameters.port_weights[edge[1]], parameters.port_weights[edge[0]])**current
                        assert not ratios or ratios == {expected}
                        modular_ratios.update(ratios)
                assert effects == [F(1)]*len(effects)
            assert tuple(output) == pi
            p = seam_kernel(parameters, edge)
            assert all(pi[i]*p[i][j] == pi[j]*p[j][i] for i in range(66) for j in range(66))
        quantum[law] = {'binary_raw_slots': raw_slots, 'binary_nonzero_slots': nonzero,
                        'binary_columns': checked_columns, 'H78_and_H66_complete': True,
                        'faithful_binary_reference_fixed_by_every_seam': True,
                        'record_modular_weight_ratios': sorted(modular_ratios),
                        'all_reference_extensions_admitted': True}
    p = kernel(parameters)
    assert all(sum(row) == 1 for row in p)
    assert all(sum(pi[i]*p[i][j] for i in range(66)) == pi[j] for j in range(66))
    # Strict positive edge support connects every two-load configuration; holds exist.
    visited = {0}; todo = [0]
    while todo:
        i = todo.pop()
        for j, value in enumerate(p[i]):
            if value and j not in visited:
                visited.add(j); todo.append(j)
    assert len(visited) == 66 and all(p[i][i] > 0 for i in range(66))
    pair_weights = {tuple(i for i, n in enumerate(x) if n): p
                    for x, p in zip(basis(), pi, strict=True)}
    factor = common_caps.factor_characterization(12, pair_weights)
    assert factor['all_bare_two_port_algebras_invariant']
    return {'parameters': parameters.description(), 'reference': pi,
            'quantum_instruments': quantum, 'classical_stationary_law_unique': True,
            'full_quantum_stationary_uniqueness_claimed': False,
            'regional_factorization': factor,
            'source_independent_scheduler_required': True,
            'modular_covariance': 'sigma^it K_b sigma^-it = (sigma_out/sigma_in)^it K_b; branch scalar cancels in the CP map',
            'physical_clock_or_energy_selected': False}


def old_reference_obstructions():
    """One two-spectator detailed-balance contradiction for every core."""
    rows = []
    for core in combinations(range(12), 2):
        a, b = next((a, b) if a in core else (b, a) for a, b in SEAMS if (a in core) != (b in core))
        inside = next(c for c in core if c != a)
        outside = next(c for c in range(12) if c not in (*core, b))
        # For a->b, detailed balance requires r_forward/r_reverse=pi_out/pi_in.
        ratios = (OLD_WEIGHTS[1]/OLD_WEIGHTS[2], OLD_WEIGHTS[0]/OLD_WEIGHTS[1])
        assert ratios == (F(5, 8), F(8, 25)) and ratios[0] != ratios[1]
        rows.append({'core': core, 'directed_seam': (a, b),
                     'spectators': (inside, outside), 'required_forward_reverse_ratios': ratios})
    return {'target_weights_by_core_occupation': OLD_WEIGHTS,
            'nonproduct_ratio_w1_squared_over_w0_w2': OLD_WEIGHTS[1]**2/(OLD_WEIGHTS[0]*OLD_WEIGHTS[2]),
            'all_66_core_witnesses': rows,
            'excluded': 'positive reversible endpoint-only jump rates with a fixed core and no spectator input',
            'not_excluded': ['irreversible stationary currents', 'larger local neighborhoods or interactions',
                             'record/controller-correlated references', 'co-moving reference metadata',
                             'different references', 'refining carriers']}


def coherence_control(parameters):
    edge = (0, 1)
    x = tuple(int(i in (1, 2)) for i in range(12))
    y = tuple(int(i in (2, 3)) for i in range(12))
    rows = {}
    for law in ('C', 'E'):
        outputs = {}
        for coin in (0, 1):
            xx, rx = records.recorded_step(law, x, edge, coin)
            yy, ry = records.recorded_step(law, y, edge, coin)
            if rx == ry:
                outputs[xx, yy] = outputs.get((xx, yy), F())+coin_probability(parameters, edge, coin)
        rows[law] = sum(outputs.values(), F())
    assert rows['C'] == coin_probability(parameters, edge, 0) and rows['E'] == 0
    return {'input_matrix_unit_pairs': [(1, 2), (2, 3)], 'seam': edge,
            'unread_output_matrix_unit_coefficient': rows,
            'population_agreement_does_not_identify_quantum_record_laws': True}


def source_bindings(paths):
    result = {}
    for name in paths:
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('local source path required')
        result[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(sorted(result.items()))


def contract(path=CONTRACT):
    value = strict_json.load(path)
    if (value['schema'] != 'oph.local_repair_family.contract.v1'
            or value['campaign_id'] != CAMPAIGN or value['boundary_status'] != 'FROZEN_FINITE_LOCAL_CLASS'
            or value['scope'] != 'FINITE_COMPONENT_FAMILY'
            or value['whole_program_frozen'] is not False or value['CP00_accepted'] is not False
            or type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0
            or not value['full_issue_acceptance_missing']):
        raise ValueError('finite class cannot close the complete CP00 gate')
    expected_boundary = {
            'port_weights': {'count': 12, 'integer_min': 1, 'integer_max': 32,
                             'primitive_gcd': 1, 'max_to_min_ratio': 4},
            'seam_weights': {'count': 30, 'integer_min': 1, 'integer_max': 16},
            'probability_max_denominator': 256, 'repair_probability_min': '1/16',
            'law_tags': ['C', 'E'], 'source_and_record_independent_selection': True}
    if json.dumps(value['parameter_boundary'], sort_keys=True) != json.dumps(expected_boundary, sort_keys=True):
        raise ValueError('declared parameter boundary differs from executable family')
    if source_bindings(value['source_sha256']) != value['source_sha256']:
        raise ValueError('predecessor sources changed')
    return value


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(x) for x in value]
    return value


def report():
    frozen = contract()
    biased = Parameters(port_weights=(25, 25)+(8,)*10)
    uniform = Parameters()
    nonuniform_selector = Parameters(port_weights=biased.port_weights,
        seam_weights=tuple(1+i%3 for i in range(30)), repair_probability=F(3, 4), c_weight=F(1, 3))
    pi = reference(biased)
    old = tuple(OLD_WEIGHTS[x[0]+x[1]] for x in basis())
    assert sum(old) == 1
    paths = ['model/local_repair_family/local_bias.py', 'model/local_repair_family/test_local_bias.py',
             'model/local_repair_family/README.md', 'model/local_repair_family/campaign.json',
             'model/completion/defect_algebra.py', 'model/completion/defect_algebra_results.json',
             'tools/strict_json.py', *frozen['source_sha256']]
    return encode({'schema': 'oph.local_repair_family.v1',
        'campaign_id': CAMPAIGN, 'contract_sha256': hashlib.sha256(CONTRACT.read_bytes()).hexdigest(),
        'source_sha256': source_bindings(paths),
        'uniform_control': exact_certificate(uniform),
        'biased_control': exact_certificate(biased),
        'held_nonuniform_scheduler_control': exact_certificate(nonuniform_selector),
        'old_reference_boundary': old_reference_obstructions(),
        'biased_reference_sector_masses': [sum(pi[i] for i, x in enumerate(basis()) if x[0]+x[1] == n) for n in range(3)],
        'biased_reference_trace_distance_from_old_target': sum(abs(a-b) for a, b in zip(pi, old, strict=True))/2,
        'quantum_distinction': coherence_control(biased),
        'theorem_reuse': [
            {'path': 'model/native_modular_reference/common_caps.py',
             'result': 'factor_characterization',
             'application': 'Actual product-pair weights satisfy every four-port identity, so all bare two-port number-preserving algebras are modular invariant.',
             'new_result': 'The explicitly changed local placement coin fixes that same nonuniform reference and is modular covariant, removing the unchanged symmetric-law restriction for this finite control.'},
            {'path': 'Lean/Thermodynamics/FiniteConditionalRepair.lean',
             'theorems': ['heatBath_row_sum', 'heatBath_detailedBalance', 'heatBath_stationary', 'mixture_stationary'],
             'application': 'Population seam repair resamples its one-load endpoint conditional from the independently declared product reference; binary0/2 endpoint totals are fixed.',
             'boundary': 'No identification of the C/E quantum instrument with an arbitrary conditional expectation or physical equilibrium is made.'},
            {'path': 'model/completion/defect_algebra.py',
             'result': 'complete C/E branch star algebras M78 and M66',
             'application': 'All admitted seam and coin probabilities remain strictly positive. Each changed Kraus map is a nonzero scalar multiple of the same native partial injection; the generated branch star algebra and scalar commutant are unchanged. Each law endpoint has this property; positive tagged mixtures and identity holds preserve it.',
             'boundary': 'No nontrivial Hermitian bare-source charge commuting with every active branch is gained. This does not exclude nonreducing invariant subspaces, arbitrary expectation-fixed observables, history/controller charges, metastability, active recovery or refining carriers.'}],
        'scope': {'CP00_accepted': False, 'whole_program_frozen': False,
                  'new_native_placement_law_declared': True,
                  'weights_selected_from_native_law': False, 'new_physical_comparisons': 0,
                  'modular_action_executed': False, 'particle_or_stress_derived': False},
        'remaining_full_issue_requirements': frozen['full_issue_acceptance_missing']})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = report()
    if args.write:
        OUTPUT.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')
    if args.check:
        actual = strict_json.load(OUTPUT)
        if json.dumps(actual, sort_keys=True) != json.dumps(value, sort_keys=True):
            raise ValueError('local-family receipt differs from exact replay')
    print('Local C/E family: complete instruments, fixed faithful reference, covariance and scope verified.')


if __name__ == '__main__':
    main()
