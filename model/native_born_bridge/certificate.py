"""Derive record weights from endpoint execution, before any matrix semantics.

Inputs are an exact classical initial law, a finite public-record policy and
independent fair placement coins. Quantum probability calculations are used
only by the separate comparison/conditional-tomography certificates.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from itertools import product
import json
from pathlib import Path

from tools import strict_json
from ..observer import native

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASIS = native.states()
INDEX = {state: i for i, state in enumerate(BASIS)}
MAX_PATHS = 100_000


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rational(value):
    require(type(value) in (int, Q), 'exact integer or Fraction required; no floats or Booleans')
    return Q(value)


@lru_cache(maxsize=1)
def carrier_seams():
    # Read the existing source geometry; no channel or trace weights are used.
    from ..observer.quantum_mixing import seams
    return seams()


def seam(value):
    require(type(value) is tuple and len(value) == 2
            and all(type(i) is int for i in value)
            and 0 <= value[0] < value[1] < 12,
            'two increasing source port indices required')
    require(value in carrier_seams(), 'selected pair is not a committed carrier seam')
    return value


def endpoint(state, edge, coin):
    """Execute the source and check its independent floor/ceiling characterization."""
    seam(edge)
    require(type(state) is tuple and all(type(x) is int for x in state)
            and state in INDEX and type(coin) is int and coin in (0, 1),
            'total-two source configuration and binary placement coin required')
    a, b = edge
    actual, record = native.native(state, edge, coin)
    total = state[a]+state[b]
    expected = list(state)
    expected[a] = (total+coin)//2
    expected[b] = total-expected[a]
    require(actual == tuple(expected) and record == (coin, state[a]-expected[a]),
            'source endpoint/record differs from independent balanced repair')
    require(actual in INDEX, 'repair left the complete total-two carrier')
    return actual, (a, b, coin, record[1])


@dataclass(frozen=True)
class PathRow:
    initial: int
    final: int
    record: tuple[tuple[int, int, int, int], ...]

    @property
    def conditional_weight(self):
        return Q(1, 2**len(self.record))


@dataclass(frozen=True)
class Tree:
    steps: int
    rules: tuple
    paths: tuple[PathRow, ...]

    @property
    def records(self):
        return tuple(sorted({p.record for p in self.paths}))


def policy(steps, selector):
    """Snapshot a declared public-prefix rule, once per prefix, on all 78 inputs.

    The returned finite table is the actual policy. The selector receives no
    hidden load or initial-law argument. It may return None to halt early.
    """
    require(type(steps) is int and 0 <= steps <= 10, 'history budget is zero through ten attempts')
    rules = {}
    frontier = {(state, ()) for state in BASIS}
    for _ in range(steps):
        following = set()
        for state, prefix in sorted(frontier):
            if prefix not in rules:
                edge = selector(prefix)
                rules[prefix] = None if edge is None else seam(edge)
            edge = rules[prefix]
            if edge is not None:
                for coin in (0, 1):
                    final, record = endpoint(state, edge, coin)
                    following.add((final, prefix+(record,)))
        frontier = following
    return rules


def compile_tree(steps, rules):
    """All fair-coin leaves for every input, including inputs of later zero mass."""
    require(type(steps) is int and 0 <= steps <= 10, 'history budget is zero through ten attempts')
    require(type(rules) is dict, 'an explicit finite policy table is required')
    used, rows = set(), []

    def descend(initial, state, prefix):
        if len(prefix) == steps:
            rows.append(PathRow(initial, INDEX[state], prefix))
        else:
            require(prefix in rules, 'policy omits a reachable public prefix')
            used.add(prefix)
            edge = rules[prefix]
            if edge is None:
                rows.append(PathRow(initial, INDEX[state], prefix))
            else:
                seam(edge)
                for coin in (0, 1):
                    final, record = endpoint(state, edge, coin)
                    descend(initial, final, prefix+(record,))
        require(len(rows) <= MAX_PATHS, 'complete history exceeds the explicit path budget')

    for i, state in enumerate(BASIS):
        descend(i, state, ())
    require(set(rules) == used, 'policy contains unreachable or post-terminal rules')
    tree = Tree(steps, tuple(sorted(rules.items())), tuple(rows))
    validate_tree(tree)
    return tree


def validate_tree(tree):
    """Exact normalization and reversible record custody, with no probability fit."""
    require(type(tree) is Tree and type(tree.steps) is int and 0 <= tree.steps <= 10
            and type(tree.paths) is tuple and 0 < len(tree.paths) <= MAX_PATHS
            and type(tree.rules) is tuple, 'complete native tree required')
    rules = dict(tree.rules)
    require(len(rules) == len(tree.rules), 'policy repeats a public prefix')
    columns = [Q() for _ in BASIS]
    outputs, used = set(), set()
    for path in tree.paths:
        require(type(path) is PathRow and type(path.initial) is int
                and type(path.final) is int and 0 <= path.initial < 78 and 0 <= path.final < 78,
                'invalid native path indices')
        state = BASIS[path.initial]
        prefix = ()
        require(type(path.record) is tuple and len(path.record) <= tree.steps,
                'bounded immutable record required')
        for a, b, coin, current in path.record:
            require(all(type(x) is int for x in (a, b, coin, current)), 'native record entries must be exact integers')
            require(rules.get(prefix) == (a, b), 'path does not follow its public policy')
            used.add(prefix)
            state, observed = endpoint(state, (a, b), coin)
            require(observed == (a, b, coin, current), 'path relabels a retained native record')
            prefix += (observed,)
        require(INDEX[state] == path.final and
                (len(prefix) == tree.steps or (prefix in rules and rules[prefix] is None)),
                'path terminated outside the stated horizon or halt')
        if len(prefix) < tree.steps:
            used.add(prefix)
        key = (path.final, path.record)
        require(key not in outputs, 'duplicated terminal path or loss of reversible record custody')
        outputs.add(key)
        restored = list(state)
        for a, b, _, current in reversed(path.record):
            restored[a] += current
            restored[b] -= current
        require(tuple(restored) == BASIS[path.initial], 'retained currents do not recover the input')
        columns[path.initial] += path.conditional_weight
    require(all(w == 1 for w in columns), 'missing native leaves or unnormalized conditional law')
    require(set(rules) == used, 'policy contains an unexecuted rule')
    return tree


def initial_law(weights):
    """Any rational classical probability on all 78 configurations; zeros retained."""
    require(type(weights) is tuple and len(weights) == len(BASIS),
            'one exact weight for every source configuration is required')
    result = tuple(rational(w) for w in weights)
    require(all(w >= 0 for w in result) and sum(result) == 1, 'initial law must be nonnegative and normalized')
    return result


def execute(tree, weights):
    """Complete joint (public record, terminal configuration) masses from coin counts."""
    validate_tree(tree)
    weights = initial_law(weights)
    output = {(record, i): Q() for record in tree.records for i in range(len(BASIS))}
    for path in tree.paths:
        output[path.record, path.final] += weights[path.initial]*path.conditional_weight
    require(sum(output.values()) == 1, 'complete classical execution lost mass')
    return output


def record_law(tree, weights):
    joint = execute(tree, weights)
    return {r: sum(joint[r, i] for i in range(len(BASIS))) for r in tree.records}


def values_on(records, values):
    require(isinstance(values, Mapping) and set(values) == set(records),
            'every public outcome, including zero-mass outcomes, is required')
    result = {r: rational(values[r]) for r in records}
    require(all(0 <= x <= 1 for x in result.values()), 'classical effects lie in the unit interval')
    return result


def valuation(distribution, values):
    """The finite-sum positive valuation derived from native record probabilities."""
    require(isinstance(distribution, Mapping) and distribution, 'complete record distribution required')
    mass = {r: rational(w) for r, w in distribution.items()}
    require(all(w >= 0 for w in mass.values()) and sum(mass.values()) == 1,
            'record weights must be a normalized nonnegative law')
    values = values_on(mass, values)
    return sum((mass[r]*values[r] for r in mass), Q())


def event_effect(tree, selected):
    """Input diagonal response, computed by counting selected leaves at each input."""
    validate_tree(tree)
    require(type(selected) is frozenset and selected <= set(tree.records), 'event must be a set of public records')
    effect = [Q() for _ in BASIS]
    for path in tree.paths:
        if path.record in selected:
            effect[path.initial] += path.conditional_weight
    return tuple(effect)


def refinement_certificate(short, long):
    """Every initial state has the same old-record law after an allowed continuation."""
    validate_tree(short)
    validate_tree(long)
    require(long.steps >= short.steps, 'future refinement cannot shorten the history')
    later = dict(long.rules)
    require(all(prefix in later and later[prefix] == edge for prefix, edge in short.rules),
            'continuation changes an earlier policy decision or resumes an explicit halt')
    old_records = set(short.records)
    require(all(r[:short.steps] in old_records for r in long.records),
            'continuation loses an earlier record')
    for record in short.records:
        before = event_effect(short, frozenset((record,)))
        after = event_effect(long, frozenset(r for r in long.records if r[:short.steps] == record))
        require(before == after, 'future pushforward changed an earlier event probability')
    return {'old_public_outcomes': len(short.records), 'new_public_outcomes': len(long.records),
            'input_columns_checked': 78, 'all_initial_laws_and_classical_references': True}


def dyadic_partition(tree, weights, counts, bits, edge=native.E):
    """Realize disjoint dyadic record effects in one complete native suffix.

    Each record supplies a tuple of nonnegative integer interval lengths,
    summing to 2^bits. Every suffix record and every outcome label is retained.
    """
    validate_tree(tree)
    weights = initial_law(weights)
    require(type(bits) is int and 0 <= bits <= 10, 'zero through ten additional coin records required')
    require(isinstance(counts, Mapping) and set(counts) == set(tree.records),
            'every original public record needs its complete interval partition')
    require(all(type(row) is tuple and row and all(type(k) is int and k >= 0 for k in row)
                and sum(row) == 2**bits for row in counts.values()),
            'nonnegative integer interval lengths must partition the whole coin domain')
    sizes = {len(row) for row in counts.values()}
    require(len(sizes) == 1, 'every record must retain the same complete outcome labels')
    labels = tuple(range(next(iter(sizes))))
    require(len(tree.paths)*2**bits <= MAX_PATHS, 'all acceptance suffixes must fit the path budget')
    seam(edge)
    masses = {(record, label): Q() for record in tree.records for label in labels}
    raw_masses = defaultdict(Q)
    suffix_count = 0
    for path in tree.paths:
        for coins in product((0, 1), repeat=bits):
            state, number, suffix = BASIS[path.final], 0, ()
            for coin in coins:
                state, record = endpoint(state, edge, coin)
                suffix += (record,)
                number = 2*number+coin
            cumulative = 0
            for label, count in enumerate(counts[path.record]):
                cumulative += count
                if number < cumulative:
                    break
            mass = weights[path.initial]*path.conditional_weight/Q(2**bits)
            masses[path.record, label] += mass
            raw_masses[path.record, suffix, label] += mass
            suffix_count += 1
    require(sum(masses.values()) == 1, 'dyadic acceptance discarded an outcome')
    distribution = record_law(tree, weights)
    outcome_masses = tuple(sum(masses[r, label] for r in tree.records) for label in labels)
    require(all(masses[r, label] == distribution[r]*Q(counts[r][label], 2**bits)
                for r in tree.records for label in labels),
            'native partition differs from the derived effect valuation')
    return {'outcome_masses': outcome_masses, 'record_masses': masses,
            'raw_record_masses': dict(raw_masses),
            'executed_suffixes': suffix_count}


def dyadic_acceptance(tree, weights, numerators, bits, edge=native.E):
    """Two-label special case, keeping both accept and reject."""
    require(type(bits) is int and 0 <= bits <= 10, 'zero through ten additional coin records required')
    require(isinstance(numerators, Mapping) and all(type(k) is int for k in numerators.values()),
            'exact integer acceptance thresholds required')
    result = dyadic_partition(tree, weights,
                             {r: (k, 2**bits-k) for r, k in numerators.items()}, bits, edge)
    return {'accept': result['outcome_masses'][0], 'reject': result['outcome_masses'][1],
            'record_masses': {(r, label == 0): w for (r, label), w in result['record_masses'].items()},
            'raw_record_masses': result['raw_record_masses'],
            'executed_suffixes': result['executed_suffixes']}


def dyadic_enclosure(distribution, values, bits):
    """Exact rational instance of the uniform order-sandwich extension to real effects."""
    require(type(bits) is int and 0 <= bits <= 64, 'a bounded dyadic precision is required')
    values = values_on(distribution, values)
    denominator = 2**bits
    lower = {r: Q((x*denominator).__floor__(), denominator) for r, x in values.items()}
    upper = {r: Q((x*denominator).__ceil__(), denominator) for r, x in values.items()}
    lo, hi = valuation(distribution, lower), valuation(distribution, upper)
    value = valuation(distribution, values)
    require(lo <= value <= hi and hi-lo <= Q(1, denominator), 'dyadic order enclosure failed')
    return {'lower': lo, 'value': value, 'upper': hi, 'width_upper': Q(1, denominator)}


def source_tree():
    """The existing two-attempt E/F feedback controller, with no initial-state selection."""
    return compile_tree(2, policy(2, lambda p: native.F if p and p[0][3] > 0 else native.E))


def quantum_diagonal_comparison():
    """Only now compare the derived classical rule with both existing quantum lifts."""
    from ..observer import quantum_mixing as quantum
    require(BASIS == quantum.STATES, 'quantum comparison uses a different input carrier')
    averaged = [[Q() for _ in BASIS] for _ in BASIS]
    checks = 0
    for edge in quantum.seams():
        for i, state in enumerate(BASIS):
            for coin in (0, 1):
                final, _ = endpoint(state, edge, coin)
                averaged[INDEX[final]][i] += Q(1, 60)
                checks += 1
    for candidate in ('C', 'E'):
        block = quantum.native_block(candidate, (0,)*12, binary=False)
        require(tuple(map(tuple, averaged)) == block, 'quantum diagonal restriction differs from native counting')
    return {'endpoint_coin_branches_checked': checks, 'configurations': 78,
            'candidate_lifts': ['C', 'E'], 'exact_diagonal_blocks_equal': True,
            'matrix_comparison_used_to_define_classical_weights': False,
            'coherent_extension_selected_by_diagonal_equality': False}


def tomography_attachment():
    """Retain the real implemented complex-tomography result and its declared services."""
    from ..native_interface import tomography
    path = ROOT/'model/native_interface/tomography_results.json'
    old = strict_json.load(path)
    for name, expected in old['source_sha256'].items():
        require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected,
                f'stale conditional tomography source: {name}')
    design = tomography.design_certificate()
    require(design == old['certificate']['exact_design'], 'tomography design differs from its source receipt')
    require(old['all_pass'] and old['certificate']['all_failures_retained'], 'incomplete tomography receipt')
    return {'fresh_exact_design': design, 'source_pins_verified': len(old['source_sha256']),
            'receipt': 'model/native_interface/tomography_results.json',
            'real_hermitian_dimension': 4356, 'native_diagonal_dimension': 66,
            'coherent_code_services_and_quantum_semantics_supplied': True,
            'full_compiler_and_runtime_replayed_by_this_package': False,
            'laboratory_data_or_universal_operational_additivity': False}


def report():
    tree = source_tree()
    uniform = tuple(Q(1, 78) for _ in BASIS)
    binary = tuple(Q(1, 66) if max(s) == 1 else Q() for s in BASIS)
    sparse = tuple(Q(2, 7) if i == 0 else Q(5, 7) if i == 17 else Q() for i in range(78))
    cases = []
    for name, law in [('all78', uniform), ('binary66', binary), ('sparse', sparse)]:
        dist = record_law(tree, law)
        numerators = {r: i % 5 for i, r in enumerate(tree.records)}
        accept = dyadic_acceptance(tree, law, numerators, 2)
        joint = dyadic_partition(tree, law,
                                {r: (i % 3, 1, 3-i % 3) for i, r in enumerate(tree.records)}, 2)
        effect = {r: Q((i % 7)+1, 9) for i, r in enumerate(tree.records)}
        enclosure = dyadic_enclosure(dist, effect, 8)
        # These events are disjoint and exhaust the original public outcome set.
        parts = [frozenset(r for r in tree.records if r[0][2] == c) for c in (0, 1)]
        pulled = [event_effect(tree, event) for event in parts]
        require(all(a+b == 1 for a, b in zip(*pulled)), 'coin partition did not pull back to identity')
        for event, diagonal in zip(parts, pulled):
            require(sum(dist[r] for r in event) == sum(law[i]*diagonal[i] for i in range(78)),
                    'classical pullback square failed')
        cases.append({'initial_law': name, 'total_weight': '1',
                      'zero_weight_public_outcomes_retained': sum(w == 0 for w in dist.values()),
                      'record_law': [{'record': r, 'weight': str(w)} for r, w in dist.items()],
                      'native_dyadic_acceptance': str(accept['accept']),
                      'native_dyadic_rejection': str(accept['reject']),
                      'acceptance_suffixes_executed': accept['executed_suffixes'],
                      'shared_coexistent_refinement_A_B_reject': list(map(str, joint['outcome_masses'])),
                      'shared_refinement_raw_record_count': len(joint['raw_record_masses']),
                      'dyadic_order_enclosure': {k: str(v) for k, v in enclosure.items()}})
    paths = [[p.initial, p.final, p.record, str(p.conditional_weight)] for p in tree.paths]
    prefix = compile_tree(1, {(): native.E})
    refinement = refinement_certificate(prefix, tree)
    inherited = {'path': 'Lean/EventAlgebra/OperationalAdditivityBoundary.lean',
                 'theorems': ['producedCubicValuation_not_additive', 'producedCubicValuation_not_born'],
                 'one_witness_sum_value': '35/64', 'sum_of_witness_values': '143/256',
                 'exact_additivity_defect': str(Q(35, 64)-Q(143, 256)),
                 'scope': 'Inherited exact countermodel on its committed tomographically complete qubit fixture; not a new native readout law.'}
    source = (ROOT/inherited['path']).read_text()
    require(all('theorem '+name in source for name in inherited['theorems']), 'missing inherited countermodel theorem')
    require(inherited['exact_additivity_defect'] == '-3/256', 'countermodel arithmetic changed')
    sources = ['model/native_born_bridge/certificate.py', 'model/native_born_bridge/test_certificate.py',
               'model/native_born_bridge/PROOF.md', 'model/observer/native.py',
               'model/observer/quantum_mixing.py', 'model/native_interface/tomography.py',
               'model/native_interface/tomography_results.json', inherited['path'],
               'Lean/EventAlgebra/FiniteBuschGleason.lean', 'Lean/EventAlgebra/PublicRecordAlgebra.lean']
    return {'schema': 'oph.native_born_bridge.v1', 'all_pass': True,
            'classical_tree': {'input_configurations': 78, 'complete_leaves': len(paths),
                               'public_outcomes': len(tree.records),
                               'path_table_sha256': hashlib.sha256(json.dumps(paths, sort_keys=True).encode()).hexdigest(),
                               'public_record_adaptive_policy': True, 'cases': cases},
            'quantum_diagonal_comparison': quantum_diagonal_comparison(),
            'future_record_pushforward': refinement,
            'conditional_complex_tomography': tomography_attachment(),
            'finite_web_countermodel': inherited,
            'scope': {'native_classical_record_measure_derived': True,
                      'diagonal_density_representation_derived': True,
                      'dyadic_record_effects_executed': True,
                      'arbitrary_finite_history_theorem': 'Analytic proof; executable limited to 78 configurations and 100000 leaves.',
                      'initial_law_selected': False, 'independent_fair_coin_measure_selected': False,
                      'unique_complex_density_from_diagonal_records': False,
                      'complex_effect_operational_equivalence_derived': False,
                      'full_quantum_theory_derived_from_scalar_consensus': False},
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = report()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    if args.check:
        require(strict_json.load(path) == json.loads(json.dumps(result)), 'native Born bridge receipt is stale')
    print(json.dumps({'all_pass': result['all_pass'], 'complete_leaves': result['classical_tree']['complete_leaves'],
                      'public_outcomes': result['classical_tree']['public_outcomes'],
                      'complex_additivity_derived': False}))


if __name__ == '__main__':
    main()
