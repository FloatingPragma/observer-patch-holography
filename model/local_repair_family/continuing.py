"""Complete fixed-selector repair histories and their accessible FIFO reference.

The source is binary H66. A committed raw record enters a bounded classical
FIFO; an overwritten record is exported, never erased. The stationary object
is the source and a full FIFO after earlier exports are traced out. The growing
archive is neither stationary nor independent of an arbitrary input reference.
Sparse zero blocks are specified by the complete 721-slot branch alphabet.
Integer and Fraction operator coefficients retain exact branch weights.
Float/complex coefficients are numerical exports: a history is refused if
positive branch weighting erases a nonzero real or imaginary component.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

from . import local_bias as local
from ..observer.source_discrimination import circuit

ROOT = Path(__file__).resolve().parents[2]
STATES = local.basis()
IDENTITY = tuple((i, i) for i in range(66))


def _fraction(value):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.denominator <= 0
            or math.gcd(value.numerator, value.denominator) != 1):
        raise ValueError('an original canonical Fraction required')
    return F(value.numerator, value.denominator)


def _parameters(value):
    if type(value) is not local.Parameters:
        raise ValueError('the declared local repair Parameters required')
    return local.Parameters(value.port_weights, value.seam_weights,
                            _fraction(value.c_weight), _fraction(value.repair_probability))


def _capacity(value):
    if type(value) is not int or not 1 <= value <= 256:
        raise ValueError('an integer FIFO capacity in1..256 required')
    return value


@dataclass(frozen=True)
class Record:
    kind: str
    seam: tuple[int, int] | None = None
    law: str | None = None
    coin: int | None = None
    raw: int | None = None
    label: int | None = None
    current: int | None = None


@dataclass(frozen=True)
class Branch:
    record: Record
    weight: F
    pairs: tuple[tuple[int, int], ...]
    modular_ratio: F


def _record(value):
    if type(value) is not Record or type(value.kind) is not str:
        raise ValueError('an original complete raw record required')
    if value.kind == 'hold':
        if value != Record('hold'):
            raise ValueError('hold has no hidden native fields')
        return Record('hold')
    if (value.kind != 'repair' or type(value.seam) is not tuple
            or value.seam not in local.SEAMS or any(type(i) is not int for i in value.seam)
            or type(value.law) is not str or value.law not in ('C', 'E')
            or type(value.coin) is not int or value.coin not in (0, 1)
            or type(value.raw) is not int or not 0 <= value.raw < 1 << circuit.RECORD_BITS[value.law]):
        raise ValueError('complete law/seam/coin/raw record required')
    label = circuit.record_values(value.law, value.coin).get(value.raw)
    current = (None if label is None else local.records.record_current(
        value.law, (*value.seam, value.coin, label)))
    if (value.label != label or value.current != current
            or (label is not None and (type(value.label) is not int or type(value.current) is not int))):
        raise ValueError('raw label and decoded current disagree')
    return Record(value.kind, value.seam, value.law, value.coin, value.raw,
                  value.label, value.current)


def _fifo(value, capacity):
    if type(value) is not tuple or len(value) > capacity:
        raise ValueError('a bounded original FIFO tuple required')
    return tuple(_record(r) for r in value)


def branches(parameters):
    """All721 native/hold slots, including zero laws and invalid raw slots."""
    p = _parameters(parameters)
    result = [Branch(Record('hold'), 1-p.repair_probability, IDENTITY, F(1))]
    total = sum(p.seam_weights)
    for edge, selection in zip(local.SEAMS, p.seam_weights, strict=True):
        for law, mix in (('C', p.c_weight), ('E', 1-p.c_weight)):
            semantic = {(coin, label): pairs for coin, label, pairs in local.raw_maps(law, edge)}
            for coin in (0, 1):
                weight = p.repair_probability*F(selection, total)*mix*local.coin_probability(p, edge, coin)
                values = circuit.record_values(law, coin)
                for raw in range(1 << circuit.RECORD_BITS[law]):
                    label = values.get(raw)
                    current = None if label is None else local.records.record_current(law, (*edge, coin, label))
                    pairs = () if label is None else semantic[coin, label]
                    ratio = F(1) if current is None else F(p.port_weights[edge[1]], p.port_weights[edge[0]])**current
                    result.append(Branch(Record('repair', edge, law, coin, raw, label, current),
                                         weight, pairs, ratio))
    return tuple(result)


def _density(value):
    """Validate finite sparse operators, including arbitrary reference columns.

    Positivity/normalization is a preparation premise. Matrix units and zero
    operators remain valid linear-map inputs; actual outcomes are not normalized.
    """
    if type(value) is not dict:
        raise ValueError('an original sparse source/reference operator required')
    result = {}
    for key, coefficient in value.items():
        if (type(key) is not tuple or len(key) != 4 or any(type(i) is not int or i < 0 for i in key)
                or key[0] >= 66 or key[2] >= 66):
            raise ValueError('density keys are source,reference,source,reference indices')
        if type(coefficient) is F:
            coefficient = _fraction(coefficient)
        elif type(coefficient) not in (int, float, complex):
            raise ValueError('original finite scalar coefficients required')
        elif type(coefficient) in (float, complex) and not (math.isfinite(coefficient.real) and math.isfinite(coefficient.imag)):
            raise ValueError('finite density coefficients required')
        if coefficient:
            result[key] = coefficient
    return result


def _apply(density, branch):
    if not branch.weight or not branch.pairs:
        return {}
    mapping = dict(branch.pairs)
    result = {}
    for (i, a, j, b), value in density.items():
        if not value or i not in mapping or j not in mapping:
            continue
        weighted = branch.weight*value
        if type(value) in (float, complex):
            if not (math.isfinite(weighted.real) and math.isfinite(weighted.imag)):
                raise ArithmeticError('nonfinite floating branch coefficient')
            if ((value.real != 0 and weighted.real == 0)
                    or (value.imag != 0 and weighted.imag == 0)):
                raise ArithmeticError('floating branch coefficient underflow; exact rational real coefficients remain available')
        result[mapping[i], a, mapping[j], b] = weighted
    return result


def apply_branch(density, parameters, record):
    record = _record(record)
    branch = next(b for b in branches(parameters) if b.record == record)
    return _apply(_density(density), branch)


def step_blocks(density, fifo, parameters, *, capacity=1):
    """Snapshot inputs, then yield all raw slots with detached output records.

    Each slot is (new FIFO, exported record or None, output operator). Mutating
    caller inputs or a previously yielded envelope cannot alter later slots.
    """
    capacity = _capacity(capacity)
    fifo, density = _fifo(fifo, capacity), _density(density)
    catalogue = branches(parameters)

    def slots():
        for branch in catalogue:
            exported = _record(fifo[0]) if len(fifo) == capacity else None
            suffix = fifo[1:] if exported is not None else fifo
            after = (*(_record(r) for r in suffix), _record(branch.record))
            yield after, exported, _apply(density, branch)

    return slots()


def history(density, tape, parameters, *, capacity=1):
    """Replay an arbitrary complete word, retaining zero weight and all exports."""
    capacity, p = _capacity(capacity), _parameters(parameters)
    if type(tape) is not tuple:
        raise ValueError('an original complete tape tuple required')
    tape = tuple(_record(record) for record in tape)
    catalogue = {b.record: b for b in branches(p)}
    current, fifo, exported = _density(density), (), []
    for record in tape:
        current = _apply(current, catalogue[record])
        if len(fifo) == capacity:
            exported.append(fifo[0]); fifo = fifo[1:]
        fifo = (*fifo, record)
    return {'density': current, 'fifo': fifo, 'exported': tuple(exported)}


def source_sum(blocks):
    if type(blocks) is not dict:
        raise ValueError('an original classical block/operator mapping required')
    result = defaultdict(F)
    for density in blocks.values():
        for key, value in _density(density).items():
            result[key] += value
    return {key: value for key, value in result.items() if value}


def advance_marginal(blocks, parameters, *, capacity=1):
    """Trace exported records explicitly; retain source/reference and FIFO.

    Zero blocks are sparse zeros, with their slots defined by branches().
    Fixed source/record-independent selection permits grouping by retained suffix.
    This is not the complete archive channel and cannot certify its factorization.
    """
    capacity = _capacity(capacity)
    if type(blocks) is not dict:
        raise ValueError('a sparse classical FIFO/operator mapping required')
    groups = {}
    for fifo, density in blocks.items():
        fifo = _fifo(fifo, capacity)
        suffix = fifo[1:] if len(fifo) == capacity else fifo
        target = groups.setdefault(suffix, defaultdict(F))
        for key, value in _density(density).items():
            target[key] += value
    catalogue = branches(parameters)
    result = {}
    for suffix, density in groups.items():
        for branch in catalogue:
            output = _apply(density, branch)
            if output:
                target = result.setdefault((*suffix, branch.record), defaultdict(F))
                for key, value in output.items():
                    target[key] += value
    return {key: dict(value) for key, value in result.items()}


def stationary(parameters):
    """Exact one-slot joint cq state Ω=J(σ), supported and generally correlated."""
    p = _parameters(parameters)
    sigma = {(i, 0, i, 0): value for i, value in enumerate(local.reference(p))}
    return advance_marginal({(): sigma}, p)


def stationary_word(parameters, tape):
    """One block of Ω_m=J_m(σ); no exponential full-FIFO state is materialized."""
    if type(tape) is not tuple or not 1 <= len(tape) <= 256:
        raise ValueError('a nonempty stationary FIFO word of length at most256 required')
    p = _parameters(parameters)
    sigma = {(i, 0, i, 0): value for i, value in enumerate(local.reference(p))}
    return history(sigma, tape, p, capacity=len(tape))['density']


def word_reference(parameters, tape):
    """Actual supported stationary block and its scalar modular restriction.

    Ω_h=(w_h/kappa_h) P_range σ P_range. For a nonzero block its modular
    conjugation on B(P_range H66) is the restriction of that of σ; the scalar
    cancels. This statement is on the committed direct-sum record algebra.
    """
    p = _parameters(parameters)
    if type(tape) is not tuple or not 1 <= len(tape) <= 256:
        raise ValueError('a finite nonempty original record tuple required')
    catalogue = {b.record: b for b in branches(p)}
    mapping, weight, ratio = dict(IDENTITY), F(1), F(1)
    for record in tape:
        b = catalogue[_record(record)]
        following = dict(b.pairs)
        mapping = {i: following[j] for i, j in mapping.items() if j in following}
        weight *= b.weight; ratio *= b.modular_ratio
    if not mapping or not weight:
        raise ValueError('zero branches have no normalized reference or modular corner')
    sigma = local.reference(p)
    support = tuple(sorted(mapping.values()))
    values = tuple(weight/ratio*sigma[j] for j in support)
    if any(sigma[j] != ratio*sigma[i] for i, j in mapping.items()):
        raise ArithmeticError('native word lost its exact modular character')
    return {'pairs': tuple(mapping.items()), 'support': support, 'weight': weight,
            'modular_ratio': ratio, 'reference_scalar': weight/ratio,
            'diagonal': values, 'mass': sum(values, F())}


def preparation_bound(parameters, attempts, *, capacity=1):
    from . import mixing
    capacity = _capacity(capacity)
    if type(attempts) is not int or attempts < capacity:
        raise ValueError('the full FIFO must be filled before comparison to Ω_m')
    return mixing.mixing_bound(_parameters(parameters), attempts-capacity)


def resource_counts(attempts, *, capacity=1):
    capacity = _capacity(capacity)
    if type(attempts) is not int or attempts < 0:
        raise ValueError('a nonnegative integer model-attempt horizon required')
    return {'source_dimension': 66, 'native_source_address_bits': 9,
            'complete_raw_slots_per_attempt': 721, 'record_encoding_bits': 10,
            'fifo_capacity': capacity, 'fifo_record_bits_upper': 10*capacity,
            'record_writes_per_attempt': 1, 'exports_per_attempt_after_FIFO_full': 1,
            'retained_records': min(attempts, capacity), 'exported_records': max(0, attempts-capacity),
            'total_committed_records': attempts, 'complete_archive_bits_upper': 10*attempts,
            'record_dependent_selection': False, 'hidden_record_erasure': False,
            'physical_memory_or_selector_hardware_supplied': False}


def certificate(parameters):
    p = _parameters(parameters)
    sigma, catalogue = local.reference(p), branches(p)
    columns = [F()]*66
    for branch in catalogue:
        for i, j in branch.pairs:
            columns[i] += branch.weight
            if sigma[j] != branch.modular_ratio*sigma[i]:
                raise ArithmeticError('branch reference ratio depends on hidden source input')
    if columns != [F(1)]*66:
        raise ArithmeticError('complete native/hold instrument is not trace preserving')
    omega = stationary(p)
    initial = {(i, 0, i, 0): value for i, value in enumerate(sigma)}
    if source_sum(omega) != initial or advance_marginal(omega, p) != omega:
        raise ArithmeticError('one-slot source/FIFO stationary reference failed')
    for fifo, block in omega.items():
        r = word_reference(p, fifo)
        expected = {(j, 0, j, 0): value for j, value in zip(r['support'], r['diagonal'], strict=True)}
        if block != expected:
            raise ArithmeticError('stationary cq reference is not its actual branch restriction')
    return {'parameters': p.description(), 'complete_raw_slots': len(catalogue),
            'positive_stationary_record_blocks': len(omega), 'all66_input_completeness': True,
            'exact_joint_source_one_slot_stationarity': True,
            'source_marginal_is_product_pair_reference': True,
            'stationary_record_blocks_are_supported_not_full_rank': True,
            'stationary_word_formula': 'Omega_h=(w_h/kappa_h) P_range sigma P_range',
            'modular_domain': 'direct sum of actual supported classical-record corners',
            'arbitrary_untouched_reference_extension': True,
            'finite_memory_identity': 'accessible_n=J_m E^(n-m), n>=m',
            'preparation_comparison': 'Omega_m tensor rho_R after exports are traced; bound inherited at n-m',
            'full_archive_stationarity_or_reference_factorization_claimed': False}


def archive_counterexample():
    """Exact first-record information persists in every later exported prefix."""
    p = local.Parameters(port_weights=(25, 25)+(8,)*10, c_weight=F(1), repair_probability=F(1))
    x = STATES.index(tuple(int(i in (0, 2)) for i in range(12)))
    y = STATES.index(tuple(int(i in (2, 3)) for i in range(12)))
    record = Record('repair', (0, 1), 'C', 0, 1, 1, 1)
    # Classically correlated input is enough; the complete theorem also permits entanglement.
    initial = {(x, 0, x, 0): F(1, 2), (y, 1, y, 1): F(1, 2)}
    selected = apply_branch(initial, p, record)
    mass0 = sum(v for (i, a, j, b), v in selected.items() if i == j and a == b == 0)
    mass1 = sum(v for (i, a, j, b), v in selected.items() if i == j and a == b == 1)
    if (mass0, mass1) != (F(1, 120), F()):
        raise ArithmeticError('persistent first-record/reference witness changed')
    omega = stationary(p)
    block = omega[(record,)]
    event_mass = sum(block.values(), F())
    sigma = local.reference(p)
    compatible_mass = sum(w for x, w in zip(STATES, sigma) if x[0] == 0 and x[1] == 1)
    if compatible_mass != F(400, 1501) or event_mass != F(20, 4503):
        raise ArithmeticError('stationary source/FIFO correlation witness changed')
    return {'input_pairs': [[0, 2], [2, 3]], 'initial_reference_probabilities': ['1/2', '1/2'],
            'first_record': ['C', [0, 1], 0, 1, 1],
            'event_and_reference_joint_probabilities': ['1/120', '0'],
            'all_later_full_archive_reference_independence_distance_lower': '1/120',
            'witness_is_unchanged_by_later_source_channels': True,
            'stationary_record_probability': str(event_mass),
            'stationary_conditional_source_event_probability': '1',
            'stationary_source_marginal_event_probability': str(compatible_mass),
            'source_times_FIFO_product_claim_is_false': True}


def correctable_clock_certificate():
    """Symbolic generator obstruction, with a finite-time countercontrol.

    Moving effects have supports n_a(1-n_b), n_b(1-n_a); their difference is
    n_a-n_b. Positive biases rescale these effects without changing their span.
    Connectivity and sum n_i=2 imply scalar compression of every n_i in an
    exactly flagged-correctable code. This is a tangent statement. Only a code
    invariant for every real modular parameter has a trivial logical group.
    """
    supports = (((0, 2), (1, 3)), ((0, 3), (1, 2)))
    indices = tuple(tuple(STATES.index(tuple(int(i in pair) for i in range(12)))
                          for pair in column) for column in supports)
    for law in ('C', 'E'):
        for a, b in local.SEAMS:
            raw = local.raw_maps(law, (a, b))
            maps = {(coin, label): {i for i, _ in pairs} for coin, label, pairs in raw}
            plus = maps[0, 1]
            minus = maps[1, -1]
            if any(int(i in plus)-int(i in minus) != x[a]-x[b] for i, x in enumerate(STATES)):
                raise ArithmeticError('actual moving effects do not generate occupation differences')
            for _, _, pairs in raw:
                mapping = dict(pairs)
                images = [{mapping[i] for i in column if i in mapping} for column in indices]
                if len(images[0]) != len(images[1]) or images[0] & images[1]:
                    raise ArithmeticError('finite tick control is not an exact all-seam flagged code')
    connected = {0}
    while True:
        grown = connected | {a for a, b in local.SEAMS if b in connected} | {b for a, b in local.SEAMS if a in connected}
        if grown == connected: break
        connected = grown
    if len(connected) != 12:
        raise ArithmeticError('native seam graph is disconnected')
    # Existing four-support code; powers of2 determine the phases exactly.
    theta = (2, 1, 2)+(1,)*9
    exponents = tuple(tuple(sum(STATES[i][port] for port, weight in enumerate(theta) if weight == 2)
                            for i in column) for column in indices)
    if [sum(row) for row in exponents] != [2, 2] or [sum((-1)**n for n in row) for row in exponents] != [2, -2]:
        raise ArithmeticError('finite tick counterexample changed')
    return {'generator': 'H=-sum_i log(theta_i) n_i (scalar logZ omitted)',
            'exact_correctability_implies': 'P H P=h P',
            'full_one_parameter_invariant_code_implies': 'exp(-itH)P=exp(-ith)P for all real t',
            'finite_compressed_phase_scalar_claimed_without_invariance': False,
            'compression_bound': '||P exp(-itH) P-exp(-ith)P|| <= t^2 L^2/2; L=||(I-P)HP||',
            'leakage_probability_upper': 'min(1,t^2 L^2)',
            'finite_tick_countercontrol': {'theta': list(theta),
                'code_supports': [[[0, 2], [1, 3]], [[0, 3], [1, 2]]],
                'H_compression': '-log(2) I_2', 'parameter': 'pi/log(2)',
                'finite_logical_action': 'diag(1,-1)', 'invariant_at_intermediate_parameters': False}}


def source_bindings():
    paths = ('model/local_repair_family/continuing.py', 'model/local_repair_family/test_continuing.py',
             'model/local_repair_family/local_bias.py', 'model/local_repair_family/mixing.py', 'model/observer/native.py',
             'model/observer/record_recovery.py', 'model/observer/source_discrimination/circuit.py')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def report():
    p = local.Parameters(port_weights=(25, 25)+(8,)*10, c_weight=F(1), repair_probability=F(1))
    return {'schema': 'oph.local_repair_family.continuing.v1',
            'reference': certificate(p), 'archive_countercontrol': archive_counterexample(),
            'memory_clock_boundary': correctable_clock_certificate(),
            'resources_at64_attempts': resource_counts(64), 'source_sha256': source_bindings(),
            'scope': {'fixed_source_and_record_independent_selector': True,
                      'complete_instrument_is_ideal_not_new_source_compilation': True,
                      'bounded_FIFO_is_not_bounded_total_archive': True,
                      'physical_time_energy_or_integrated_observer_acceptance': False}}
