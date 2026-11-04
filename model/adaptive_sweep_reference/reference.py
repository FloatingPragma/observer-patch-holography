"""Exact adaptive preparation on the entire native total-two carrier.

All proof calculations use integer branch counts. Numerical execution is an
independent diagnostic of the archive-traced channel, not its certificate.
The continuing archive and the deterministic sweep phase are retained; only
the field and two-bit memory mix at a common sweep boundary.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from tools import strict_json
from ..adaptive_reference import controller
from ..observer import native, quantum_mixing
from ..observer.scalar_response.experiment import canonical, compare_receipt
from ..repair_scheduler import coins


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
STATES = native.states(binary=True)
SEAMS = quantum_mixing.seams()
D = 66
MEMORIES = 4
DENOMINATOR = 1 << 30
BLOCK_TICKS = 4
SOURCE_PATHS = (
    'model/adaptive_sweep_reference/__init__.py',
    'model/adaptive_sweep_reference/__main__.py',
    'model/adaptive_sweep_reference/reference.py',
    'model/adaptive_sweep_reference/test_reference.py',
    'model/adaptive_sweep_reference/README.md',
    'model/adaptive_sweep_reference/PROOF.md',
    'model/adaptive_sweep_reference/contract.json',
    'model/adaptive_reference/controller.py',
    'model/observer/native.py', 'model/observer/quantum_mixing.py',
    'model/repair_scheduler/coins.py',
    'Lean/QFT/AdaptiveReference.lean',
)


def _selector(value):
    if type(value) is not int or value not in (0, 1):
        raise ValueError('a selector bit is required')
    return value


def _count(value):
    if type(value) is not int or not 0 <= value <= 100000:
        raise ValueError('an integer sweep count between zero and 100000 is required')
    return value


def order(selector):
    return SEAMS if _selector(selector) == 0 else SEAMS[::-1]


@dataclass(frozen=True)
class Checkpoint:
    """Executable classical controller; its archive contains every raw event."""

    memory: tuple[int, int]
    phase: int = 0
    history: tuple[tuple, ...] = ()

    def __post_init__(self):
        controller.bits(self.memory)
        if type(self.phase) is not int or not 0 <= self.phase < 30:
            raise ValueError('a classical sweep phase from zero to 29 is required')
        if type(self.history) is not tuple:
            raise ValueError('an immutable complete raw archive is required')
        for record in self.history:
            controller.label(record)

    def selected_event(self):
        selector = self.memory[0 if self.phase == 0 else 1]
        return 'E' if selector else 'C', order(selector)[self.phase]

    def commit(self, record):
        record = controller.label(record)
        memory = controller.advance(self.memory, record) if self.phase == 0 else self.memory
        return Checkpoint(memory, (self.phase+1) % 30, self.history+(record,))

    def recover_initial_control(self):
        memory, phase = self.memory, self.phase
        for record in reversed(self.history):
            phase = (phase-1) % 30
            if phase == 0:
                memory = controller.retreat(memory, record)
        return memory, phase

    def export(self):
        return {'memory': list(self.memory), 'phase': self.phase,
                'history': [list(record) for record in self.history]}

    @staticmethod
    def restore(value):
        if type(value) is not dict or set(value) != {'memory', 'phase', 'history'}:
            raise ValueError('an exact controller checkpoint is required')
        return Checkpoint(tuple(value['memory']), value['phase'],
                          tuple(tuple(record) for record in value['history']))


@lru_cache(maxsize=120)
def source_maps(selector, edge, full=False):
    """Seventeen integer partial injections, including every zero record.

    The source Kraus operator is this map divided by sqrt(2). The outside
    record is zero on the native carrier. No projection occurs at execution.
    """
    _selector(selector)
    if type(full) is not bool or type(edge) is not tuple or edge not in SEAMS:
        raise ValueError('a native seam and Boolean carrier flag are required')
    basis = native.states(binary=not full)
    index = {state: i for i, state in enumerate(basis)}
    maps = {r: [] for r in controller.LABELS}
    for i, state in enumerate(basis):
        for coin in (0, 1):
            out, (_, current) = native.native(state, edge, coin)
            raw = {0: 0, 1: 1, -1: 2}[current]
            if selector:
                raw += 4*((state[edge[0]]+state[edge[1]]) % 2)
            maps['native', coin, raw].append((i, index[out]))
    result = tuple(tuple(maps[r]) for r in controller.LABELS)
    for pairs in result:
        if (len({i for i, _ in pairs}) != len(pairs)
                or len({j for _, j in pairs}) != len(pairs)):
            raise ArithmeticError('native recorded map is not a partial injection')
    if len(basis) != (78 if full else D):
        raise ArithmeticError('native total-two carrier changed')
    return result


def structural_certificate():
    """Attach labels to the independent current/contrast instrument."""
    checked = 0
    for full in (False, True):
        dimension = 78 if full else D
        for selector in (0, 1):
            for edge in SEAMS:
                maps = source_maps(selector, edge, full)
                incoming = [0]*dimension
                outgoing = [0]*dimension
                for pairs in maps:
                    for i, j in pairs:
                        incoming[i] += 1
                        outgoing[j] += 1
                if incoming != [2]*dimension or (not full and outgoing != [2]*dimension):
                    raise ArithmeticError('native completeness or binary unitality failed')
                expected = {}
                for coin, label, pairs in coins.coin_maps('E' if selector else 'C', edge, not full):
                    current = ((label-(2*coin-1)*(abs(label) % 2))//2
                               if selector else label)
                    raw = {0: 0, 1: 1, -1: 2}[current]
                    if selector:
                        raw += 4*(abs(label) % 2)
                    expected['native', coin, raw] = tuple(pairs)
                for record, pairs in zip(controller.LABELS, maps, strict=True):
                    if pairs != expected.get(record, ()):
                        raise ArithmeticError('raw label disagrees with independent current/contrast source')
                checked += len(maps)
    # A positive matrix-unit path cannot cancel another path. If no distinct
    # pair survives the full E word, each fine-record word has rank <= 1.
    surviving = {(i, j) for i in range(D) for j in range(D) if i != j}
    survival_counts = [len(surviving)]
    for edge in order(1):
        next_pairs = set()
        for pairs in source_maps(1, edge):
            image = dict(pairs)
            for i, j in surviving:
                if i in image and j in image:
                    next_pairs.add((image[i], image[j]))
        surviving = next_pairs
        survival_counts.append(len(surviving))
    if surviving:
        raise ArithmeticError('a full E sweep has a non-rank-one record path')
    # Full-sector absorption is support propagation, without conditioning on
    # outcomes. Both selected words visit every source seam.
    absorption = []
    full_basis = native.states()
    binary_set = set(STATES)
    for selector in (0, 1):
        support = set(range(78))
        for edge in order(selector):
            support = {j for pairs in source_maps(selector, edge, True)
                       for i, j in pairs if i in support}
        if {full_basis[i] for i in support} != binary_set:
            raise ArithmeticError('complete sweep did not absorb the full carrier into its binary sector')
        absorption.append(len(support))
    return {
        'source_dimensions': [78, D], 'checked_raw_maps': checked,
        'all_native_records_partial_injections': True,
        'all_native_instruments_complete': True,
        'binary_instruments_unital': True,
        'independent_current_contrast_labels_agree': True,
        'E_distinct_pair_survivors_after_each_attempt': survival_counts,
        'every_fine_E_sweep_record_has_rank_at_most_one': True,
        'full_carrier_support_after_C_and_E_sweeps': absorption,
        'absorption_is_complete_output_support_without_postselection': True,
    }


@lru_cache(maxsize=1)
def population_channel():
    """Exact numerator on 4*66 populations; denominator is 2^30."""
    channel = np.zeros((MEMORIES*D, MEMORIES*D), dtype=object)
    for selector in (0, 1):
        words = order(selector)
        conditional = np.zeros((2, D, D), dtype=object)
        for record, pairs in zip(controller.LABELS, source_maps(selector, words[0]), strict=True):
            bit = controller.feedback(record)
            for i, j in pairs:
                conditional[bit, j, i] += 1
        for edge in words[1:]:
            updated = np.zeros_like(conditional)
            for pairs in source_maps(selector, edge):
                for i, j in pairs:
                    updated[:, j, :] += conditional[:, i, :]
            conditional = updated
        for second in (0, 1):
            before = 2*selector+second
            for feedback in (0, 1):
                after = 2*(second ^ feedback)+selector
                channel[D*after:D*(after+1), D*before:D*(before+1)] += conditional[feedback]
    if (np.any(channel.sum(axis=0) != DENOMINATOR)
            or np.any(channel.sum(axis=1) != DENOMINATOR)):
        raise ArithmeticError('adaptive source population channel is not doubly stochastic')
    channel.flags.writeable = False
    return channel


def _digest(matrix):
    return hashlib.sha256(json.dumps([[int(x) for x in row] for row in matrix],
                                    separators=(',', ':')).encode()).hexdigest()


@lru_cache(maxsize=1)
def minorization():
    """Count only complete histories containing an E sweep.

    Each counted fine history has a single matrix entry, since its E factor
    has rank at most one and every other factor is a partial injection.
    Positive classical entries therefore give a *quantum* CP minorization.
    Classical positivity alone would not prove that implication.
    """
    population = population_channel()
    no_e = population.copy()
    no_e[:, 2*D:] = 0
    total = np.eye(MEMORIES*D, dtype=object)
    only_c = total.copy()
    minimums = []
    for _ in range(BLOCK_TICKS):
        total = total@population
        only_c = only_c@no_e
        minimums.append(int((total-only_c).min()))
    measured = total-only_c
    minimum = minimums[-1]
    alpha = F(MEMORIES*D*minimum, DENOMINATOR**BLOCK_TICKS)
    if minimum <= 0 or not 0 < alpha < 1 or any(x < 0 for x in minimums):
        raise ArithmeticError('four adaptive sweeps have no strict quantum minorization')
    return {
        'population_dimension': MEMORIES*D,
        'CQ_algebra_complex_dimension': MEMORIES*D**2,
        'native_attempts_per_sweep': 30,
        'sweeps_per_minorization': BLOCK_TICKS,
        'population_denominator_per_sweep': DENOMINATOR,
        'population_numerator_sha256': _digest(population),
        'rank_one_path_count_sha256': _digest(measured),
        'minimum_measured_path_counts_by_sweep': minimums,
        'rank_one_path_denominator': DENOMINATOR**BLOCK_TICKS,
        'depolarizing_weight': str(alpha),
        'depolarizing_weight_lower': '3/10',
        'depolarizing_weight_upper': '1/3',
        'CP_minorization': 'Phi^4-alpha D is completely positive; D(X)=Tr(X) I264/264 on direct_sum_4 M66',
        'unique_stationary_boundary_reference': 'I264/264 on direct_sum_4 M66',
        'nontrivial_boundary_spectrum_radius_upper': '(1-alpha)^(1/4)',
        'arbitrary_reference_half_trace_bound': '(1-alpha)^floor(n/4) times the initial half trace distance',
        'stationary_elementary_reference': 'I7920/7920 on direct_sum_120 M66, including thirty classical sweep phases',
        'elementary_phase_period': 30,
        'arbitrary_phase_inputs_mix_to_uniform_phase': False,
        'archive_traced_in_mixing_comparison': True,
        'stationary_reference_modular_automorphism': 'identity',
    }


def preparation_bound(sweeps, full=False):
    """A worst-case bound, also for memory/source/reference correlations."""
    sweeps = _count(sweeps)
    if type(full) is not bool:
        raise ValueError('Boolean full-carrier flag required')
    effective = max(0, sweeps-int(full))
    alpha = F(minorization()['depolarizing_weight'])
    return (1-alpha)**(effective//BLOCK_TICKS)


def preparation_budget(tolerance, full=False):
    if type(tolerance) not in (int, F) or not 0 < tolerance <= 1:
        raise ValueError('an exact rational tolerance in (0,1] is required')
    if type(full) is not bool:
        raise ValueError('Boolean full-carrier flag required')
    for count in range(100001):
        bound = preparation_bound(count, full)
        if bound <= tolerance:
            return {'sweeps': count, 'native_attempts': 30*count,
                    'full_78_state_input': full, 'half_trace_bound': str(bound),
                    'requested_tolerance': str(F(tolerance))}
    raise ValueError('requested accuracy is outside the declared finite budget')


def _apply(value, maps, dimension, reference):
    output = np.zeros_like(value)
    extra = np.arange(reference)
    for pairs in maps:
        if not pairs:
            continue
        sources = (np.array([i for i, _ in pairs])[:, None]*reference+extra).ravel()
        targets = (np.array([j for _, j in pairs])[:, None]*reference+extra).ravel()
        output[np.ix_(targets, targets)] += value[np.ix_(sources, sources)]/2
    return output


def archive_traced_tick(blocks, reference_dimension=1, full=False):
    """Numerical source execution; all fine outcomes contribute their mass."""
    if type(full) is not bool or type(reference_dimension) is not int or reference_dimension < 1:
        raise ValueError('Boolean carrier flag and positive reference dimension required')
    d = 78 if full else D
    dimension = d*reference_dimension
    if (type(blocks) is not np.ndarray or blocks.dtype.kind not in 'fc'
            or blocks.shape != (4, dimension, dimension) or not np.isfinite(blocks).all()):
        raise ValueError('four finite floating/complex density blocks are required')
    mass = sum(np.trace(block) for block in blocks)
    if (abs(mass.imag) > 1e-10 or not -1e-10 <= mass.real <= 1+1e-10
            or any(np.max(abs(block-block.conj().T)) > 1e-10
                   or np.linalg.eigvalsh(block).min() < -1e-10 for block in blocks)):
        raise ValueError('positive Hermitian CQ input with total trace at most one required')
    output = np.zeros_like(blocks, dtype=complex)
    for before in range(4):
        selector, second = divmod(before, 2)
        words = order(selector)
        routed = [np.zeros((dimension, dimension), dtype=complex) for _ in (0, 1)]
        for record, pairs in zip(controller.LABELS, source_maps(selector, words[0], full), strict=True):
            routed[controller.feedback(record)] += _apply(blocks[before], (pairs,), d, reference_dimension)
        for feedback, branch in enumerate(routed):
            for edge in words[1:]:
                branch = _apply(branch, source_maps(selector, edge, full), d, reference_dimension)
            output[2*(second ^ feedback)+selector] += branch
    if not np.isfinite(output).all() or abs(sum(np.trace(b) for b in output)-mass) > 2e-10:
        raise ArithmeticError('complete adaptive sweep lost finite probability mass')
    return output


def _negative_control():
    """The first-record feedback changes the classical next-selector law."""
    edge = order(0)[0]
    responses = []
    for state in STATES:
        weights = [0, 0]
        for coin in (0, 1):
            _, (_, current) = native.native(state, edge, coin)
            raw = {0: 0, 1: 1, -1: 2}[current]
            weights[controller.feedback(('native', coin, raw))] += 1
        responses.append(tuple(weights))
    if len(set(responses)) < 2:
        raise ArithmeticError('declared field-record feedback has become field independent')
    return {'first_record_feedback_coin_numerators_over2': [list(row) for row in sorted(set(responses))],
            'next_selector_depends_on_actual_source_record': True,
            'dropping_E_histories_leaves_zero_minorization': True,
            'initial_phase_zero_stays_phase_zero_at_sweep_boundaries': True,
            'global_archive_stationarity_claimed': False,
            'nontrivial_modular_or_physical_clock_claimed': False}


def run():
    structure = structural_certificate()
    contraction = minorization()
    alpha = F(contraction['depolarizing_weight'])
    if not F(3, 10) < alpha < F(1, 3):
        raise ArithmeticError('exact minorization left its declared coarse interval')
    return {
        'schema': 'oph.adaptive_sweep_reference.v1', 'issues': [1138],
        'gates': ['CP-02'], 'all_pass': True,
        'contract': strict_json.loads(Path(__file__).with_name('contract.json').read_text()),
        'source_attachment': structure,
        'memory': controller.permutation_certificate(),
        'full_carrier_mixing': contraction,
        'finite_preparation_budgets': [preparation_budget(F(1, 10**power), full)
                                     for power in (3, 6) for full in (False, True)],
        'controls': _negative_control(),
        'canonical_physical_premises_discharged': [],
        'integrated_gate_accepted': False,
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                          for p in SOURCE_PATHS},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = run()
    if args.write:
        OUTPUT.write_text(canonical(result))
    if args.check:
        compare_receipt(strict_json.loads(args.check.read_text()), result)
    print(json.dumps({'all_pass': True,
                      'dimension': MEMORIES*D,
                      'four_sweep_depolarizing_weight': result['full_carrier_mixing']['depolarizing_weight'],
                      'full_input_microaccuracy_attempts': result['finite_preparation_budgets'][-1]['native_attempts'],
                      'integrated_gate_accepted': False}, indent=2))
