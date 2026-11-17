"""Modular transport through the complete, committed native repair history.

The implemented carrier is the binary total-two shell, in native.states(True)
order (66 coordinates). A bounded classical policy selects an oriented native
seam or halts, using only the complete earlier seam/coin/C-current or E-contrast
tape. Each branch is K_h=sqrt(w_h) U_h P_h, w_h=2**(-len(h)); U_h is an
injective translation on the diagonal domain P_h. No branch is selected away.

For ANY faithful diagonal sigma and q_h=Tr(P_h sigma), the actual reference
block tau_h=K_h sigma K_h* is faithful on range(K_h), of mass w_h q_h.
On that support, for real t,

    tau_h**(it) K_h = w_h**(it) K_h sigma**(it),
    (tau_h/(w_h*q_h))**(it) K_h = q_h**(-it) K_h sigma**(it).

Indeed [P_h,sigma]=0 and tau_h(U_h e_i)=w_h sigma_i U_h e_i.
The scalar cancels in conjugation: the complete recorded channel intertwines
the input modular group and the direct sum of its ACTUAL output-corner
modular groups. Entrywise proof on matrix units tensors with any untouched
reference. It does not require [K_h,sigma]=0. Zero corners have no normalized
state or logarithm and are refused rather than padded with a fictitious state.

Binary occupancy gives a canonical faithful full-shell extension. A zero
current leaves both endpoint bits fixed; a nonzero current exchanges them.
Compose these recorded port transpositions to obtain pi_h. Then
K_h=sqrt(w_h) U_pi_h P_h, and sigma_h=U_pi_h sigma U_pi_h* satisfies
sigma_h**(it) K_h=K_h sigma**(it). Its compression to the actual range is
tau_h/w_h. The full reference remains faithful with unit mass, while its
weighted compression retains the actual branch mass. This extension composes by
recorded swaps; it is not the full-shell posterior or an executed routing.
Port transpositions need not preserve native graph edges or cap locality.
The extension fails for double occupations and is restricted to binary66.

For two diagonal references sigma,omega the transported raw relative cocycle
has ratio omega_i/sigma_i on U_h e_i. Normalization only multiplies that ratio
by q_sigma/q_omega, a scalar that again cancels in conjugation. Successive
updates multiply the coin weights and restrict the actual injective map;
there is no reset. Unread repair, a fixed endpoint reference, and coherent
interventions need not satisfy these identities. This theorem supplies no
fixed bare-cap locality, autonomous unitary waiting, or physical clock units.

The finite implementation accepts at most 12 attempts and 65,536 reachable
decision nodes. The proof is uniform in any finite bounded native policy.
No dense quantum channel, logarithm, or exponentially long JSON is emitted.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import math
from pathlib import Path

from ..observer import native, record_recovery as records

DIMENSION = 66
MAX_ATTEMPTS = 12
MAX_NODES = 65536
ROOT = Path(__file__).resolve().parents[2]
Branch = records.Branch


def _fraction(value, *, positive=False):
    """Reject subclasses and corrupt writable Fraction slots before arithmetic."""
    if type(value) is int:
        result = F(value)
    elif type(value) is F:
        n, d = value.numerator, value.denominator
        if (type(n) is not int or type(d) is not int or d <= 0
                or math.gcd(abs(n), d) != 1):
            raise ValueError('canonical original rational required')
        result = F(n, d)
    else:
        raise ValueError('original integer or Fraction required')
    if positive and result <= 0:
        raise ValueError('strictly positive rational required')
    return result


def _reference(weights):
    if type(weights) is not tuple or len(weights) != DIMENSION:
        raise ValueError('66 exact positive diagonal weights in native order required')
    values = tuple(_fraction(w, positive=True) for w in weights)
    if sum(values, F()) != 1:
        raise ValueError('reference must have exact unit mass')
    return values


def reference(kind='local'):
    """Fresh exact weights from existing native preparation, not fitted Gibbs data."""
    if type(kind) is not str or kind not in ('local', 'palm', 'common', 'tracial'):
        raise ValueError('local, palm, common, or tracial reference required')
    from . import local_preparation, common_caps, preparation
    if kind == 'local':
        values = local_preparation.reference()
    elif kind == 'palm':
        values = preparation.anchor()[1]
    elif kind == 'common':
        values = common_caps.reference()
    else:
        values = (F(1, DIMENSION),)*DIMENSION
    result = _reference(tuple(values))
    # Existing producer caches contain mutable Fraction instances. Recheck the
    # source-selected core law rather than accepting a poisoned cached tuple.
    if kind in ('local', 'common'):
        by_charge = ((F(2, 225), F(1, 36), F(2, 45)) if kind == 'local'
                     else (F(68, 4500), F(69, 4500), F(60, 4500)))
        expected = tuple(by_charge[x[0]+x[1]] for x in native.states(binary=True))
        if result != expected:
            raise ArithmeticError('cached reference differs from its source-selected core law')
    return result


def _edge(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(x) is not int for x in value)):
        raise ValueError('only an original native seam tuple or halt is supported')
    return records.seam(value)


def _history(law, history):
    if type(history) is not tuple or len(history) > MAX_ATTEMPTS:
        raise ValueError('bounded original complete-history tuple required')
    for record in history:
        if (type(record) is not tuple or len(record) != 4
                or any(type(x) is not int for x in record)):
            raise ValueError('complete original seam/coin/value record required')
        _edge(record[:2])
        records.record_current(law, record)
    return history


def _branch(law, branch):
    """Replay ALL compatible inputs; caller-provided endpoint tables are not authority."""
    if type(branch) is not Branch:
        raise ValueError('original native Branch required')
    history = _history(law, branch.history)
    weight = _fraction(branch.weight, positive=True)
    if weight != F(1, 2**len(history)):
        raise ValueError('branch weight must retain every fair native coin')
    for indices in (branch.inputs, branch.outputs):
        if (type(indices) is not tuple or not indices
                or any(type(i) is not int or not 0 <= i < DIMENSION for i in indices)
                or len(set(indices)) != len(indices)):
            raise ValueError('nonempty injective binary-shell support required')
    basis = native.states(binary=True)
    positions = {state: i for i, state in enumerate(basis)}
    pairs = []
    for i, state in enumerate(basis):
        current = state
        for record in history:
            current, actual = records.recorded_step(law, current, record[:2], record[2])
            if actual != record:
                break
        else:
            if records.restore_load(law, current, history) != state:
                raise ArithmeticError('complete native tape failed exact inversion')
            pairs.append((i, positions[current]))
    if not pairs:
        raise ValueError('zero branch has no supported reference or modular logarithm')
    inputs, outputs = tuple(zip(*pairs, strict=True))
    if branch.inputs != inputs or branch.outputs != outputs:
        raise ValueError('branch must include every actual compatible native input and endpoint')
    return Branch(history, inputs, outputs, weight)


def validate_instrument(law, branches):
    """Detach and replay a complete bounded, record-only native instrument."""
    law = records.candidate(law)
    if type(branches) is not tuple or not branches or len(branches) > MAX_NODES:
        raise ValueError('nonempty bounded tuple of complete native branches required')
    result = tuple(_branch(law, b) for b in branches)
    decisions, seen = {}, set()
    for branch in result:
        if branch.history in seen:
            raise ValueError('duplicate complete record history')
        seen.add(branch.history)
        for n in range(len(branch.history)+1):
            prefix = branch.history[:n]
            action = None if n == len(branch.history) else branch.history[n][:2]
            if prefix in decisions and decisions[prefix] != action:
                raise ValueError('control or stopping must depend only on previous records')
            decisions[prefix] = action
    if len(decisions) > MAX_NODES:
        raise ValueError('finite controller resource ceiling exceeded')
    totals = [F() for _ in range(DIMENSION)]
    for branch in result:
        for i in branch.inputs:
            totals[i] += branch.weight
    if any(w != 1 for w in totals):
        raise ValueError('every coin, record, stopped branch and input column must be retained')
    return result


def _expand(law, choose):
    law = records.candidate(law)
    basis = native.states(binary=True)
    positions = {state: i for i, state in enumerate(basis)}
    stack = [((), tuple((i, i) for i in range(DIMENSION)))]
    terminal, visited = [], set()
    while stack:
        history, pairs = stack.pop()
        visited.add(history)
        if len(visited) > MAX_NODES:
            raise ValueError('finite controller resource ceiling exceeded')
        edge = choose(history)
        if edge is None:
            inputs, outputs = tuple(zip(*pairs, strict=True))
            terminal.append(Branch(history, inputs, outputs, F(1, 2**len(history))))
            continue
        edge = _edge(edge)
        if len(history) == MAX_ATTEMPTS:
            raise ValueError('policy must halt within the implemented attempt ceiling')
        groups = defaultdict(list)
        for coin in (0, 1):
            for initial, current in pairs:
                output, record = records.recorded_step(law, basis[current], edge, coin)
                groups[record].append((initial, positions[output]))
        for record, images in sorted(groups.items(), reverse=True):
            stack.append(((*history, record), tuple(images)))
    result = validate_instrument(law, tuple(sorted(terminal, key=lambda b: b.history)))
    return result, visited


def word_instrument(law, word):
    if type(word) is not tuple or len(word) > MAX_ATTEMPTS:
        raise ValueError('bounded original tuple of native seams required')
    word = tuple(_edge(edge) for edge in word)
    return _expand(law, lambda h: word[len(h)] if len(h) < len(word) else None)[0]


def instrument(law, policy):
    """Explicit finite history->seam/None table; no callable or quantum controls."""
    law = records.candidate(law)
    if type(policy) is not dict or not policy or len(policy) > MAX_NODES:
        raise ValueError('bounded plain record-only decision dictionary required')
    decisions = {_history(law, key): None if action is None else _edge(action)
                 for key, action in policy.items()}

    def choose(history):
        if history not in decisions:
            raise ValueError('policy omits a reachable complete record prefix')
        return decisions[history]

    result, visited = _expand(law, choose)
    if visited != set(decisions):
        raise ValueError('policy contains unreachable record decisions')
    return result


@dataclass(frozen=True)
class ReferenceBlock:
    branch: Branch
    support: tuple[int, ...]  # branch output order, paired with input order
    unnormalized: tuple[F, ...]
    normalized: tuple[F, ...]
    mass: F
    domain_mass: F


def _reference_block(branch, sigma):
    q = sum((sigma[i] for i in branch.inputs), F())
    if not q:
        raise ValueError('zero branch cannot be normalized')
    raw = tuple(branch.weight*sigma[i] for i in branch.inputs)
    return ReferenceBlock(branch, branch.outputs, raw,
                          tuple(sigma[i]/q for i in branch.inputs), branch.weight*q, q)


def branch_reference(law, branch, weights):
    """One actual nonzero branch, with raw and normalized support-corner state."""
    return _reference_block(_branch(records.candidate(law), branch), _reference(weights))


def port_permutation(law, branch):
    """A record-derived unitary extension, not a symmetry of the native graph.

    On binary loads a nonzero current swaps the two endpoint bits, while a
    zero current fixes every bit. Their chronological port transpositions
    extend the actual partial injection to all 66 binary configurations.
    The extension is only a clock compiler choice: its off-domain action is
    not a repair branch, an executed routing service or a prepared posterior.
    """
    law = records.candidate(law)
    branch = _branch(law, branch)
    permutation = tuple(range(12))
    for record in branch.history:
        a, b = record[:2]
        if records.record_current(law, record):
            permutation = tuple(b if p == a else a if p == b else p for p in permutation)
    basis = native.states(binary=True)
    for i, j in zip(branch.inputs, branch.outputs, strict=True):
        output = [0]*12
        for source, target in enumerate(permutation):
            output[target] = basis[i][source]
        if tuple(output) != basis[j]:
            raise ArithmeticError('record permutation fails the actual native branch')
    return permutation


def extended_reference(law, branch, weights):
    """Faithful full-shell extension whose compression is tau_h / w_h."""
    sigma = _reference(weights)
    permutation = port_permutation(law, branch)
    basis = native.states(binary=True)
    index = {state: i for i, state in enumerate(basis)}
    result = [F()]*DIMENSION
    for state, weight in zip(basis, sigma, strict=True):
        output = [0]*12
        for source, target in enumerate(permutation):
            output[target] = state[source]
        result[index[tuple(output)]] = weight
    return tuple(result)


def transport(law, branches, weights):
    sigma = _reference(weights)
    return tuple(_reference_block(b, sigma) for b in validate_instrument(law, branches))


def phase_ratio_certificate(law, branches, weights, other=None):
    """Exact all-matrix-unit certificate, including relative cocycle scalars."""
    sigma = _reference(weights)
    omega = sigma if other is None else _reference(other)
    blocks = transport(law, branches, sigma)
    checked = 0
    for block in blocks:
        second = _reference_block(block.branch, omega)
        for a, i in enumerate(block.branch.inputs):
            if block.unnormalized[a]/sigma[i] != block.branch.weight:
                raise ArithmeticError('raw modular intertwining scalar changed')
            if block.normalized[a]/sigma[i] != 1/block.domain_mass:
                raise ArithmeticError('normalized modular intertwining scalar changed')
            if second.unnormalized[a]/block.unnormalized[a] != omega[i]/sigma[i]:
                raise ArithmeticError('raw relative cocycle changed')
            expected = omega[i]/sigma[i]*block.domain_mass/second.domain_mass
            if second.normalized[a]/block.normalized[a] != expected:
                raise ArithmeticError('normalized relative cocycle scalar changed')
            for b, j in enumerate(block.branch.inputs):
                if (block.unnormalized[a]/block.unnormalized[b] != sigma[i]/sigma[j]
                        or block.normalized[a]/block.normalized[b] != sigma[i]/sigma[j]):
                    raise ArithmeticError('transported matrix-unit modular frequency changed')
                checked += 1
    if sum((b.mass for b in blocks), F()) != 1:
        raise ArithmeticError('complete transported reference lost mass')
    return {'branches': len(blocks), 'matrix_units_checked': checked,
            'support_rank_counts': {str(k): v for k, v in sorted(Counter(len(b.support) for b in blocks).items())},
            'complete_reference_mass': '1', 'actual_support_faithfulness': True,
            'raw_and_normalized_phase_ratios_exact': True,
            'relative_cocycle_ratios_exact': True,
            'arbitrary_untouched_reference_identity': True}


def sandwich_bound(epsilon, parameter):
    """History-uniform modular errors from an unnormalized CP sandwich.

    Assume (1-epsilon)sigma <= rho <= (1+epsilon)sigma; rho need not be
    diagonal or normalized. Every branch has the SAME sandwich on its actual
    range, however rare it is. Operator monotonicity of log bounds the
    logarithm difference between log(1-epsilon)I and log(1+epsilon)I, hence
    its norm by epsilon/(1-epsilon). Normalization subtracts a scalar log,
    which cancels in modular conjugation. Duhamel bounds a unitary channel's
    half-diamond distance by |s| times this generator error, capped at one.

    The complete output is a fixed direct-sum channel, not the nonlinear
    posterior normalization map. Its norm is the maximum over blocks, so no
    rare-record denominator or branch-count factor occurs. A triangle bound
    on input and output modular channels gives the approximate intertwiner.
    The modular maps act only on work; an arbitrary correlated JOINT
    reference's own modular group is not claimed to intertwine.
    """
    epsilon, parameter = _fraction(epsilon), _fraction(parameter)
    if not 0 <= epsilon < 1:
        raise ValueError('faithful sandwich requires exact 0 <= epsilon < 1')
    generator = epsilon/(1-epsilon)
    return {'sandwich_radius': epsilon, 'modular_parameter': parameter,
            'generator_error_up_to_scalar': generator,
            'input_modular_half_diamond_upper': min(F(1), abs(parameter)*generator),
            'output_direct_sum_modular_half_diamond_upper': min(F(1), abs(parameter)*generator),
            'complete_intertwining_half_diamond_upper': min(F(1), 2*abs(parameter)*generator),
            'uniform_over_nonzero_history_probabilities': True,
            'normalized_posterior_is_a_linear_channel': False}


def preparation_certificate(parameter=F(1), precision_bits=24):
    """Consume the existing complete native preparation contract, with its flags.

    Its successful UNNORMALIZED state divided by the ideal herald mass has
    the sandwich radius below. The unknown actual success normalization is
    only a modular scalar. This does not expose the private preparation tape
    or treat preparation success as certain; the four flags and their full
    channel budget remain the predecessor contract.
    """
    from . import local_preparation
    budget = local_preparation.budget(precision_bits=precision_bits)
    bounds = sandwich_bound(budget['operator_sandwich_radius'], parameter)
    return {**{k: str(v) if type(v) is F else v for k, v in bounds.items()},
            'preparation_complete_flagged_halftrace_upper': str(budget['complete_flagged_reference_halftrace_upper']),
            'preparation_success_probability_interval': list(map(str, budget['success_probability_interval'])),
            'preparation_flags': list(budget['terminal_flags_in_comparison']),
            'preparation_fine_tapes_accessible_to_later_controller': False,
            'total_preparation_native_attempts_upper': budget['total_native_attempts_upper'],
            'history_ticks_and_modular_parameter_identified': False}


def covariance_witness(law, branches, weights, *, forget_records=False):
    """First exact nonzero matrix-unit obstruction for fixed or unread reference.

    With records retained, compare to the unchanged source sigma at endpoints.
    With records forgotten, compare to the ACTUAL unread output Phi(sigma).
    A changed positive eigenvalue ratio proves failure for the full real
    modular group; isolated times may still agree by phase periodicity.
    """
    if type(forget_records) is not bool:
        raise ValueError('original boolean record-forgetting choice required')
    sigma = _reference(weights)
    branches = validate_instrument(law, branches)
    target = list(sigma) if not forget_records else [F()]*DIMENSION
    if forget_records:
        for branch in branches:
            for i, out in zip(branch.inputs, branch.outputs, strict=True):
                target[out] += branch.weight*sigma[i]
    for branch in branches:
        for i, out in zip(branch.inputs, branch.outputs, strict=True):
            for j, other_out in zip(branch.inputs, branch.outputs, strict=True):
                left, right = sigma[i]/sigma[j], target[out]/target[other_out]
                if left != right:
                    coefficient = branch.weight
                    if forget_records:
                        coefficient = sum((b.weight for b in branches
                            if (i, out) in tuple(zip(b.inputs, b.outputs))
                            and (j, other_out) in tuple(zip(b.inputs, b.outputs))), F())
                    return {'history': [list(r) for r in branch.history],
                            'input_matrix_unit': [i, j], 'output_matrix_unit': [out, other_out],
                            'input_ratio': str(left), 'comparison_output_ratio': str(right),
                            'ratio_mismatch': str(left/right),
                            'nonzero_channel_coefficient': str(coefficient),
                            'comparison': 'actual_unread_reference' if forget_records else 'fixed_source_reference'}
    return None


def adaptive_policy(law):
    """Two attempts crossing the protected core, followed by record-only halt."""
    law = records.candidate(law)
    policy = {(): (0, 5)}
    basis = native.states(binary=True)
    for first in word_instrument(law, ((0, 5),)):
        edge = (0, 1) if records.record_current(law, first.history[-1]) else (1, 8)
        policy[first.history] = edge
        for output in first.outputs:
            for coin in (0, 1):
                _, record = records.recorded_step(law, basis[output], edge, coin)
                policy[(*first.history, record)] = None
    return policy


def lost_coin_control():
    """Exact smaller erasure: retain C current zero, but hide its fair coin."""
    sigma = reference()
    basis = native.states(binary=True)
    configuration = lambda pair: tuple(int(i in pair) for i in range(12))
    u, v = basis.index(configuration((0, 1))), basis.index(configuration((1, 2)))
    branches = tuple(b for b in word_instrument('C', ((0, 5),)) if b.history[-1][-1] == 0)
    diagonal = [F()]*DIMENSION
    coefficient = F()
    for branch in branches:
        pairs = dict(zip(branch.inputs, branch.outputs))
        for i, output in pairs.items():
            diagonal[output] += branch.weight*sigma[i]
        if pairs.get(u) == u and pairs.get(v) == v:
            coefficient += branch.weight
    before, after = sigma[u]/sigma[v], diagonal[u]/diagonal[v]
    if (before, after, coefficient) != (F(8, 5), F(4, 5), F(1, 2)):
        raise ArithmeticError('coin-erasure native counterexample changed')
    return {'retained_record': 'C current zero at seam (0,5); coin erased',
            'input_and_output_matrix_unit': [u, v], 'input_ratio': str(before),
            'coarse_output_ratio': str(after), 'surviving_coefficient': str(coefficient),
            'modular_frequency_difference': 'log(2)',
            'matrix_unit_entry_defect_at_s_pi_over_log2': '1'}


def noncommuting_reference_control():
    """A faithful, isospectral Hadamard rotation violates the required domain condition.

    Rotate u inside P and v outside P. At t=pi/log(a/b), the two eigenphases
    are opposite, so rho^(it) swaps these coordinates up to a scalar. K kills
    v whereas K rho^(it) sends it to the image of u with modulus sqrt(1/2).
    The defective two-column row has two entries of that modulus and norm1.
    This is a negative control, not an admitted coherent policy operation.
    """
    sigma = reference()
    basis = native.states(binary=True)
    configuration = lambda pair: tuple(int(i in pair) for i in range(12))
    u, v = basis.index(configuration((0, 1))), basis.index(configuration((1, 5)))
    branch = next(b for b in word_instrument('C', ((0, 5),)) if b.history[-1][-1] == 1)
    if u not in branch.inputs or v in branch.inputs or sigma[u]/sigma[v] != F(8, 5):
        raise ArithmeticError('noncommuting domain control changed')
    return {'history': [list(r) for r in branch.history], 'rotated_coordinates': [u, v],
            'rotation': 'Hadamard on two coordinates, identity on their complement',
            'positive_eigenvalues_preserved': True, 'domain_commutator_nonzero': True,
            'test_parameter': 'pi/log(8/5)', 'outside_column_defect_squared': '1/2',
            'raw_intertwiner_operator_defect_squared': '1',
            'admitted_as_exact_diagonal_reference': False}


def exact_certificate():
    sigma, omega = reference('local'), reference('palm')
    rows = {}
    for law in ('C', 'E'):
        branches = instrument(law, adaptive_policy(law))
        simple = word_instrument(law, ((0, 5),))
        rows[law] = phase_ratio_certificate(law, branches, sigma, omega)
        rows[law]['fixed_reference_counterexample'] = covariance_witness(law, simple, sigma)
        rows[law]['forgotten_record_counterexample'] = covariance_witness(law, simple, sigma, forget_records=True)
        if any(rows[law][key] is None for key in ('fixed_reference_counterexample', 'forgotten_record_counterexample')):
            raise ArithmeticError('negative controls stopped distinguishing transported references')
    return {'schema': 'oph.native_modular_reference.history_transport.v1',
            'carrier_dimension': DIMENSION, 'carrier': 'binary total-two native shell',
            'reference': {'origin': 'existing locally incident-prepared common-cap reference',
                          'weight_multiplicities': {str(k): v for k, v in sorted(Counter(sigma).items())},
                          'minimum_eigenvalue': str(min(sigma))},
            'laws': rows,
            'finite_preparation': preparation_certificate(),
            'lost_coin_control': lost_coin_control(),
            'noncommuting_reference_control': noncommuting_reference_control(),
            'domain': {'fixed_law_C_or_E': True, 'record_adaptive_seams_and_stopping': True,
                       'maximum_implemented_attempts': MAX_ATTEMPTS, 'maximum_decision_nodes': MAX_NODES,
                       'coherent_interventions_accepted': False, 'full_78_carrier_claimed': False,
                       'full_66_faithfulness_of_each_record_block': False,
                       'fixed_bare_cap_locality_claimed': False, 'autonomous_modular_waiting_claimed': False,
                       'physical_duration_or_energy_normalization_derived': False}}


def source_bindings():
    from . import local_preparation, common_caps
    paths = set(local_preparation.SOURCE_PATHS) | set(common_caps.SOURCE_PATHS)
    paths.update(('model/native_modular_reference/history_transport.py',
                  'model/native_modular_reference/test_history_transport.py',
                  'model/observer/record_recovery.py', 'model/observer/native.py',
                  'model/common_source/slow_geometry.py'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}
