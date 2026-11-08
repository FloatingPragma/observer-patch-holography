"""Operator covariance certified on actual native partial-injective columns."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F

from model.common_source import slow_geometry as geometry
from model.observer import native
from model.protected_spin_limit.algebra import protected_columns, spin_mask


@dataclass(frozen=True)
class Branch:
    seam: tuple
    coin: int
    current: int
    displacement: tuple
    pairs: tuple


def _seam(edge):
    if (type(edge) is not tuple or len(edge) != 2
            or any(type(i) is not int for i in edge)
            or edge not in geometry.source().seams()):
        raise ValueError('an ordered native seam in the source chart is required')
    return edge


def branches(edge):
    """All six complete (coin,current) slots, K=S/sqrt(2)."""
    edge = _seam(edge)
    states = native.states()
    index = {state: i for i, state in enumerate(states)}
    rows = []
    for coin in (0, 1):
        slots = {j: [] for j in (-1, 0, 1)}
        for i, state in enumerate(states):
            output, record = native.native(state, edge, coin)
            if record[0] != coin:
                raise ArithmeticError('the native record lost its placement coin')
            slots[record[1]].append((i, index[output]))
        for current, pairs in slots.items():
            if len({b for _, b in pairs}) != len(pairs):
                raise ArithmeticError('the committed native branch lost injectivity')
            rows.append(Branch(edge, coin, current,
                               geometry.current_displacement(edge, current), tuple(pairs)))
    return tuple(rows)


def branch_certificate():
    """Every entry of X_i S-S X_i=d_i S and of sum K* K=I."""
    c = geometry.source()
    states = native.states()
    coordinates = tuple(geometry.load_readout(state) for state in states)
    count = 0
    support_sizes = Counter()
    for edge in c.seams():
        incoming = Counter()
        for branch in branches(edge):
            support_sizes[len(branch.pairs)] += 1
            for source, target in branch.pairs:
                difference = geometry._subtract(coordinates[target], coordinates[source])
                if difference != branch.displacement:
                    raise ArithmeticError('a native branch fails operator position covariance')
                incoming[source] += 1
                count += 1
        if incoming != Counter({i: 2 for i in range(len(states))}):
            raise ArithmeticError('complete native branches lost normalization')
    if count != 4680:
        raise ArithmeticError('the whole native total-two carrier is required')
    return {'states': len(states), 'seams': len(c.seams()),
            'branch_slots': 6*len(c.seams()), 'source_columns': count,
            'coordinate_entry_identities': 3*count,
            'branch_support_size_histogram': {str(k): v for k, v in sorted(support_sizes.items())},
            'kraus_square_amplitude': '1/2',
            'matrix_identity': 'X_i K_r-K_r X_i=d_(r,i) K_r',
            'weyl_identity': 'exp(i q.X) K_r exp(-i q.X)=exp(i q.d_r) K_r',
            'complete_instrument_identity': '(X_out-D_record) V=V X_in',
            'arbitrary_coherent_input_and_untouched_reference': True}


def complete_observer_certificate():
    """Actual record-adaptive two-event observer, with all four coin words."""
    checked = 0
    terminal_labels = set()
    for initial in native.states():
        for _, first in native.tick(native.State(initial)):
            edge = native.selected_seam(first)
            for _, final in native.tick(first):
                displacement = geometry._add(
                    geometry.current_displacement(native.E, first.prefix[0][1]),
                    geometry.current_displacement(edge, final.prefix[1][1]))
                if geometry._subtract(geometry.load_readout(final.loads),
                                      geometry.load_readout(initial)) != displacement:
                    raise ArithmeticError('the complete observer lost its displacement sum')
                label = (final.loads, final.prefix)
                if label in terminal_labels:
                    raise ArithmeticError('complete observer columns ceased to be disjoint')
                terminal_labels.add(label)
                checked += 1
    if checked != 312:
        raise ArithmeticError('the complete adaptive observer has missing outcomes')
    return {'initial_columns': 78, 'complete_history_columns': checked,
            'each_coin_word_probability': '1/4',
            'actual_policy': 'first E; positive committed current chooses F, otherwise E',
            'history_displacement': 'sum_r d_r=X_terminal-X_initial',
            'complete_isometry_and_reference_extension': True}


def mean_certificate():
    """Uniform-seam/fair-coin pullback is contraction, not a translation."""
    c = geometry.source()
    alpha = c.q5(F(55, 60), F(1, 60))
    for state in native.states():
        value = (c.q5(0),)*3
        for edge in c.seams():
            for coin in (0, 1):
                output, _ = native.native(state, edge, coin)
                value = geometry._add(value, geometry.load_readout(output))
        value = geometry._scale(F(1, 60), value)
        wanted = tuple(c.q5_mul(alpha, x) for x in geometry.load_readout(state))
        if value != wanted:
            raise ArithmeticError('uniform native clock lost its slow-band contraction')
    return {'diagonal_operator_entries_checked': 234,
            'slow_laplacian_eigenvalue': '5-sqrt(5)',
            'uniform_channel_pullback': 'Phi*(X_i)=(55+sqrt(5))/60 X_i',
            'uniform_calibration_is_the_adaptive_observer_schedule': False,
            'unitary_translation_or_momentum_generator_identified': False}


def projected_load_certificate():
    """Actual protected cube: all total onsite occupations equal one."""
    sites = 8
    for word in range(1 << sites):
        mask = spin_mask(word, sites)
        for site in range(sites):
            total = ((mask >> (2*site)) & 1)+((mask >> (2*site+1)) & 1)
            if total != 1:
                raise ArithmeticError('a protected source spin has nonunit total occupation')
    original = protected_columns()
    return {'protected_spin_columns': 256, 'onsite_occupation_entries_checked': 2048,
            'original_repair_certificate': original,
            'exact_identity_for_any_real_weights': '(sum_x f_x N_x) P=(sum_x f_x) P',
            'centered_load_on_every_protected_spin_state': '0',
            'dressed_bound': 'norm(W*X_f W-c_f I)<=4 sum_x abs(f_x) (tau/U)^2',
            'fixed_local_weight_bound_on_source_trajectory': '4 sum_x abs(f_x) / E^6',
            'same_total_load_detector_can_locate_protected_spin_excitations': False,
            'twelve_port_native_chart_identified_with_cubic_site_chart': False}


def dressed_load_bound(weight_l1, edges):
    if (type(weight_l1) is not F or weight_l1 < 0
            or type(edges) is not int or edges < 2):
        raise ValueError('nonnegative rational total absolute weight and at least two edges required')
    return 4*weight_l1/F(edges**6)


def finite_translation_residual(shift):
    """Dimension-independent lower bound from the trace of a finite matrix."""
    if type(shift) is not F:
        raise ValueError('an exact rational coordinate shift is required')
    return abs(shift)
