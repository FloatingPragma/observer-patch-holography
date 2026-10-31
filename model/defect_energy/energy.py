"""Exact clock tests on the same native defect as defect_protection.

No Hamiltonian is supplied.  The finite instrument fixes a history-dependent
memory gate, which can be tested against a proposed count-only energy clock.
All arithmetic in the certificates is integer or rational.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from ..defect_protection import protection
from ..defect_protection import service
from ..exchange_sector import statistics as exchange


INITIAL_PAIR = (0, 11)
ODD_WORD = ((0, 5), (0, 11), (5, 11))
EVEN_WORD = ((0, 1), (1, 5), (0, 5))
MATCHED_EVEN_WORD = ((5, 11), (1, 5), (0, 5), (5, 11), (0, 11), (1, 5), (5, 11))
MATCHED_ODD_WORD = ((0, 5), (0, 11), (5, 11), (5, 11), (1, 5), (1, 5), (5, 11))


def _law(law):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("native law C or E required")
    return law


@lru_cache(maxsize=2)
def _rows(law):
    return protection.branches(_law(law))


def word_instrument(law, initial_pair, word):
    """All coin/record outcomes on a supplied source basis preparation.

    The initial memory may be arbitrary, including entanglement with an
    untouched reference.  Each output applies I or X to that memory, with
    its full native record retained.  No successful branch is renormalized.
    The length-eight bound caps enumeration at 256 outcomes.
    """
    _law(law)
    if (type(initial_pair) is not tuple or len(initial_pair) != 2
            or any(type(p) is not int for p in initial_pair)
            or initial_pair not in exchange.INDEX):
        raise ValueError("a sorted binary-source pair required")
    if (type(word) is not tuple or len(word) > 8
            or any(type(edge) is not tuple or len(edge) != 2
                   or any(type(p) is not int for p in edge)
                   or edge not in exchange.SEAMS for edge in word)):
        raise ValueError("at most eight canonical native seams required")
    rows = _rows(law)
    result = []
    for coins in product((0, 1), repeat=len(word)):
        index, parity, movements, toggles = exchange.INDEX[initial_pair], 0, 0, 0
        path, history = [initial_pair], []
        for seam, coin in zip(word, coins, strict=True):
            matches = [row for row in rows if row.seam == seam
                       and row.coin == coin and index in row.inputs]
            if len(matches) != 1:
                raise ArithmeticError("complete raw menu lost or duplicated a basis outcome")
            row = matches[0]
            position = row.inputs.index(index)
            target, sign = row.outputs[position], row.signs[position]
            lifted = dict(protection.lifted_pairs(row))
            flip = int(sign == -1)
            if any(lifted[r*66+index] != (r ^ flip)*66+target for r in (0, 1)):
                raise ArithmeticError("history differs from the declared memory lift")
            if row.abort:
                raise ArithmeticError("an abort acquired positive support on the binary source")
            parity ^= flip
            toggles += flip
            movements += index != target
            index = target
            history.append((law, seam, coin, row.raw))
            path.append(exchange.PAIRS[index])
        result.append({"coins": coins, "history": tuple(history), "path": tuple(path),
                       "probability": F(1, 2**len(word)), "movements": movements,
                       "orientation_toggles": toggles,
                       "memory_flip": parity, "final_pair": exchange.PAIRS[index]})
    if sum((row["probability"] for row in result), F()) != 1:
        raise ArithmeticError("complete word outcomes are not normalized")
    return tuple(result)


@lru_cache(maxsize=2)
def transfer_counts(law):
    """All 132 basis transitions, including every zero-current hold.

    A uniform independent seam and fair coin give denominator 60.  Empty
    raw and abort slots are visited in the supplied complete branch menu.
    """
    rows = [Counter() for _ in range(132)]
    movement = [0]*132
    for branch in _rows(law):
        for source, target in protection.lifted_pairs(branch):
            rows[source][target] += 1
            movement[source] += source % 66 != target % 66
    if any(sum(row.values()) != 60 for row in rows):
        raise ArithmeticError("uniform native transfer lost probability")
    if any(sum(row.get(j, 0) for row in rows) != 60 for j in range(132)):
        raise ArithmeticError("ordered-cover reference is not stationary")
    return tuple(tuple(sorted(row.items())) for row in rows), tuple(movement)


def memory_channel(law, attempts):
    """Reduced memory channel for I66/66 tensor rho, without source resets.

    The result is (1-p) rho + p X rho X, also for rho with a reference.
    This exact transfer computation retains source-memory correlations;
    iterating the one-step reduced memory channel would discard them.
    """
    _law(law)
    if type(attempts) is not int or not 0 <= attempts <= 64:
        raise ValueError("an integer attempt count from zero through 64 required")
    rows, movement = transfer_counts(law)
    distribution = [int(i < 66) for i in range(132)]
    work = F()
    for n in range(attempts):
        work += F(sum(p*count for p, count in zip(distribution, movement, strict=True)),
                  66*60**(n+1))
        following = [0]*132
        for probability, row in zip(distribution, rows, strict=True):
            for target, count in row:
                following[target] += probability*count
        distribution = following
    denominator = 66*60**attempts
    if sum(distribution) != denominator:
        raise ArithmeticError("continuing-carrier memory law is not normalized")
    flip = F(sum(distribution[66:]), denominator)
    return {"attempts": attempts, "memory_flip_probability": flip,
            "X_eigenvalue": F(1), "Y_Z_eigenvalue": 1-2*flip,
            "expected_native_movements": work}


def _serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item) for item in value]
    return value


def clock_certificate(law):
    """An equal-count closed-loop obstruction on the same initial source."""
    even = word_instrument(law, INITIAL_PAIR, EVEN_WORD)
    odd = word_instrument(law, INITIAL_PAIR, ODD_WORD)
    witnesses = []
    for name, instrument, expected in (("even", even, 0), ("odd", odd, 1)):
        moving = [row for row in instrument if row["movements"] == 3]
        if (len(moving) != 1 or moving[0]["final_pair"] != INITIAL_PAIR
                or moving[0]["memory_flip"] != expected):
            raise ArithmeticError("equal-count closed-loop witness changed")
        witnesses.append({"name": name, **moving[0]})
    # Gate identity and X are distinct even after arbitrary branch-global
    # phases: acting on memory |0> gives orthogonal terminal readouts.
    return _serialize({"initial_pair": INITIAL_PAIR,
                       "complete_even_word": even, "complete_odd_word": odd,
                       "selected_nonrenormalized_witnesses": witnesses,
                       "same_attempts": 3, "same_movements": 3,
                       "same_returned_source": True,
                       "memory_Z_readout_on_initial_zero": [0, 1],
                       "witness_probability_each_fixed_word": F(1, 8),
                       "witness_probability_each_uniform_seam_word": F(1, 60**3),
                       "count_only_common_memory_unitary_exists": False,
                       "closed_loop_parity_invariant_under_vertex_rechart": True})


@lru_cache(maxsize=60)
def _resource_rows(law, seam):
    return service._resource_rows(law, seam)


def history_resources(history):
    """Count the actual declared compiler menu along one retained history.

    These are word inventories, not another physical replay of the expanded
    pulses.  Preparation, storage duration and terminal read are separate.
    """
    totals = Counter()
    for law, seam, coin, raw in history:
        matches = [row for row in _resource_rows(law, seam)
                   if row["coin"] == coin and row["raw"] == raw]
        if len(matches) != 1:
            raise ValueError("history contains no unique complete resource row")
        row = matches[0]
        for key, value in row.items():
            if key not in ("coin", "raw", "abort"):
                totals[key] += value
        totals["native_attempts"] += 1
        totals["fresh_coin_bits"] += 1
        totals["retained_native_tape_bits"] += 3 if law == "C" else 4
    return dict(sorted(totals.items()))


def matched_cost_certificate(law):
    """Equal complete branch multisets, different ordered memory action.

    Every scalar sum c(law,seam,coin,raw) is therefore equal, without
    selecting or fitting c.  The actual compiler counts are one example.
    The memory's relative action is instead pi times its state-dependent
    orientation-toggle count, modulo 2*pi.
    """
    instruments = [word_instrument(law, INITIAL_PAIR, word)
                   for word in (MATCHED_EVEN_WORD, MATCHED_ODD_WORD)]
    chosen = []
    complete = []
    for parity, instrument in enumerate(instruments):
        moving = [row for row in instrument if row["movements"] == 7]
        if (len(moving) != 1 or moving[0]["final_pair"] != INITIAL_PAIR
                or moving[0]["memory_flip"] != parity):
            raise ArithmeticError("matched-cost closed history changed")
        chosen.append(moving[0])
        # Every outcome remains present, including all zero-current holds.
        complete.append([{**row, "resources": history_resources(row["history"])}
                         for row in instrument])
    labels = [Counter(row["history"]) for row in chosen]
    costs = [history_resources(row["history"]) for row in chosen]
    if labels[0] != labels[1] or costs[0] != costs[1]:
        raise ArithmeticError("matched native branch resources differ")
    if [row["orientation_toggles"] for row in chosen] != [2, 1]:
        raise ArithmeticError("state-dependent memory action changed")
    return _serialize({"complete_word_instruments": complete,
                       "selected_nonrenormalized_witnesses": chosen,
                       "complete_branch_multisets_equal": True,
                       "all_additive_raw_branch_costs_equal": True,
                       "common_compiled_resources": costs[0],
                       "state_dependent_orientation_toggle_counts": [2, 1],
                       "relative_memory_action_in_pi_units": [2, 1],
                       "relative_memory_action_modulo_2pi_in_pi_units": [0, 1],
                       "single_scalar_raw_branch_cost_determines_action": False,
                       "initial_source_preparation_counted": False,
                       "persistent_storage_duration_counted": False,
                       "terminal_computational_memory_read_bits": 1,
                       "terminal_Z_read_basis_change_pulses": 0,
                       "expanded_pulse_word_executed_in_this_packet": False})


def action_certificate(law):
    """Check the constructive repair-action identity on every raw column.

    P_minus=(I-X)/2 and exp(-i*pi*p*P_minus)=I-2*p*P_minus
    for p in {0,1}.  The latter has only exact integer entries.  Applying
    it after the source partial injection gives the supplied lifted Kraus
    map, still with amplitude 1/sqrt(2), on every input and reference.
    """
    checked, nontrivial = 0, 0
    for row in _rows(law):
        actual = dict(protection.lifted_pairs(row))
        for source, target, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
            parity = int(sign == -1)
            # Exact exponential of the projector: [[1-p,p],[p,1-p]].
            exponential = ((1-parity, parity), (parity, 1-parity))
            for r in (0, 1):
                support = [s*66+target for s in (0, 1) if exponential[s][r]]
                if support != [actual[r*66+source]]:
                    raise ArithmeticError("operator-valued repair action differs from the complete map")
            checked += 1
            nontrivial += parity
    return {"raw_slots_including_empty_and_abort": len(_rows(law)),
            "source_columns_checked": checked, "nonzero_parity_columns": nontrivial,
            "memory_projector": "P_minus=(I-X)/2",
            "on_each_source_input_column": "exp(-i*pi*p(x,h)*P_minus)=X^p(x,h)",
            "complete_branch": "K_h=exp(-i*A_h)*(M_h tensor I)/sqrt(2)",
            "output_supported_action": "A_h=pi*sum_x p(x,h)|f_h(x)><f_h(x)| tensor P_minus",
            "M_h": "unnormalized source partial injection, zero outside its input support",
            "source_superpositions_and_arbitrary_reference_included": True,
            "history_relative_action": "pi*sum_t p(x_t,h_t), modulo 2*pi",
            "energy_from_action_requires_duration_and_logarithm_choice": True}


def report():
    """Finite exact consequences, with scope boundaries kept in the receipt."""
    laws = {}
    for law in ("C", "E"):
        channels = [memory_channel(law, n) for n in (0, 1, 2, 3, 6, 12)]
        one, two = channels[1], channels[2]
        if any(row["expected_native_movements"] != F(5*row["attempts"], 33)
               for row in channels):
            raise ArithmeticError("stationary native movement rate changed")
        flip = one["memory_flip_probability"]
        if flip != F(26, 495):
            raise ArithmeticError("one-step coherent memory response changed")
        laws[law] = {"closed_clock_witness": clock_certificate(law),
                     "complete_branch_action_identity": action_certificate(law),
                     "matched_resource_action_witness": matched_cost_certificate(law),
                     "continuing_carrier_memory_channels": _serialize(channels),
                     "one_step_distance_from_every_unitary_lower_bound": str(flip),
                     "two_step_Z_multiplier_minus_one_step_multiplier_squared": str(
                         two["Y_Z_eigenvalue"]-one["Y_Z_eigenvalue"]**2),
                     "expected_native_movements_per_attempt": "5/33",
                     "arbitrary_memory_and_untouched_reference": True}
    return {"schema": "oph.defect_energy.results.v1", "laws": laws,
            "same_carrier": "H66 tensor one persistent coherent orientation qubit",
            "same_memory_law": "model.defect_protection.protection.lifted_pairs",
            "preparation": "fixed basis source for word instruments; I66/66 for reduced channels",
            "scope": "FINITE_COMPONENT_DIAGNOSTICS",
            "conclusions": {
                "persistent_observable": "I66 tensor X; not a localized source defect",
                "action_selected_by_specification": "ordered-path memory holonomy X^parity",
                "cost_to_relative_action_identity": "A_minus-A_plus = pi times logical orientation toggles, modulo 2*pi",
                "scalar_branch_cost_obstruction": "equal complete raw-branch multisets can have different relative action",
                "maintenance_count_selected": "native movements, not joules or stored energy",
                "one_step_memory_dynamics": "nonunitary after the continuing source is traced out",
                "memory_alone_is_a_discrete_semigroup": False,
                "count_only_rest_gap_excluded_for_full_history_law": True,
                "history_dependent_or_enlarged_dynamics_excluded": False,
                "Hamiltonian_part_plus_noise_excluded": False,
                "physical_energy_mass_or_time_selected": False},
            "physical_comparisons_added": 0, "integrated_gate_accepted": False}
