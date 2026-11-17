"""Exact exchange-character lifts of the complete binary native instrument.

The positive lift is the existing C/E instrument.  The negative lift is a
different, explicitly supplied operation: its moving columns are the
number-conserving CAR bilinears on Lambda^2(C^12).  Both retain every native
record.  Their classical histories agree, while coherent source outputs
can distinguish them.  Neither construction selects physical statistics.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product

from ..clock_defect_energy import defect_history as topology
from ..common_source import slow_geometry
from ..observer import adaptive_record_recovery as adaptive
from ..observer import native, record_recovery as records
from ..observer.source_discrimination import circuit


STATES = native.states(binary=True)
PAIRS = tuple(tuple(i for i, n in enumerate(state) if n) for state in STATES)
INDEX = {pair: i for i, pair in enumerate(PAIRS)}
SEAMS = records.seams()
WITNESS_SEAM = (0, 5)
WITNESS_INPUT_PAIRS = ((0, 1), (0, 6))


def _character(value):
    if type(value) is not int or value not in (-1, 1):
        raise ValueError("an integer exchange character +1 or -1 required")
    return value


def _port(value):
    if type(value) is not int or not 0 <= value < 12:
        raise ValueError("an integer native port from 0 through 11 required")
    return value


def _pair(value):
    if (type(value) is not tuple or len(value) != 2
            or any(type(p) is not int for p in value) or value not in INDEX):
        raise ValueError("a sorted pair of distinct integer native ports required")
    return value


def _seam(value):
    _pair(value)
    if value not in SEAMS:
        raise ValueError("an actual native seam required")
    return value


def _label(law, seam, coin, raw):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("native law C or E required")
    _seam(seam)
    if type(coin) is not int or coin not in (0, 1):
        raise ValueError("an integer native coin required")
    if type(raw) is not int or not 0 <= raw < (4 if law == "C" else 8):
        raise ValueError("a complete native raw record address required")
    return law, seam, coin, raw


def sign(source_pair, target_pair, character=-1):
    """Character of the actual ordered one-particle lift, with identity holds."""
    source_pair, target_pair = _pair(source_pair), _pair(target_pair)
    character = _character(character)
    if source_pair == target_pair:
        return 1
    removed, added = set(source_pair)-set(target_pair), set(target_pair)-set(source_pair)
    if len(removed) != 1 or len(added) != 1:
        raise ValueError("one collision-free native movement required")
    a, b = next(iter(removed)), next(iter(added))
    _seam(tuple(sorted((a, b))))
    lifted = tuple(b if p == a else p for p in source_pair)
    return character if lifted != target_pair else 1


def phase_on_output(seam, current, pair, character=-1):
    """Output phase extension; identity on holds and outside the moving image.

    This is a full H66 diagonal extension for a record-selected source word.
    Padding and collision addresses require their separately declared identity
    extension in the full-register compiler.
    """
    a, b = _seam(seam)
    pair, character = _pair(pair), _character(character)
    if type(current) is not int or current not in (-1, 0, 1):
        raise ValueError("a binary current -1, 0, or 1 required")
    if current == 0:
        return 1
    origin, destination = (a, b) if current == 1 else (b, a)
    if destination not in pair or origin in pair:
        return 1
    spectator = next(p for p in pair if p != destination)
    return character if a < spectator < b else 1


@dataclass(frozen=True, slots=True)
class SignedBranch:
    law: str
    seam: tuple[int, int]
    coin: int
    raw: int
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]
    abort: bool


def branch(law, seam, coin, raw, character=-1):
    """One full raw slot; its displayed integer map has amplitude 1/sqrt(2)."""
    _label(law, seam, coin, raw)
    character = _character(character)
    rows = []
    a, b = seam
    for i, state in enumerate(STATES):
        out, (_, current) = native.native(state, seam, coin)
        record = {0: 0, 1: 1, -1: 2}[current]
        if law == "E":
            record += 4*((state[a]+state[b]) % 2)
        if record == raw:
            pair = tuple(p for p, n in enumerate(out) if n)
            rows.append((i, INDEX[pair], sign(PAIRS[i], pair, character)))
    inputs, outputs, phases = tuple(zip(*rows)) if rows else ((), (), ())
    return SignedBranch(law, seam, coin, raw, inputs, outputs, phases,
                        raw not in circuit.record_values(law, coin))


def branches(law, character=-1):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("native law C or E required")
    _character(character)
    return tuple(branch(law, seam, coin, raw, character)
                 for seam in SEAMS for coin in (0, 1)
                 for raw in range(4 if law == "C" else 8))


@dataclass(frozen=True, slots=True)
class SignedHistory:
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    signs: tuple[int, ...]
    squared_amplitude: F


def history_map(history, character=-1):
    """Exact composition without summing amplitudes across committed flags."""
    if type(history) is not tuple:
        raise ValueError("a tuple of complete native record labels required")
    _character(character)
    image = {i: (i, 1) for i in range(66)}
    for label in history:
        if type(label) is not tuple or len(label) != 4:
            raise ValueError("each history entry needs law, seam, coin and raw")
        step = branch(*label, character)
        lookup = {i: (j, s) for i, j, s in zip(step.inputs, step.outputs, step.signs, strict=True)}
        image = {i: (lookup[j][0], s*lookup[j][1])
                 for i, (j, s) in image.items() if j in lookup}
    inputs = tuple(image)
    return SignedHistory(inputs, tuple(image[i][0] for i in inputs),
                         tuple(image[i][1] for i in inputs), F(1, 2**len(history)))


def car_bilinear(destination, origin):
    """Restriction of c_destination^dagger c_origin to Lambda^2(C^12).

    Creation or annihilation separately leaves H66.  Their signed product
    below is number conserving and acts on this sector without such a claim.
    """
    destination, origin = _port(destination), _port(origin)
    output = []
    for pair in PAIRS:
        if origin not in pair:
            output.append((-1, 0))
            continue
        place = pair.index(origin)
        remaining = tuple(p for p in pair if p != origin)
        if destination in remaining:
            output.append((-1, 0))
            continue
        phase = (-1)**(place+sum(p < destination for p in remaining))
        output.append((INDEX[tuple(sorted((*remaining, destination)))], phase))
    return tuple(output)


def car_certificate():
    """Exact u(12) bilinear algebra and the actual directed native currents."""
    maps = {(a, b): car_bilinear(a, b) for a in range(12) for b in range(12)}
    checks = 0
    for a, b, c, d in product(range(12), repeat=4):
        ab, cd = maps[a, b], maps[c, d]
        for i in range(66):
            left, right = Counter(), Counter()
            j, s = cd[i]
            if s and ab[j][1]:
                left[ab[j][0]] += s*ab[j][1]
            j, s = ab[i]
            if s and cd[j][1]:
                left[cd[j][0]] -= s*cd[j][1]
            if b == c:
                j, s = maps[a, d][i]
                right[j] += s
            if a == d:
                j, s = maps[c, b][i]
                right[j] -= s
            if any(left[k] != right[k] for k in left.keys() | right.keys()):
                raise ArithmeticError("exterior-square bilinear commutator failed")
            checks += 1
    moving = 0
    for law in ("C", "E"):
        for row in branches(law):
            a, b = row.seam
            for i, j, phase in zip(row.inputs, row.outputs, row.signs, strict=True):
                current = STATES[i][a]-STATES[j][a]
                if current:
                    target, expected = maps[b, a][i] if current == 1 else maps[a, b][i]
                    if (j, phase) != (target, expected):
                        raise ArithmeticError("native movement is not the claimed CAR bilinear")
                    moving += 1
                elif phase != 1:
                    raise ArithmeticError("the declared no-movement operator is not identity")
                if any(STATES[j][p]-STATES[i][p] != current*(int(p == b)-int(p == a))
                       for p in range(12)):
                    raise ArithmeticError("port-number current continuity failed")
    return {"one_particle_modes": 12, "fixed_number": 2, "sector_dimension": 66,
            "bilinear_commutator_columns": checks, "moving_CAR_columns": moving,
            "identity": "[E_ab,E_cd]=delta_bc E_ad-delta_ad E_cb",
            "port_number_current_continuity": True,
            "standalone_creation_annihilation_CAR_on_H66": False,
            "physical_fermion_statistics_selected": False}


def complete_certificate():
    """Bind both character families to the native maps and check exact effects.

    Counts below describe one complete raw grid.  Mixed-sign counts refer
    to the negative character; the positive family is separately enumerated
    and must have exactly +1 in every supported column.
    """
    originals = {(adaptive.actions()[r.action][0], adaptive.actions()[r.action][1], r.coin, r.raw): r
                 for r in adaptive.branch_maps(binary=True)}
    result = {}
    for law in ("C", "E"):
        rows = branches(law)
        expected_labels = {(law, seam, coin, raw) for seam in SEAMS for coin in (0, 1)
                           for raw in range(4 if law == "C" else 8)}
        labels = [(row.law, row.seam, row.coin, row.raw) for row in rows]
        if len(labels) != len(expected_labels) or set(labels) != expected_labels:
            raise ArithmeticError("complete unique raw branch grid required, including zero and abort slots")
        endpoint = {label: [] for label in expected_labels}
        for seam in SEAMS:
            a, b = seam
            for coin in (0, 1):
                for i, state in enumerate(STATES):
                    total = state[a]+state[b]
                    target = list(state)
                    target[a], target[b] = (total+coin)//2, (total+1-coin)//2
                    current = state[a]-target[a]
                    raw = {0: 0, 1: 1, -1: 2}[current]+(4*(total % 2) if law == "E" else 0)
                    endpoint[law, seam, coin, raw].append((i, INDEX[tuple(p for p, n in enumerate(target) if n)]))
        counts = {seam: [0]*66 for seam in SEAMS}
        mixed, zeros, aborts, columns = 0, 0, 0, 0
        for row in rows:
            _label(row.law, row.seam, row.coin, row.raw)
            admitted = ({0, 1, 2} if law == "C" else
                        {0, 1, 2, 4, 5} if row.coin == 0 else {0, 1, 2, 4, 6})
            if type(row.abort) is not bool or row.abort != (row.raw not in admitted):
                raise ArithmeticError("native abort metadata disagrees with the complete raw menu")
            original = originals[law, row.seam, row.coin, row.raw]
            if ((row.inputs, row.outputs) != (original.inputs, original.outputs)
                    or tuple(zip(row.inputs, row.outputs)) != tuple(endpoint[law, row.seam, row.coin, row.raw])):
                raise ArithmeticError("signed lift changed the native branch support")
            if (len(row.inputs) != len(row.outputs) or len(row.inputs) != len(row.signs)
                    or len(set(row.outputs)) != len(row.outputs)
                    or any(type(i) is not int or not 0 <= i < 66 for i in (*row.inputs, *row.outputs))
                    or any(type(s) is not int or s not in (-1, 1) for s in row.signs)):
                raise ArithmeticError("signed branch is not a unit-modulus partial injection")
            if row.abort and row.inputs:
                raise ArithmeticError("native abort address acquired support")
            for i, j, phase in zip(row.inputs, row.outputs, row.signs, strict=True):
                counts[row.seam][i] += 1
                current = STATES[i][row.seam[0]]-STATES[j][row.seam[0]]
                if phase != phase_on_output(row.seam, current, PAIRS[j]):
                    raise ArithmeticError("compiled output phase extension disagrees with the lift")
            mixed += len(set(row.signs)) > 1
            zeros += not row.inputs
            aborts += row.abort
            columns += len(row.inputs)
        if any(value != 2 for values in counts.values() for value in values):
            raise ArithmeticError("complete fair-coin instrument is not trace preserving")
        positive = branches(law, 1)
        positive_labels = [(row.law, row.seam, row.coin, row.raw) for row in positive]
        if len(positive_labels) != len(expected_labels) or set(positive_labels) != expected_labels:
            raise ArithmeticError("positive lift requires the complete unique native raw branch grid")
        for row in positive:
            _label(row.law, row.seam, row.coin, row.raw)
            admitted = ({0, 1, 2} if law == "C" else
                        {0, 1, 2, 4, 5} if row.coin == 0 else {0, 1, 2, 4, 6})
            original = originals[law, row.seam, row.coin, row.raw]
            if (type(row.abort) is not bool or row.abort != (row.raw not in admitted)
                    or (row.inputs, row.outputs) != (original.inputs, original.outputs)
                    or tuple(zip(row.inputs, row.outputs)) != tuple(endpoint[law, row.seam, row.coin, row.raw])
                    or len(row.signs) != len(row.inputs)
                    or any(type(i) is not int or not 0 <= i < 66 for i in (*row.inputs, *row.outputs))
                    or any(type(phase) is not int or phase != 1 for phase in row.signs)):
                raise ArithmeticError("positive lift differs from the complete native instrument")
        result[law] = {"raw_slots": len(rows), "zero_slots": zeros, "abort_slots": aborts,
                       "input_columns": columns, "mixed_sign_slots": mixed,
                       "counts_are_per_character": True,
                       "mixed_sign_count_character": -1,
                       "characters_checked": (1, -1),
                       "positive_lift_matches_complete_native_instrument": True,
                       "every_KstarK_unchanged": True, "complete_trace_preservation": True}
    return result


def rotation_action(permutation, character=-1):
    """Ordinary unordered/exterior-square port action in the sorted-pair chart."""
    if (type(permutation) is not tuple or len(permutation) != 12
            or any(type(p) is not int for p in permutation)
            or set(permutation) != set(range(12))):
        raise ValueError("a permutation of the twelve integer ports required")
    _character(character)
    images, signs = [], []
    for a, b in PAIRS:
        lifted = (permutation[a], permutation[b])
        images.append(INDEX[tuple(sorted(lifted))])
        signs.append(character if lifted[0] > lifted[1] else 1)
    return tuple(images), tuple(signs)


def covariance_certificate():
    rotations = tuple(tuple(map(int, p)) for p in slow_geometry.source().rotations())
    rows = tuple(r for law in ("C", "E") for r in branches(law))
    lookup = {(r.law, r.seam, r.coin, r.raw): r for r in rows}
    checks = 0
    actions = {p: rotation_action(p) for p in rotations}
    for p in rotations:
        images, phases = actions[p]
        for row in rows:
            endpoints = tuple(p[a] for a in row.seam)
            flipped = endpoints[0] > endpoints[1]
            raw = row.raw
            if flipped:
                raw = 4*(raw//4)+{0: 0, 1: 2, 2: 1, 3: 3}[raw % 4]
            target = lookup[row.law, tuple(sorted(endpoints)), 1-row.coin if flipped else row.coin, raw]
            if target.abort != row.abort:
                raise ArithmeticError("proper rotation changed native abort classification")
            translated = {i: (j, s) for i, j, s in zip(target.inputs, target.outputs, target.signs, strict=True)}
            if len(translated) != len(row.inputs):
                raise ArithmeticError("proper rotation changed a complete branch domain")
            for i, j, s in zip(row.inputs, row.outputs, row.signs, strict=True):
                mapped, target_sign = translated[images[i]]
                if mapped != images[j] or s*phases[j] != phases[i]*target_sign:
                    raise ArithmeticError("exterior-square rotation covariance failed")
                checks += 1
    for p, q in product(rotations, repeat=2):
        combined = tuple(p[q[i]] for i in range(12))
        ci, cs = actions[combined]
        pi, ps = actions[p]
        qi, qs = actions[q]
        if any((pi[qi[i]], ps[qi[i]]*qs[i]) != (ci[i], cs[i]) for i in range(66)):
            raise ArithmeticError("rotation lift ceased to be an ordinary representation")
    pairs, edges, cells, _, _, _ = topology.geometry()
    if pairs != PAIRS:
        raise ArithmeticError("native topology and quantum source charts differ")
    for _, vertices in cells:
        phase = 1
        for i, j in zip(vertices, vertices[1:]+vertices[:1], strict=True):
            phase *= sign(PAIRS[i], PAIRS[j])
        if phase != 1:
            raise ArithmeticError("negative lift is not flat on an admitted local cell")
    _, loop = topology.exchange_loop()
    phase = 1
    for i, j in zip(loop, loop[1:]):
        phase *= sign(PAIRS[i], PAIRS[j])
    if phase != -1:
        raise ArithmeticError("exchange character lost its negative loop")
    return {"proper_rotations": len(rotations), "branch_covariance_columns": checks,
            "ordinary_group_product_columns": len(rotations)**2*66,
            "flat_local_cells": len(cells), "movement_edges": len(edges),
            "negative_exchange_loop": [list(PAIRS[i]) for i in loop],
            "negative_exchange_holonomy": -1, "positive_exchange_holonomy": 1,
            "chart_rule": "s'_xy=g_y s_xy g_x; K'_h=G K_h G*",
            "rotation_lift_is_projective_spin_half": False,
            "occupation_preserving_diagonal_gauge_equivalence": False,
            "gauge_obstruction": "A diagonal endpoint gauge has unit product around every loop; the exchange product is -1."}


def gauge_obstruction_certificate():
    """Exclude even diagonal equivalence up to an arbitrary scalar per flag.

    A spectator loop and an exchange loop use the same three directed moving
    branch labels, in different orders.  Endpoint phases cancel on each loop;
    scalar branch phases have equal products.  Their lift signs are opposite.
    Thus record-dependent scalar Kraus gauges cannot remove this obstruction.
    """
    a, b, c = topology.exchange_loop()[0]
    d = next(i for i in range(12) if i not in (a, b, c))
    single = tuple(tuple(sorted(p)) for p in ((a, d), (b, d), (c, d), (a, d)))
    exchange = tuple(tuple(sorted(p)) for p in ((a, b), (a, c), (b, c), (a, b)))
    loops = []
    for pairs in (single, exchange):
        phase, labels = 1, {"C": [], "E": []}
        for source, target in zip(pairs, pairs[1:]):
            origin, destination = next(iter(set(source)-set(target))), next(iter(set(target)-set(source)))
            seam = tuple(sorted((origin, destination)))
            coin = int(origin > destination)
            raw = 1 if origin < destination else 2
            phase *= sign(source, target)
            for law in ("C", "E"):
                labels[law].append((seam, coin, raw+(4 if law == "E" else 0)))
        loops.append({"pairs": [list(p) for p in pairs], "sign": phase,
                      "directed_record_labels": labels})
    if (loops[0]["sign"], loops[1]["sign"]) != (1, -1):
        raise ArithmeticError("spectator and exchange loops lost their opposite characters")
    for law in ("C", "E"):
        if Counter(loops[0]["directed_record_labels"][law]) != Counter(loops[1]["directed_record_labels"][law]):
            raise ArithmeticError("the scalar-branch gauge witness uses different labels")
    return {"witnesses": loops, "same_multiset_of_directed_branch_labels": True,
            "diagonal_occupation_preserving_unitary_plus_scalar_branch_gauge_possible": False,
            "arbitrary_unitary_ignoring_port_occupations_excluded": False}


def channel_distance_certificate():
    """Sharp half diamond norm for the *committed complete* one-seam channel.

    Only moving blocks differ. Their total effect is half the projector onto
    exactly one occupied endpoint. For any source/reference state, triangle
    bounds each block difference by twice its probability, so half the full
    trace norm is at most 1/2. Two opposite-sign spectator columns attain it.
    Scalar block signs cancel and give zero distance instead.
    """
    rows = []
    for seam in SEAMS:
        a, b = seam
        inside = tuple(p for p in range(12) if a < p < b)
        outside = tuple(p for p in range(12) if p not in seam and p not in inside)
        distance = F(1, 2) if inside and outside else F()
        laws = {}
        for law in ("C", "E"):
            moving_effect = [0]*66
            mixed = 0
            for coin in (0, 1):
                for raw in range(4 if law == "C" else 8):
                    row = branch(law, seam, coin, raw)
                    mixed += len(set(row.signs)) > 1
                    for i, j in zip(row.inputs, row.outputs, strict=True):
                        moving_effect[i] += i != j
            if any(count != int(STATES[i][a]+STATES[i][b] == 1)
                   for i, count in enumerate(moving_effect)):
                raise ArithmeticError("moving-branch effect is not the claimed half projector")
            if (mixed != 0) != bool(distance):
                raise ArithmeticError("mixed signs do not match the sharp distance witness")
            laws[law] = {"half_diamond_distance": str(distance), "mixed_sign_branches": mixed}
        rows.append({"seam": list(seam), "inside_spectators": len(inside),
                     "outside_spectators": len(outside), "laws": laws,
                     "attaining_input_pairs": [list(tuple(sorted((a, p)))) for p in (inside[0], outside[0])]
                     if distance else None})
    return {"metric": "one half diamond norm of complete source-plus-classical-record channels",
            "identification": "fixed sorted-pair amplitude chart, with the same input state for both lifts",
            "relative_chart_changes_require_transporting_preparation_and_analyzer": True,
            "arbitrary_untouched_reference_in_upper_bound": True,
            "moving_effect": "sum_moving K_h* K_h = P_(n_a+n_b=1)/2",
            "readout_requires_coherent_source_access": True, "rows": rows}


def witness():
    inputs = tuple(INDEX[p] for p in WITNESS_INPUT_PAIRS)
    rows = []
    for law in ("C", "E"):
        for coin in (0, 1):
            for raw in range(4 if law == "C" else 8):
                negative = branch(law, WITNESS_SEAM, coin, raw)
                lookup = {i: (j, s) for i, j, s in zip(negative.inputs, negative.outputs, negative.signs, strict=True)}
                selected = [lookup[i] for i in inputs if i in lookup]
                if selected and len(selected) != 2:
                    raise ArithmeticError("witness components obtained different native flags")
                probability = F(len(selected), 4)
                plus = F(sum(s for _, s in selected)**2, 8) if selected else F()
                minus = probability-plus
                rows.append({"law": law, "coin": coin, "raw": raw,
                             "abort": negative.abort, "native_probability": str(probability),
                             "output_pairs": [list(PAIRS[j]) for j, _ in selected],
                             "negative_amplitude_signs": [s for _, s in selected],
                             "positive_read_plus_probability": str(probability),
                             "positive_read_minus_probability": "0",
                             "negative_read_plus_probability": str(plus),
                             "negative_read_minus_probability": str(minus)})
    for law in ("C", "E"):
        selected = [row for row in rows if row["law"] == law]
        if (sum(F(row["native_probability"]) for row in selected) != 1
                or sum(F(row["negative_read_minus_probability"]) for row in selected) != F(1, 2)):
            raise ArithmeticError("complete witness readout changed")
    return {"seam": list(WITNESS_SEAM), "input_pairs": [list(p) for p in WITNESS_INPUT_PAIRS],
            "input": "equal positive coherent superposition of the two displayed pairs",
            "read": "record-selected inverse movement, then X on the two input addresses",
            "all_raw_outputs": rows, "complete_readout_total_variation": "1/2",
            "equal_prior_single_attempt_discrimination_success": "3/4",
            "postselection_used": False,
            "negative_controls": {
                "dephase_input": "both lifts give plus and minus probability 1/2",
                "scalar_sign_per_committed_branch": "all complete quantum outputs unchanged",
                "read_only_native_records": "both laws and both lifts have identical distributions",
                "restore_ordered_pair_chart_in_preparation_and_read": "all probabilities unchanged"}}


def certificate():
    return {"scope": "EXACT_FINITE_EXCHANGE_CHARACTER_COMPARISON",
            "positive_lift_is_existing_native_instrument": True,
            "negative_lift_changes_quantum_operation": True,
            "complete_native_maps": complete_certificate(),
            "number_conserving_CAR": car_certificate(),
            "rotation_and_chart_covariance": covariance_certificate(),
            "occupation_and_record_gauge_obstruction": gauge_obstruction_certificate(),
            "complete_channel_distances": channel_distance_certificate(),
            "one_attempt_witness": witness(),
            "all_finite_history_identity": {
                "map": "K_h=sqrt(w_h) sum_(x in D_h) chi_h(x)|f_h(x)><x|, with injective f_h",
                "same_history_effect": "K_h* K_h=w_h P_Dh for both characters",
                "same_classical_history_and_terminal_population_laws": True,
                "same_history_plus_untouched_reference_marginal": True,
                "arbitrary_source_reference_inputs": True,
                "record_driven_classical_feedback_allowed": True,
                "intermediate_coherent_quantum_controls_in_this_identity": False,
                "complete_quantum_source_outputs_identical": False},
            "physical_spin_or_species_identified": False,
            "native_action_selects_negative_character": False,
            "mass_energy_or_physical_clock_identified": False}


report = certificate
