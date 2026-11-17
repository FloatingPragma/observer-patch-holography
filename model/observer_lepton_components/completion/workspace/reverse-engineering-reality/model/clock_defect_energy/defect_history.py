"""A native exchange-history class, with an explicit finite deformation domain.

The binary configuration complex has H_1(Z)=Z/2. Its sign is derived by
lifting native paths to ordered trajectories, not inserted as a quantum
phase. Native repair does not thereby carry a protected matter charge.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

import sympy as sp

from . import joins_contract as contract
from ..observer import native, native_loop_readout, record_recovery as records
from ..observer import adaptive_record_recovery as adaptive
from ..observer.source_discrimination import circuit
from ..native_modular_reference import history_transport

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "defect_history_results.json"
# An integral certificate, not a physics parameter. These actual cell columns
# give a determinant-two minor in the spanning-tree cycle coordinates.
INTEGRAL_MINOR = (393, *range(1, 185), 187, 188, 189, 194, 201, 202, 203, 204,
                  205, 206, 207, 208, 211, 222, 223, 224, 225, 226, 227, 228,
                  229, 232, 245, 246, 247, 248, 249, 250, 253, 266, 267, 268,
                  269, 270, 271, 321, 322, 323, 324, 328, 338, 339, 340, 341,
                  355, 360, 367, 368, 379, 380)


@lru_cache(maxsize=1)
def configurations():
    return tuple(tuple(i for i, value in enumerate(state) if value)
                 for state in native.states(binary=True))


@lru_cache(maxsize=1)
def geometry():
    """Compatibility/deformations defined from incidence before native updates."""
    pairs = configurations()
    index = {pair: i for i, pair in enumerate(pairs)}
    seams, faces = native_loop_readout.carrier()
    if tuple(seams) != records.seams():
        raise ArithmeticError("source graph and face complex disagree")
    edges = set()
    for i, pair in enumerate(pairs):
        for a, b in seams:
            if (a in pair) != (b in pair):
                out = tuple(sorted(b if p == a else a if p == b else p for p in pair))
                edges.add(tuple(sorted((i, index[out]))))
    edges = tuple(sorted(edges))
    def at(a, b):
        return index[tuple(sorted((a, b)))]
    cells = []
    for face in faces:
        for spectator in range(12):
            if spectator not in face:
                cells.append(("one_particle_face", tuple(at(spectator, p) for p in face)))
    for (a, b), (c, d) in combinations(seams, 2):
        if len({a, b, c, d}) == 4:
            cells.append(("disjoint_square", (at(a, c), at(b, c), at(b, d), at(a, d))))
    reached, tree = {0}, []
    while len(reached) < len(pairs):
        previous = len(reached)
        for i, (a, b) in enumerate(edges):
            if (a in reached) != (b in reached):
                reached.update((a, b))
                tree.append(i)
        if len(reached) == previous:
            raise ArithmeticError("binary configuration graph disconnected")
    chords = tuple(i for i in range(len(edges)) if i not in tree)
    if (len(pairs), len(edges), len(cells), len(tree), len(chords)) != (66, 300, 495, 65, 235):
        raise ArithmeticError("native configuration complex census changed")
    return pairs, edges, tuple(cells), tuple(tree), chords, tuple(faces)


def _vertex(value):
    if type(value) is not int or not 0 <= value < 66:
        raise ValueError("an original configuration index from0 through65 required")
    return value


def edge_parity(source, target):
    """Change of sorted pair chart for the actual one-particle lift, in Z/2."""
    source, target = _vertex(source), _vertex(target)
    if source == target:
        return 0
    pairs = configurations()
    left, right = pairs[source], pairs[target]
    removed, added = set(left) - set(right), set(right) - set(left)
    if len(removed) != 1 or len(added) != 1:
        raise ValueError("one native particle move required")
    a, b = next(iter(removed)), next(iter(added))
    if tuple(sorted((a, b))) not in records.seams():
        raise ValueError("particle move must use a native seam")
    lifted = tuple(b if p == a else p for p in left)
    return int(lifted != right)


def transport(path):
    if type(path) is not tuple or not path:
        raise ValueError("a nonempty tuple of native configuration indices required")
    for vertex in path:
        _vertex(vertex)
    parity = 0
    for a, b in zip(path, path[1:]):
        parity ^= edge_parity(a, b)
    return parity


def boundary(vertices, edges):
    """Integer oriented chain of a closed polygon; canonical edges point up."""
    index = {edge: i for i, edge in enumerate(edges)}
    result = [0] * len(edges)
    for a, b in zip(vertices, vertices[1:] + vertices[:1], strict=True):
        if a != b:
            result[index[tuple(sorted((a, b)))]] += 1 if a < b else -1
    return tuple(result)


def gf2_rank(rows):
    pivots = {}
    for row in rows:
        value = sum((entry & 1) << i for i, entry in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def exchange_loop():
    """Three native moves exchange two particles around one native triangle."""
    pairs, _, _, _, _, faces = geometry()
    a, b, c = faces[0]
    vertices = tuple(pairs.index(tuple(sorted(pair))) for pair in ((a, b), (b, c), (a, c)))
    return faces[0], vertices + vertices[:1]


def topology_certificate():
    pairs, edges, cells, tree, chords, _ = geometry()
    chains = tuple(boundary(vertices, edges) for _, vertices in cells)
    signs = tuple(edge_parity(a, b) for a, b in edges)
    if any(sum(bit * coefficient for bit, coefficient in zip(signs, chain)) % 2 for chain in chains):
        raise ArithmeticError("path-lift parity is not flat on an admitted local cell")
    if any(edge_parity(a, b) != edge_parity(b, a) for a, b in edges):
        raise ArithmeticError("path reversal is not group inverse")
    face, loop = exchange_loop()
    loop_chain = boundary(loop[:-1], edges)
    if transport(loop) != 1:
        raise ArithmeticError("local exchange lost odd holonomy")
    rank = gf2_rank(chains)
    if rank != 234 or len(INTEGRAL_MINOR) != 235 or len(set(INTEGRAL_MINOR)) != 235:
        raise ArithmeticError("local-cell rank or integral certificate changed")
    minor = sp.Matrix([[chains[c][edge] for c in INTEGRAL_MINOR] for edge in chords])
    determinant = int(minor.det(method="domain-ge"))
    if abs(determinant) != 2:
        raise ArithmeticError("integral boundary minor no longer has index two")
    # This constructive two-chain certifies twice the exchange is a boundary.
    coefficients = minor.inv(method="DM") * sp.Matrix([2 * loop_chain[i] for i in chords])
    if any(value.q != 1 for value in coefficients):
        raise ArithmeticError("twice the exchange has no integral filling certificate")
    filling = tuple((cell, int(value)) for cell, value in zip(INTEGRAL_MINOR, coefficients) if value)
    if any(sum(value * chains[cell][i] for cell, value in filling) != 2 * loop_chain[i]
           for i in range(len(edges))):
        raise ArithmeticError("integral exchange filling fails on actual native edges")
    augmented = gf2_rank((*chains, loop_chain))
    if augmented != 235:
        raise ArithmeticError("erasing the exchange did not change the deformation class")
    return {
        "vertices": len(pairs), "movement_edges": len(edges),
        "local_cells": dict(Counter(kind for kind, _ in cells)),
        "spanning_tree_edges": [list(edges[i]) for i in tree],
        "cycle_lattice_rank": len(chords), "local_boundary_rank_mod2": rank,
        "integral_minor_cell_indices": INTEGRAL_MINOR, "integral_minor_determinant": determinant,
        "integral_first_homology": "Z/2", "first_homology_mod2_dimension": 1,
        "double_exchange_integral_filling": filling,
        "all_local_cell_holonomies": 1,
        "exchange_face": face, "exchange_vertices": loop,
        "exchange_holonomy": -1,
        "loop_occupied_pairs": [pairs[i] for i in loop],
        "negative_control_filling_exchange_as_extra_cell": {
            "boundary_rank_mod2": augmented, "homology_mod2_dimension": 0,
            "changes_declared_collision_free_deformation_contract": True},
        "homotopy_group_computed": False,
        "chain_data_sha256": hashlib.sha256(json.dumps([pairs, edges, cells, chains], separators=(",", ":")).encode()).hexdigest(),
    }


def native_branch_certificate():
    """Both complete quantum grammars have precisely the declared move graph."""
    states = native.states(binary=True)
    _, edges, _, _, _, _ = geometry()
    outputs = {}
    for law in ("C", "E"):
        reached, normalization = set(), [[0] * 66 for _ in records.seams()]
        selected = [b for b in adaptive.branch_maps(binary=True) if adaptive.actions()[b.action][0] == law]
        for branch in selected:
            _, seam = adaptive.actions()[branch.action]
            slot = records.seams().index(seam)
            for i, j in zip(branch.inputs, branch.outputs, strict=True):
                normalization[slot][i] += 1
                value = circuit.record_values(law, branch.coin)[branch.raw]
                current = records.record_current(law, (*seam, branch.coin, value))
                if current != states[i][seam[0]] - states[j][seam[0]]:
                    raise ArithmeticError("native record does not reconstruct this move")
                if i != j:
                    reached.add(tuple(sorted((i, j))))
                    edge_parity(i, j)
        if reached != set(edges) or any(n != 2 for row in normalization for n in row):
            raise ArithmeticError("complete native maps do not instantiate the configuration graph")
        outputs[law] = {"raw_maps": len(selected), "zero_raw_maps": sum(not b.inputs for b in selected),
                        "input_columns": sum(len(b.inputs) for b in selected),
                        "movement_edges": len(reached), "normalization_per_seam_input": 2}
    return outputs


def _lift_step(ordered, seam, current):
    if current == 0:
        return ordered
    if current not in (-1, 1):
        raise ValueError("binary current required")
    source, target = seam if current == 1 else seam[::-1]
    if source not in ordered or target in ordered:
        raise ValueError("current is incompatible with the lifted binary path")
    return tuple(target if p == source else p for p in ordered)


def local_protocol(law):
    """Every coin word, complete record kernel and a phase-erasure control."""
    law = records.candidate(law)
    face, loop = exchange_loop()
    pairs = configurations()
    initial = pairs[loop[0]]
    schedule = tuple(tuple(sorted(set(pairs[a]) ^ set(pairs[b]))) for a, b in zip(loop, loop[1:]))
    rows = []
    for coins in product((0, 1), repeat=3):
        state = tuple(int(i in initial) for i in range(12))
        ordered, tape, path = initial, [], [loop[0]]
        for seam, coin in zip(schedule, coins, strict=True):
            state, record = records.recorded_step(law, state, seam, coin)
            current = records.record_current(law, record)
            ordered = _lift_step(ordered, seam, current)
            pair = tuple(i for i, value in enumerate(state) if value)
            if tuple(sorted(ordered)) != pair:
                raise ArithmeticError("record-only ordered lift disagrees with native source")
            path.append(pairs.index(pair))
            tape.append(record)
        closed = path[-1] == path[0]
        parity = transport(tuple(path))
        if parity != int(ordered != tuple(sorted(ordered))):
            raise ArithmeticError("edge-chart transport and direct worldline lift disagree")
        rows.append({"coins": coins, "complete_records": tape, "vertices": path,
                     "ordered_endpoint": ordered, "closed": closed,
                     "closed_exchange_parity": parity if closed else None,
                     "probability": F(1, 8)})
    branches = history_transport.word_instrument(law, schedule)
    if any(sum(b.weight for b in branches if i in b.inputs) != 1 for i in range(66)):
        raise ArithmeticError("three-attempt complete quantum instrument lost mass")
    # Multiplication of one committed Kraus block by a scalar sign cancels on
    # every matrix unit. This is a countercontrol, not an admitted phase gate.
    odd_tapes = {tuple(tuple(r) for r in row["complete_records"]) for row in rows
                 if row["closed_exchange_parity"] == 1}
    phased_blocks = [b for b in branches if b.history in odd_tapes]
    if not phased_blocks:
        raise ArithmeticError("positive odd history absent from actual complete instrument")
    return {"initial_pair": initial, "supporting_native_face": face, "schedule": schedule,
            "coin_words": rows,
            "closed_even_mass": sum(row["probability"] for row in rows if row["closed_exchange_parity"] == 0),
            "closed_odd_mass": sum(row["probability"] for row in rows if row["closed_exchange_parity"] == 1),
            "open_mass": sum(row["probability"] for row in rows if not row["closed"]),
            "all_paths_remain_inside_face": all(set(pairs[i]) <= set(face) for row in rows for i in row["vertices"]),
            "complete_raw_record_slots": (2 * (1 << circuit.RECORD_BITS[law])) ** 3,
            "nonzero_full_input_quantum_blocks": len(branches),
            "matrix_unit_branch_checks": sum(len(b.inputs) ** 2 for b in branches),
            "native_branch_amplitudes_have_exchange_sign": False,
            "record_conditioned_scalar_sign_changes_committed_instrument": False,
            "signed_block_control": {"blocks_with_sign_minus_one": len(phased_blocks),
                                     "identity": "(-K_h) rho (-K_h)* = K_h rho K_h* on every input/reference"},
            "fixed_local_schedule_is_persistently_exploring": False}


def exploring_protocol(law, epsilon=F(1, 16)):
    """All seams retained; exact three-attempt diagonal history marginal."""
    law = records.candidate(law)
    epsilon = adaptive._fraction(epsilon, positive=True)
    local = local_protocol(law)
    initial, face, schedule = tuple(local["initial_pair"]), set(local["supporting_native_face"]), local["schedule"]
    paths = {(initial, False): F(1)}
    for preferred in schedule:
        following = defaultdict(F)
        for (ordered, left_face), weight in paths.items():
            state = tuple(int(p in ordered) for p in range(12))
            for seam in records.seams():
                probability = epsilon / 30 + (1 - epsilon if seam == preferred else 0)
                for coin in (0, 1):
                    _, record = records.recorded_step(law, state, seam, coin)
                    after = _lift_step(ordered, seam, records.record_current(law, record))
                    following[after, left_face or not set(after) <= face] += weight * probability / 2
        paths = following
    masses = defaultdict(F)
    for (ordered, left_face), weight in paths.items():
        if set(ordered) == set(initial):
            label = "closed_even" if ordered == initial else "closed_odd"
        else:
            label = "open"
        masses[label] += weight
        if left_face:
            masses["ever_left_face"] += weight
        elif label == "closed_odd":
            masses["localized_closed_odd"] += weight
    if sum(paths.values()) != 1 or sum(masses[k] for k in ("closed_even", "closed_odd", "open")) != 1:
        raise ArithmeticError("exploring loop experiment dropped branch mass")
    selected_word_lower = ((epsilon / 30 + 1 - epsilon) / 2) ** 3
    if masses["localized_closed_odd"] < selected_word_lower:
        raise ArithmeticError("retained local exchange word missing from exploring history")
    return {"epsilon": epsilon, "attempts": 3, "masses": dict(masses),
            "specified_three_move_word_probability": selected_word_lower,
            "positive_compressed_states": len(paths),
            "terminal_ordered_pair_and_exit_flag_counts_are_offline": True,
            "all_raw_native_outcomes_defined_by_complete_instrument": True,
            "quantum_histories_not_replaced_by_diagonal_marginal": True}


def signed_gram():
    pairs, edges, _, _, _, _ = geometry()
    incidence = sp.zeros(len(edges), len(pairs))
    for e, (a, b) in enumerate(edges):
        incidence[e, a] = -1
        incidence[e, b] = 1 if edge_parity(a, b) == 0 else -1
    return incidence.T * incidence


def positive_leading_minors(matrix):
    """Fraction-free Sylvester certificate, using exact integer divisions."""
    work = [list(map(int, row)) for row in matrix.tolist()]
    size, previous, minors = len(work), 1, []
    for k in range(size):
        pivot = work[k][k]
        if pivot <= 0:
            raise ArithmeticError("matrix is not certified positive definite")
        minors.append(pivot)
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                numerator = pivot * work[i][j] - work[i][k] * work[k][j]
                if numerator % previous:
                    raise ArithmeticError("fraction-free determinant division failed")
                work[i][j] = numerator // previous
        for i in range(k + 1, size):
            work[i][k] = 0
        previous = pivot
    return tuple(minors)


def signed_compatibility_control():
    """Instantiate SignedGap on derived labels, without assigning a native action."""
    pairs, edges, _, _, _, _ = geometry()
    gram = signed_gram()
    if any(gram[i, i] > 10 or sum(abs(gram[i, j]) for j in range(66)) > 20 for i in range(66)):
        raise ArithmeticError("signed configuration degree bound changed")
    determinant = int(gram.det(method="domain-ge"))
    if determinant <= 0:
        raise ArithmeticError("negative-cycle signed Gram lost injectivity")
    unsigned = sp.zeros(len(edges), len(pairs))
    for e, (a, b) in enumerate(edges):
        unsigned[e, a], unsigned[e, b] = -1, 1
    if unsigned * sp.ones(66, 1) != sp.zeros(len(edges), 1):
        raise ArithmeticError("unsigned compatibility control lost its constant section")
    minors = positive_leading_minors(gram - 2 * sp.eye(66))
    return {"dimension": 66, "signed_gram_integer_determinant": determinant,
            "kernel_dimension": 0, "maximum_configuration_degree": 10,
            "signed_gram_minus_2I_positive_leading_minors": minors,
            "certified_eigenvalue_interval": "2 < lambda <= 20",
            "unsigned_native_population_form_has_constant_zero_mode": True,
            "signed_gram_is_the_native_repair_generator": False,
            "native_dynamics_minimizes_this_signed_cost": False,
            "energy_or_mass_identification": False}


def parity_transfer_counts(law):
    """Unsigned endpoint and signed record-character kernels, both /60.

    The signed kernel weights each classical trajectory by its ordered-lift
    chart parity. It is not a replacement Kraus map on native quantum states.
    """
    law = records.candidate(law)
    pairs, states = configurations(), native.states(binary=True)
    index = {pair: i for i, pair in enumerate(pairs)}
    unsigned, signed = [[0] * 66 for _ in pairs], [[0] * 66 for _ in pairs]
    for i, (pair, state) in enumerate(zip(pairs, states, strict=True)):
        for seam in records.seams():
            for coin in (0, 1):
                out, record = records.recorded_step(law, state, seam, coin)
                after = _lift_step(pair, seam, records.record_current(law, record))
                target = tuple(p for p, value in enumerate(out) if value)
                if tuple(sorted(after)) != target:
                    raise ArithmeticError("ordered lift left its native endpoint")
                j = index[target]
                unsigned[i][j] += 1
                signed[i][j] += 1 if after == target else -1
    return tuple(map(tuple, unsigned)), tuple(map(tuple, signed))


def particle_readiness():
    unsigned, signed = parity_transfer_counts("C")
    if (unsigned, signed) != parity_transfer_counts("E"):
        raise ArithmeticError("C/E population character kernels differ")
    gram = signed_gram()
    if sp.Matrix(signed) != 60 * sp.eye(66) - gram:
        raise ArithmeticError("native record character does not match signed incidence")
    cluster = {i for i, pair in enumerate(configurations()) if pair in records.seams()}
    exit_counts = {sum(unsigned[i][j] for j in range(66) if j not in cluster) for i in cluster}
    if exit_counts != {4}:
        raise ArithmeticError("native localized pair escape changed")
    # The positive gap and degree bound imply ||I-L_sign/60||_2 <=29/30.
    positive_leading_minors(gram - 2 * sp.eye(66))
    _, loop = exchange_loop()
    start = loop[0]
    endpoint, character = [0] * 66, [0] * 66
    endpoint[start] = character[start] = 1
    denominator, rows = 1, []
    horizons = (0, 1, 2, 3, 8, 16, 32, 64)
    for n in range(max(horizons) + 1):
        if n in horizons:
            even = F(endpoint[start] + character[start], 2 * denominator)
            odd = F(endpoint[start] - character[start], 2 * denominator)
            squared_norm = F(sum(v * v for v in character), denominator * denominator)
            if (min(even, odd) < 0 or sum(endpoint) != denominator
                    or squared_norm > F(29, 30) ** (2 * n)):
                raise ArithmeticError("complete parity-return law or contraction failed")
            rows.append({"attempts": n, "closed_even_mass": even, "closed_odd_mass": odd,
                         "open_mass": 1 - even - odd,
                         "closed_return_character": F(character[start], denominator),
                         "endpoint_character_l2_squared": squared_norm,
                         "l2_squared_upper": F(29, 30) ** (2 * n)})
        if n < max(horizons):
            endpoint = [sum(endpoint[i] * unsigned[i][j] for i in range(66)) for j in range(66)]
            character = [sum(character[i] * signed[i][j] for i in range(66)) for j in range(66)]
            denominator *= 60
    return {"uniform_complete_record_return_laws": rows,
            "signed_transfer_identity": "S=I-L_sign/60; S^n[i,j]=P(j,even|i)-P(j,odd|i)",
            "signed_transfer_norm_upper": F(29, 30),
            "parameter": "uniform native repair count; no physical duration",
            "signed_transfer_is_native_quantum_source_channel": False,
            "source_endpoint_determines_history_parity": False,
            "stored_complete_classical_history_is_erased_by_this_bound": False,
            "native_Kraus_exchange_word_amplitude_in_fixed_record_chart": "+sqrt(1/8)",
            "relative_phase_between_committed_record_blocks_is_observable": False,
            "localized_cluster_mean_first_exit_uniform_attempts": F(60, next(iter(exit_counts))),
            "missing_particle_join": [
                "a live localized sector and reference-like exterior, not only a completed path class",
                "a native phase-sensitive representation and complete readout if exchange statistics are claimed",
                "protection and transport under the full admitted dynamics, with any active services counted",
                "same-process rest action, kinetic response and stress; neither signed gap nor exit rate is mass"],
            "actionable_control": "compare unchanged positive C/E Kraus branches with any proposed signed-amplitude transport on the same complete quantum records; record-dependent scalar phase alone cannot change an instrument"}


def theorem_maps():
    return [
        {"source": "Lean/Screen/HolonomyInterference.lean",
         "theorems": ["transport_append", "transport_reverse", "transportRatio_eq_closedLoopHolonomy",
                      "holonomy_rechart_invariant", "relativeCharacterPhase_eq_closedLoopPhase"],
         "actual_objects": "V=66 unordered configurations; G={+1,-1}; edge labels from ordered native path lift, extended by+1 on nonedges for the theorem's total edgeLabel; physical paths use native edges and stutters only",
         "checked_premises": "complete native movement graph; reversal-compatible signs; all495 admitted cells flat; explicit odd loop",
         "conclusion": "composable, reversal-compatible, rechart-invariant exchange holonomy and its mathematical character",
         "remaining": "no coherent character service, interference law, spatial defect charge or fermionic statistics is supplied"},
        {"source": "Lean/OphGap/SignedGap.lean",
         "theorems": ["kernel_trivial_of_negative_cycle", "gram_posDef_of_negative_cycle"],
         "actual_objects": "300-by66 signed incidence on the native configuration graph",
         "checked_premises": "loop-free edges; signs plus/minus one;65-edge spanning tree; three-edge negative cycle",
         "conclusion": "the separately defined signed compatibility Gram is positive definite",
         "remaining": "native C/E amplitudes are positive partial injections; no minimization law, clock, action or mass for this signed form is derived"},
        {"source": "model/observer/record_recovery.py",
         "theorems": ["recorded_step", "record_current", "restore_load"],
         "actual_objects": "both complete C/E raw records, fixed-face loop and positive-exploration histories",
         "checked_premises": "binary total-two preparation, all coins/outcomes, known initial configuration and retained tape",
         "conclusion": "native current records reconstruct the ordered path and its closed exchange parity",
         "remaining": "reconstruction is a bounded offline calculation, not an extra hidden particle label or quantum phase memory"},
    ]


def _json(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def report():
    topology = topology_certificate()
    local = {law: local_protocol(law) for law in ("C", "E")}
    exploring = {law: exploring_protocol(law) for law in ("C", "E")}
    if exploring["C"]["masses"] != exploring["E"]["masses"]:
        raise ArithmeticError("C/E diagonal loop statistics differ")
    paths = ["model/clock_defect_energy/defect_history.py", "model/clock_defect_energy/test_defect_history.py",
             "model/clock_defect_energy/DEFECT_HISTORY.md",
             "model/clock_defect_energy/joins_contract.py",
             "model/clock_defect_energy/contract.py", "tools/strict_json.py",
             "model/observer/native.py", "model/observer/record_recovery.py",
             "model/observer/adaptive_record_recovery.py", "model/observer/native_loop_readout.py",
             "model/observer/source_discrimination/circuit.py", "model/native_modular_reference/history_transport.py",
             "Lean/Screen/HolonomyInterference.lean", "Lean/OphGap/SignedGap.lean",
             "model/common_source/slow_geometry.py", "evidence/native_geometric_source_20260925/spec.json"]
    return _json({"schema": "oph.finite_theorem_joins.defect_history.v1", "campaign": contract.binding(),
                  "topology": topology, "native_instruments": native_branch_certificate(),
                  "local_protocols": local, "exploring_protocols": exploring,
                  "signed_compatibility_control": signed_compatibility_control(),
                  "particle_readiness": particle_readiness(),
                  "theorem_maps": theorem_maps(),
                  "resources": {"native_preparation": "declared known binary pair; no new preparation service executed",
                                "local_protocol_attempts": 3,
                                "local_positive_coin_histories_per_law": 8,
                                "uniform_record_character_horizon": 64,
                                "offline_ordered_pair_addresses": 132,
                                "offline_address_bits_upper": 8,
                                "extra_physical_registers": 0,
                                "online_path_tracker_implemented": False,
                                "complete_tape_not_replaced_by_a_quantum_parity_qubit": True},
                  "scope": {"positive_native_history_class": True,
                            "full_fundamental_group_computed": False,
                            "native_autonomous_source_charge": False,
                            "new_physical_label_register": False,
                            "native_fermionic_phase_or_statistics": False,
                            "protected_localized_matter_or_mass": False,
                            "integrated_CP03_accepted": False,
                            "physical_comparison_count_added": 0},
                  "source_sha256": contract.source_bindings(paths)})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", nargs="?", const=str(OUTPUT))
    args = parser.parse_args()
    result = report()
    if args.check:
        contract.check_receipt(Path(args.check), result)
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print("native history topology: exact Z/2 and complete-record checks passed")


if __name__ == "__main__":
    main()
