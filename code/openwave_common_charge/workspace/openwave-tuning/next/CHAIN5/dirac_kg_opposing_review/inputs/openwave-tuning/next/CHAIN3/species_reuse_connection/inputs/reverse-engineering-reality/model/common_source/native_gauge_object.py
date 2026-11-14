"""A protected native object with its actual icosahedral current carrier.

The thirty native edges are forbidden supports of any nontrivial code that
corrects every complete recorded C/E branch. The remaining thirty-six pairs
admit two icosahedrally equivariant perfect-matching colorings by the six
antipodal axes. Each gives a protected six-state object with actual port
action 1+5. Its sum-zero five-state part is protected as well. Two such
objects contain a canonical geometric 3+3' sector under simultaneous native
port rotations. No Standard Model current matrix is inserted into this code.

Encoding, record-dependent recovery and coherent native services remain
implementation choices. A geometric carrier does not itself select its
continuous response grammar, spin lift, matter exterior algebra or masses.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path

import sympy as s

from tools import strict_json
from . import slow_geometry
from ..observer import native
from ..observer.source_discrimination import circuit
from ..repair_scheduler import coins


ROOT = Path(__file__).resolve().parents[2]


def _clean(matrix):
    matrix = s.SparseMatrix(matrix)
    entries = {key: s.expand(s.radsimp(value)) for key, value in matrix.todok().items()}
    return s.ImmutableSparseMatrix(matrix.rows, matrix.cols,
        {key: value for key, value in entries.items() if value != 0})


def _zero(matrix):
    return not _clean(matrix).todok()


def _entries(matrix):
    matrix = _clean(matrix)
    return {"shape": list(matrix.shape),
            "entries": [[i, j, str(value)] for (i, j), value in sorted(matrix.todok().items())]}


def permutation_order(permutation):
    p = tuple(permutation)
    if sorted(p) != list(range(len(p))):
        raise ValueError("a complete permutation is required")
    seen, order = set(), 1
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j); length += 1; j = p[j]
        order = int(s.ilcm(order, length))
    return order


@lru_cache(maxsize=1)
def geometry():
    source = slow_geometry.source()
    frame = s.ImmutableMatrix([[s.Rational(a)+s.Rational(b)*s.sqrt(5) for a, b in row]
                               for row in slow_geometry.raw_frame()])
    pairs = tuple(combinations(range(12), 2))
    seams = tuple(sorted(tuple(sorted(edge)) for edge in source.seams()))
    axes = tuple(pair for pair in pairs if _zero(frame[pair[0], :]+frame[pair[1], :]))
    nonedges = tuple(pair for pair in pairs if pair not in seams)
    far = tuple(pair for pair in nonedges if pair not in axes)
    rotations = tuple(tuple(p) for p in source.rotations())
    if (len(pairs), len(seams), len(axes), len(far), len(rotations)) != (66, 30, 6, 30, 60):
        raise ArithmeticError("native pair or proper-rotation geometry changed")
    states = tuple(native.states(binary=True))
    if {tuple(i for i, n in enumerate(state) if n) for state in states} != set(pairs):
        raise ArithmeticError("native binary load-two state space differs from the pair carrier")
    return {"frame": frame, "pairs": pairs, "seams": seams, "axes": axes,
            "nonedges": nonedges, "far_pairs": far, "rotations": rotations}


def axis_action(permutation):
    g = tuple(permutation)
    if g not in geometry()["rotations"]:
        raise ValueError("a native proper rotation is required")
    index = {pair: i for i, pair in enumerate(geometry()["axes"])}
    return tuple(index[tuple(sorted((g[i], g[j])))] for i, j in geometry()["axes"])


def pair_action(permutation):
    g = tuple(permutation)
    if g not in geometry()["rotations"]:
        raise ValueError("a native proper rotation is required")
    index = {pair: i for i, pair in enumerate(geometry()["pairs"])}
    return tuple(index[tuple(sorted((g[i], g[j])))] for i, j in geometry()["pairs"])


def permutation_matrix(permutation):
    p = tuple(permutation)
    if sorted(p) != list(range(len(p))):
        raise ValueError("a complete permutation is required")
    return s.ImmutableSparseMatrix(len(p), len(p), {(p[j], j): 1 for j in range(len(p))})


@lru_cache(maxsize=2)
def equivariant_matchings(choice=0):
    """The two equivariant source colorings, ordered before any physics target.

    An equivariant map on the transitive thirty-pair orbit is determined by
    the color of its first pair. That pair's order-two stabilizer fixes exactly
    two antipodal axes. The antipodal pairs use their own axis as color.
    Both choices obey the complete perfect-matching and covariance checks.
    """
    if type(choice) is not int or choice not in (0, 1):
        raise ValueError("one of the two exact native coloring choices is required")
    data = geometry()
    base = data["far_pairs"][0]
    stabilizer = [g for g in data["rotations"] if tuple(sorted(g[i] for i in base)) == base]
    fixed = [i for i in range(6) if all(axis_action(g)[i] == i for g in stabilizer)]
    if len(stabilizer) != 2 or len(fixed) != 2:
        raise ArithmeticError("native pair stabilizer does not give two axis colors")
    colors = {pair: i for i, pair in enumerate(data["axes"])}
    for g in data["rotations"]:
        pair = tuple(sorted(g[i] for i in base))
        color = axis_action(g)[fixed[choice]]
        if pair in colors and colors[pair] != color:
            raise ArithmeticError("equivariant coloring is not well-defined")
        colors[pair] = color
    matchings = tuple(tuple(sorted(pair for pair, color in colors.items() if color == i)) for i in range(6))
    if set(colors) != set(data["nonedges"]):
        raise ArithmeticError("coloring lost a native nonedge")
    for matching in matchings:
        if len(matching) != 6 or Counter(i for pair in matching for i in pair) != Counter(range(12)):
            raise ArithmeticError("axis color is not a perfect native matching")
    for g in data["rotations"]:
        action = axis_action(g)
        if any(colors[tuple(sorted((g[i], g[j])))] != action[color]
               for (i, j), color in colors.items()):
            raise ArithmeticError("native port action does not intertwine the color action")
    return matchings


def color_matrix(choice=0):
    """Integer 66-by-6 source incidence; normalized encoder is W/sqrt(6)."""
    index = {pair: i for i, pair in enumerate(geometry()["pairs"])}
    return s.ImmutableSparseMatrix(66, 6, {(index[pair], color): 1
        for color, matching in enumerate(equivariant_matchings(choice)) for pair in matching})


def matching_isometry(choice=0):
    return _clean(color_matrix(choice)/s.sqrt(6))


def five_projector():
    return s.ImmutableMatrix(s.eye(6)-s.ones(6)/6)


@lru_cache(maxsize=1)
def geometric_projectors():
    """Actual two-carrier geometric sectors in the 36 color-product labels.

    The two three-dimensional irreps each occur once in 5 tensor5. Their
    central projectors are (3/60)sum_g character(g)*[R6(g) tensor R6(g)].
    They are derived from actual simultaneous port rotations; their sum lies
    in the antisymmetric, sum-zero sector, with no supplied current matrices.
    """
    data = geometry(); B = data["frame"]
    norm = s.expand((B.T*B)[0, 0])
    p3, p3p = s.zeros(36), s.zeros(36)
    rows = []
    for g in data["rotations"]:
        inverse = tuple(g.index(i) for i in range(12))
        R = _clean(B.T*s.Matrix([list(B[i, :]) for i in inverse])/norm)
        char = s.expand(s.radsimp(s.trace(R))); conjugate = s.expand(char.xreplace({s.sqrt(5): -s.sqrt(5)}))
        A = permutation_matrix(axis_action(g))
        U = s.kronecker_product(A, A)
        p3 += char*U/20; p3p += conjugate*U/20
        rows.append({"port_action": list(g), "axis_action": list(axis_action(g)),
                     "order": permutation_order(g), "character3": str(char), "character3prime": str(conjugate),
                     "rotation3": _entries(R)})
    p3, p3p = _clean(p3), _clean(p3p)
    total = _clean(p3+p3p)
    swap = s.ImmutableSparseMatrix(36, 36, {(6*j+i, 6*i+j): 1 for i in range(6) for j in range(6)})
    q = s.kronecker_product(five_projector(), five_projector())
    if (not _zero(p3*p3-p3) or not _zero(p3p*p3p-p3p) or not _zero(p3*p3p)
            or s.trace(p3) != 3 or s.trace(p3p) != 3 or not _zero(q*total-total)
            or not _zero(swap*total+total)):
        raise ArithmeticError("native two-object geometric projector failed")
    return {"three": p3, "three_prime": p3p, "six": total, "rotation_rows": tuple(rows)}


@lru_cache(maxsize=1)
def geometric_intertwiners():
    """Explicit exact metric frames from actual geometric3 and3' to code6²."""
    rows = geometric_projectors()["rotation_rows"]
    result = {}
    for name, conjugate in (("three", False), ("three_prime", True)):
        actions = []
        for row in rows:
            record = row["rotation3"]
            R = s.SparseMatrix(*record["shape"], {(i, j): s.sympify(v) for i, j, v in record["entries"]})
            if conjugate:
                R = R.xreplace({s.sqrt(5): -s.sqrt(5)})
            axis = row["axis_action"]
            action = tuple(6*axis[i]+axis[j] for i in range(6) for j in range(6))
            actions.append((action, R))
        frame = None
        for index in range(36):
            for component in range(3):
                average = s.zeros(36, 3)
                # For X=e_index e_component^T, U X R^T has just one
                # nonzero row. Accumulate its three exact entries directly.
                for action, R in actions:
                    for j in range(3):
                        average[action[index], j] += R[j, component]/60
                average = _clean(average)
                if not _zero(average):
                    frame = average
                    break
            if frame is not None:
                break
        if frame is None:
            raise ArithmeticError("native geometric intertwiner vanished")
        gram = _clean(frame.T*frame); metric = gram[0, 0]
        if metric <= 0 or not _zero(gram-metric*s.eye(3)):
            raise ArithmeticError("native geometric frame has wrong metric")
        projector = geometric_projectors()[name]
        if not _zero(frame*frame.T/metric-projector):
            raise ArithmeticError("native geometric frame does not span its selected sector")
        for action, R in actions:
            inverse = tuple(action.index(i) for i in range(36))
            moved = frame.extract(inverse, range(3))
            if not _zero(moved-frame*R):
                raise ArithmeticError("actual simultaneous native rotation does not intertwine the geometric frame")
        result[name] = {"frame": frame, "metric": metric}
    if not _zero(result["three"]["frame"].T*result["three_prime"]["frame"]):
        raise ArithmeticError("native geometric sectors are not orthogonal")
    return result


def coloring_relation():
    W0, W1 = color_matrix(0), color_matrix(1)
    overlap = _clean(W0.T*W1/6)
    Q = five_projector()
    data = geometry()
    antipode = tuple(next(j for j in range(12) if _zero(data["frame"][i, :]+data["frame"][j, :]))
                     for i in range(12))
    pair_index = {pair: i for i, pair in enumerate(data["pairs"])}
    actual = tuple(pair_index[tuple(sorted((antipode[i], antipode[j])))] for i, j in data["pairs"])
    inversion = permutation_matrix(actual)
    if overlap != s.ones(6)/6 or not _zero(Q*overlap*Q) or not _zero(inversion*W0-W1):
        raise ArithmeticError("the two native orientation choices lost their actual relation")
    return {"encoder_overlap": _entries(overlap), "common_uniform_rank": 1,
            "protected_five_subspaces_are_orthogonal": True,
            "native_antipodal_inversion": list(antipode),
            "actual_inversion_exchanges_colorings_with_identity_axis_map": True,
            "orientation_exchange_is_a_physical_parity_or_CP_claim": False,
            "two_choices_exhaust_all_possible_correctable_codes": False,
            "two_choices_exhaust_equivariant_axis_colorings_with_own_axis_antipodal_color": True}


def _branch_maps(law, seam):
    """Every raw source record, including identically zero abort labels."""
    records = 1 << circuit.RECORD_BITS[law]
    result = {(coin, raw): {} for coin in (0, 1) for raw in range(records)}
    for coin, label, pairs in coins.coin_maps(law, seam):
        reverse = {value: raw for raw, value in circuit.record_values(law, coin).items()}
        if label not in reverse:
            raise ArithmeticError("the native source emitted an unrecognized record")
        result[coin, reverse[label]] = dict(pairs)
    return result


def recorded_protection():
    """Exact source K_r†K_r, with the recorded fair-coin factor retained."""
    states = tuple(native.states(binary=True))
    state_index = {tuple(i for i, value in enumerate(state) if value): j for j, state in enumerate(states)}
    profiles = [tuple(tuple(state_index[pair] for pair in matching) for matching in equivariant_matchings(choice))
                for choice in (0, 1)]
    branches, empty, addresses = 0, 0, 0
    probability_counts = {"C": Counter(), "E": Counter()}
    digest = hashlib.sha256()
    for law in ("C", "E"):
        for seam in geometry()["seams"]:
            normalization = F()
            for (coin, raw), mapping in _branch_maps(law, seam).items():
                branches += 1; empty += not mapping; addresses += len(mapping)
                if len(set(mapping.values())) != len(mapping):
                    raise ArithmeticError("a native branch is not an injection on its retained inputs")
                probabilities = []
                for profile in profiles:
                    counts = [sum(i in mapping for i in support) for support in profile]
                    if len(set(counts)) != 1:
                        raise ArithmeticError("recorded native branch distinguishes the new logical colors")
                    probabilities.append(F(counts[0], 12))
                if probabilities[0] != probabilities[1]:
                    raise ArithmeticError("orientation choice changes the complete logical record law")
                probability_counts[law][str(probabilities[0])] += 1
                normalization += probabilities[0]
                digest.update(json.dumps([law, list(seam), coin, raw, sorted(mapping.items())], separators=(",", ":")).encode())
            if normalization != 1:
                raise ArithmeticError("complete native instrument lost probability mass")
    return {"native_binary_rows": 66, "logical_dimension": 6, "colorings": 2,
            "seams": 30, "laws": ["C", "E"], "complete_raw_branches": branches,
            "identically_zero_raw_branches": empty, "input_coin_addresses": addresses,
            "complete_branch_map_sha256": digest.hexdigest(),
            "probability_counts_by_law": {law: dict(sorted(counts.items())) for law, counts in probability_counts.items()},
            "all_original_branch_Gram_compressions_are_scalar": True,
            "five_subspace_and_two_object_geometric_six_inherit_record_correctability": True,
            "correction_requires_retained_record_and_coherent_recovery": True,
            "unread_native_work_channel_preserves_the_code": False,
            "new_source_pulse_execution_claimed_by_this_exact_certificate": False}


def effect_span_certificate():
    """Explicit native-branch combinations yielding n_i and each edge projector.

    A_ab=2K_(coin0,current+1)^dagger K=n_a(1-n_b),
    B_ab=2K_(coin1,current-1)^dagger K=n_b(1-n_a).
    A-B=n_a-n_b; connectedness and sum n_i=2 reconstruct all n_i.
    (n_a+n_b-A-B)/2 is the rank-one projector onto occupied native edge ab.
    On a correctable code of rank>=2 this projector must compress to zero.
    """
    data = geometry(); seams, pairs = data["seams"], data["pairs"]
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    states = tuple(native.states(binary=True))
    native_to_pair = {i: pair_index[tuple(j for j, value in enumerate(state) if value)]
                      for i, state in enumerate(states)}
    basis = [s.ones(66, 1)]; labels = ["identity"]
    adjacency = [[] for _ in range(12)]
    for index, (a, b) in enumerate(seams):
        maps = _branch_maps("C", (a, b))
        first, second = [s.zeros(66, 1) for _ in range(2)]
        for i in maps[0, 1]: first[native_to_pair[i]] = 1
        for i in maps[1, 2]: second[native_to_pair[i]] = 1
        expected_a = s.Matrix([int(a in pair and b not in pair) for pair in pairs])
        expected_b = s.Matrix([int(b in pair and a not in pair) for pair in pairs])
        if first != expected_a or second != expected_b:
            raise ArithmeticError("native recorded effects differ from the endpoint occupancy identities")
        basis += [first, second]
        labels += [f"2*C({a},{b};coin0;raw1)^dagger*C", f"2*C({a},{b};coin1;raw2)^dagger*C"]
        adjacency[a].append((b, index)); adjacency[b].append((a, index))
    differences = {0: s.zeros(61, 1)}; queue = deque([0])
    while queue:
        a = queue.popleft()
        for b, index in sorted(adjacency[a]):
            if b in differences:
                continue
            delta = s.zeros(61, 1)
            sign = 1 if b == seams[index][0] else -1
            delta[1+2*index], delta[2+2*index] = sign, -sign
            differences[b] = differences[a]+delta
            queue.append(b)
    if len(differences) != 12:
        raise ArithmeticError("native seam graph is not connected")
    unit = s.zeros(61, 1); unit[0] = 2
    base = (unit-sum(differences.values(), s.zeros(61, 1)))/12
    ports = [base+differences[i] for i in range(12)]
    table = s.Matrix.hstack(*basis)
    for i, coefficient in enumerate(ports):
        if table*coefficient != s.Matrix([int(i in pair) for pair in pairs]):
            raise ArithmeticError("native one-port effect span witness failed")
    edge_coefficients = []
    for index, (a, b) in enumerate(seams):
        coefficient = (ports[a]+ports[b])/2
        coefficient[1+2*index] -= s.Rational(1, 2)
        coefficient[2+2*index] -= s.Rational(1, 2)
        expected = s.zeros(66, 1); expected[pair_index[a, b]] = 1
        if table*coefficient != expected:
            raise ArithmeticError("native edge-projector effect span witness failed")
        edge_coefficients.append(coefficient)
    return {"basis_labels": labels, "one_port_coefficients": _entries(s.Matrix.hstack(*ports)),
            "native_edge_projector_coefficients": _entries(s.Matrix.hstack(*edge_coefficients)),
            "all_12_port_and_30_edge_identities_exact": True,
            "rank_one_native_edge_compression_vanishes_for_correctable_rank_at_least_two": True,
            "domain": "binary native total load two, every retained C/E seam branch, rank>=2"}


def single_carrier_obstruction():
    """Why a second native object is required for this geometric6 route."""
    data = geometry(); rows = geometric_projectors()["rotation_rows"]
    B = data["frame"]
    pairs = data["nonedges"]
    witnesses = {}
    for name, conjugate in (("three", False), ("three_prime", True)):
        frame = B.xreplace({s.sqrt(5): -s.sqrt(5)}) if conjugate else B
        lifted = _clean(s.Matrix([list(frame[i, :]+frame[j, :]) for i, j in pairs]))
        gram = _clean(lifted.T*lifted)
        if gram != 40*s.eye(3):
            raise ArithmeticError("native nonedge vector-frame metric differs")
        C = _clean(lifted.T*s.diag(*[int(0 in pair) for pair in pairs])*lifted/40)
        defect = s.expand(s.trace((C-s.eye(3)/6)**2))
        expected = s.Rational(7, 240)+(1 if conjugate else -1)*s.sqrt(5)/80
        if defect != expected or s.simplify(defect <= 0) is not s.false:
            raise ArithmeticError("native vector irrep lacks its exact occupation obstruction")
        multiplicity = 0
        for row in rows:
            g = row["port_action"]
            fixed = sum(tuple(sorted((g[i], g[j]))) == (i, j) for i, j in pairs)
            multiplicity += fixed*s.sympify(row["character3prime" if conjugate else "character3"])/60
        multiplicity = s.simplify(multiplicity)
        if multiplicity != 1:
            raise ArithmeticError("native nonedge vector irrep multiplicity differs")
        witnesses[name] = {"unique_nonedge_copy": True, "raw_frame_metric": "40",
            "normalized_port0_compression": _entries(C), "positive_scalar_compression_defect_squared": str(defect)}
    return {"effect_span": effect_span_certificate(), "vector_irreps": witnesses,
            "single_carrier_A5_invariant_correctable_code_can_contain_geometric3_or3prime": False,
            "multiple_corrected_carriers_excluded": False,
            "different_error_menus_load_sectors_or_projective_actions_excluded": False}


def finite_response_comparison():
    """Compare the actual protected representation with the existing theorem input.

    Only the character of the retained response witness is read here. Its
    infinitesimal generators and response law are not copied into the code.
    """
    path = ROOT/"code/source_selection_model/response.json"
    packet = strict_json.load(path)
    if packet.get("schema") != "oph.source-selection-response.v1":
        raise ArithmeticError("retained finite gauge response schema differs")
    observed = Counter()
    for row in packet["closed_paths"]:
        value = 0
        for block in row["observed_response"]:
            for i in range(3):
                a, b, c, d = map(s.Rational, block[i][i])
                value += a+b*s.sqrt(5)+s.I*(c+d*s.sqrt(5))
        value = s.simplify(value)
        observed[permutation_order(row["port_action"]), str(value)] += 1
    expected = Counter()
    for row in geometric_projectors()["rotation_rows"]:
        expected[row["order"], str(s.simplify(s.sympify(row["character3"])+s.sympify(row["character3prime"])))] += 1
    if observed != expected:
        raise ArithmeticError("actual native geometric carrier differs from the finite response M6 type")
    return {"existing_response_source": "code/source_selection_model/response.json",
            "actual_geometric_character_matches_retained_M6": True,
            "order_character_counts": [{"order": order, "character": char, "count": count}
                                       for (order, char), count in sorted(observed.items())],
            "retained_infinitesimal_response_generators_inserted_into_native_code": False,
            "same_continuous_port_current_law_derived_here": False,
            "spin_half_lift_or_3_plus_2_internal_module_derived_here": False,
            "exterior_M15_matter_or_particle_masses_derived_here": False}


def report():
    data = geometry(); projectors = geometric_projectors(); frames = geometric_intertwiners()
    W = color_matrix(0); Q = five_projector()
    if W.T*W != 6*s.eye(6):
        raise ArithmeticError("native matching encoder is not an exact isometry")
    for i in range(12):
        occupation = s.diag(*[int(i in pair) for pair in data["pairs"]])
        if W.T*occupation*W != s.eye(6):
            raise ArithmeticError("native code has a logical one-port signal")
    for g in data["rotations"]:
        if not _zero(permutation_matrix(pair_action(g))*W-W*permutation_matrix(axis_action(g))):
            raise ArithmeticError("actual native rotation does not implement the protected axis action")
    return {"schema": "oph.native_gauge_object.v1", "source": {
        "binary_total_load": 2, "port_count": 12, "native_basis_pairs": [list(p) for p in data["pairs"]],
        "native_edges": [list(e) for e in data["seams"]], "antipodal_axes": [list(a) for a in data["axes"]],
        "empirical_mass_or_mixing_inputs": []},
        "protected_object": {"kind": "active record-correctable native matching code",
            "equivariant_matchings": [[[list(pair) for pair in matching] for matching in equivariant_matchings(choice)]
                                     for choice in (0, 1)],
            "encoder_formula": "V_choice[pair,axis]=1/sqrt(6) if pair has that axis color, else0",
            "full_logical_dimension": 6, "actual_rotation_representation": "1+5",
            "five_projector": _entries(Q), "all_12_occupation_compressions": "I6/6",
            "all_60_actual_port_rotation_intertwiners_exact": True,
            "orientation_relation": coloring_relation(), "record_protection": recorded_protection()},
        "two_object_geometry": {"protected_product_dimension": 36,
            "five_product_dimension": 25, "selected_geometric_dimension": 6,
            "projector_three": _entries(projectors["three"]),
            "projector_three_prime": _entries(projectors["three_prime"]),
            "projector_six": _entries(projectors["six"]),
            "intertwiner_frames": {name: {"frame": _entries(value["frame"]), "metric": str(value["metric"])}
                                   for name, value in frames.items()},
            "selected_six_is_antisymmetric_and_sum_zero_on_each_color_factor": True,
            "single_carrier_rank_six_limit_applies_to_this_two_carrier_subcode": False,
            "independent_recorded_errors_have_scalar_Gram_on_selected_six": True,
            "rotation_rows": list(projectors["rotation_rows"])},
        "one_carrier_boundary": single_carrier_obstruction(),
        "finite_gauge_theorem_attachment": finite_response_comparison(),
        "scope": {"same_native_persistent_object_and_geometric_M6_action_constructed": True,
            "source_choices": ["one of two native equivariant colorings", "coherent encoder and retained-record recovery",
                               "two native carriers", "joint geometric six-sector preparation"],
            "autonomous_native_defect_creation_or_topological_protection_derived": False,
            "same_object_continuous_u3_plus_so3_response_selected": False,
            "physical_spin_chirality_gauge_dynamics_or_mass_spectrum_recovered": False,
            "antisymmetric_carrier_tensor_sector_establishes_particle_statistics": False,
            "geometric_code5_is_the_spin_lifted_3_plus_2_by_dimension": False,
            "coherent_source_service_accuracy_or_physical_clock_certified_here": False},
        "nearest_existing_results": {"record_protection": "model/common_source/record_code.py",
            "finite_same_response_geometric_hypothesis": "code/source_selection_model/DERIVATION.md",
            "new_result": "native equivariant protected object and actual two-carrier3+3prime intertwiner; current law remains independent"}}


def source_bindings():
    paths = {"model/common_source/native_gauge_object.py", "model/common_source/test_native_gauge_object.py",
             "model/common_source/record_code.py", "model/common_source/slow_geometry.py",
             "model/observer/native.py", "model/observer/quantum_mixing.py", "model/repair_scheduler/coins.py",
             "model/observer/source_discrimination/circuit.py", "model/common_source/native_routing.py",
             "tools/strict_json.py",
             "code/source_selection_model/response.py", "code/source_selection_model/response.json",
             "code/source_selection_model/DERIVATION.md", "code/m1_source_realization/model.py",
             "Lean/Screen/A5ResponseWordAlgebra.lean", "Lean/Screen/A5PortModule.lean",
             slow_geometry.SPEC.relative_to(ROOT).as_posix(),
             (slow_geometry.ARCHIVE/"oph_exact/carrier.py").relative_to(ROOT).as_posix(),
             (slow_geometry.ARCHIVE/"oph_fpe/dynamics/self_readback_repair_closure.py").relative_to(ROOT).as_posix()}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def check(path):
    packet = strict_json.load(Path(path))
    if packet.get("schema") != "oph.native_gauge_object.v1":
        raise ValueError("wrong native gauge-object schema")
    pins = source_bindings()
    if packet.get("source_bindings") != pins:
        raise ValueError("native gauge-object source bindings differ")
    expected = report()
    actual = {key: value for key, value in packet.items() if key != "source_bindings"}
    _same(expected, actual)
    if source_bindings() != pins:
        raise ValueError("native source changed during replay")
    return actual


def _same(expected, actual, path="receipt"):
    if type(expected) is not type(actual):
        raise ValueError("exact native protected-object JSON type differs: "+path)
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise ValueError("exact native protected-object keys differ: "+path)
        for key in expected:
            _same(expected[key], actual[key], path+"."+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise ValueError("exact native protected-object length differs: "+path)
        for i, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _same(left, right, path+"["+str(i)+"]")
    elif expected != actual:
        raise ValueError("exact native protected-object replay differs: "+path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--output", type=Path)
    choice.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print("native equivariant protected object and geometric M6 attachment: PASS")
        return
    if args.output.exists():
        parser.error("a fresh exact receipt path is required")
    pins = source_bindings()
    packet = report()
    if source_bindings() != pins:
        raise RuntimeError("source changed during exact construction")
    packet["source_bindings"] = pins
    args.output.write_text(json.dumps(packet, sort_keys=True, indent=2, allow_nan=False)+"\n")
    print(packet["two_object_geometry"]["selected_geometric_dimension"])


if __name__ == "__main__":
    main()
