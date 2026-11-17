"""Complete finite seam-phase reconstruction, with explicit clock-phase nuisance.

The tree is the tree of SeamU1HolonomyClassification, transported through the
already checked native/Lean port permutation. In tree gauge a connection is
nineteen chord phases x. Oriented face phases are Bx in R/Z. Deleting ANY one
face gives an integer matrix with integer inverse: the exact identities
B^-1 B = B B^-1 = I prove reconstruction for every U(1) connection, not only
the rational-turn examples replayed here. The omitted face is minus the sum
of the others. This instantiates the existing seam gauge-class theorem and
extends it from chord-loop labels to a complete local face measurement menu.

For coherent two-route triangle comparisons suppose each selected moving
event contributes the SAME additional phase mu, and the idle contributes nu.
The two-versus-one moving-event comparison gives F_f = h_f + delta; repeating
the oppositely oriented comparison gives R_f = -h_f + delta, delta=mu-nu.
Thus F_f-R_f=2h_f cancels the common dynamics exactly, but does not determine
h_f. All twenty pairs and sphere closure leave precisely TWO solutions,
(h_f,delta) and (h_f+1/2,delta+1/2). One independently calibrated face phase
selects the solution. If each seam instead has its own unknown even dynamic
phase, the squared connection has 2^19 gauge-inequivalent square roots.

These are conditional operational identifiability results. They neither
select coherent transport from the committed C/E instrument nor establish
phase calibration, an observer clock, a protected defect or a particle mass.
All phases are dimensionless turns; no physical time parameter is assigned.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re

from model.common_source import slow_geometry
from model.observer import native, native_loop_readout


HERE = Path(__file__).resolve().parent
DEFAULT_RECEIPT = HERE / "results" / "holonomy_reconstruction.json"
# Literal tree edges of OPH.SeamU1HolonomyClassification.treeSeams.
LEAN_TREE = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 6), (1, 5),
             (1, 7), (2, 8), (3, 9), (4, 10), (5, 11))


def integer_inverse(matrix):
    """Return an integer inverse, checking both products exactly over Z."""
    n = len(matrix)
    if not n or any(len(row) != n for row in matrix):
        raise ValueError("nonempty square matrix required")
    augmented = [[F(value) for value in row] + [F(i == j) for j in range(n)]
                 for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if augmented[i][j]), None)
        if pivot is None:
            raise ValueError("face menu is not independent")
        augmented[j], augmented[pivot] = augmented[pivot], augmented[j]
        factor = augmented[j][j]
        augmented[j] = [value / factor for value in augmented[j]]
        for i in range(n):
            if i != j:
                factor = augmented[i][j]
                augmented[i] = [a - factor * b for a, b in
                                zip(augmented[i], augmented[j], strict=True)]
    if any(value.denominator != 1 for row in augmented for value in row[n:]):
        raise ValueError("face inverse is not integral; phase roots are ambiguous")
    inverse = tuple(tuple(int(value) for value in row[n:]) for row in augmented)
    for left, right in ((matrix, inverse), (inverse, matrix)):
        if any(sum(left[i][k] * right[k][j] for k in range(n)) != int(i == j)
               for i in range(n) for j in range(n)):
            raise ArithmeticError("integer inverse identity failed")
    return inverse


@lru_cache(maxsize=1)
def incidence():
    slow_geometry.registered_binding()
    lean = HERE.parents[1] / "Lean" / "Screen"
    index_text = (lean / "SeamU1HolonomyClassification.lean").read_text(encoding="utf-8")
    endpoints_text = (lean / "SeamCurrentCarrierQuotient.lean").read_text(encoding="utf-8")
    tree_match = re.search(r"def treeSeams[^=]*:=\s*\{([^}]+)\}", index_text)
    endpoints = [re.search(r"def " + name + r"[^=]*:=\s*!\[([^]]+)\]", endpoints_text)
                 for name in ("seamLeft", "seamRight")]
    if tree_match is None or any(match is None for match in endpoints):
        raise ArithmeticError("registered Lean seam/tree definitions unavailable")
    indices = tuple(map(int, tree_match.group(1).split(",")))
    left, right = [tuple(map(int, match.group(1).split(","))) for match in endpoints]
    if tuple((left[i], right[i]) for i in indices) != LEAN_TREE:
        raise ArithmeticError("native reconstruction tree differs from the Lean classifier")
    seams, faces = native_loop_readout.carrier()
    inverse_ports = {v: i for i, v in enumerate(slow_geometry.SOURCE_TO_RER)}
    tree = tuple(sorted(tuple(sorted((inverse_ports[a], inverse_ports[b])))
                        for a, b in LEAN_TREE))
    if len(set(tree)) != 11 or not set(tree) <= set(seams):
        raise ArithmeticError("registered Lean tree left the native carrier")
    chords = tuple(edge for edge in seams if edge not in tree)
    matrix = tuple(tuple(sum((1 if a < b else -1)
                             for a, b in zip(face, face[1:] + face[:1], strict=True)
                             if tuple(sorted((a, b))) == edge)
                         for edge in chords) for face in faces)
    if len(chords) != 19 or any(sum(row[j] for row in matrix) for j in range(19)):
        raise ArithmeticError("closed oriented sphere incidence failed")
    return seams, faces, tree, chords, matrix


def _turns(values):
    return tuple(F(value) % 1 for value in values)


def face_phases(connection):
    seams, faces, _, _, _ = incidence()
    if set(connection) != set(seams):
        raise ValueError("one phase for every canonical native seam is required")
    return tuple(sum((F(connection[tuple(sorted((a, b)))]) if a < b else
                      -F(connection[tuple(sorted((a, b)))]))
                     for a, b in zip(face, face[1:] + face[:1], strict=True)) % 1
                 for face in faces)


def tree_representative(connection):
    seams, _, tree, chords, _ = incidence()
    if set(connection) != set(seams):
        raise ValueError("one phase for every canonical native seam is required")
    root = slow_geometry.SOURCE_TO_RER.index(0)
    potentials = {root: F(0)}
    queue = [root]
    for a in queue:
        for edge in tree:
            if a in edge:
                b = edge[1] if a == edge[0] else edge[0]
                if b not in potentials:
                    potentials[b] = potentials[a] + (F(connection[edge]) if a < b
                                                     else -F(connection[edge]))
                    queue.append(b)
    if len(potentials) != 12:
        raise ArithmeticError("declared tree is disconnected")
    return tuple((F(connection[a, b]) + potentials[a] - potentials[b]) % 1
                 for a, b in chords)


def reconstruct_faces(readouts, omitted=19):
    """Recover all chord phases from any 19 ordered oriented face phases."""
    if not isinstance(omitted, int) or not 0 <= omitted < 20:
        raise ValueError("omitted face index must be in 0..19")
    phases = _turns(readouts)
    if len(phases) != 19:
        raise ValueError("exactly nineteen face phases required")
    matrix = incidence()[4]
    reduced = tuple(row for i, row in enumerate(matrix) if i != omitted)
    inverse = integer_inverse(reduced)
    return tuple(sum(coefficient * phase for coefficient, phase in
                     zip(row, phases, strict=True)) % 1 for row in inverse)


def paired_orientations(phases, offset):
    """Conditional triangle quadrature phases; offset = move phase-idle phase."""
    h = _turns(phases)
    if len(h) != 20 or sum(h) % 1:
        raise ValueError("twenty face phases satisfying sphere closure required")
    return (_turns(phase + F(offset) for phase in h),
            _turns(-phase + F(offset) for phase in h))


def directed_phases(connection, dynamics):
    """Total phase = orientation-even seam dynamics + oriented link phase."""
    seams = incidence()[0]
    if set(connection) != set(seams) or set(dynamics) != set(seams):
        raise ValueError("one link and one even dynamic phase per seam required")
    return {directed: (F(dynamics[a, b]) + sign * F(connection[a, b])) % 1
            for a, b in seams for directed, sign in (((a, b), 1), ((b, a), -1))}


def recover_common_offset(forward, reverse, anchor=None):
    """Exactly two solutions, or one with an independently calibrated face."""
    forward, reverse = _turns(forward), _turns(reverse)
    if len(forward) != 20 or len(reverse) != 20:
        raise ValueError("twenty forward and reverse phase readouts required")
    sums = tuple((a + b) % 1 for a, b in zip(forward, reverse, strict=True))
    if len(set(sums)) != 1:
        raise ValueError("a common move/idle dynamical offset is inconsistent")
    if anchor is not None and (len(anchor) != 2 or type(anchor[0]) is not int
                               or not 0 <= anchor[0] < 20):
        raise ValueError("anchor must name a face and its calibrated phase")
    candidates = []
    for offset in (sums[0] / 2, sums[0] / 2 + F(1, 2)):
        phases = _turns(value - offset for value in forward)
        if sum(phases) % 1:
            continue
        if anchor is not None and phases[anchor[0]] != F(anchor[1]) % 1:
            continue
        candidates.append((phases, offset))
    if not candidates:
        raise ValueError("face closure or independent calibration is inconsistent")
    return tuple(candidates)


def quadrature_character(x_pair, y_pair):
    """Decode a balanced pure two-route character from complete success bins.

    The successful comparison mass is measured, not discarded. A mixed or
    zero-visibility comparison cannot be promoted to a unit transport phase.
    Values here are exact probabilities, not finite-sample estimates.
    """
    if len(x_pair) != 2 or len(y_pair) != 2:
        raise ValueError("plus/minus outcome pairs required")
    x, y = tuple(map(F, x_pair)), tuple(map(F, y_pair))
    weight = sum(x)
    if any(v < 0 for v in x + y) or not 0 < weight <= 1 or sum(y) != weight:
        raise ValueError("both quadratures need the same positive success mass")
    character = ((x[0] - x[1]) / weight, (y[0] - y[1]) / weight)
    if sum(value * value for value in character) != 1:
        raise ValueError("contrast does not identify a pure unit transport phase")
    return character


def native_face_menu():
    """Compile all forty oriented comparisons through the actual native rule.

    A fresh two-load preparation and controllable route/coin schedule are
    supplied. Both arms have two attempts; the short arm has one empty-seam
    idle. Retain all eight equiprobable route/coin histories. Select one of
    the two idle outcomes for the coherent comparison, without discarding
    the other outcome from the probability law.
    """
    seams, faces, _, _, _ = incidence()
    protocols = []
    for face_index, face in enumerate(faces):
        spectator = next(port for port in range(12) if port not in face)
        idle = next(edge for edge in seams if not set(edge) & (set(face) | {spectator}))
        for orientation, (a, b, c) in ((1, face), (-1, (face[0], face[2], face[1]))):
            initial = tuple(int(port in (a, spectator)) for port in range(12))
            terminal = tuple(int(port in (c, spectator)) for port in range(12))
            routes = ((tuple(sorted((a, b))), tuple(sorted((b, c)))),
                      (tuple(sorted((a, c))), idle))
            histories = []
            for route, (first_edge, second_edge) in enumerate(routes):
                for first_coin in (0, 1):
                    middle, first = native.native(initial, first_edge, first_coin)
                    for second_coin in (0, 1):
                        final, second = native.native(middle, second_edge, second_coin)
                        histories.append({"route": route, "coins": [first_coin, second_coin],
                                          "currents": [first[1], second[1]],
                                          "final": list(final)})
            long = [row for row in histories if row["route"] == 0 and tuple(row["final"]) == terminal]
            short = [row for row in histories if row["route"] == 1 and tuple(row["final"]) == terminal]
            if (len(histories) != 8 or len(long) != 1 or len(short) != 2
                    or any(row["currents"][1] != 0 for row in short)
                    or any(sum(value != 0 for value in row["currents"]) != 2 for row in long)):
                raise ArithmeticError("native complete face menu lost its two/one moving arms")
            protocols.append({"face": face_index, "orientation": orientation,
                              "routes": [[a, b, c], [a, c]], "spectator": spectator,
                              "idle_seam": list(idle), "histories": histories,
                              "selected_coins": [long[0]["coins"], short[0]["coins"]]})
    return protocols


def run():
    seams, faces, tree, chords, matrix = incidence()
    connection = {edge: F(i * i + 3, 101) for i, edge in enumerate(seams)}
    phases = face_phases(connection)
    expected = tree_representative(connection)
    inverse_hashes = []
    for omitted in range(20):
        reduced = tuple(row for i, row in enumerate(matrix) if i != omitted)
        inverse = integer_inverse(reduced)
        if reconstruct_faces([h for i, h in enumerate(phases) if i != omitted], omitted) != expected:
            raise ArithmeticError("face reconstruction lost a seam gauge class")
        inverse_hashes.append(hashlib.sha256(canonical(inverse).encode()).hexdigest())
    # All 19 square-root toggles independently alter the gauge class without
    # changing any squared face holonomy. Their products give all 2^19 roots.
    dynamics = {edge: F(i + 2, 73) for i, edge in enumerate(seams)}
    full_directed = directed_phases(connection, dynamics)
    for edge in chords:
        changed = dict(connection)
        changed[edge] += F(1, 2)
        changed_dynamics = dict(dynamics)
        changed_dynamics[edge] += F(1, 2)
        if (_turns(2 * h for h in face_phases(changed)) != _turns(2 * h for h in phases)
                or tree_representative(changed) == expected
                or directed_phases(changed, changed_dynamics) != full_directed):
            raise ArithmeticError("squared-holonomy kernel witness failed")
    forward, reverse = paired_orientations(phases, F(7, 23))
    candidates = recover_common_offset(forward, reverse)
    if (len(candidates) != 2 or candidates[0][0] == candidates[1][0]
            or recover_common_offset(forward, reverse, (0, phases[0])) != ((phases, F(7, 23)),)):
        raise ArithmeticError("common-phase calibration result failed")
    shifted = {edge: phase + F(1, 2) for edge, phase in connection.items()}
    if paired_orientations(face_phases(shifted), F(7, 23) + F(1, 2)) != (forward, reverse):
        raise ArithmeticError("uniform-minus coherent completion was distinguished without calibration")
    # Real source-pulse comparator checks the probability -> character bridge.
    triangle = native_loop_readout.registered_triangle_control()
    protocols = native_face_menu()
    characters = {}
    for name, row in triangle["readout_cases"].items():
        q = row["quadratures"]
        character = quadrature_character(q["X"]["exact_selected_plus_minus"],
                                          q["Y"]["exact_selected_plus_minus"])
        expected_character = (F(1), F(0)) if name == "flat" else (F(-1), F(0))
        if character != expected_character:
            raise ArithmeticError("source-executed comparison character changed")
        characters[name] = list(map(str, character))
    return {
        "schema": "oph.holonomy_reconstruction.v1",
        "status": "COMPLETE_CONDITIONAL_FACE_RECONSTRUCTION_WITH_CALIBRATION_KERNEL",
        "seams": len(seams), "faces": len(faces), "tree_seams": len(tree),
        "gauge_class_phases": len(chords),
        "source_to_registered_lean_port": list(slow_geometry.SOURCE_TO_RER),
        "all_twenty_omitted_face_menus_have_verified_integer_two_sided_inverse": True,
        "integer_inverse_sha256_by_omitted_face": inverse_hashes,
        "squared_holonomy_unresolved_gauge_classes": 2**19,
        "all_sixty_directed_total_phases_preserved_by_each_independent_sign_dynamic_shift": True,
        "common_move_idle_offset_remaining_solutions": len(candidates),
        "one_independently_calibrated_face_resolves_common_offset": True,
        "common_offset_testable_by_all_forward_reverse_phase_sums": True,
        "source_executed_triangle_characters": characters,
        "native_oriented_face_protocols_checked": len(protocols),
        "native_face_route_coin_histories_checked": sum(len(row["histories"]) for row in protocols),
        "complete_native_face_protocols_sha256": hashlib.sha256(canonical(protocols).encode()).hexdigest(),
        "each_face_comparison_mass": "1/4",
        "each_face_unselected_same_endpoint_mass": "1/8",
        "source_selected_coherent_transport": False,
        "independent_move_idle_phase_calibration_provided_by_native_records": False,
        "protected_defect_or_particle_identification": False,
        "finite_theorems": ["OPH.HolonomyInterference.transportRatio_eq_closedLoopHolonomy",
                            "OPH.SeamU1HolonomyClassification.gaugeEquiv_iff_chordHolonomy_eq"],
        "scope": "Exact finite conditional U(1) reconstruction on the registered carrier. "
                 "The 19-class graph connection does not assert 19 protected defects on the filled sphere. "
                 "Complete coherent face access and a stable common move/idle phase are additional "
                 "protocol premises, tested by the paired menu but not selected by native C/E repair.",
        "all_pass": True,
    }


def canonical(payload):
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", type=Path, default=DEFAULT_RECEIPT)
    group.add_argument("--write", type=Path)
    args = parser.parse_args()
    encoded = canonical(run())
    if args.write is not None:
        with args.write.open("x", encoding="utf-8") as stream:
            stream.write(encoded)
    elif args.check.read_text(encoding="utf-8") != encoded:
        raise SystemExit("Holonomy reconstruction receipt mismatch")
    print(json.dumps({"all_pass": True, "sha256": hashlib.sha256(encoded.encode()).hexdigest()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
