"""Exact source-charge obstruction for the bounded C/E completion campaign.

The proof is a constructive matrix-algebra certificate, not a numerical rank
calculation. It covers source-only reducing sectors, including arbitrary
quantum basis changes. It does not classify history or controller charges,
nonreducing invariant subspaces, or metastable/refining defects.
"""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
from itertools import combinations
import json
from pathlib import Path

from tools import strict_json
from model.common_source import native_defect_decision as native_decision
from model.observer import native
from .family import binding
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]


def _canonical_digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _branch_label(row):
    return {key: row[key] for key in ("law", "seam", "coin", "raw")}


def certify_maps(rows, indices):
    """Generate and verify an exact diagonal-separation/transition certificate.

    ``rows`` contain integer partial injections M; positive Kraus scalar
    factors do not change the generated unital *-algebra. ``indices`` gives
    a closed basis domain in the original 78-state index convention.
    """
    rows = tuple(rows)
    indices = tuple(indices)
    domain = set(indices)
    if not rows or not indices or len(domain) != len(indices):
        raise ValueError("nonempty branch list and distinct state indices required")
    if any(type(i) is not int or not 0 <= i < 78 for i in indices):
        raise ValueError("native 78-state integer indices required")
    maps, supports = [], []
    for row in rows:
        pairs = tuple(tuple(pair) for pair in row["pairs"] if pair[0] in domain)
        if any(b not in domain for _, b in pairs):
            raise ValueError("the restricted domain is not invariant")
        if (len({a for a, _ in pairs}) != len(pairs)
                or len({b for _, b in pairs}) != len(pairs)):
            raise ValueError("integer partial injections required")
        maps.append(pairs)
        supports.append(frozenset(a for a, _ in pairs))

    # M* M is exactly the diagonal support projection. Greedily retain a
    # small subset separating every basis pair; all arithmetic is discrete.
    unresolved = list(combinations(indices, 2))
    selected = []
    while unresolved:
        scores = [sum((a in support) != (b in support) for a, b in unresolved)
                  for support in supports]
        choice = max(range(len(rows)), key=lambda i: scores[i])
        if scores[choice] == 0:
            raise ValueError("record effects do not separate this basis; no certificate")
        selected.append(choice)
        support = supports[choice]
        unresolved = [(a, b) for a, b in unresolved
                      if (a in support) == (b in support)]
    signatures = {i: sum((i in supports[j]) << bit
                         for bit, j in enumerate(selected)) for i in indices}
    if len(set(signatures.values())) != len(indices):
        raise ArithmeticError("selected support projectors do not separate the domain")

    # P_b M P_a = |b><a| for each native transition. Adjoint closure allows
    # both directions. A spanning forest constructs all component matrix
    # units by multiplying along paths; no spectral tolerance enters.
    neighbors = {i: {} for i in indices}
    for branch, pairs in enumerate(maps):
        for source, target in pairs:
            if source != target:
                edge = {"branch": branch, "source": source, "target": target}
                neighbors[source].setdefault(target, edge)
                neighbors[target].setdefault(source, edge)
    remaining, components, forest = set(indices), [], []
    while remaining:
        first = min(remaining)
        queue, reached = deque([first]), {first}
        while queue:
            parent = queue.popleft()
            for child, edge in sorted(neighbors[parent].items()):
                if child in reached:
                    continue
                reached.add(child)
                queue.append(child)
                forest.append({"parent": parent, "child": child,
                               "branch": _branch_label(rows[edge["branch"]]),
                               "native_source": edge["source"],
                               "native_target": edge["target"],
                               "adjoint_for_parent_to_child": parent != edge["source"]})
        components.append(sorted(reached))
        remaining.difference_update(reached)
    return {
        "dimension": len(indices),
        "raw_branches_including_zero": len(rows),
        "zero_branches_on_domain": sum(not pairs for pairs in maps),
        "partial_injection_addresses": sum(map(len, maps)),
        "separated_basis_pairs": len(indices)*(len(indices)-1)//2,
        "separating_effects": [_branch_label(rows[j]) for j in selected],
        "effect_signatures": [[i, str(signatures[i])] for i in indices],
        "components": components,
        "spanning_forest": forest,
        "generated_star_algebra_complex_dimension": sum(len(c)**2 for c in components),
        "commutant_complex_dimension": len(components),
        "algebra": "direct sum of full matrix algebras on the listed components",
        "primitive_projector_formula": "P_i=product_j(D_j if signature_ij=1 else I-D_j), D_j=M_j* M_j",
        "matrix_unit_formula": "P_target M_branch P_source=|target><source|; multiply along forest paths, using adjoints when needed",
    }


def _binary_invariant_control(rows, states):
    binary = {i for i, state in enumerate(states) if max(state) <= 1}
    out = [(row, a, b) for row in rows for a, b in row["pairs"]
           if a in binary and b not in binary]
    incoming = [(row, a, b) for row in rows for a, b in row["pairs"]
                if a not in binary and b in binary]
    if out or not incoming:
        raise ArithmeticError("the binary invariant/nonreducing control changed")
    row, a, b = incoming[0]
    return {"binary_dimension": len(binary), "binary_to_collision_addresses": len(out),
            "collision_to_binary_addresses": len(incoming),
            "binary_subspace_invariant": True, "binary_projector_commutes_with_every_branch": False,
            "commutator_witness": {"branch": _branch_label(row), "source": a, "target": b,
                                   "P_M_minus_M_P_integer_entry": 1}}


def report():
    states = native.states()
    rows = native_decision.complete_record_maps()
    binary = tuple(i for i, state in enumerate(states) if max(state) <= 1)
    reports = {}
    for law in ("C", "E"):
        branches = [row for row in rows if row["law"] == law]
        full = certify_maps(branches, tuple(range(78)))
        shell = certify_maps(branches, binary)
        if (full["commutant_complex_dimension"], shell["commutant_complex_dimension"]) != (1, 1):
            raise ArithmeticError("the native full-matrix-algebra conclusion changed")
        held = [row for row in branches if 0 not in row["seam"]]
        held_certificate = certify_maps(held, tuple(range(78)))
        if sorted(map(len, held_certificate["components"])) != [1, 11, 66]:
            raise ArithmeticError("held-port control changed")
        if any(states[a][0] != states[b][0] for row in held for a, b in row["pairs"]):
            raise ArithmeticError("held-port occupation fails to commute")
        reports[law] = {"full78": full, "binary66": shell,
                        "invariant_not_reducing_control": _binary_invariant_control(branches, states),
                        "held_port_control": {"held_port": 0, "deleted_seams": 5,
                            "outside_campaign_positive_exploration": True,
                            "conserved_source_charge": "load at port 0",
                            "certificate": held_certificate}}

    # E refines C's raw current label. This exact identity also covers a
    # different supported C/E choice on each seam: every such choice's
    # *-algebra contains the complete C branch family after rescaling.
    by_key = {(row["law"], tuple(row["seam"]), row["coin"], row["raw"]): row
              for row in rows}
    refinement_checks = 0
    for row in rows:
        if row["law"] != "C":
            continue
        prefix = ("E", tuple(row["seam"]), row["coin"])
        merged = (by_key[(*prefix, row["raw"])]["pairs"]
                  + by_key[(*prefix, row["raw"]+4)]["pairs"])
        if len(merged) != len({tuple(pair) for pair in merged}) or sorted(merged) != row["pairs"]:
            raise ArithmeticError("E branch refinement of C changed")
        refinement_checks += 1

    paths = [Path(__file__), Path(native.__file__), Path(native_decision.__file__),
             Path(native_decision.coins.__file__), Path(native_decision.circuit.__file__),
             Path(native_decision.slow_geometry.__file__),
             ROOT/"evidence/native_geometric_source_20260925/spec.json",
             ROOT/"evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
             ROOT/"evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py"]
    return {
        "schema": "oph.completion.defect_algebra.v1",
        "gate": "CP-03", "issue": 1139,
        "campaign": binding(),
        "decision": "No nontrivial Hermitian source charge commuting with every admitted active branch, or common reducing sector, on bare H78 or H66",
        "canonical_native_maps_sha256": _canonical_digest(rows),
        "C_E": reports,
        "E_refines_C": {"identities_checked": refinement_checks,
                         "integer_identity": "M_C(raw)=M_E(raw)+M_E(raw+4)"},
        "family_extension": {
            "C_weight": "0<=lambda<=1; both endpoints included",
            "repair_and_exploration": "strictly positive repair probability and probability for every native seam at each nonhalted controller state",
            "campaign_lower_bounds": "repair>=1/16; uniform exploration weight>=1/16",
            "identity_hold": "adding a positive-weight identity branch does not change the generated algebra",
            "controller": "any classical history/controller selecting only these instruments; fixed source-only charge Q tensor I_memory",
            "tag": "C/E selection is a committed classical label, not a coherent superposition of laws",
            "halt_boundary": "the theorem concerns admitted active branches; a permanent full-buffer halt freezes the source and is a termination/failure control, not continued autonomous protection",
            "argument": "positive rescaling preserves each supported Kraus generator; E contains C by the checked refinement identity; all-seam C algebra is full",
            "requires_IID_scheduler": False,
            "new_adaptive_mixing_rate_claimed": False,
        },
        "scope": {
            "arbitrary_quantum_basis_for_source_charge": True,
            "exact_source_only_reducing_sectors_excluded": True,
            "arbitrary_Heisenberg_fixed_observable_classification_from_this_certificate": False,
            "all_source_invariant_subspaces_excluded": False,
            "controller_or_history_encoded_charges_excluded": False,
            "metastability_or_controlled_lifetimes_excluded": False,
            "active_recovery_or_coherent_feedback_excluded": False,
            "refined_or_federated_defects_excluded": False,
            "physical_energy_mass_or_gravity_derived": False,
            "full_CP03_gate_accepted": False,
            "whole_campaign_or_program_falsified": False,
            "time_parameter": "none needed: exact branch-support statement, not a lifetime or physical-clock claim",
        },
        "nearest_existing_results": {
            "native_defect_decision": "classical invariant and collision-subspace exclusion; reused complete C/E Kraus maps",
            "quantum_mixing_and_burn_in": "source forgetting for exogenous schedules with discarded records; not extended to adaptive retained-record dynamics",
            "new_content": "explicit full source branch *-algebra certificate and support-wise extension to the bounded classical-controller campaign",
        },
        "source_sha256": bound_pins({path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in paths}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--output", type=Path)
    modes.add_argument("--verify", "--check", dest="verify", type=Path)
    args = parser.parse_args()
    result = report()
    if args.verify:
        if strict_json.load(args.verify) != result:
            raise ValueError("defect algebra certificate, scope or source custody changed")
        print("Verified exact C/E source-algebra obstruction and negative controls")
    else:
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
