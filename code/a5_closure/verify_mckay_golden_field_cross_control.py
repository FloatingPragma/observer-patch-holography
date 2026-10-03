#!/usr/bin/env python3
"""Fail-closed crossover check against the pinned raw OPH producer output.

No OPH code is imported. Structural facts are extracted from the preserved raw
JSON, while both Arithmon receipts are byte-pinned to commit c6f69b5.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIXTURE_PATH = ROOT / "code/a5_closure/receipts/arithmon_mckay_golden_field_cross_control.json"
PINNED_ARITHMON_RECEIPTS = {
    "producer": {
        "path": "code/a5_closure/receipts/arithmon_mckay_golden_field_certificate.json",
        "git_blob_sha": "62e5d622f00cb9df920503014b8d558227e94412",
        "sha256": "179702705753810b70fa0551c1f9847bfa54594fcfb9c02cf46e8041d546b67e",
    },
    "independent_reference": {
        "path": "code/a5_closure/receipts/arithmon_mckay_golden_field_reference.json",
        "git_blob_sha": "847117b3e509af40b4e7b815ba0ab6178810bade",
        "sha256": "f2c1f5741c52a6cd3495e933d5ae542a98c0ea078e49ac733ae5e9135ed62cee",
    },
}
PINNED_OPH_RECEIPT_SHA256 = "6b2daf8fe387cff880aa50af5edc150fd41cfdb0eb6f03ae3ace082ea79a3b55"
PINNED_OPH_COMMIT = "4ae2148a26ce15591adaac78ff408fb2cc32d3a2"
PINNED_OPH_PRODUCER_BLOB = "2938db085477737cab6eb09a3a87f7330c37f473"
E8_NODES = {f"e{i}" for i in range(9)}
E8_EDGES = {
    frozenset(edge) for edge in (
        ("e0", "e1"), ("e1", "e2"), ("e2", "e3"), ("e3", "e4"),
        ("e4", "e5"), ("e5", "e6"), ("e6", "e7"), ("e5", "e8"),
    )
}

def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def parse_q5(text: str) -> tuple[Fraction, Fraction]:
    text = text.strip()
    marker = "*sqrt(5)"
    if marker not in text:
        return Fraction(text), Fraction(0)
    if " + " in text:
        a, b = text.split(" + ", 1)
        if not b.endswith(marker):
            raise ValueError(f"unsupported exact Q(sqrt(5)) value: {text}")
        return Fraction(a), Fraction(b[:-len(marker)])
    if text.endswith(marker):
        return Fraction(0), Fraction(text[:-len(marker)])
    raise ValueError(f"unsupported exact Q(sqrt(5)) value: {text}")

def is_affine_e8(nodes, edges) -> bool:
    nodes = set(nodes)
    if len(nodes) != 9 or len(edges) != 8:
        return False
    adjacency = {node: set() for node in nodes}
    normalized = set()
    for edge in edges:
        if len(edge) != 2:
            return False
        a, b = edge
        if a == b or a not in nodes or b not in nodes:
            return False
        e = frozenset((a, b))
        if e in normalized:
            return False
        normalized.add(e)
        adjacency[a].add(b)
        adjacency[b].add(a)
    seen = set()
    stack = [next(iter(nodes))]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(adjacency[node] - seen)
    if seen != nodes or len(normalized) != len(nodes) - 1:
        return False
    centers = [node for node in nodes if len(adjacency[node]) == 3]
    if len(centers) != 1 or any(len(adjacency[node]) > 3 for node in nodes):
        return False
    center = centers[0]
    arm_lengths = []
    for neighbor in adjacency[center]:
        previous, current, length = center, neighbor, 1
        while len(adjacency[current]) == 2:
            following = next(x for x in adjacency[current] if x != previous)
            previous, current = current, following
            length += 1
        if len(adjacency[current]) != 1:
            return False
        arm_lengths.append(length)
    return sorted(arm_lengths) == [1, 2, 5]

def explicit_e8_map_valid(mapping, nodes, edges) -> bool:
    if set(mapping) != set(nodes) or set(mapping.values()) != E8_NODES or len(mapping) != 9:
        return False
    mapped = {frozenset((mapping[a], mapping[b])) for a, b in edges}
    return mapped == E8_EDGES

def exact_character_norm(table, label: str) -> Fraction:
    total = Fraction(0)
    order = sum(row["size"] for row in table)
    for row in table:
        a, b = parse_q5(row["characters"][label])
        total += row["size"] * (a*a + 5*b*b)
    return total / order

def trace_kernel_is_trivial(table, label: str) -> bool:
    trace_two = [row for row in table if parse_q5(row["characters"][label]) == (Fraction(2), Fraction(0))]
    return len(trace_two) == 1 and trace_two[0]["element_order"] == 1 and trace_two[0]["size"] == 1

def verify() -> dict:
    x = json.loads(FIXTURE_PATH.read_text())
    if x["schema"] != "arithmon.mckay_golden_field.cross_control.v2":
        raise ValueError("unexpected crossover schema")
    if x["arithmon_frozen_commit"] != "c6f69b502d30f465868e0127c2f2e1d30f5791ed":
        raise ValueError("Arithmon freeze commit changed")
    if x["oph"]["commit"] != PINNED_OPH_COMMIT:
        raise ValueError("pinned OPH comparison commit changed")
    if x["oph"]["producer_blob_sha"] != PINNED_OPH_PRODUCER_BLOB:
        raise ValueError("pinned OPH producer blob changed")
    if x["arithmon_receipts"] != PINNED_ARITHMON_RECEIPTS:
        raise ValueError("Arithmon receipt pin set changed")
    if x["oph"]["raw_receipt_sha256"] != PINNED_OPH_RECEIPT_SHA256:
        raise ValueError("pinned OPH raw receipt digest changed")
    if x.get("upstream_current_main") != {
        "commit": "d2fb7da33aa1056dd3feebec50937497bfdb5da9",
        "producer_blob_sha": PINNED_OPH_PRODUCER_BLOB,
        "verified": True,
    }:
        raise ValueError("reported current OPH main/blob note changed")
    if x["oph"]["exact_producer_executed"] is not True:
        raise ValueError("pinned OPH producer execution is not recorded")
    if x["claim_boundary"] != {
        "no_explicit_isomorphism_between_carriers": True,
        "oph_used_as_derivation_input": False,
        "koide_used": False,
        "experimental_data_used": False,
        "observable_changed": False,
        "physical_identification": False,
    }:
        raise ValueError("cross-control claim boundary changed")

    pinned_bytes = {}
    for role, record in x["arithmon_receipts"].items():
        data = (ROOT / record["path"]).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        blob = git_blob_sha(data)
        if digest != record["sha256"] or blob != record["git_blob_sha"]:
            raise ValueError(f"{role} receipt is not the frozen c6f69b5 object")
        pinned_bytes[role] = data
    raw = json.loads(pinned_bytes["producer"])
    ref = json.loads(pinned_bytes["independent_reference"])

    oph_path = ROOT / x["oph"]["raw_receipt_path"]
    oph_bytes = oph_path.read_bytes()
    if hashlib.sha256(oph_bytes).hexdigest() != x["oph"]["raw_receipt_sha256"]:
        raise ValueError("pinned raw OPH receipt SHA256 mismatch")
    oph = json.loads(oph_bytes)
    if oph.get("schema") != x["oph"]["output_schema"] or oph.get("schema") != "oph.sl2f5_mckay_e8.v1":
        raise ValueError("pinned raw OPH receipt schema mismatch")

    oph_ir = oph["irreducible_recovery"]
    oph_dims = oph_ir["dimensions"]
    oph_table = oph["character_table"]
    oph_mckay = oph["mckay"]
    oph_galois = oph["galois_control"]
    oph_spin = oph_mckay["tensor_doublet"]
    oph_conj = oph_galois["conjugate_doublet"]

    # Recheck exact character norms and kernel witnesses from OPH's raw table,
    # rather than accepting booleans copied into the crossover fixture.
    oph_spin_norm = exact_character_norm(oph_table, oph_spin)
    oph_conj_norm = exact_character_norm(oph_table, oph_conj)
    oph_spin_faithful = trace_kernel_is_trivial(oph_table, oph_spin)
    oph_conj_faithful = trace_kernel_is_trivial(oph_table, oph_conj)
    oph_golden_values = [
        parse_q5(row["characters"][oph_spin])
        for row in oph_table
    ]
    oph_field_sensitive = any(b != 0 for _, b in oph_golden_values)
    oph_galois_matches = all(
        parse_q5(row["characters"][oph_conj])
        == (parse_q5(row["characters"][oph_spin])[0], -parse_q5(row["characters"][oph_spin])[1])
        for row in oph_table
    )
    oph_base_map_valid = explicit_e8_map_valid(
        oph_mckay["affine_e8_isomorphism"], oph_dims.keys(), oph_mckay["edges"])
    oph_base_e8 = is_affine_e8(oph_dims.keys(), oph_mckay["edges"]) and oph_base_map_valid
    oph_galois_e8 = is_affine_e8(oph_dims.keys(), oph_galois["edges"])
    oph_galois_recomputed = (
        len(oph_galois["fusion_matrix"]) == len(oph_dims)
        and oph_galois["fusion_matrix"] != oph_mckay["fusion_matrix"]
        and oph_galois["same_labeled_graph"] is False
    )

    ar_raw = json.loads(pinned_bytes["producer"])
    ar_ref = json.loads(pinned_bytes["independent_reference"])
    ar_dims = ar_raw["irreducibles"]["dimensions"]
    ar_edges = ar_ref["fusion"]["edges"]
    ar_sig_edges = ar_ref["galois"]["conjugate_graph"]["edges"]
    ar_e8 = is_affine_e8(range(len(ar_dims)), ar_edges) and ar_ref["affine_e8"]["isomorphic"]
    ar_sig_e8 = is_affine_e8(range(len(ar_dims)), ar_sig_edges) and ar_ref["galois"]["same_affine_e8_graph_type"]
    ar_spin = ar_raw["faithful_doublet"]
    ar_galois = ar_raw["galois_conjugate_fusion"]
    ar_witness = ar_raw["golden_character_field"]["nonrational_trace_witness"]["trace"]

    comparisons = {
        "group_order": (ar_raw["source_group"]["order"], oph["source"]["order"]),
        "conjugacy_class_count": (ar_raw["conjugacy_classes"]["count"], oph["source"]["conjugacy_class_count"]),
        "irreducible_dimension_multiset": (sorted(ar_dims), sorted(oph_dims.values())),
        "sum_squared_dimensions": (ar_raw["irreducibles"]["sum_squared_dimensions"], oph_ir["sum_squared_dimensions"]),
        "faithful_doublet_dimension": (ar_spin["dimension"], oph_dims[oph_spin]),
        "faithful_irreducible_doublet": (ar_spin["faithful"] and ar_spin["irreducible"], oph_spin_norm == 1 and oph_spin_faithful),
        "fusion_edge_count": (len(ar_edges), len(oph_mckay["edges"])),
        "affine_e8": (ar_e8, oph_base_e8),
        "golden_character_field_sensitive": (
            ar_raw["golden_character_field"]["field_exactly"] == "Q(sqrt(5))" and "sqrt(5)" in ar_witness,
            oph_field_sensitive and all(isinstance(v, tuple) for v in oph_golden_values)),
        "galois_doublet_distinct_faithful_irreducible": (
            ar_galois["doublet_distinct"] and ar_galois["doublet_faithful"] and ar_galois["doublet_irreducible"],
            oph_conj != oph_spin and oph_conj_norm == 1 and oph_conj_faithful and oph_galois_matches),
        "galois_fusion_recomputed": (ar_galois["fusion_recomputed"], oph_galois_recomputed),
        "conjugate_affine_e8": (ar_sig_e8, oph_galois_e8),
    }
    mismatches = {name: pair for name, pair in comparisons.items() if pair[0] != pair[1]}
    if mismatches:
        raise ValueError("structural disagreement in raw pinned receipts: " + repr(mismatches))
    if any(value[0] is not True for value in comparisons.values() if isinstance(value[0], bool)):
        raise ValueError("an exact structural gate is false")
    verdict = "INDEPENDENT_EXACT_AGREEMENT"
    if verdict != x["expected_comparison_result"]:
        raise ValueError("computed crossover verdict differs from fixture declaration")
    return {
        "schema": x["schema"],
        "comparison_result": verdict,
        "oph_receipt_sha256": hashlib.sha256(oph_bytes).hexdigest(),
        "arithmon_receipt_blobs": {key: x["arithmon_receipts"][key]["git_blob_sha"] for key in x["arithmon_receipts"]},
        "compared_invariants": {name: {"arithmon": pair[0], "oph": pair[1], "match": True}
                                for name, pair in comparisons.items()},
    }

if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
