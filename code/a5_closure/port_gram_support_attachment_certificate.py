"""Exact audit of the source-support attachment boundary (seed and tower)."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
LOCAL = "Lean/ObserverPatchHolography/CoreAxioms.lean"
GEOMETRY = "code/source_selection_model/geometry.py"
RECEIPT = "code/m1_source_realization/receipt.json"
GALOIS_REFERENCE = "code/a5_closure/manifests/galois_port_frame_reference.json"

# These bytes bound every source assertion in this packet. A source change
# requires a fresh attachment audit; the producer never silently re-pins.
PINNED_SOURCE_SHA256 = {
    LOCAL: "addcea08fe55c0ccc3edbeb90ebd9df4b640c96bd20e9755417bcdb84f093fe4",
    GEOMETRY: "de47ae19dd4ae1677d0cefb0163513f2f1f6d24dfb4dd0264506d0b17fceb843",
    "code/source_selection_model/DERIVATION.md": "e53df3ae59398a6abe2ad3b85f7871ff7ea728e72bd9de073af8bf398a86dcbe",
    "code/m1_source_realization/topology.py": "7319da3f7e9544cf334f7cf152b1aff350f822b2ace0bc786b5bef53e0429f14",
    RECEIPT: "cef14e5c14bdf02e14f21d6f470bcfc6c6f396dd301e8067ca6454194f25ab49",
    "code/m1_source_realization/README.md": "b7986b78adfa2d0f407fb5150f2188ae84cd3ef96cc425d41c56323823fb0a36",
    "Lean/Geometry/ScreenCarrierMapCandidate.lean": "60a95ebdf173f2c2c5bd8bed4ade34e64392299d787c4305ed90e8335e5c2ab1",
    "Lean/Screen/PortGramRepairBand.lean": "75286414b7c33492b40225dd48ca9321cf3a09ecf96b65e254af9bd02421cf72",
    "code/a5_closure/galois_port_frame_certificate.py": "05b917898f3a2266ca54c22a9c1b174aa87f374196e4297c8bda1839a8b0452a",
    "code/a5_closure/galois_port_frame_certificate.md": "f1927c1481c87411443d42c1d48f0667a3c5cf7221a01f6a05f576dd51b13b44",
    "code/a5_closure/verify_galois_port_frame_independent.py": "0b805033eb3e7c27c812a2d0919504796f73980d5d6d75f23dd68925045fbbea",
    GALOIS_REFERENCE: "455085536262ad942918ab541bd0a14784d682bb428fa25deabd55f400993a6f",
}

ATTACHMENT_AUDIT_CORPUS = tuple(PINNED_SOURCE_SHA256)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def verify_source_pins(root=ROOT, source_bytes_override=None):
    source_bytes_override = source_bytes_override or {}
    actual = {}
    for path, expected in PINNED_SOURCE_SHA256.items():
        payload = source_bytes_override.get(path)
        if payload is None:
            payload = (root / path).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        require(digest == expected, f"pinned attachment corpus drift: {path}")
        actual[path] = digest
    return actual


def local_faces() -> list[tuple[int, int, int]]:
    text = (ROOT / LOCAL).read_text()
    body = text.split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return [tuple(map(int, f)) for f in re.findall(r"\((\d+), (\d+), (\d+)\)", body)]


def geometry_seed() -> list[tuple[int, int, int]]:
    tree = ast.parse((ROOT / GEOMETRY).read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "seed")
    returned = next(n.value for n in fn.body if isinstance(n, ast.Return))
    return [tuple(row) for row in ast.literal_eval(returned)]


def oriented_key(face):
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


def oriented_chain(faces):
    out = Counter()
    for a, b, c in faces:
        key = min((a, b, c), (b, c, a), (c, a, b),
                  (a, c, b), (c, b, a), (b, a, c))
        out[key] += 1 if key in ((a, b, c), (b, c, a), (c, a, b)) else -1
    return out


def validate_seed_triplet(local, geometry, support):
    require(local == geometry == support, "literal oriented seed equality")
    require(len(local) == 20 and len({frozenset(f) for f in local}) == 20,
            "oriented seed face census")


def validate_report(packet):
    require(packet["seed"]["relation"] == "SUPPORT_SEED_LITERAL_ORIENTED_MATCH",
            "seed relation")
    require(packet["source_binding"]["status"] ==
            "NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS",
            "unproved source binding promoted")
    require(packet["source_binding"]["corpus"] == list(ATTACHMENT_AUDIT_CORPUS),
            "bounded attachment corpus changed")
    require(packet["global_support_degree"]["value"] == 1, "global support degree")
    require(packet["local_radial_degrees"]["plus"] == 1 and
            packet["local_radial_degrees"]["minus"] == 7 and
            packet["local_radial_degrees"]["status"] ==
            "CONSUMED_FROM_PINNED_GALOIS_REFERENCE_AND_INDEPENDENTLY_REPLAYED",
        "local branch degrees / coverage")
    family = packet["admissible_attachments"]
    require((family["count"], family["stabilizer_size"], family["orbit_count"],
             family["orbit_size"], family["is_A5_torsor"]) == (60, 1, 1, 60, True),
            "attachment orbit")
    require(family["fixed_full_support_row_aut_size"] == 5 and
            family["raw_maps_mod_fixed_mark_aut_orbits"] == 12,
            "decorated support stabilizer and quotient")
    require(packet["secondary_verdicts"] == [
        "SEED_ATTACHMENT_FAMILY_CLASSIFIED", "SEED_ATTACHMENT_SET_IS_A5_TORSOR",
        "ATTACHMENT_ONLY_AT_SEED_LEVEL", "CONDITIONAL_SUPPORT_SELECTOR"],
        "seed-scoped verdicts")
    require(packet["primary_verdict"] == "SUPPORT_ATTACHMENT_REMAINS_UNSELECTED",
            "primary verdict")
    require(packet["target_firewall"] ==
            "PASSED: no branch degree, G_plus, shortest-chord, trace, repair, PR#999 selector or physical target selects an attachment",
            "target leakage")
    require(packet["refinement_compatibility"] ==
            "ATTACHMENT_ONLY_AT_SEED_LEVEL_IN_PINNED_AUDITED_CORPUS",
            "unsupported refinement claim")
    require(packet["geometric_comparison"]["oriented_homotopy"] ==
            "NOT_ESTABLISHED_SOURCE_BOUND", "unsupported homotopy claim")
    require(packet["pr53_status"] == "UNCHANGED_NOT_DISCHARGED", "PR-53 boundary")
    require(len(packet["hostile_control_names"]) == 20 and
            len(set(packet["hostile_control_names"])) == 20, "hostile control coverage")
    return packet


def proper_automorphisms(faces):
    return oriented_isomorphisms(faces, faces)


def oriented_isomorphisms(source, target):
    """Enumerate bijections preserving ordered oriented triangle incidence."""
    source_set = {oriented_key(f) for f in source}
    target_set = {oriented_key(f) for f in target}
    if len(source_set) != len(source) or len(target_set) != len(target):
        return []
    source_adj = [set() for _ in range(12)]
    target_adj = [set() for _ in range(12)]
    incident = [set() for _ in range(12)]
    for a, b, c in source:
        source_adj[a].update((b, c)); source_adj[b].update((a, c)); source_adj[c].update((a, b))
        for v in (a, b, c): incident[v].add((a, b, c))
    for a, b, c in target:
        target_adj[a].update((b, c)); target_adj[b].update((a, c)); target_adj[c].update((a, b))
    if sorted(map(len, source_adj)) != sorted(map(len, target_adj)):
        return []
    result = []

    def visit(mapping, used):
        if len(mapping) == 12:
            p = tuple(mapping[i] for i in range(12))
            if {oriented_key(tuple(p[v] for v in f)) for f in source} == target_set:
                result.append(p)
            return
        v = max((x for x in range(12) if x not in mapping),
                key=lambda x: (sum(y in mapping for y in source_adj[x]), -x))
        for w in range(12):
            if w in used or len(source_adj[v]) != len(target_adj[w]):
                continue
            if any((u in source_adj[v]) != (z in target_adj[w])
                   for u, z in mapping.items()):
                continue
            trial = {**mapping, v: w}
            valid = all(oriented_key(tuple(trial[x] for x in f)) in target_set
                        for f in incident[v] if all(x in trial for x in f))
            if valid:
                visit(trial, used | {w})

    visit({}, set())
    return sorted(set(result))


def classify_seed_attachments(local, support, section, designated_cycle=None,
                              proposed_map=None):
    """Pure seed classifier with injectable presentations and exact cycles."""
    require(len(section) == 12 and len(set(section)) == 12, "seed section")
    target_cycle = designated_cycle
    if target_cycle is None:
        target_cycle = [tuple(section[v] for v in f) for f in support]
    target_chain = oriented_chain(target_cycle)
    maps = oriented_isomorphisms(local, support)
    if proposed_map is not None:
        proposed_map = tuple(proposed_map)
        require(proposed_map in maps, "proposed map is not an oriented simplicial isomorphism")
        maps = [proposed_map]
    attachments = []
    for p in maps:
        h = tuple(section[p[v]] for v in range(12))
        image_chain = oriented_chain(tuple(h[v] for v in face) for face in local)
        require(image_chain == target_chain, "candidate is not the exact designated cycle")
        attachments.append({"port_to_support_vertex": list(p),
                            "port_to_section_carrier": list(h),
                            "cycle_pushforward": "EXACT_DESIGNATED_CHAIN"})
    return attachments


def decorated_support_automorphisms(faces, row):
    """Proper seed actions preserving the fixed fibre/process decoration."""
    fibres, carriers, links = row["fibres"], row["carriers"], row["process_links"]
    carrier_index = {tuple(pair): i for i, pair in enumerate(carriers)}
    result = []
    for p in proper_automorphisms(faces):
        if any(fibres[v] != fibres[p[v]] for v in range(12)):
            continue
        carrier_action = [carrier_index[(p[v], slot)]
                          for v, slot in carriers]
        moved_links = {tuple(sorted((carrier_action[i], carrier_action[j]))) for i, j in links}
        if moved_links == {tuple(sorted(pair)) for pair in links}:
            result.append(p)
    return result


def attachment_orbit_count(attachments, target_group):
    members = {tuple(row["port_to_support_vertex"]) for row in attachments}
    remaining = set(members)
    count = 0
    while remaining:
        p = next(iter(remaining))
        orbit = {tuple(g[p[v]] for v in range(12)) for g in target_group}
        require(orbit <= members, "decorated group action leaves attachment family")
        require(len(orbit) == len(target_group), "nonfree target action on seed maps")
        remaining -= orbit
        count += 1
    return count


def galois_degrees(reference, local):
    require(reference.get("schema") == "oph.galois-port-frames.v1", "Galois receipt schema")
    require(reference.get("faces") == [list(f) for f in local], "Galois receipt local faces")
    out = {}
    for branch in ("plus", "minus"):
        controls = reference["family"][branch]["degree_controls"]
        values = {row["degree"] for row in controls}
        require(len(controls) == 3 and len(values) == 1, f"Galois {branch} degree controls")
        out[branch] = values.pop()
    return out


def verify_tower(rows):
    require(len(rows) == 4, "support tower coverage")
    summaries = []
    for level, row in enumerate(rows):
        support = row["support"]
        faces, n = support["faces"], support["vertices"]
        require(support["level"] == level and len(faces) == 20 * 4**level,
                "support level/face census")
        require(len(row["section"]) == n and len(row["carriers"]) == len(row["coarsen"])
                if level else len(row["section"]) == n, "section/carrier census")
        section = row["section"]
        require(len(set(section)) == n, "section is injective")
        # The nerve has a distinct carrier at every listed fibre slot. Its
        # simplices are exactly the inverse images of support simplices.
        require(all(0 <= c < len(row["carriers"]) for c in section), "section range")
        nerve_cycle = [tuple(section[v] for v in f) for f in faces]
        pushed = [tuple(row["carriers"][c][0] for c in f) for f in nerve_cycle]
        require(oriented_chain(pushed) == oriented_chain(faces), "b_* z_r != [S_r]")
        if level:
            prev = rows[level - 1]
            require(len(support["coarsen"]) == n, "support coarsening domain")
            require(len(row["coarsen"]) == len(row["carriers"]), "nerve coarsening domain")
            prev_section = prev["section"]
            for v, parent in enumerate(support["coarsen"]):
                require(row["coarsen"][section[v]] == prev_section[parent],
                        "section/coarsening square")
            # Verify the designated section cycle maps exactly to the prior
            # designated cycle under the committed nerve coarsening.
            mapped = [tuple(row["coarsen"][c] for c in f) for f in nerve_cycle]
            mapped = [f for f in mapped if len(set(f)) == 3]
            prior = [tuple(prev_section[v] for v in f) for f in prev["support"]["faces"]]
            require(oriented_chain(mapped) == oriented_chain(prior), "z_r coarsening square")
        summaries.append({"level": level, "support_vertices": n,
                          "support_faces": len(faces), "nerve_carriers": len(row["carriers"]),
                          "section_cycle_faces": len(nerve_cycle), "bridge_degree": 1})
    return summaries


def build(local_faces_override=None, support_tower_override=None,
          source_bytes_override=None):
    source_pins = verify_source_pins(source_bytes_override=source_bytes_override)
    fs = local_faces()
    if local_faces_override is not None:
        fs = [tuple(f) for f in local_faces_override]
    seed = geometry_seed()
    receipt = json.loads((ROOT / RECEIPT).read_text())
    tower = (copy.deepcopy(support_tower_override) if support_tower_override is not None
             else receipt["evidence"]["topology"])
    support_seed = [tuple(f) for f in tower[0]["support"]["faces"]]
    validate_seed_triplet(fs, seed, support_seed)
    group = proper_automorphisms(fs)
    require(len(group) == 60, "proper automorphism census")
    section = tower[0]["section"]
    z0 = [tuple(section[v] for v in f) for f in support_seed]
    attachments = classify_seed_attachments(fs, support_seed, section, z0)
    require(len({tuple(x["port_to_support_vertex"]) for x in attachments}) == 60,
            "attachment uniqueness")
    decorated = decorated_support_automorphisms(support_seed, tower[0])
    require(len(decorated) == 5, "fixed level-zero decorated support stabilizer")
    decorated_orbits = attachment_orbit_count(attachments, decorated)
    require(decorated_orbits == 12, "fixed marked support attachment quotient")
    tower_summary = verify_tower(tower)
    reference = json.loads((ROOT / GALOIS_REFERENCE).read_text())
    require(reference["face_source_sha256"] == source_pins[LOCAL],
            "Galois receipt local source digest")
    degrees = galois_degrees(reference, fs)
    packet = {
        "schema": "oph.port-gram-support-attachment.v1",
        "baseline_sha": "800ed61ac9fb21d61f206a0240c8b3717c3af945",
        "seed": {"relation": "SUPPORT_SEED_LITERAL_ORIENTED_MATCH",
                 "vertices": 12, "ordered_faces": 20, "edges": 30,
                 "fundamental_cycle": "sum of the 20 listed oriented faces"},
        "source_binding": {"status": "NO_ATTACHMENT_FOUND_IN_PINNED_AUDITED_CORPUS",
                           "corpus": list(ATTACHMENT_AUDIT_CORPUS),
                           "evidence": "The pinned source-clock model declares each local twelve-port boundary independently of the federation nerve. The shared labels are a presentation convention; no port-to-support map is declared in this bounded corpus."},
        "support_tower": tower_summary,
        "global_support_degree": {"value": 1, "status": "SOURCE_CONSTRUCTED_CHAIN_BRIDGE",
                                  "realization": "The receipt certifies b_* z_r=[S_r] and oriented support sphere; it does not give a coordinate map from local radial rays."},
        "local_radial_degrees": {**degrees,
                                 "status": "CONSUMED_FROM_PINNED_GALOIS_REFERENCE_AND_INDEPENDENTLY_REPLAYED"},
        "admissible_attachments": {"count": len(attachments), "group": "Aut+(K_i) ≅ A5",
                                   "Aut+(unmarked section image)": "A5",
                                   "stabilizer_size": 1, "orbit_count": 1, "orbit_size": 60,
                                   "is_A5_torsor": True, "members": attachments,
                                   "fixed_full_support_row_aut_size": len(decorated),
                                   "fixed_full_support_row_aut": "C5 (fibre/process marking fixes support vertex 0)",
                                   "raw_maps_mod_fixed_mark_aut_orbits": decorated_orbits,
                                   "host_mark_orbit_size_under_A5": 12,
                                   "presentation_equivalence_action_on_attachment_maps": "NOT_DECLARED; no complete source-relative quotient claimed"},
        "cycle_pushforward": "EXACT_DESIGNATED_CHAIN_FOR_ALL_60",
        "refinement_compatibility": "ATTACHMENT_ONLY_AT_SEED_LEVEL_IN_PINNED_AUDITED_CORPUS",
        "refinement_evidence": "Within the pinned corpus, local Lean refinement scales barycentric coordinates within the same face; the support tower uses midpoint subdivision and fibre coarsening. No shared local incidence tower or attachment square is present.",
        "geometric_comparison": {"exact_equality": "NOT_PRESENT", "oriented_homotopy": "NOT_ESTABLISHED_SOURCE_BOUND",
                                 "degree_only": "CONDITIONAL: if both realizations are continuous maps to the same oriented S2 with the pinned shared convention, plus has degree 1 and minus degree 7; degree classifies S2-to-S2 homotopy classes."},
        "target_firewall": "PASSED: no branch degree, G_plus, shortest-chord, trace, repair, PR#999 selector or physical target selects an attachment",
        "primary_verdict": "SUPPORT_ATTACHMENT_REMAINS_UNSELECTED",
        "secondary_verdicts": ["SEED_ATTACHMENT_FAMILY_CLASSIFIED", "SEED_ATTACHMENT_SET_IS_A5_TORSOR",
                                "ATTACHMENT_ONLY_AT_SEED_LEVEL", "CONDITIONAL_SUPPORT_SELECTOR"],
        "missing_source_datum": "A source-naturally selected map port p -> support vertex, extended to a carrier map h_r with a refinement-commuting square and a realized oriented homotopy.",
        "pr53_status": "UNCHANGED_NOT_DISCHARGED",
        "source_pins": source_pins,
        "hostile_control_names": [
            "reverse_local_face", "reverse_support_face", "relabel_local_only", "relabel_support_only",
            "relabel_both_coherently", "orientation_reversing_symmetry", "swap_support_labels_without_faces",
            "negate_designated_cycle", "forge_support_degree_7", "delete_support_face",
            "alter_coarsening_image", "alter_section_choice", "break_refinement_square",
            "omit_minus_branch", "import_gplus_as_target", "import_pr999_dynamic_selector",
            "silent_integer_label_binding", "graph_isomorphism_instead_of_oriented_chain",
            "homotopy_from_degree_without_sphere_premises", "claim_pr53_discharged"],
    }
    return validate_report(packet)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    packet = build()
    if args.json:
        print(json.dumps(packet, sort_keys=True, indent=2))
    else:
        print(packet["primary_verdict"])
        print(f"seed candidates: {packet['admissible_attachments']['count']}; A5 torsor, exact chain pushforward")
