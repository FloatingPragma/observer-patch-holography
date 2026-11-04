"""Independent reconstruction of seed attachments and support-chain data.

This module deliberately does not import the certificate producer and never
reads its candidate list. Its permutation search is driven by partial triangle
incidence checks; chain pushforwards are recomputed from the source faces.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FACES_PATH = ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean"
GEOMETRY_PATH = ROOT / "code/source_selection_model/geometry.py"
RECEIPT_PATH = ROOT / "code/m1_source_realization/receipt.json"
DEGREE_VERIFIER = ROOT / "code/a5_closure/verify_galois_port_frame_independent.py"
PINNED_SOURCE_SHA256 = {
    "Lean/ObserverPatchHolography/CoreAxioms.lean": "addcea08fe55c0ccc3edbeb90ebd9df4b640c96bd20e9755417bcdb84f093fe4",
    "code/source_selection_model/geometry.py": "de47ae19dd4ae1677d0cefb0163513f2f1f6d24dfb4dd0264506d0b17fceb843",
    "code/source_selection_model/DERIVATION.md": "e53df3ae59398a6abe2ad3b85f7871ff7ea728e72bd9de073af8bf398a86dcbe",
    "code/m1_source_realization/topology.py": "7319da3f7e9544cf334f7cf152b1aff350f822b2ace0bc786b5bef53e0429f14",
    "code/m1_source_realization/receipt.json": "cef14e5c14bdf02e14f21d6f470bcfc6c6f396dd301e8067ca6454194f25ab49",
    "code/m1_source_realization/README.md": "b7986b78adfa2d0f407fb5150f2188ae84cd3ef96cc425d41c56323823fb0a36",
    "Lean/Geometry/ScreenCarrierMapCandidate.lean": "60a95ebdf173f2c2c5bd8bed4ade34e64392299d787c4305ed90e8335e5c2ab1",
    "Lean/Screen/PortGramRepairBand.lean": "75286414b7c33492b40225dd48ca9321cf3a09ecf96b65e254af9bd02421cf72",
    "code/a5_closure/galois_port_frame_certificate.py": "05b917898f3a2266ca54c22a9c1b174aa87f374196e4297c8bda1839a8b0452a",
    "code/a5_closure/galois_port_frame_certificate.md": "f1927c1481c87411443d42c1d48f0667a3c5cf7221a01f6a05f576dd51b13b44",
    "code/a5_closure/verify_galois_port_frame_independent.py": "0b805033eb3e7c27c812a2d0919504796f73980d5d6d75f23dd68925045fbbea",
    "code/a5_closure/manifests/galois_port_frame_reference.json": "455085536262ad942918ab541bd0a14784d682bb428fa25deabd55f400993a6f",
}


def need(value, label):
    if not value:
        raise ValueError(label)


def cyclic(face):
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


def chain(faces):
    out = Counter()
    for face in faces:
        fwd = (face, face[1:] + face[:1], face[2:] + face[:2])
        rev = ((face[0], face[2], face[1]),
               (face[2], face[1], face[0]),
               (face[1], face[0], face[2]))
        canonical = min(fwd + rev)
        out[canonical] += 1 if canonical in fwd else -1
    return Counter({key: value for key, value in out.items() if value})


def read_local():
    text = FACES_PATH.read_text()
    block = text.split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return [tuple(map(int, x)) for x in re.findall(r"\((\d+), (\d+), (\d+)\)", block)]


def read_geometry_seed():
    text = GEOMETRY_PATH.read_text()
    block = text.split("def seed():", 1)[1].split("\n\n\ndef subdivide", 1)[0]
    return [tuple(map(int, x)) for x in re.findall(r"\((\d+), (\d+), (\d+)\)", block)]


def automorphisms(faces):
    # Each mapped oriented face must still be a source face whenever its
    # three vertices have been assigned. This is independent of the producer's
    # distance-table refinement strategy.
    face_set = {cyclic(f) for f in faces}
    incident = [set() for _ in range(12)]
    adjacent = [set() for _ in range(12)]
    for f in faces:
        for v in f:
            incident[v].add(f)
            adjacent[v].update(set(f) - {v})
    found = []

    def extend(mapping, used):
        if len(mapping) == 12:
            p = tuple(mapping[i] for i in range(12))
            if {cyclic(tuple(p[v] for v in f)) for f in faces} == face_set:
                found.append(p)
            return
        # Choose the vertex constrained by the most already mapped neighbors.
        v = max((x for x in range(12) if x not in mapping),
                key=lambda x: (sum(n in mapping for n in adjacent[x]), len(incident[x]), -x))
        for w in range(12):
            if w in used or len(adjacent[v]) != len(adjacent[w]):
                continue
            if any((u in adjacent[v]) != (image in adjacent[w])
                   for u, image in mapping.items()):
                continue
            trial = {**mapping, v: w}
            valid = True
            for f in incident[v]:
                if all(x in trial for x in f):
                    if cyclic(tuple(trial[x] for x in f)) not in face_set:
                        valid = False
                        break
            if valid:
                extend(trial, used | {w})

    extend({}, set())
    return sorted(set(found))


def verify_support(tower):
    need(len(tower) == 4, "tower row count")
    levels = []
    for level, row in enumerate(tower):
        s = row["support"]
        faces, section = [tuple(f) for f in s["faces"]], row["section"]
        need(s["level"] == level and len(faces) == 20 * 4**level, "support census")
        need(len(section) == s["vertices"] and len(set(section)) == len(section), "section")
        # The section cycle is pushed by b to the entire oriented support cycle.
        section_chain = [tuple(section[v] for v in f) for f in faces]
        vertex_projection = [tuple(row["carriers"][c][0] for c in f) for f in section_chain]
        need(chain(vertex_projection) == chain(faces), "b_* z != [S]")
        if level:
            prev = tower[level - 1]
            support_parent = s["coarsen"]
            nerve_parent = row["coarsen"]
            need(len(support_parent) == s["vertices"], "support parent map")
            need(all(nerve_parent[section[v]] == prev["section"][support_parent[v]]
                     for v in range(s["vertices"])), "section naturality")
            pushed = [tuple(nerve_parent[c] for c in f) for f in section_chain]
            pushed = [f for f in pushed if len(set(f)) == 3]
            prev_chain = [tuple(prev["section"][v] for v in f)
                          for f in prev["support"]["faces"]]
            need(chain(pushed) == chain(prev_chain), "section cycle refinement square")
        levels.append({"level": level, "vertices": s["vertices"],
                       "faces": len(faces), "carriers": len(row["carriers"]),
                       "bridge_degree": 1})
    return levels


def verify():
    source_hashes = {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
                     for path in PINNED_SOURCE_SHA256}
    need(source_hashes == PINNED_SOURCE_SHA256, "pinned attachment corpus drift")
    local = read_local()
    seed = read_geometry_seed()
    packet = json.loads(RECEIPT_PATH.read_text())
    tower = packet["evidence"]["topology"]
    support = [tuple(f) for f in tower[0]["support"]["faces"]]
    need(local == seed == support, "literal ordered seed equality")
    need(len(local) == 20 and len({frozenset(f) for f in local}) == 20, "local faces")
    group = automorphisms(local)
    need(len(group) == 60, "Aut+ census")
    section = tower[0]["section"]
    z0 = [tuple(section[v] for v in f) for f in support]
    target = chain(z0)
    for perm in group:
        h = [section[perm[v]] for v in range(12)]
        image = [tuple(h[v] for v in f) for f in local]
        need(chain(image) == target, "exact integer-chain image")
    # Composition is a free transitive action on the 60 bijective seed maps.
    sample = group[0]
    images = {tuple(g[sample[v]] for v in range(12)) for g in group}
    need(len(images) == 60, "free transitive group action")
    # Independently inspect the fixed level-zero source row. Its larger fibre
    # and sole process link mark support vertex zero, reducing its proper
    # rechart stabilizer from A5 to C5.
    row = tower[0]
    fibres, carriers, links = row["fibres"], row["carriers"], row["process_links"]
    carrier_index = {tuple(pair): i for i, pair in enumerate(carriers)}
    decorated = []
    for p in group:
        if any(fibres[v] != fibres[p[v]] for v in range(12)):
            continue
        action = [carrier_index[(p[v], slot)] for v, slot in carriers]
        moved = {tuple(sorted((action[i], action[j]))) for i, j in links}
        if moved == {tuple(sorted(pair)) for pair in links}:
            decorated.append(p)
    need(len(decorated) == 5, "fixed full-row proper automorphism stabilizer")
    maps = {tuple(p) for p in group}
    remaining, quotient_orbits = set(maps), 0
    while remaining:
        p = next(iter(remaining))
        orbit = {tuple(g[p[v]] for v in range(12)) for g in decorated}
        need(orbit <= maps and len(orbit) == 5, "marked-support action on seed maps")
        remaining -= orbit
        quotient_orbits += 1
    need(quotient_orbits == 12, "fixed-mark quotient orbit count")
    support_levels = verify_support(tower)

    # Degree data comes from the separately maintained independent Galois
    # verifier, pinned here by source bytes. It does not supply attachments.
    verifier_sha = hashlib.sha256(DEGREE_VERIFIER.read_bytes()).hexdigest()
    expected_sha = "0b805033eb3e7c27c812a2d0919504796f73980d5d6d75f23dd68925045fbbea"
    need(verifier_sha == expected_sha, "pinned degree verifier changed")
    run = subprocess.run([sys.executable, "-B", str(DEGREE_VERIFIER)], cwd=ROOT,
                         check=True, capture_output=True, text=True)
    degree_receipt = json.loads(run.stdout)
    degrees = {name: degree_receipt["verified"][name]["degree"] for name in ("plus", "minus")}
    need(degrees == {"plus": 1, "minus": 7}, "independent local degrees")
    return {"verified": True, "seed_relation": "literal oriented equality",
            "candidate_count": 60, "exact_chain_pushforwards": 60,
            "stabilizer_size": 1, "orbit_count": 1, "support_levels": support_levels,
            "fixed_full_support_row_aut_size": len(decorated),
            "raw_maps_mod_fixed_mark_aut_orbits": quotient_orbits,
            "host_mark_orbit_size_under_A5": 12,
            "degree_receipt_sha256": verifier_sha, "local_degrees": degrees,
            "producer_candidate_list_consumed": False}


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True))
