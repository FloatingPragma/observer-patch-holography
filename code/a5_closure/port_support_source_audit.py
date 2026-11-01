"""Live-source relative seed/refinement audit; no golden target selects a map.

The only executed scientific source is geometry.py and topology.py. The
attachment search consumes oriented incidence and the source decoration.
This is a bounded audit of those objects, not an absence theorem for OPH.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys
import types

ROOT = Path(__file__).resolve().parents[2]
BASE = "83113a9cc33d82a56c84a8634445c8edd1e599a2"
LOCAL = "Lean/ObserverPatchHolography/CoreAxioms.lean"
GEOMETRY = "code/source_selection_model/geometry.py"
TOPOLOGY = "code/m1_source_realization/topology.py"
SOURCES = (
    LOCAL, GEOMETRY, TOPOLOGY,
    "code/source_selection_model/__init__.py",
    "code/source_selection_model/DERIVATION.md",
    "code/m1_source_realization/README.md",
    "code/m1_source_realization/receipt.json",
    "Lean/Geometry/ScreenCarrierMapCandidate.lean",
    "Lean/Screen/PortFrameGram.lean",
    "Lean/Screen/PortGramRepairBand.lean",
    "code/a5_closure/port_gram_support_attachment_certificate.py",
    "code/a5_closure/verify_port_gram_support_attachment_independent.py",
    "code/a5_closure/galois_source_response_control.py",
)
RECEIPT = ROOT / "code/a5_closure/receipts/port_support_source_audit.json"


def need(ok, message):
    if not ok:
        raise ValueError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def chain(faces, mapping=None):
    result = Counter()
    for face in faces:
        a, b, c = face if mapping is None else tuple(mapping[x] for x in face)
        if len({a, b, c}) < 3:
            continue
        cyclic = ((a, b, c), (b, c, a), (c, a, b))
        reverse = ((a, c, b), (c, b, a), (b, a, c))
        key = min(cyclic + reverse)
        result[key] += 1 if key in cyclic else -1
    return +result, -result


def compose(p, q):
    return tuple(p[x] for x in q)


def inverse(p):
    out = [0] * len(p)
    for i, v in enumerate(p):
        out[v] = i
    return tuple(out)


def order(p):
    value = tuple(range(len(p)))
    for n in range(1, 61):
        value = compose(p, value)
        if value == tuple(range(len(p))):
            return n
    raise ValueError("finite permutation order")


def normalizer_audit(group, faces):
    """Finite assumptions for the separate continuous proof candidate."""
    edges = sorted({tuple(sorted(e)) for a, b, c in faces for e in ((a, b), (b, c), (c, a))})
    rows = []
    for prime in (2, 3, 5):
        g = next(p for p in group if order(p) == prime)
        # Generate directly, without relying on permutation sort order.
        power = tuple(range(12))
        powers = set()
        for _ in range(prime):
            powers.add(power)
            power = compose(g, power)
        normalizer = [p for p in group if compose(compose(p, g), inverse(p)) in powers]
        objects = [(v,) for v in range(12)] if prime == 5 else edges if prime == 2 else faces
        fixed = sorted({tuple(sorted(o)) for o in objects if {g[x] for x in o} == set(o)})
        need(len(fixed) == 2 and len(normalizer) == {2: 4, 3: 6, 5: 10}[prime], "axial normalizer")
        swap = next(p for p in normalizer if order(p) == 2 and {p[x] for x in fixed[0]} == set(fixed[1]))
        rows.append(dict(prime=prime, generator=list(g), normalizer_size=len(normalizer),
                         domain_axis_poles=fixed, pole_swap_involution=list(swap)))
    return rows


def local_faces():
    text = (ROOT / LOCAL).read_text().split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return ast.literal_eval(text[text.index("["):].strip())


def automorphisms(faces):
    """Propagate each of the 60 possible oriented images of one triangle."""
    third = {}
    for a, b, c in faces:
        for x, y, z in ((a, b, c), (b, c, a), (c, a, b)):
            need((x, y) not in third, "unique directed triangle")
            third[x, y] = z
    start = faces[0]
    found = []
    for a, b, c in faces:
        for image in ((a, b, c), (b, c, a), (c, a, b)):
            p = dict(zip(start, image))
            good = True
            while good:
                size = len(p)
                for (x, y), z in third.items():
                    if x not in p or y not in p:
                        continue
                    w = third.get((p[x], p[y]))
                    if w is None or (z in p and p[z] != w) or (z not in p and w in p.values()):
                        good = False
                        break
                    p[z] = w
                if len(p) == size:
                    break
            if good and len(p) == 12:
                found.append(tuple(p[i] for i in range(12)))
    return sorted(set(found))


def lift(p, count, faces):
    edges = sorted({tuple(sorted((a, b))) for f in faces
                    for a, b in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0]))})
    middle = {edge: count + i for i, edge in enumerate(edges)}
    return tuple(p) + tuple(middle[tuple(sorted(p[x] for x in edge))] for edge in edges)


def refinement_audit(group, rows):
    need([r["level"] for r in rows] == [0, 1, 2, 3], "consecutive support levels")
    lifts = {p: p for p in group}
    cumulative = set(group)
    summaries = []
    for old, new in zip(rows, rows[1:]):
        coarse = new["coarsen"]
        compatible = []
        counts = []
        witness = None
        for p in group:
            before = lifts[p]
            after = lift(before, old["vertices"], old["faces"])
            need(chain(new["faces"], after) == chain(new["faces"]), "natural oriented subdivision")
            failures = [i for i in range(new["vertices"])
                        if coarse[after[i]] != before[coarse[i]]]
            if not failures:
                compatible.append(p)
            counts.append(len(failures))
            if failures and witness is None:
                v = failures[0]
                witness = dict(seed_permutation=p, fine_vertex=v,
                               relabeled_fine_vertex=after[v],
                               coarsen_after_relabel=coarse[after[v]],
                               relabel_after_coarsen=before[coarse[v]])
            lifts[p] = after
        cumulative &= set(compatible)
        need(chain(new["faces"], coarse) == chain(old["faces"]), "oriented coarsening degree")
        summaries.append(dict(from_level=old["level"], to_level=new["level"],
                              vertices=new["vertices"], faces=len(new["faces"]),
                              subdivision_chain_symmetries=len(group),
                              coarsening_symmetries=len(compatible),
                              cumulative_coarsening_symmetries=len(cumulative),
                              failure_count_census=sorted(Counter(counts).items()),
                              first_counterexample=witness))
    return summaries


def source_rows():
    # topology imports the live geometry module; run in the same isolated repo.
    name = "source_selection_model.geometry"
    saved = sys.modules.get(name)
    module = types.ModuleType(name)
    exec(compile((ROOT / GEOMETRY).read_bytes(), str(ROOT / GEOMETRY), "exec"), module.__dict__)
    sys.modules[name] = module
    sys.path.insert(0, str(ROOT / "code"))
    try:
        return runpy.run_path(str(ROOT / TOPOLOGY))["candidate"]()
    finally:
        sys.path.pop(0)
        if saved is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = saved


def seed_audit(faces, rows):
    group = automorphisms(faces)
    identity = tuple(range(12))
    need(len(group) == 60 and identity in group, "proper seed group")
    need(all(compose(p, q) in group for p in group for q in group), "group closure")
    base = rows[0]
    need(base["support"]["faces"] == [list(f) for f in faces], "literal live source seed")
    fibres = base["fibres"]
    fixed = [p for p in group if all(fibres[p[v]] == fibres[v] for v in range(12))]
    # Process links are checked on the carrier set, including the two slots.
    links = {tuple(sorted(e)) for e in base["process_links"]}
    for p in fixed:
        carrier_action = [base["carriers"].index([p[v], slot]) for v, slot in base["carriers"]]
        need({tuple(sorted(carrier_action[i] for i in e)) for e in links} == links,
             "full process-link decoration")
    need(len(fixed) == 5, "decorated C5")
    need(all(chain(faces, p) == chain(faces) for p in group), "integral oriented seed maps")
    # h is port -> support. Independent changes act by b h a^{-1}.
    left_orbits = {frozenset(compose(b, h) for b in fixed) for h in group}
    diagonal = {frozenset(compose(compose(a, h), inverse(a)) for a in group) for h in group}
    # Product action with fixed support mark is transitive, stabilizer size 5.
    product_orbit = {compose(b, inverse(a)) for a in group for b in fixed}
    stabilizer = sum(compose(b, inverse(a)) == identity for a in group for b in fixed)
    same_mark = [p for p in group if p[0] == 0]
    edges = sorted({tuple(sorted(e)) for a, b, c in faces for e in ((a, b), (b, c), (c, a))})
    edge_swaps = []
    for a, b in edges:
        swaps = [p for p in group if p[a] == b and p[b] == a]
        need(len(swaps) == 1 and compose(swaps[0], swaps[0]) == identity,
             "endpoint swap involution")
        edge_swaps.append(dict(edge=[a, b], swap=list(swaps[0])))
    return dict(raw_maps=[list(p) for p in group], oriented_chain_equalities=len(group),
                support_fibres=fibres, process_links=base["process_links"],
                fixed_support_stabilizer=[list(p) for p in fixed],
                fixed_support_left_orbit_sizes=sorted(map(len, left_orbits)),
                diagonal_conjugacy_orbit_sizes=sorted(map(len, diagonal)),
                fixed_support_independent_local_action=dict(group_size=300,
                    orbit_size=len(product_orbit), stabilizer_size=stabilizer),
                preimage_of_support_mark_census=sorted(Counter(p.index(0) for p in group).items()),
                indistinguishable_with_added_local_mark=[list(p) for p in same_mark[:2]],
                endpoint_choice_obstruction=dict(candidate_class="A5-equivariant endpoint-valued midpoint rounding with the standard lifted seed action",
                    edge_swap_witnesses=edge_swaps),
                cyclic_axis_normalizers=normalizer_audit(group, faces),
                source_scope="oriented incidence, fixed fibre/slot section and process links",
                quotient_warning="Product/diagonal actions are mathematical controls; no source law identifies independent local and support chart changes.")


def tower_audit(rows):
    summaries = []
    for level, row in enumerate(rows):
        support = row["support"]
        section = row["section"]
        bridge = [v for v, _ in row["carriers"]]
        need([bridge[x] for x in section] == list(range(support["vertices"])), "section right inverse")
        cycle = [tuple(section[x] for x in face) for face in support["faces"]]
        need(chain(cycle, bridge) == chain(support["faces"]), "bridge integral cycle")
        if level:
            prev = rows[level-1]
            v = support["coarsen"]
            q = row["coarsen"]
            need([prev["carriers"][q[x]][0] for x in range(len(bridge))] ==
                 [v[x] for x in bridge], "nerve bridge square")
            need([q[x] for x in section] == [prev["section"][x] for x in v], "section square")
            need(chain(cycle, q) == chain([tuple(prev["section"][x] for x in f)
                                        for f in prev["support"]["faces"]]), "nerve cycle coarsening")
        summaries.append(dict(level=level, support_vertices=support["vertices"],
                              support_faces=len(support["faces"]), nerve_carriers=len(bridge),
                              bridge_chain_degree=1))
    return summaries


def custody():
    result = {}
    for path in SOURCES:
        payload = (ROOT / path).read_bytes()
        result[path] = dict(sha256=sha(payload), git_blob=git("hash-object", path),
                           base_blob=git("rev-parse", BASE + ":" + path))
    return result


def historical_drift():
    namespace = runpy.run_path(str(ROOT / "code/a5_closure/port_gram_support_attachment_certificate.py"))
    return {path: dict(frozen_sha256=expected, live_sha256=sha((ROOT/path).read_bytes()),
                       classification="BASELINE_PIN_DRIFT")
            for path, expected in namespace["PINNED_SOURCE_SHA256"].items()
            if sha((ROOT/path).read_bytes()) != expected}


def build():
    faces = local_faces()
    rows = source_rows()
    committed = json.loads((ROOT / "code/m1_source_realization/receipt.json").read_text())
    need(rows == committed["evidence"]["topology"], "live topology equals committed source rows")
    implementation = {p: sha((ROOT/p).read_bytes()) for p in (
        "code/a5_closure/port_support_source_audit.py",
        "code/a5_closure/verify_port_support_source_audit.py")}
    return dict(schema="oph.port_support_source_audit.v1", base=BASE, implementation_sha256=implementation,
                base_tree=git("rev-parse", BASE + "^{tree}"), sources=custody(),
                historical_drift=historical_drift(), seed=seed_audit(faces, rows),
                support_tower=tower_audit(rows),
                native_support_refinement=refinement_audit(automorphisms(faces), [r["support"] for r in rows]),
                verdicts=["SEED_TORSOR_REMAINS_UNSELECTED", "DECLARED_COARSENING_NONNATURAL", "ATTACHMENT_REMAINS_OPEN"],
                local_refinement="No local incidence tower supplied by ScreenCarrierMapCandidate; support automorphism audit only.",
                common_sphere="REALIZED_SPHERE_COMPARISON_MISSING", physical_selection=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the scoped receipt")
    parser.add_argument("--check", action="store_true", help="byte-exact receipt replay")
    args = parser.parse_args()
    payload = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.write:
        RECEIPT.write_text(payload)
    elif args.check:
        need(RECEIPT.read_text() == payload, "receipt replay")
        print("PORT_SUPPORT_SOURCE_AUDIT_REPLAY_VALID")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
