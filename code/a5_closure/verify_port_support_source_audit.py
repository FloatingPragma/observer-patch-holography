"""Independent adjacency-search and incidence replay; never imports the producer."""

import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
BASE = "31e1506c30114e37b45f6858d0029a8ff0b89db5"
DEFAULT = ROOT / "code/a5_closure/receipts/port_support_source_audit.json"
SOURCE_PATHS = {
    "Lean/ObserverPatchHolography/CoreAxioms.lean",
    "code/source_selection_model/geometry.py",
    "code/source_selection_model/__init__.py",
    "code/m1_source_realization/topology.py",
    "code/source_selection_model/DERIVATION.md",
    "code/m1_source_realization/README.md",
    "code/m1_source_realization/receipt.json",
    "Lean/Geometry/ScreenCarrierMapCandidate.lean",
    "Lean/Screen/PortFrameGram.lean",
    "Lean/Screen/PortGramRepairBand.lean",
    "code/a5_closure/port_gram_support_attachment_certificate.py",
    "code/a5_closure/verify_port_gram_support_attachment_independent.py",
    "code/a5_closure/galois_source_response_control.py",
}


def require(ok, why):
    if not ok:
        raise ValueError(why)


def integral_chain(faces, p=None):
    out = Counter()
    for f in faces:
        xs = list(f) if p is None else [p[v] for v in f]
        if len(set(xs)) < 3:
            continue
        parity = sum(xs[i] > xs[j] for i in range(3) for j in range(i+1, 3)) % 2
        out[tuple(sorted(xs))] += (-1) ** parity
    return {f: n for f, n in out.items() if n}


def search(faces):
    adjacency = [set() for _ in range(12)]
    cyclic = {tuple(f[k:]+f[:k]) for f in faces for k in range(3)}
    for a, b, c in faces:
        adjacency[a].update((b, c))
        adjacency[b].update((a, c))
        adjacency[c].update((a, b))
    found = []

    def visit(p, used):
        if len(p) == 12:
            found.append(tuple(p[i] for i in range(12)))
            return
        v = max((i for i in range(12) if i not in p), key=lambda i: (len(adjacency[i] & p.keys()), -i))
        for w in range(12):
            if w in used or any((x in adjacency[v]) != (y in adjacency[w]) for x, y in p.items()):
                continue
            trial = dict(p)
            trial[v] = w
            if all(tuple(trial[x] for x in f) in cyclic for f in faces if set(f) <= trial.keys()):
                visit(trial, used | {w})

    visit({}, set())
    return sorted(found)


def cp(a, b):
    return tuple(a[x] for x in b)


def inv(a):
    return tuple(a.index(i) for i in range(len(a)))


def midpoint_identities(old, fine):
    adjacency = {v: set() for v in range(fine["vertices"])}
    for f in fine["faces"]:
        for v in f:
            adjacency[v].update(set(f)-{v})
    n = old["vertices"]
    endpoints = {v: frozenset(x for x in adjacency[v] if x < n)
                 for v in range(n, fine["vertices"])}
    require(all(len(e) == 2 for e in endpoints.values()), "midpoint endpoint types")
    require(len(set(endpoints.values())) == len(endpoints), "midpoint injectivity")
    return endpoints


def check_refinements(group, tower, submitted):
    require(len(tower) == 4 and len(submitted) == 3, "no skipped levels")
    actions = {p: p for p in group}
    cumulative = set(group)
    for r in range(1, 4):
        old, fine = tower[r-1], tower[r]
        endpoints = midpoint_identities(old, fine)
        middle = {e: v for v, e in endpoints.items()}
        c = fine["coarsen"]
        require(c == list(range(old["vertices"])) + [max(endpoints[v])
                for v in range(old["vertices"], fine["vertices"])], "declared max coarsening")
        require(integral_chain(fine["faces"], c) == integral_chain(old["faces"]), "signed coarse cycle")
        failures = {}
        first = None
        for p in group:
            q = actions[p]
            lifted = q + tuple(middle[frozenset(q[x] for x in endpoints[v])]
                               for v in range(old["vertices"], fine["vertices"]))
            require(integral_chain(fine["faces"], lifted) == integral_chain(fine["faces"]), "subdivision chain symmetry")
            bad = [i for i in range(len(lifted)) if c[lifted[i]] != q[c[i]]]
            failures[p] = len(bad)
            if bad and first is None:
                v = bad[0]
                first = dict(seed_permutation=list(p), fine_vertex=v, relabeled_fine_vertex=lifted[v],
                             coarsen_after_relabel=c[lifted[v]], relabel_after_coarsen=q[c[v]])
            actions[p] = lifted
        compatible = {p for p, n in failures.items() if not n}
        cumulative &= compatible
        expected = dict(from_level=r-1, to_level=r, vertices=fine["vertices"], faces=len(fine["faces"]),
                        subdivision_chain_symmetries=len(group), coarsening_symmetries=len(compatible),
                        cumulative_coarsening_symmetries=len(cumulative),
                        failure_count_census=[list(x) for x in sorted(Counter(failures.values()).items())],
                        first_counterexample=first)
        require(submitted[r-1] == expected, "exact refinement audit")


def verify(packet=None):
    p = json.loads(DEFAULT.read_text()) if packet is None else packet
    require(p["schema"] == "oph.port_support_source_audit.v1" and p["base"] == BASE, "schema/base")
    require(p["implementation_sha256"] == {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in ("code/a5_closure/port_support_source_audit.py",
                         "code/a5_closure/verify_port_support_source_audit.py")}, "implementation custody")
    git = lambda *args: subprocess.check_output(["git", *args], cwd=ROOT).decode().strip()
    require(p["base_tree"] == git("rev-parse", BASE + "^{tree}"), "base tree")
    require(set(p["sources"]) == SOURCE_PATHS, "explicit source closure")
    for path, pin in p["sources"].items():
        data = (ROOT / path).read_bytes()
        expected = dict(sha256=hashlib.sha256(data).hexdigest(), git_blob=git("hash-object", path),
                        base_blob=git("rev-parse", BASE + ":" + path))
        require(pin == expected and pin["git_blob"] == pin["base_blob"], "source custody: " + path)
    # Read the historical pins as data, resolving just the named path keys.
    text = (ROOT / "code/a5_closure/port_gram_support_attachment_certificate.py").read_text()
    tree = ast.parse(text)
    names = {}
    for n in tree.body:
        if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name):
            if isinstance(n.value, ast.Constant) and isinstance(n.value.value, str):
                names[n.targets[0].id] = n.value.value
            if n.targets[0].id == "PINNED_SOURCE_SHA256":
                pins = {(k.value if isinstance(k, ast.Constant) else names[k.id]): ast.literal_eval(v)
                        for k, v in zip(n.value.keys, n.value.values)}
    drift = {}
    for path, frozen in pins.items():
        live = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        if frozen != live:
            drift[path] = dict(frozen_sha256=frozen, live_sha256=live, classification="BASELINE_PIN_DRIFT")
    require(p["historical_drift"] == drift, "historical pin drift disclosure")
    body = (ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean").read_text().split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    faces = [list(f) for f in ast.literal_eval(body[body.index("["):].strip())]
    geometry = runpy.run_path(str(ROOT / "code/source_selection_model/geometry.py"))
    require(faces == [list(f) for f in geometry["seed"]()], "live seed equality")
    group = search(faces)
    require(len(group) == 60, "independent search census")
    s = p["seed"]
    require(s["raw_maps"] == [list(a) for a in group], "complete oriented map family")
    require(all(integral_chain(faces, a) == integral_chain(faces) for a in group)
            and s["oriented_chain_equalities"] == 60, "integral seed pushforward")
    require(s["support_fibres"] == [2]+[1]*11 and s["process_links"] == [[0, 1]], "fixed decoration")
    h = [a for a in group if a[0] == 0]
    require(s["fixed_support_stabilizer"] == [list(a) for a in h] and len(h) == 5, "full C5 stabilizer")
    left = {frozenset(cp(b, a) for b in h) for a in group}
    diag = {frozenset(cp(cp(b, a), inv(b)) for b in group) for a in group}
    require(s["fixed_support_left_orbit_sizes"] == sorted(map(len, left)), "left quotient")
    require(s["diagonal_conjugacy_orbit_sizes"] == sorted(map(len, diag)), "diagonal conjugation control")
    product = {cp(b, inv(a)) for a in group for b in h}
    stabilizer = sum(cp(b, inv(a)) == tuple(range(12)) for a in group for b in h)
    require(s["fixed_support_independent_local_action"] ==
            dict(group_size=300, orbit_size=len(product), stabilizer_size=stabilizer), "product control")
    require(s["preimage_of_support_mark_census"] == [[i, 5] for i in range(12)], "mark preimages")
    require(s["indistinguishable_with_added_local_mark"] == [list(a) for a in h[:2]], "explicit two-map witness")
    edges = sorted({tuple(sorted((a, b))) for f in faces for a in f for b in f if a != b})
    swaps = [dict(edge=[a, b], swap=list(next(p for p in group if p[a] == b and p[b] == a)))
             for a, b in edges]
    require(len(edges) == 30 and all(cp(tuple(w["swap"]), tuple(w["swap"])) == tuple(range(12)) for w in swaps),
            "all edge-swap involutions")
    require(s["endpoint_choice_obstruction"] == dict(candidate_class="A5-equivariant endpoint-valued midpoint rounding with the standard lifted seed action", edge_swap_witnesses=swaps), "endpoint rounding obstruction")
    def permutation_order(a):
        unused, lengths = set(range(12)), []
        while unused:
            start = unused.pop()
            v, length = a[start], 1
            while v != start:
                unused.remove(v)
                v, length = a[v], length+1
            lengths.append(length)
        from math import lcm
        return lcm(*lengths)
    require(len(s["cyclic_axis_normalizers"]) == 3, "normalizer coverage")
    for row, prime in zip(s["cyclic_axis_normalizers"], (2, 3, 5)):
        g = next(a for a in group if permutation_order(a) == prime)
        # A cyclic subgroup's generator conjugates to a nonidentity power.
        powers, value = set(), tuple(range(12))
        for _ in range(prime):
            powers.add(value)
            value = cp(g, value)
        normalizer = [a for a in group if cp(cp(a, g), inv(a)) in powers]
        objects = [(v,) for v in range(12)] if prime == 5 else edges if prime == 2 else faces
        fixed = sorted({tuple(sorted(o)) for o in objects if {g[v] for v in o} == set(o)})
        require(len(fixed) == 2, "two axial poles")
        swap = next(a for a in normalizer if permutation_order(a) == 2 and {a[v] for v in fixed[0]} == set(fixed[1]))
        require(row == dict(prime=prime, generator=list(g), normalizer_size=len(normalizer),
                            domain_axis_poles=[list(x) for x in fixed], pole_swap_involution=list(swap)), "finite normalizer assumptions")
        require(len(normalizer) == {2: 4, 3: 6, 5: 10}[prime], "normalizer census")
    require(s["source_scope"] == "oriented incidence, fixed fibre/slot section and process links", "bounded scope")
    require(s["quotient_warning"] == "Product/diagonal actions are mathematical controls; no source law identifies independent local and support chart changes.", "quotient boundary")
    tower = geometry["tower"]()
    check_refinements(group, tower, p["native_support_refinement"])
    # Independent nerve/section construction from the live support stages.
    summaries, previous = [], None
    committed_rows = json.loads((ROOT / "code/m1_source_realization/receipt.json").read_text())["evidence"]["topology"]
    require(len(committed_rows) == 4, "committed topology coverage")
    for r, row in enumerate(tower):
        n = row["vertices"]
        carriers = [(v, slot) for v in range(n) for slot in range(r+2 if v == 0 else 1)]
        section = [carriers.index((v, 0)) for v in range(n)]
        bridge = [v for v, slot in carriers]
        z = [[section[v] for v in f] for f in row["faces"]]
        require(integral_chain(z, bridge) == integral_chain(row["faces"]), "bridge signed cycle")
        if previous is not None:
            old_carriers, old_section, old_z = previous
            c = row["coarsen"]
            q = [old_carriers.index((c[v], slot) if (v, slot) in old_carriers else (c[v], 0))
                 for v, slot in carriers]
            require([old_carriers[q[i]][0] for i in range(len(carriers))] == [c[v] for v in bridge], "bridge square")
            require([q[i] for i in section] == [old_section[v] for v in c], "section square")
            require(integral_chain(z, q) == integral_chain(old_z), "nerve coarse cycle")
        else:
            q = None
        expected_row = dict(support=json.loads(json.dumps(row)), fibres=[r+2]+[1]*(n-1),
                            carriers=[list(v) for v in carriers], coarsen=q, section=section,
                            process_links=[[i, j] for i, (v, _) in enumerate(carriers)
                                           for j, (w, _) in enumerate(carriers) if i < j and v == w])
        require(committed_rows[r] == expected_row, "independent live topology/receipt comparison")
        summaries.append(dict(level=r, support_vertices=n, support_faces=len(row["faces"]),
                              nerve_carriers=len(carriers), bridge_chain_degree=1))
        previous = carriers, section, z
    require(p["support_tower"] == summaries, "all support levels")
    require(p["verdicts"] == ["SEED_TORSOR_REMAINS_UNSELECTED", "DECLARED_COARSENING_NONNATURAL", "ATTACHMENT_REMAINS_OPEN"], "unselected verdict")
    require(p["local_refinement"] == "No local incidence tower supplied by ScreenCarrierMapCandidate; support automorphism audit only.", "no invented local tower")
    require(p["common_sphere"] == "REALIZED_SPHERE_COMPARISON_MISSING" and p["physical_selection"] is False, "no physical promotion")
    require(set(p) == {"schema", "base", "base_tree", "implementation_sha256", "sources", "historical_drift", "seed", "support_tower", "native_support_refinement", "verdicts", "local_refinement", "common_sphere", "physical_selection"}, "closed schema")
    return dict(raw_maps=len(group), fixed_support_stabilizer=len(h), left_orbits=len(left),
                subdivision_symmetries=[60]*3, coarsening_symmetries=[x["coarsening_symmetries"] for x in p["native_support_refinement"]],
                verdicts=p["verdicts"])


if __name__ == "__main__":
    require(len(sys.argv) == 1, "no CLI arguments; verifies the committed live-source receipt")
    print(json.dumps(verify(), sort_keys=True))
