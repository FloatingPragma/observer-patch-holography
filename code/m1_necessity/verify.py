"""Producer-free verification of every claimed reference graph, route and cut."""

import argparse
import hashlib
from itertools import product
from pathlib import Path

from . import check, balanced_check
from .audit import consumer_audit, downstream_audit
import pin_scope

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def pins():
    names = [p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
             if p.suffix in (".py",".md") or p.name in ("consumers.json","downstream.json")]
    names += ["Lean/Geometry/M1Necessity.lean","Lean/Geometry/M1NecessityAxiomAudit.lean",
              "paper/tex_fragments/CROSSING_READ_AREA_LAW.tex",
              "paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex",
              "paper/tex_fragments/SOURCE_METRIC_SCALAR_CONTINUUM.tex",
              "paper/tex_fragments/M1_NECESSITY.tex",
              "Lean/Geometry/FlatDiamondNormalization.lean"]
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(pin_scope.bound(names))}


def route_targets():
    # The complete fixed catalog is independent of the producer's emitted rows.
    targets = []
    for n in range(1,9):
        r,q = check.level(n)
        targets.extend((n,d) for d in [(0,0,0),(q-1,0,0),(q//2,-q//3,q//5),
                                      (-q+1,q//2,-1),(1,1,1),(r-1,r+1,-2)])
    targets.extend((n,d) for n in (1,2) for d in product(range(-2,3),repeat=3))
    return targets


def verify_routes(rows):
    check.require(type(rows) is list and len(rows) == len(route_targets()), "incomplete route catalog")
    sets = {n:check.expected_sparse(n) for n in range(1,9)}
    total = 0
    for row,(n,d) in zip(rows,route_targets()):
        check.require(type(row) is dict and set(row) == {"n","displacement","instructions","steps","points_sha256"},
                      "route schema")
        check.require(check.canonical([row["n"],row["displacement"]]) == check.canonical([n,list(d)]),
                      "substituted route target")
        check.require(type(row["instructions"]) is list, "route instructions required")
        expanded = []
        for instruction in row["instructions"]:
            check.require(type(instruction) is list and len(instruction) == 4 and
                          all(type(x) is int for x in instruction), "route instruction schema")
            count,*v = instruction
            check.require(1 <= count <= 10000 and len(expanded)+count <= 10000,
                          "invalid repetition count")
            expanded.extend([tuple(v)]*count)
        actual = check.check_route(n,d,expanded,sets[n])
        check.require(check.canonical(actual) == check.canonical({k:row[k] for k in actual}),
                      "route ledger mismatch")
        total += actual["steps"]
    return total


def summarize(sites,outgoing,states):
    return {"states":[check.digest(s) for s in states],"events":len(sites)*len(states),
            "reads":3*sum(map(len,outgoing))}


def intervention_replay(sites,outgoing,initial,layer,site):
    states = [list(initial)]
    index = sites.index(tuple(site))
    for j in range(4):
        if j:
            new = [1]*len(sites)
            for i,value in enumerate(states[-1]):
                for target in outgoing[i]:
                    new[target] += value
            states.append(new)
        if j == layer:
            states[-1][index] += 1
    return states


def graph_checks(n,family):
    sites,outgoing,baseline = check.graph_replay(n,family)
    _,q = check.level(n)
    result = {"baseline":summarize(sites,outgoing,baseline)}
    for name,layer in (("initial",0),("intermediate",1)):
        changed = intervention_replay(sites,outgoing,baseline[0],layer,[q//2]*3)
        reached = {sites.index((q//2,)*3)}
        for j in range(4):
            different = {i for i,(a,b) in enumerate(zip(baseline[j],changed[j])) if a != b}
            check.require(different == (reached if j >= layer else set()), "writer intervention disagrees with ancestry")
            if j >= layer:
                reached = {target for i in reached for target in outgoing[i]}
        result[name] = summarize(sites,outgoing,changed)
    counts = {}
    for name,diagonal in (("coordinate",False),("diagonal",True)):
        threshold = 2*q-1 if diagonal else q-1
        projection = [2*(x[0]+x[1] if diagonal else x[0]) for x in sites]
        counts[name] = sum(1 for i,targets in enumerate(outgoing) if projection[i] < threshold
                           for target in targets if projection[target] > threshold)
    return result,counts


def verify(packet):
    check.require(type(packet) is dict and set(packet) == {"schema","sources","consumers","cuts","routes","graphs","balanced","downstream"},
                  "exact receipt schema required")
    check.require(packet["schema"] == "oph-m1-necessity-v1", "schema version")
    check.require(check.canonical(pin_scope.bound_pins(packet["sources"])) == check.canonical(pins()), "scientific source custody changed")
    check.require(check.canonical(packet["consumers"]) == check.canonical(consumer_audit()), "consumer evidence changed")
    downstream=downstream_audit()
    check.require(check.canonical(packet['downstream'])==check.canonical(downstream),'downstream contract evidence changed')
    levels = (1,2,3,4,6,8,12,16,20)
    check.require(type(packet["cuts"]) is dict and set(packet["cuts"]) == set(map(str,levels)), "cut catalog changed")
    for n in levels:
        submitted = packet["cuts"][str(n)]
        vectors = sorted(check.expected_sparse(n))
        cuts = check.expected_cuts(n,vectors)
        families = ("sparse","axial")+(("dense",) if n<=4 else ())
        check.require(type(submitted) is dict and set(submitted) == {"stencil_sha256","cuts","actions"}, "cut entry schema")
        check.require(submitted["stencil_sha256"] == check.digest(vectors), "stencil differs from independent definition")
        check.require(check.canonical(submitted["cuts"]) == check.canonical(cuts), "false crossing count")
        check.require(type(submitted["actions"]) is dict and set(submitted["actions"]) == set(families), "action catalog")
        for family in families:
            check.check_action(n,family,submitted["actions"][family])
    total_steps = verify_routes(packet["routes"])
    expected_keys = {f"{n}:{f}" for n in (1,2) for f in ("sparse","axial","dense")}
    check.require(type(packet["graphs"]) is dict and set(packet["graphs"]) == expected_keys, "graph catalog changed")
    total_events,total_reads = 0,0
    for n in (1,2):
        for family in ("sparse","axial","dense"):
            graph,cuts = graph_checks(n,family)
            check.require(check.canonical(packet["graphs"][f"{n}:{family}"]) == check.canonical(graph), "false graph execution")
            for cut,value in cuts.items():
                check.require(value == packet["cuts"][str(n)]["cuts"][f"{family}_{cut}"], "full pair census disagrees with cut formula")
            total_events += sum(item["events"] for item in graph.values())
            total_reads += sum(item["reads"] for item in graph.values())
    return {"routes":len(packet["routes"]),"route_steps":total_steps,
            "executed_events":total_events,"executed_reads":total_reads,"consumer_claims":8,
            "downstream_claims":len(downstream['claims']),**balanced_check.verify(packet['balanced'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt",type=Path,default=HERE/"receipt.json")
    args = parser.parse_args()
    print(check.canonical(verify(check.strict_load(args.receipt))).decode())


if __name__ == "__main__":
    main()
