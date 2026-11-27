"""Replay source formation, all refusals and all terminal address probabilities."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from model.executable_model.carrier import HERE as CARRIER_DIR
from tools.strict_json import load
from .graph import DONOR_HASHES, History, donor_selection, donor_waits
from .local import check_local, fuel_probability, require
from .source import attempt, bare_floor, branch_law, exact_floor, install, triangles


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE/"results.json"
STRENGTHS = (0, .5, -.5, 8, -8)
TIMES = (0, 1, 7, 37, 150)
CAP = 96
BOARD_SHA = (
    "99870398e3d516f899c600e9bf25ba53c7f513142ba8facbe1cedd94b23f915b",
    "c7d4ad87b09d9158a3d9be60c537d1aaecc69e951e6a9979e505e0162e134d4b",
    "e7273bbd0e814abf49f3a333ecfe58a542262aee7d2399db677d15c19fcac673",
)


def terminal(tree, choices):
    E = len(tree.board.seams)
    en,u = eigh(tree.twice_scaled)
    ls,v = eigh(tree.lap)
    read = np.asarray(tree.read_ports)
    ports = sorted(set(tree.read_ports))
    records = []
    for item in choices:
        x,z = [tree.index[tuple(path)] for path in item["paths"]]
        a,b = [int(tree.sites[i]) for i in (x,z)]
        parent = tuple(item["paths"][0][:-1])
        require(tree.paths[z][:-1] == parent and tree.board.degree[parent[-1]]==5,
                "Unavailable common parent vacancy meter")
        wanted = np.eye(len(read))[[x,z]]
        minus, plus = wanted[0]-wanted[1], wanted[0]+wanted[1]
        np.testing.assert_array_equal(tree.twice_scaled @ minus, minus)
        np.testing.assert_array_equal(tree.lap @ minus, minus)
        require((tree.twice_scaled @ plus)[tree.index[parent]] == -2, "No transmitting exit for bright control")
        depth = len(parent)-1
        floors = []
        for n in (1,8,16,32,CAP):
            floors.append({"retry_cap":n,"matched":float(exact_floor(E,depth,n)),
                           "bare_abs_k_at_most_128":float(bare_floor(E,depth,n))})
        require(min(floors[-1]["matched"],floors[-1]["bare_abs_k_at_most_128"])>.98,
                "Source fails historical spectral precision")
        appointments = []
        for k in STRENGTHS:
            ps,fail = branch_law(E,depth,CAP,-1,.5)
            pp,fail_plus = branch_law(E,depth,CAP,1,.5)
            p_bare,bare_fail = branch_law(E,depth,CAP,-1,fuel_probability(k,matched=False))
            for t in TIMES:
                U = (u*np.exp(-.5j*k*t*en)) @ u.T
                channels = U @ wanted.T
                phase = np.exp(-.5j*k*t)
                np.testing.assert_allclose(channels[:,0]-channels[:,1],phase*minus,atol=4e-12,rtol=0)
                incoherent = (abs(channels[:,0])**2+abs(channels[:,1])**2)/2
                swaps = v @ (np.exp(-abs(k)*t*ls/2)*(v.T @ ((minus**2)/2)))
                arms = {}
                for label,admission,prob,refusals in (
                        ("minus",ps,minus**2/2,fail), ("bare_fuel_minus",p_bare,minus**2/2,bare_fail),
                        ("plus",pp,abs(channels[:,0]+channels[:,1])**2/2,fail_plus),
                        ("dephase_at_release",ps,incoherent,fail), ("random_swaps",ps,swaps,fail)):
                    require(np.isfinite(prob).all() and min(prob)>-2e-12 and abs(sum(prob)-1)<6e-12,
                            f"Invalid {label} conditional law")
                    dist = np.bincount(read,weights=admission*prob,minlength=tree.board.n)
                    require(abs(dist.sum()+sum(refusals.values())-1)<6e-12,"Invalid original-start denominator")
                    pair_read = float(dist[a]+dist[b])
                    if label=="minus":
                        require(pair_read+1e-12>=float(exact_floor(E,depth,CAP)), "All-time floor violated")
                    if label in ("plus","random_swaps","dephase_at_release") and k==8 and t==37:
                        require(pair_read<.94, f"Control {label} does not discriminate at registered appointment")
                    # Store 12-decimal address probabilities and separate (unrounded)
                    # source refusals. Unreachable physical ports have exactly zero weight.
                    arms[label] = {"ports":np.round(dist[ports],12).tolist(),
                                   "home":round(pair_read,12),"refused":refusals}
                appointments.append({"strength":k,"exposure":t,"arms":arms})
        records.append({"donor_column":item["analyzer"],"histories":item["paths"],"prefix_hops":depth,
                        "floor_by_cap":floors,"appointments":appointments})
    return ports,records


def order_control(board, choices, order):
    ranks = {p:i for i,p in enumerate(order)}
    heads = {e:min(e,key=ranks.__getitem__) for e in board.seams}
    # Count all monotone paths exactly without building a dense alternative Hamiltonian.
    paths = [0]*board.n
    paths[0] = 1
    for p in order:
        for q in board.neighbours[p]:
            if ranks[q]>ranks[p]:
                paths[q]+=paths[p]
    return {"reachable_configurations":sum(paths),
            "donor_words_reachable":[all(heads[tuple(sorted((p,q)))]==p
                                         for h in c["paths"] for p,q in zip(h,h[1:])) for c in choices]}


def run():
    check_local()
    boards = []
    selection = None
    for level in (0,1,2):
        fixture = CARRIER_DIR/"carriers"/f"FEDERATION_LEVEL{level}.npz"
        require(sha256(fixture.read_bytes()).hexdigest()==BOARD_SHA[level], "Different spherical carrier")
        board = federation(level)
        tri = triangles(board)
        heads,events = install(board,tri[0])
        for alt in tri[1:]:
            require(install(board,alt)[0]==heads, "Work route changes formed material")
        tree = History(board,heads)
        if level==0:
            selection = donor_selection(tree)
            old = [{"analyzer":c["analyzer"],"settings":[{"k":k,**donor_waits(tree,c,k)}
                                                         for k in (0,.5,1,2,4,8,16,64,128)]}
                   for c in selection["choices"]]
        ports,records = terminal(tree,selection["choices"])
        boards.append({"cells":board.cells,"seams":len(heads),"parking_triangle_orders":len(tri),
                       "builder_actions_on_accept":len(events),"complete_history_dimension":len(tree.paths),
                       "physical_port_order":ports,"all_unlisted_ports":0,
                       "alternative_orders":{
                           "reverse_local":order_control(board,selection["choices"],sorted(range(board.n),key=lambda p:(p//12,-(p%12)%12))),
                           "descending":order_control(board,selection["choices"],list(reversed(range(board.n))))},
                       "records":records})
    return {"schema":"oph.interval_echo.spectral_release.v1",
            "donor":"44f7680df55b4e68b89c301896df3e10237f046d","donor_sha256":DONOR_HASHES,
            "design_sha256":sha256((HERE/"DESIGN.md").read_bytes().replace(b"\r\n",b"\n")).hexdigest(),
            "selection":selection,"archived_waits":old,
            "first_wrong_then_parent_left_right":{str(s):[[str(x) for x in row] for row in attempt(s)[2]] for s in (-1,1)},
            "boards":boards}


def compare(a,b,path="root"):
    if type(a) is not type(b):
        raise ValueError(f"Type differs: {path}")
    if isinstance(a,dict):
        if a.keys()!=b.keys():
            raise ValueError(f"Fields differ: {path}")
        for k,v in a.items():
            compare(v,b[k],path+"/"+k)
    elif isinstance(a,list):
        for i,(x,y) in enumerate(zip(a,b,strict=True)):
            compare(x,y,path+"/"+str(i))
    elif isinstance(a,float):
        # Small source-failure probabilities are real records, not numerical zero.
        if not (np.isfinite(a) and np.isfinite(b)):
            raise ValueError(f"Nonfinite: {path}")
        refusal = "/refused/" in path or path.endswith("/source_refused")
        probability = (refusal or "/arms/" in path or "/fidelity/" in path
                       or "/outcomes/" in path or "/floor_by_cap/" in path
                       or "/returns_37/" in path or path.rsplit("/",1)[-1] in
                       ("dark_weight","two_port_return_37","minimum_bound"))
        if probability and not (0<=a<=1 and 0<=b<=1):
            raise ValueError(f"Out-of-range probability: {path}")
        if refusal and (a==0)!=(b==0):
            raise ValueError(f"Missing source outcome: {path}")
        tol = max(1e-320,2e-12*abs(a)) if refusal else 5e-11
        if abs(a-b)>tol:
            raise ValueError(f"Numerical replay differs: {path}")
    elif a!=b:
        raise ValueError(f"Exact record differs: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record",type=Path,help="create a first receipt; refuses existing file")
    args = parser.parse_args()
    if args.record and args.record.exists():
        parser.error("record exists")
    # JSON canonicalization only normalizes tuple/list for comparison.
    result = json.loads(json.dumps(run(),allow_nan=False))
    if args.record:
        with args.record.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,separators=(",",":"),allow_nan=False)
            f.write("\n")
    else:
        compare(result,load(ARCHIVE))
    for b in result["boards"]:
        print(b["cells"],"cells:",b["complete_history_dimension"],"histories; worst inclusive lower bounds:",
              {k:min(r["floor_by_cap"][-1][k] for r in b["records"]) for k in ("matched","bare_abs_k_at_most_128")})


if __name__=="__main__":
    main()
