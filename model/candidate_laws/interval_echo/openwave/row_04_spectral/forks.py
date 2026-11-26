"""Does the source itself enforce persistence? Compile both signs at every eligible fork."""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from .graph import History, donor_selection
from .local import check_local, require
from .source import edge, install, triangles, attempt, exact_floor, bare_floor


def all_forks(tree, choices):
    # These are physical seam queries acting on every history, not matrices
    # addressed to a particular hidden path index.
    n = len(tree.paths)
    transfers = defaultdict(list)
    children = defaultdict(list)
    for j,w in enumerate(tree.paths[1:],1):
        i = tree.index[w[:-1]]
        children[i].append(j)
        transfers[edge(w[-1],w[-2])].append((i,j))
    def measure(e, state, selected):
        out = np.zeros_like(state)
        for i,j in transfers[e]:
            out[i] = out[j] = (state[i]+state[j])/2
        return out if selected else state-out

    vals,u = eigh(tree.twice_scaled)
    forward = (u*np.exp(-4j*37*vals)) @ u.T
    read = np.array(tree.read_ports)
    enrolled = [p["paths"] for p in choices]
    rows = []
    for i, sons in children.items():
        p = tree.paths[i][-1]
        if tree.board.degree[p] != 5:
            continue  # The native vacant-reference center meter must exist.
        for a,b in combinations(sons,2):
            tips = tree.sites[[a,b]]
            e1,e2 = [edge(p,int(port)) for port in tips]
            receiver = np.isin(read,tips)
            result = []
            for sign in (-1,1):
                state = np.zeros(n)
                state[i] = 1
                loss = []
                for contact,answer in ((e1,sign==-1),(e2,False),(None,False),(e1,False),(None,False)):
                    kept = (state*(read!=p)) if contact is None else measure(contact,state,answer)
                    other = state-kept
                    local = other[[i,a,b]]**2
                    require(local.sum()==other @ other, "Failure leaves the three source histories")
                    loss.append(local.tolist())
                    state = kept
                success,_,fail = attempt(sign)
                require(np.array_equal(np.asarray(fail,dtype=float),loss) and state@state==success,
                        "Full physical queries disagree with source recurrence")
                v = np.zeros(n,dtype=int)
                v[a],v[b] = 1,sign
                require(np.array_equal(state, v*state[a]),"Compiled a different wave")
                # Column iff, allowing any eigenvalue, retaining every possible destination.
                f = tree.twice_scaled @ v
                eigen = np.array_equal(f,f[a]*v)
                return_p = float(sum(abs((forward @ v)[receiver])**2)/2)
                result.append((eigen,return_p))
            rows.append({"histories":[tree.paths[a],tree.paths[b]],"donor":[tree.paths[a],tree.paths[b]] in enrolled,
                         "minus_eigenray":result[0][0],"minus_conditional_home":result[0][1],
                         "plus_conditional_home":result[1][1]})
    return rows


def summary():
    check_local()
    choices = None
    out = []
    for l in range(3):
        board = federation(l)
        tree = History(board,install(board,triangles(board)[0])[0])
        if choices is None:
            choices = donor_selection(tree)["choices"]
        rows = all_forks(tree,choices)
        other = [r for r in rows if not r["donor"]]
        require(sum(r["donor"] for r in rows)==2 and min(r["minus_conditional_home"] for r in other)<.94,
                "The source did not provide a negative on a different fork")
        out.append({"cells":board.cells,"compiled_both_signs":len(rows),
                    "minus_exact_eigenrays":sum(r["minus_eigenray"] for r in rows),
                    "other_forks_below_94":sum(r["minus_conditional_home"]<.94 for r in other),
                    "lowest_other":min(other,key=lambda r:r["minus_conditional_home"]),
                    "minimum_caps_from_uniform_bounds":{str(g):[next(n for n in range(1,97) if bound(len(board.seams),2,n)>=g)
                                                         for bound in (exact_floor,bare_floor)]
                                                       for g in (Fraction(47,50),Fraction(49,50))}})
    return out


if __name__=="__main__":
    import json
    print(json.dumps(summary(),indent=2,allow_nan=False))
