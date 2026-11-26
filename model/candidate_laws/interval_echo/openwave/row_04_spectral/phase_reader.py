"""Distinguish both signed material histories by a destructive terminal read on fresh copies."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.executable_model import federation
from tools.strict_json import load
from .experiment import HERE, CAP, compare
from .graph import History, donor_selection
from .local import check_local, require
from .source import branch_law, edge, install, triangles


def effects(tree,paths,sign):
    """Effect of every terminal record; only the prefix head, local flags and positions are read."""
    parent = tuple(paths[0][:-1])
    p,a,b = parent[-1],paths[0][-1],paths[1][-1]
    require(parent==tuple(paths[1][:-1]) and tree.board.degree[p]==5,"No common native center meter")
    flag_edge=(0,1)  # The donor's paths arrive directly or after flipping this first seam.
    flipped=lambda w: any(edge(v,q)==flag_edge for v,q in zip(w,w[1:]))
    same = np.array([flipped(w)==flipped(parent) for w in tree.paths])
    require([w for w,v in zip(tree.paths,same) if v and w[-1]==p]==[parent],
            "The observed prefix and occupied-center record is not unique")
    dim=len(tree.paths)
    ident = np.eye(dim)
    select=np.diag(same.astype(int))
    effects=[ident-select]  # wrong prefix: neither leaf has been measured
    K=select

    def symmetric(h):
        # The projector is a direct sum over ALL admissible local hops on this
        # physical edge, including any other parent with identical physical address.
        op=np.zeros_like(ident)
        for j,w in enumerate(tree.paths[1:],1):
            if edge(w[-1],w[-2])==h:
                i=tree.index[w[:-1]]
                op[np.ix_([i,j],[i,j])]=.5
        return op

    center=np.diag((tree.sites==p).astype(int))
    ca,cb=symmetric(edge(p,a)),symmetric(edge(p,b))
    for keep in (ident-center,ident-ca,ident-center,ident-cb,ca if sign==-1 else ident-ca):
        wrong=(ident-keep) @ K
        effects.append(wrong.T @ wrong)
        K=keep @ K
    effects.extend([K.T @ (ident-center) @ K,K.T @ center @ K])
    # Dyadic operations in binary64 are exact here (five factors of at most 1/2).
    v=ident[:,tree.index[tuple(paths[0])]]+sign*ident[:,tree.index[tuple(paths[1])]]
    np.testing.assert_array_equal(16*effects[-1],np.outer(v,v))
    np.testing.assert_array_equal(sum(effects),ident)
    return np.array(effects)


def run():
    check_local()
    labels=["wrong_prefix","reverse_refusal_1","reverse_refusal_2","reverse_refusal_3",
            "reverse_refusal_4","reverse_refusal_5","parent_vacant","signed_ray_click"]
    boards=[]
    choices=None
    for l in range(3):
        board=federation(l)
        tree=History(board,install(board,triangles(board)[0])[0])
        if choices is None:
            choices=donor_selection(tree)["choices"]
        records=[]
        for read in choices:
            for read_sign in (-1,1):
                Es=effects(tree,read["paths"],read_sign)
                for release in choices:
                    target=[tree.index[tuple(w)] for w in release["paths"]]
                    en,u=eigh(tree.twice_scaled)
                    U=(u*np.exp(-4j*37*en)) @ u.T
                    cols=U[:,target]
                    for name,sign in (("minus",-1),("plus",1),("dephased",-1)):
                        if name=="dephased":
                            dist=np.einsum("gi,kgh,hi->k",cols.conj(),Es,cols).real/2
                        else:
                            v=(cols[:,0]+sign*cols[:,1])/np.sqrt(2)
                            dist=np.einsum("g,kgh,h->k",v.conj(),Es,v).real
                        mass,fail=branch_law(len(board.seams),len(release["paths"][0])-2,CAP,sign,.5)
                        require(min(dist)>-1e-12 and abs(sum(dist)-1)<3e-12,"Incomplete phase receiver")
                        if read_sign==-1:
                            weight=(1 if name=="minus" else .5 if name=="dephased" else 0)
                            expected=weight*int(read["analyzer"]==release["analyzer"])/8
                            np.testing.assert_allclose(dist[-1],expected,rtol=0,atol=2e-12)
                        records.append({"reader":read["analyzer"],"read_sign":read_sign,
                                        "source":release["analyzer"],"arm":name,
                                        "outcomes":(mass*dist).round(12).tolist(),"source_refused":sum(fail.values())})
        boards.append({"level":l,"records":records})
    return {"strength":8,"exposure":37,"record_order":labels,"boards":boards}


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--record",type=Path)
    args=cli.parse_args()
    data=json.loads(json.dumps(run(),allow_nan=False))
    if args.record:
        with args.record.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(data,f,separators=(",",":"),allow_nan=False)
            f.write("\n")
    else:
        compare(data,load(HERE/"phase.json"))
    print("Verified both prefix and sign effects on the complete 32-, 70- and 70-word sectors.")
