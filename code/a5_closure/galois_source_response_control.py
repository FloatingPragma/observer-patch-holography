"""Bounded Galois control of the integrated source-selection response tape.

Uses the existing independent checker's scalar/matrix primitives, not either
source producer. Conjugation fixes i and exchanges both real embeddings in
the entire source law, including its quadrupole. It is not a block swap or
a positivity-preserving map on arbitrary quantum states.
"""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import hashlib
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "code/source_selection_model"
SOURCE_PINS = {
    "DERIVATION.md": "4881f7c7e41b6294d3b164e943952438f09a3aa03041bf152a382d8ff42aba85",
    "geometry.py": "de47ae19dd4ae1677d0cefb0163513f2f1f6d24dfb4dd0264506d0b17fceb843",
    "response.json": "226ab001894f357d3732d8a898fbae80e28b9e16ca41104a6f293190b441e053",
    "verify_response.py": "1aaa22195b7241aaf61f0bda60940907f2416f406c53f81ebb3418005d95900e",
}
spec = importlib.util.spec_from_file_location("galois_existing_response_checker", MODEL/"verify_response.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def conjugate_serialized(x):
    if isinstance(x,dict):
        return {k:conjugate_serialized(v) for k,v in x.items()}
    if isinstance(x,list):
        if len(x) in (2,4) and all(type(v) is str for v in x):
            return [str(-Q(v)) if i%2 else v for i,v in enumerate(x)]
        return [conjugate_serialized(v) for v in x]
    return x


def replay(packet, conjugated=False):
    """Replay fixed abstract actions with independently reconstructed D_p."""
    c.keys(packet,("schema","generators","matrix_unit_response_probes",
                   "ordered_mixed_responses","cayley_factors","closed_paths"))
    c.require(packet["schema"]=="oph.source-selection-response.v1","source schema")
    original = c.vertices()
    vs = [tuple(x.conj() for x in v) for v in original] if conjugated else original
    gs = [c.blocks(g) for g in packet["generators"]]
    c.require(len(gs)==12,"twelve source directions")
    expected = []
    norm = c.dot(vs[0],vs[0])
    for v in vs:
        blocks = []
        for branch in (0,1):
            w = [x.conj()/c.F5(2) if branch else x/c.F5(2) for x in v]
            h = [[c.Z,-w[2],w[1]],[w[2],c.Z,-w[0]],[-w[1],w[0],c.Z]]
            blocks.append([[c.C5(h[i][j],c.Z if branch else
                                v[i]*v[j]/c.F5(2)+(c.F5(Q(1,12))-norm/c.F5(6) if i==j else c.Z))
                            for j in range(3)] for i in range(3)])
        expected.append(tuple(blocks))
    c.require(gs==expected,"conjugated source formula including quadrupole")
    flat = lambda g: [part for b in g for row in b for x in row for part in (x.re,x.im)]
    c.require(c.rank([flat(g) for g in gs])==12,"response rank twelve")
    c.require(all(c.star(b)==c.neg(b) for g in gs for b in g),"skew-adjoint source")
    response_rows, restricted_rows = [],[]
    c.require(len(packet["matrix_unit_response_probes"])==12,"tomography coverage")
    for g,tape in zip(gs,packet["matrix_unit_response_probes"]):
        full = [[g[i//3][i%3][j%3] if i//3==j//3 else c.CZ for j in range(6)] for i in range(6)]
        values = [full[i][j] for i in range(6) for j in range(6) if i!=j]
        values += [full[i][i]-full[5][5] for i in range(5)]
        observed = [c.C5(c.real(x[:2]),c.real(x[2:])) for x in tape]
        c.require(observed==values,"all matrix-unit probes")
        response_rows.append([part for x in observed for part in (x.re,x.im)])
        restricted = [b[i][j] for b in g for i in range(3) for j in range(3) if i!=j]
        restricted += [b[i][i]-b[2][2] for b in g for i in range(2)]
        restricted_rows.append([part for x in restricted for part in (x.re,x.im)])
    c.require(c.rank(response_rows)==12 and c.rank(restricted_rows)==11,"complete/restricted tomography ranks")
    pairs = list(combinations(range(12),2))
    c.require(len(packet["ordered_mixed_responses"])==66,"ordered bracket coverage")
    for row,(p,q) in zip(packet["ordered_mixed_responses"],pairs):
        c.require(row["ports"]==[p,q],"ordered pair")
        bracket = tuple(c.add(c.mm(a,b),c.neg(c.mm(b,a))) for a,b in zip(gs[p],gs[q]))
        c.require(c.blocks(row["mixed_coefficient"])==bracket and
                  c.combine(gs,list(map(c.real,row["in_port_basis"])))==bracket,"ordered response closure")
    anti = [vs.index(tuple(-x for x in v)) for v in vs]
    odd = [[gs[p][b][i][j].re-gs[anti[p]][b][i][j].re
            for b in range(2) for i,j in ((0,1),(0,2),(1,2))] for p in range(12) if p<anti[p]]
    c.require(c.rank(odd)==6,"endogenous rotation logarithms")
    # The six even imaginary symmetric directions fill Sym(3). Together
    # with the six odd directions this is the full closed declared tangent.
    even = [[gs[p][0][i][j].im+gs[anti[p]][0][i][j].im
             for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))]
            for p in range(12) if p<anti[p]]
    c.require(c.rank(even)==6,"all even symmetric directions")
    phi = c.F5(Q(1,2),Q(1,2))
    arcs = {(i,j) for i in range(12) for j in range(12) if c.dot(original[i],original[j])==phi}

    def action(u,p):
        c.require([[x.re.conj() for x in row] for row in u[0]]==
                  [[x.re for x in row] for row in u[1]],"Galois pair of implementers")
        for b in u:
            c.require(all(x.im==c.Z for row in b for x in row) and
                      c.mm(c.star(b),b)==c.eye() and
                      c.determinant([[x.re for x in row] for row in b])==c.O,"proper carrier rotation")
        c.require({(p[i],p[j]) for i,j in arcs}==arcs,"fixed abstract incidence action")
        for i,v in enumerate(vs):
            c.require(tuple(sum((u[0][k][j].re*v[j] for j in range(3)),c.Z)
                            for k in range(3))==vs[p[i]],"executed chart action")

    factors = packet["cayley_factors"]
    c.require(len(factors)==20,"order-three factor coverage")
    executed, actions = [],[]
    for f in factors:
        p = c.permutation(f["port_action"])
        c.require(p!=tuple(range(12)) and tuple(p[p[p[i]]] for i in range(12))==tuple(range(12)),"factor order")
        k = c.combine(gs,list(map(c.real,f["control"])))
        u = c.blocks(f["observed_response"])
        for kb,ub in zip(k,u):
            c.require(c.star(kb)==c.neg(kb) and c.mm(ub,c.add(c.eye(),c.neg(kb)))==c.add(c.eye(),kb),"Cayley response execution")
        action(u,p)
        executed.append(u)
        actions.append(p)
    c.require(len(set(actions))==20,"distinct factors")
    table = {}
    c.require(len(packet["closed_paths"])==60,"holonomy coverage")
    for row in packet["closed_paths"]:
        p = c.permutation(row["port_action"])
        u,perm = (c.eye(),c.eye()),tuple(range(12))
        word = row["factor_indices"]
        c.require(len(word)<=2 and all(type(i) is int and 0<=i<20 for i in word),"response word")
        for k in word:
            u = c.bmul(executed[k],u)
            perm = tuple(actions[k][i] for i in perm)
        c.require(perm==p and u==c.blocks(row["observed_response"]) and p not in table,"closed-path replay")
        action(u,p)
        c.require(all(c.bmul(c.bmul(u,gs[i]),tuple(c.star(b) for b in u))==gs[p[i]]
                      for i in range(12)),"same-response naturality")
        table[p] = u
    for g,u in table.items():
        for h,v in table.items():
            c.require(c.bmul(u,v)==table[tuple(g[h[i]] for i in range(12))],"all path compositions")
    return {"response_rank":12,"observable_rank":12,"restricted_observable_rank":11,
            "even_symmetric_rank":6,"odd_rotation_rank":6,"brackets":66,
            "cayley_factors":20,"proper_actions":60,"naturality_pairs":720,"compositions":3600}


def verify():
    c.require(all(hashlib.sha256((MODEL/name).read_bytes()).hexdigest()==digest
                  for name,digest in SOURCE_PINS.items()),"audited source drift")
    packet = c.strict_load(MODEL/"response.json")
    swapped = conjugate_serialized(packet)
    c.require(conjugate_serialized(swapped)==packet,"Galois involution on full tape")
    plus,minus = replay(packet),replay(swapped,True)
    c.require(plus==minus,"same audited clause counts")
    literal_swap = [[g[1],g[0]] for g in packet["generators"]]
    c.require(literal_swap!=swapped["generators"],"block swap is not field conjugation")
    return {"scope":"source_selection_model response formula, tomography, closure, Cayley holonomy and naturality",
            "plus":plus,"minus":minus,"literal_block_swap_equals_galois":False,
            "positivity_argument":"Each replay is skew-adjoint and injective over the real embedding; -Tr(X^2)>0 for nonzero X.",
            "bounded_nonselection":True,
            "source_pins":SOURCE_PINS}


if __name__ == "__main__":
    print(json.dumps(verify(),indent=2,sort_keys=True))
