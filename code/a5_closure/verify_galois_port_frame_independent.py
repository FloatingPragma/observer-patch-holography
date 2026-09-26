"""Independent replay in Q[phi]/(phi^2-phi-1); never imports the producer.

Reconstructs incidence, both Gram tables and both coordinate realizations.
Uses Cramer's rule, whereas the producer uses Gaussian elimination.
"""
from fractions import Fraction as Q
from pathlib import Path
from collections import Counter
import ast
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "code/a5_closure/manifests/galois_port_frame_reference.json"
SOURCE = "Lean/ObserverPatchHolography/CoreAxioms.lean"
BASELINE = "4511928b1d4bee6bca5601e309b44595a76f797f"
ZERO, ONE, PHI = (Q(0),Q(0)), (Q(1),Q(0)), (Q(0),Q(1))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def add(x,y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x,y):
    return add(x,neg(y))


def mul(x,y):
    a,b = x
    c,d = y
    return a*c+b*d, a*d+b*c+b*d


def sigma(x):
    return x[0]+x[1], -x[1]


def div(x,y):
    n = mul(y,sigma(y))
    need(n[1] == 0 and n[0] != 0, "field denominator")
    return tuple(t/n[0] for t in mul(x,sigma(y)))


def scalar(n):
    return Q(n),Q(0)


def total(xs):
    out = ZERO
    for x in xs:
        out = add(out,x)
    return out


def sign(x):
    a,b = x[0]+x[1]/2, x[1]/2
    if b == 0:
        return (a>0)-(a<0)
    if a == 0 or a*b > 0:
        return 1 if b>0 else -1
    return (1 if a>0 else -1)*(1 if a*a>5*b*b else -1)


def read(x):
    need(type(x) is list and len(x)==2, "scalar shape")
    need(all(type(t) is str and len(t)<100 for t in x), "rational serialization")
    a,b = map(Q,x)
    need(list(map(str,(a,b))) == x, "canonical rationals")
    return a-b,2*b


def matrix(x, n, m):
    need(type(x) is list and len(x)==n and all(type(r) is list and len(r)==m for r in x), "matrix shape")
    return [[read(t) for t in r] for r in x]


def dot(x,y):
    return total(mul(a,b) for a,b in zip(x,y,strict=True))


def determinant(m):
    from itertools import permutations
    result = ZERO
    for p in permutations(range(3)):
        term = ONE
        for i in range(3):
            term = mul(term,m[i][p[i]])
        inversions = sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        result = add(result,neg(term) if inversions%2 else term)
    return result


def source_contract(text):
    """Check the producer's declared dependency boundary without importing it."""
    tree = ast.parse(text)
    permitted = {"dataclasses","fractions","itertools","pathlib","hashlib","json","re"}
    for node in ast.walk(tree):
        if isinstance(node,ast.Import):
            need(all(x.name in permitted for x in node.names), "producer import boundary")
        if isinstance(node,ast.ImportFrom):
            need(node.module in permitted, "producer import boundary")
        if isinstance(node,ast.Constant) and isinstance(node.value,str):
            need(not any(x in node.value.lower() for x in
                         ("portframegram", "27^", "27**", "mass_relation", "coupling_prediction")),
                 "selected Gram or downstream target in producer")
    # This is a targeted mutation gate, not a sandbox for arbitrary Python.


def verify(packet, root=ROOT):
    expected_keys = {"schema","baseline","mathematical_inputs","face_source_sha256",
                     "faces","distances","generator","orientation","family","similarity",
                     "support_attachment","scope","verdict"}
    need(type(packet) is dict and set(packet)==expected_keys, "packet schema")
    need(packet["schema"]=="oph.galois-port-frames.v1" and packet["baseline"]==BASELINE, "baseline/schema")
    need(packet["mathematical_inputs"]==[SOURCE+"#orientedFaces"], "incidence-only input")
    source = (root/SOURCE).read_bytes()
    need(hashlib.sha256(source).hexdigest()==packet["face_source_sha256"], "source drift")
    text = source.decode().split("def orientedFaces :",1)[1].split("/--",1)[0]
    faces = ast.literal_eval(text[text.index("["):].strip())
    need(packet["faces"]==[list(f) for f in faces], "committed oriented faces")
    need(len(faces)==20 and len({frozenset(f) for f in faces})==20, "face census")
    arcs = Counter(e for a,b,c in faces for e in ((a,b),(b,c),(c,a)))
    need(len(arcs)==60 and all(n==1 and arcs[b,a]==1 for (a,b),n in arcs.items()), "fundamental cycle")
    adjacency = [[int((i,j) in arcs) for j in range(12)] for i in range(12)]
    ds = []
    for i in range(12):
        row, frontier, d = [-1]*12, {i}, 0
        while frontier:
            for j in frontier:
                row[j] = d
            frontier = {k for j in frontier for k in range(12) if adjacency[j][k] and row[k]<0}
            d += 1
        need(sorted(row)==[0]+[1]*5+[2]*5+[3], "distance distribution")
        ds.append(row)
    need(packet["distances"]==ds, "abstract adjacency/distance confusion")
    need(packet["orientation"]=={"ambient_relative_to_xyz":-1,"convex_degree":1}, "orientation convention")
    need(packet["support_attachment"]=={"interface":"PORT-GRAM-SUPPORT-ATTACHMENT",
          "status":"NOT_PRESENT","source_theorem":None}, "unsupported attachment assertion")
    need(packet["verdict"]=="GALOIS_FRAME_AMBIGUITY_PERSISTS_WITHOUT_SUPPORT_ATTACHMENT", "verdict boundary")
    need(packet["scope"]=="equivalent within the two Galois-conjugate rank-three realizations", "two-frame scope")
    source_contract((root/"code/a5_closure/galois_port_frame_certificate.py").read_text())
    p = packet["generator"]
    need(type(p) is list and all(type(i) is int for i in p) and sorted(p)==list(range(12)), "generator permutation")
    need(p[0]==0 and p[1]==2, "order-five class marked by one link step")
    power = list(range(12))
    for n in range(1,6):
        power = [p[i] for i in power]
        need((power==list(range(12)))==(n==5), "exact order five")
    cyclic = {tuple(f[k:]+f[:k]) for f in map(list,faces) for k in range(3)}
    need(all(tuple(p[i] for i in f) in cyclic for f in faces), "oriented proper action")
    need(set(packet["family"])=={"plus","minus"}, "both Galois branches required")
    # Coordinates are rebuilt in the phi basis, independently of receipt rows.
    z,o,t = ZERO,ONE,PHI
    vp = [[z,o,t],[o,t,z],[t,z,o],[neg(o),t,z],[z,neg(o),t],[t,z,neg(o)],
          [neg(t),z,o],[z,o,neg(t)],[o,neg(t),z],[neg(t),z,neg(o)],
          [neg(o),neg(t),z],[z,neg(o),neg(t)]]
    summary, computed = {},{}
    branch_keys = {"gram","vectors","norm_squared","rank","trace","psd","rotation",
                   "rotation_trace","adjacency_eigenvalue","laplacian_eigenvalue",
                   "chord_squared_by_distance","chord_distance_order","degree_controls"}
    for name,sg in (("plus",1),("minus",-1)):
        row = packet["family"][name]
        need(set(row)==branch_keys, "branch schema")
        v = vp if sg==1 else [[sigma(x) for x in w] for w in vp]
        need(matrix(row["vectors"],12,3)==v, "actual field automorphism on coordinates")
        s = (Q(-sg),Q(2*sg))
        entries = [ONE,div(s,scalar(5)),neg(div(s,scalar(5))),neg(ONE)]
        g = [[entries[d] for d in r] for r in ds]
        need(matrix(row["gram"],12,12)==g, "incidence-derived signed Gram")
        need(all(g[i][j]==g[j][i] for i in range(12) for j in range(12)), "symmetric Gram")
        need(all(total(mul(g[i][k],g[k][j]) for k in range(12))==mul(scalar(4),g[i][j])
                 for i in range(12) for j in range(12)), "G squared equals 4G")
        need(total(g[i][i] for i in range(12))==scalar(12), "Gram trace")
        norm = dot(v[0],v[0])
        need(sign(norm)>0 and read(row["norm_squared"])==norm, "positive normalization")
        need(all(div(dot(v[i],v[j]),norm)==g[i][j] for i in range(12) for j in range(12)), "Gram realization")
        need(determinant(v[:3])!=ZERO, "rank three witness")
        # G/4 = V V^T/(4 norm): a real sum of squares. No ordered-field
        # positivity is transported through sigma without this new witness.
        need(row["rank"]==3 and row["trace"]==12 and row["psd"] is True, "rank/PSD report")
        need(read(row["adjacency_eigenvalue"])==s and read(row["laplacian_eigenvalue"])==sub(scalar(5),s), "Laplacian sign")
        need(all(total(v[j][k] for j in range(12) if adjacency[i][j])==mul(s,v[i][k])
                 for i in range(12) for k in range(3)), "coordinate eigenmode")
        r = matrix(row["rotation"],3,3)
        need(determinant(r)==ONE and all(dot(r[i],r[j])==(ONE if i==j else ZERO)
                 for i in range(3) for j in range(3)), "proper orthogonal rotation")
        need(all(dot(r[k],v[i])==v[p[i]][k] for i in range(12) for k in range(3)), "carrier rotation action")
        trace = total(r[i][i] for i in range(3))
        need(trace==(PHI if sg==1 else sigma(PHI)) and read(row["rotation_trace"])==trace, "order-five trace")
        chords = [mul(scalar(2),sub(ONE,x)) for x in entries]
        need(list(map(read,row["chord_squared_by_distance"]))==chords, "chord labels")
        order = [1,2,3] if sg==1 else [2,1,3]
        need(row["chord_distance_order"]==order and all(sign(sub(chords[b],chords[a]))>0
             for a,b in zip([0]+order,order)), "chord ranking")
        need(all(div(dot([sub(a,b) for a,b in zip(v[i],v[j])],
                         [sub(a,b) for a,b in zip(v[i],v[j])]),norm)==chords[ds[i][j]]
                 for i in range(12) for j in range(12)), "all pairwise squared chords")
        controls = row["degree_controls"]
        need(len(controls)==3 and [c["direction"] for c in controls]==[[1,2,4],[-3,5,2],[7,-2,3]], "regular value coverage")
        for control in controls:
            need(set(control)=={"direction","coefficients","signed_preimages","raw_degree","degree"}, "degree schema")
            y = list(map(scalar,control["direction"]))
            expected, signs = [],[]
            for face in faces:
                columns = [v[i] for i in face]
                denominator = determinant(columns)
                need(denominator!=ZERO, "face cone has no zero on simplex")
                coeffs = []
                for k in range(3):
                    replacement = [y if j==k else w for j,w in enumerate(columns)]
                    coeffs.append(div(determinant(replacement),denominator))
                need(all(x!=ZERO for x in coeffs), "no face boundary at regular value")
                expected.append(coeffs)
                signs.append(sign(denominator) if all(sign(x)>0 for x in coeffs) else 0)
            need(matrix(control["coefficients"],20,3)==expected, "Cramer preimage coefficients")
            need(control["signed_preimages"]==signs and control["raw_degree"]==sum(signs), "signed degree sum")
            need(control["degree"]==-sum(signs)==(1 if sg==1 else 7), "oriented degree pair")
        computed[name] = (g,v,r)
        summary[name] = {"rank":3,"raw_degree":controls[0]["raw_degree"],
                         "degree":controls[0]["degree"],"edge_chord_rank":order.index(1)+1}
    gp,_,rp = computed["plus"]
    gm,vm,rm = computed["minus"]
    need([[sigma(x) for x in row] for row in gp]==gm and
         [[sigma(x) for x in row] for row in rp]==rm, "Galois Gram/rotation exchange")
    similarity = packet["similarity"]
    need(set(similarity)=={"scale","coordinate_permutation","port_relabel"}, "similarity schema")
    need(read(similarity["scale"])==PHI and sign(PHI)>0 and
         similarity["coordinate_permutation"]==[1,0,2], "positive similarity")
    relabel = similarity["port_relabel"]
    need(sorted(relabel)==list(range(12)), "similarity bijection")
    need(all([mul(PHI,v[k]) for k in (1,0,2)]==vp[relabel[i]] for i,v in enumerate(vm)), "unlabeled point identity")
    need(all(ds[relabel[i]][relabel[j]]==[0,2,1,3][ds[i][j]] for i in range(12) for j in range(12)), "nearest/second-nearest exchange")
    return {"verified":summary,"scope":packet["scope"],"verdict":packet["verdict"]}


def load(path):
    def unique(pairs):
        out = {}
        for key,value in pairs:
            need(key not in out,"duplicate JSON key")
            out[key] = value
        return out
    return json.loads(Path(path).read_text(),object_pairs_hook=unique,
                      parse_constant=lambda _: need(False,"nonfinite JSON"))


if __name__ == "__main__":
    print(json.dumps(verify(load(sys.argv[1] if len(sys.argv)>1 else REFERENCE)),sort_keys=True))
