"""Two incidence-bound frames. Exact arithmetic; no selected Gram import.

Emit to stdout; verification and source-model controls are separate commands.
The sole mathematical source input is CoreAxioms.orientedFaces.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[2]
BASELINE = "4511928b1d4bee6bca5601e309b44595a76f797f"
FACE_SOURCE = "Lean/ObserverPatchHolography/CoreAxioms.lean"


@dataclass(frozen=True)
class F:
    """a + b sqrt(5), in the real embedding sqrt(5)>0."""
    a: Q = Q(0)
    b: Q = Q(0)

    def __post_init__(self):
        object.__setattr__(self, "a", Q(self.a))
        object.__setattr__(self, "b", Q(self.b))

    def __add__(self, other):
        other = other if isinstance(other, F) else F(other)
        return F(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        other = other if isinstance(other, F) else F(other)
        return F(self.a*other.a + 5*self.b*other.b,
                 self.a*other.b + self.b*other.a)

    __rmul__ = __mul__

    def sigma(self):
        return F(self.a, -self.b)

    def __truediv__(self, other):
        other = other if isinstance(other, F) else F(other)
        norm = other.a**2 - 5*other.b**2
        if not norm:
            raise ValueError("zero denominator")
        return self * F(other.a/norm, -other.b/norm)

    def sign(self):
        if not self.b:
            return (self.a > 0) - (self.a < 0)
        if not self.a or self.a*self.b > 0:
            return 1 if self.b > 0 else -1
        # Opposite signs: compare rational squares, never floating roots.
        norm = self.a**2 - 5*self.b**2
        return (1 if self.a > 0 else -1) * (1 if norm > 0 else -1)


Z, O, S = F(), F(1), F(0, 1)
PHI = F(Q(1, 2), Q(1, 2))


def encode(x):
    if isinstance(x, F):
        return [str(x.a), str(x.b)]
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    return x


def require(ok, message):
    if not ok:
        raise ValueError(message)


def transpose(a):
    return list(map(list, zip(*a)))


def dot(a, b):
    return sum((x*y for x, y in zip(a, b, strict=True)), Z)


def mm(a, b):
    return [[dot(row, col) for col in transpose(b)] for row in a]


def eye(n):
    return [[O if i == j else Z for j in range(n)] for i in range(n)]


def det(m):
    a, b, c = m
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def inverse(m):
    n = len(m)
    rows = [a+b for a, b in zip(m, eye(n))]
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j] != Z)
        rows[j], rows[pivot] = rows[pivot], rows[j]
        scale = rows[j][j]
        rows[j] = [x/scale for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x-scale*y for x, y in zip(rows[i], rows[j])]
    return [r[n:] for r in rows]


def faces():
    text = (ROOT / FACE_SOURCE).read_text()
    body = text.split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return [list(map(int, f)) for f in re.findall(r"\((\d+), (\d+), (\d+)\)", body)]


def distances(fs):
    adjacent = [{b for f in fs if a in f for b in f if b != a} for a in range(12)]
    out = []
    for a in range(12):
        row = [99]*12
        row[a] = 0
        for _ in range(3):
            row = [min(row[b], min(row[c]+1 for c in adjacent[b])) for b in range(12)]
        out.append(row)
    require(all(sorted(r) == [0]+[1]*5+[2]*5+[3] for r in out), "distance census")
    return out


def coordinates():
    p = PHI
    return [[Z,O,p],[O,p,Z],[p,Z,O],[-O,p,Z],[Z,-O,p],[p,Z,-O],
            [-p,Z,O],[Z,O,-p],[O,-p,Z],[-p,Z,-O],[-O,-p,Z],[Z,-O,-p]]


def oriented_key(face):
    a,b,c = face
    return min((a,b,c),(b,c,a),(c,a,b))


def proper_actions(ds, fs):
    """Enumerate abstract incidence automorphisms before using coordinates."""
    target_faces = {oriented_key(f) for f in fs}
    found = []

    def extend(partial):
        if len(partial) == 12:
            p = [partial[i] for i in range(12)]
            if {oriented_key([p[i] for i in f]) for f in fs} == target_faces:
                found.append(p)
            return
        i = max((i for i in range(12) if i not in partial),
                key=lambda i: (sum(ds[i][j] == 1 for j in partial), -i))
        for a in range(12):
            if a not in partial.values() and all(ds[i][j] == ds[a][b] for j,b in partial.items()):
                extend({**partial, i:a})

    extend({})
    require(len(found) == 60, "proper action census")
    return sorted(found)


def degree(v, fs, y):
    coefficients, signs = [], []
    for face in fs:
        matrix = transpose([v[i] for i in face])
        require(det(matrix) != Z, "nondegenerate face cone")
        c = [dot(row, list(map(F, y))) for row in inverse(matrix)]
        require(all(x != Z for x in c), "direction on a face plane")
        coefficients.append(c)
        signs.append(det(matrix).sign() if all(x.sign() > 0 for x in c) else 0)
    return {"direction": y, "coefficients": coefficients,
            "signed_preimages": signs, "raw_degree": sum(signs),
            "degree": -sum(signs)}


def build():
    fs = faces()
    require(len(fs) == 20 and len({tuple(sorted(f)) for f in fs}) == 20, "face census")
    ds = distances(fs)
    group = proper_actions(ds, fs)
    # One link step at fixed port 0; not the other order-five conjugacy class.
    generator = next(p for p in group if p[0] == 0 and p[1] == 2)
    power = list(range(12))
    for n in range(1, 6):
        power = [generator[i] for i in power]
        require((power == list(range(12))) == (n == 5), "generator order five")
    vp = coordinates()
    vm = [[x.sigma() for x in v] for v in vp]
    family = {}
    for name, sign, v in (("plus", 1, vp), ("minus", -1, vm)):
        table = [O, F(0, Q(sign,5)), F(0, Q(-sign,5)), -O]
        gram = [[table[d] for d in row] for row in ds]
        require(mm(gram, gram) == [[4*x for x in row] for row in gram], "G squared")
        require(sum((gram[i][i] for i in range(12)), Z) == F(12), "trace")
        norm = dot(v[0], v[0])
        require(norm.sign() > 0, "positive norm")
        require([[dot(a,b)/norm for b in v] for a in v] == gram, "coordinate Gram")
        require(det(v[:3]) != Z, "coordinate rank three")
        for i in range(12):
            for k in range(3):
                require(sum((v[j][k] for j in range(12) if ds[i][j] == 1), Z)
                        == sign*S*v[i][k], "adjacency coordinate mode")
        basis = transpose(v[:3])
        rotation = mm(transpose([v[generator[i]] for i in range(3)]), inverse(basis))
        require(mm(transpose(rotation), rotation) == eye(3) and det(rotation) == O, "proper rotation")
        require(mm(v, transpose(rotation)) == [v[generator[i]] for i in range(12)], "rotation action")
        trace = sum((rotation[i][i] for i in range(3)), Z)
        require(trace == (PHI if sign == 1 else PHI.sigma()), "trace pair")
        chords = [2*(O-x) for x in table]
        order = [1,2,3] if sign == 1 else [2,1,3]
        require(all((chords[b]-chords[a]).sign() > 0 for a,b in zip([0]+order,order)), "chord order")
        degrees = [degree(v, fs, y) for y in [[1,2,4],[-3,5,2],[7,-2,3]]]
        require(all(row["degree"] == (1 if sign == 1 else 7) for row in degrees), "degree pair")
        family[name] = {"gram": gram, "vectors": v, "norm_squared": norm,
                        "rank": 3, "trace": 12, "psd": True,
                        "rotation": rotation, "rotation_trace": trace,
                        "adjacency_eigenvalue": sign*S, "laplacian_eigenvalue": F(5)-sign*S,
                        "chord_squared_by_distance": chords, "chord_distance_order": order,
                        "degree_controls": degrees}
    require([[x.sigma() for x in row] for row in family["plus"]["gram"]] == family["minus"]["gram"], "Galois exchange")
    # Positive similarity phi times an odd coordinate permutation. This
    # identifies point sets, not the labeled oriented face complexes.
    relabel = []
    for v in vm:
        w = [PHI*v[1], PHI*v[0], PHI*v[2]]
        relabel.append(vp.index(w))
    require(sorted(relabel) == list(range(12)), "point set similarity")
    require(all(ds[relabel[i]][relabel[j]] == {0:0,1:2,2:1,3:3}[ds[i][j]]
                for i in range(12) for j in range(12)), "similarity swaps chord classes")
    return encode({"schema": "oph.galois-port-frames.v1", "baseline": BASELINE,
                   "mathematical_inputs": [FACE_SOURCE+"#orientedFaces"],
                   "face_source_sha256": hashlib.sha256((ROOT/FACE_SOURCE).read_bytes()).hexdigest(),
                   "faces": fs, "distances": ds, "generator": generator,
                   "orientation": {"ambient_relative_to_xyz": -1, "convex_degree": 1},
                   "family": family,
                   "similarity": {"scale": PHI, "coordinate_permutation": [1,0,2], "port_relabel": relabel},
                   "support_attachment": {"interface": "PORT-GRAM-SUPPORT-ATTACHMENT",
                                          "status": "NOT_PRESENT", "source_theorem": None},
                   "scope": "equivalent within the two Galois-conjugate rank-three realizations",
                   "verdict": "GALOIS_FRAME_AMBIGUITY_PERSISTS_WITHOUT_SUPPORT_ATTACHMENT"})


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
