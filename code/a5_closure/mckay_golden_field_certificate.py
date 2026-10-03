#!/usr/bin/env python3
"""Exact quaternionic construction and McKay-data producer for 2I.

This module uses only Python's standard library and exact Fraction arithmetic.
It contains no OPH imports or copied OPH representation data.  It derives the
120-unit-quaternion carrier from the standard H4 root coordinates, embeds each
unit quaternion into SL(2, Q(sqrt(5), i)), and recovers irreducible characters
from symmetric powers and their exact Galois conjugates.  The affine-E8 target
is intentionally absent here; graph identification belongs to the independent
verifier.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "code/a5_closure/receipts/arithmon_mckay_golden_field_certificate.json"

@dataclass(frozen=True, order=True)
class Q5:
    """a + b sqrt(5), with a,b rational."""
    a: F = F(0)
    b: F = F(0)
    def __add__(x, y):
        y = q5(y); return Q5(x.a + y.a, x.b + y.b)
    __radd__ = __add__
    def __neg__(x): return Q5(-x.a, -x.b)
    def __sub__(x, y): return x + (-q5(y))
    def __rsub__(x, y): return q5(y) + (-x)
    def __mul__(x, y):
        y = q5(y); return Q5(x.a*y.a + 5*x.b*y.b, x.a*y.b + x.b*y.a)
    __rmul__ = __mul__
    def __truediv__(x, y):
        y = q5(y); den = y.a*y.a - 5*y.b*y.b
        if den == 0: raise ZeroDivisionError
        return Q5((x.a*y.a - 5*x.b*y.b)/den, (x.b*y.a - x.a*y.b)/den)
    def __pow__(x, n: int):
        if n < 0: return (Q5(1)/x)**(-n)
        out, base, k = Q5(1), x, n
        while k:
            if k & 1: out = out*base
            base = base*base; k >>= 1
        return out
    def is_zero(x): return x.a == 0 and x.b == 0
    def text(x):
        if x.b == 0: return frac_text(x.a)
        if x.a == 0: return f"{frac_text(x.b)}*sqrt(5)"
        sign = "+" if x.b > 0 else "-"
        return f"{frac_text(x.a)}{sign}{frac_text(abs(x.b))}*sqrt(5)"

def q5(x) -> Q5:
    if isinstance(x, Q5): return x
    return Q5(F(x), F(0))

def frac_text(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

def quadratic_minimal_polynomial_text(value: Q5) -> str:
    """Compute the monic quadratic over Q and verify it vanishes exactly."""
    if value.b == 0:
        raise ValueError("a rational trace does not generate Q(sqrt(5))")
    linear = -2 * value.a
    constant = value.a * value.a - 5 * value.b * value.b
    if not (value * value + Q5(linear) * value + Q5(constant)).is_zero():
        raise ValueError("computed quadratic does not annihilate the exact trace")
    def term(coefficient: F, variable: str) -> str:
        if coefficient == 0:
            return ""
        sign = "+" if coefficient > 0 else "-"
        magnitude = abs(coefficient)
        body = variable if variable and magnitude == 1 else (
            f"{frac_text(magnitude)}*{variable}" if variable else frac_text(magnitude))
        return sign + body
    return "X^2" + term(linear, "X") + term(constant, "")

def sigma_q5(x: Q5) -> Q5: return Q5(x.a, -x.b)

@dataclass(frozen=True)
class C5:
    """x + i y with x,y in Q(sqrt(5)); sigma fixes i."""
    re: Q5 = Q5()
    im: Q5 = Q5()
    def __add__(x, y):
        y = c5(y); return C5(x.re+y.re, x.im+y.im)
    __radd__ = __add__
    def __neg__(x): return C5(-x.re, -x.im)
    def __sub__(x, y): return x + (-c5(y))
    def __mul__(x, y):
        y = c5(y); return C5(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re)
    __rmul__ = __mul__
    def conj(x): return C5(x.re, -x.im)
    def is_zero(x): return x.re.is_zero() and x.im.is_zero()
    def text(x):
        if x.im.is_zero(): return x.re.text()
        return f"({x.re.text()})+i*({x.im.text()})"

def c5(x) -> C5:
    if isinstance(x, C5): return x
    return C5(q5(x), Q5())

def sigma_c5(x: C5) -> C5: return C5(sigma_q5(x.re), sigma_q5(x.im))

ZERO, ONE = Q5(), Q5(1)
HALF = Q5(F(1, 2))
PHI = Q5(F(1, 2), F(1, 2))
PHI_INV = Q5(F(-1, 2), F(1, 2))
I = C5(Q5(), ONE)
Q = tuple[Q5, Q5, Q5, Q5]
IDENTITY: Q = (ONE, ZERO, ZERO, ZERO)
NEG_IDENTITY: Q = (-ONE, ZERO, ZERO, ZERO)


def permutation_sign(p: Sequence[int]) -> int:
    inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))
    return -1 if inversions & 1 else 1


def qkey(x: Q):
    return tuple((v.a.numerator, v.a.denominator, v.b.numerator, v.b.denominator) for v in x)


def enumerate_2i() -> set[Q]:
    """Enumerate the standard H4/600-cell unit-quaternion coordinates.

    The three disjoint coordinate families are the 8 signed coordinate units,
    the 16 half-coordinate units, and the 96 even permutations of
    (0, +/-1/2, +/-phi/2, +/-phi^{-1}/2).  Signs on nonzero coordinates vary
    independently.  This is a direct coordinate enumeration, not a table or
    multiplication law copied from another certificate.
    """
    elements: set[Q] = set()
    for slot in range(4):
        for sign in (-1, 1):
            q = [ZERO]*4; q[slot] = Q5(sign); elements.add(tuple(q))
    for signs in product((-1, 1), repeat=4):
        elements.add(tuple(Q5(F(sign, 2)) for sign in signs))
    base = (ZERO, HALF, PHI*HALF, PHI_INV*HALF)
    for p in permutations(range(4)):
        if permutation_sign(p) < 0: continue
        values = [base[index] for index in p]
        slots = [i for i, value in enumerate(values) if not value.is_zero()]
        for signs in product((-1, 1), repeat=3):
            signed = list(values)
            for slot, sign in zip(slots, signs):
                if sign < 0: signed[slot] = -signed[slot]
            elements.add(tuple(signed))
    return elements


def qmul(x: Q, y: Q) -> Q:
    a,b,c,d=x; e,f,g,h=y
    return (a*e-b*f-c*g-d*h, a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f, a*h+b*g-c*f+d*e)


def qconj(x: Q) -> Q:
    a,b,c,d=x; return (a,-b,-c,-d)


def qnorm(x: Q) -> Q:
    return sum((v*v for v in x), ZERO)


def exact_group_checks(elements: set[Q]) -> dict:
    require(len(elements) == 120, "2I_COORDINATE_COUNT", f"coordinate families produced {len(elements)} unique elements")
    require(IDENTITY in elements, "2I_IDENTITY", "identity missing")
    require(all(qnorm(x) == ONE for x in elements), "2I_UNIT_NORM", "a coordinate has nonunit quaternion norm")
    require(all(qmul(x,y) in elements for x in elements for y in elements), "2I_CLOSURE", "multiplication closure failed")
    require(all(qconj(x) in elements for x in elements), "2I_INVERSE", "inverse closure failed")
    center = {x for x in elements if all(qmul(x,y) == qmul(y,x) for y in elements)}
    require(center == {IDENTITY, NEG_IDENTITY}, "2I_CENTER", "center is not exactly +/-1")
    # Associativity is inherited from Hamilton's associative quaternion algebra;
    # this additionally checks the implemented multiplication on all basis triples.
    basis = [(ONE,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
             (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ONE)]
    require(all(qmul(qmul(a,b),c) == qmul(a,qmul(b,c)) for a in basis for b in basis for c in basis),
            "QUATERNION_ASSOCIATIVITY_BASIS", "quaternion multiplication implementation failed basis associativity")
    return {"order": len(elements), "center_size": len(center),
            "center": ["+1", "-1"], "closure": True, "inverse_closure": True,
            "associativity_source": "Hamilton quaternion algebra; basis-triple implementation check passed",
            "construction": "exact standard H4/600-cell unit-quaternion coordinate set",
            "canonical_element_set_sha256": hashlib.sha256(json.dumps(sorted([qkey(x) for x in elements]),separators=(",",":")).encode()).hexdigest()}


def matrix_of_quaternion(q: Q) -> tuple[tuple[C5,C5], tuple[C5,C5]]:
    a,b,c,d=q
    return ((c5(a)+c5(b)*I, c5(c)+c5(d)*I),
            (-c5(c)+c5(d)*I, c5(a)-c5(b)*I))


def matrix_mul(a, b):
    return tuple(tuple(a[i][0]*b[0][j]+a[i][1]*b[1][j] for j in range(2)) for i in range(2))


def matrix_det(m): return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def matrix_sigma(m): return tuple(tuple(sigma_c5(x) for x in row) for row in m)


def chars_from_doublet(elements: set[Q], sigma: bool = False) -> dict[Q,Q5]:
    out = {}
    for g in elements:
        m = matrix_of_quaternion(g)
        if sigma: m = matrix_sigma(m)
        tr = m[0][0]+m[1][1]
        require(tr.im.is_zero(), "DOUBLET_REAL_TRACE", "spin character acquired imaginary trace")
        out[g] = tr.re
    return out


def check_faithful_kernel(representation: dict[Q,object], identity_matrix: object) -> set[Q]:
    kernel={g for g,m in representation.items() if m==identity_matrix}
    require(kernel=={IDENTITY}, "DOUBLET_FAITHFUL", "matrix representation has nontrivial kernel")
    return kernel


def exact_representation_checks(elements: set[Q], classes: list[set[Q]]) -> dict:
    rho = {g: matrix_of_quaternion(g) for g in elements}
    imat = ((c5(ONE),c5(ZERO)),(c5(ZERO),c5(ONE)))
    require(rho[IDENTITY] == imat, "DOUBLEt_IDENTITY", "rho(1) != I")
    require(all(matrix_mul(rho[x],rho[y]) == rho[qmul(x,y)] for x in elements for y in elements),
            "DOUBLEt_HOMOMORPHISM", "rho(xy)=rho(x)rho(y) failed")
    require(all(matrix_det(m) == c5(ONE) for m in rho.values()), "DOUBLEt_DETERMINANT", "determinant is not one")
    kernel = {g for g,m in rho.items() if m == imat}
    require(kernel == {IDENTITY}, "DOUBLEt_FAITHFUL", "matrix representation has nontrivial kernel")
    minus = tuple(tuple(-c5(ONE) if i==j else c5(ZERO) for j in range(2)) for i in range(2))
    require(rho[NEG_IDENTITY] == minus, "DOUBLEt_CENTER_ACTION", "-1 does not act by -I2")
    chi = chars_from_doublet(elements)
    require(chi[IDENTITY] == Q5(2), "DOUBLEt_DIMENSION", "identity trace is not two")
    require(all(all(chi[g] == chi[next(iter(cls))] for g in cls) for cls in classes),
            "DOUBLEt_CLASS_FUNCTION", "trace is not constant on conjugacy classes")
    require(int_value(character_inner(chi,chi,elements), "DOUBLET_IRREDUCIBLE_NORM") == 1, "DOUBLET_IRREDUCIBLE", "character norm is not one")
    # Matrix coefficients generate Q(sqrt(5),i): i occurs at q=i; sqrt(5)
    # occurs in the standard golden-coordinate entries, while all entries lie
    # in the displayed biquadratic coefficient model by construction.
    sqrt5_witness = next(g for g in elements if any(v.b != 0 for row in rho[g] for z in row for v in (z.re,z.im)))
    i_witness = (ZERO,ONE,ZERO,ZERO)
    require(rho[i_witness][0][0] == I, "COEFFICIENT_I_WITNESS", "i absent from coefficient field")
    return {"dimension": 2, "faithful": True, "irreducible": True, "determinant_one": True,
            "central_minus_one_action": "-I_2", "coefficient_field": "Q(sqrt(5), i)",
            "sqrt5_coefficient_witness_quaternion": [v.text() for v in sqrt5_witness],
            "i_coefficient_witness_quaternion": [v.text() for v in i_witness],
            "character_norm": 1, "trace_field_candidate": "Q(sqrt(5))"}


def exact_galois_representation_checks(elements: set[Q]) -> dict:
    rho={g:matrix_sigma(matrix_of_quaternion(g)) for g in elements}
    imat=((c5(ONE),c5(ZERO)),(c5(ZERO),c5(ONE)))
    require(rho[IDENTITY]==imat, "GALOIS_REP_IDENTITY", "conjugate representation misses identity")
    require(all(matrix_mul(rho[x],rho[y])==rho[qmul(x,y)] for x in elements for y in elements),
            "GALOIS_REP_HOMOMORPHISM", "entrywise Galois conjugate is not a representation")
    require(all(matrix_det(m)==c5(ONE) for m in rho.values()), "GALOIS_REP_DETERMINANT", "conjugate determinant is not one")
    check_faithful_kernel(rho,imat)
    chi={g:(m[0][0]+m[1][1]).re for g,m in rho.items()}
    require(all((m[0][0]+m[1][1]).im.is_zero() for m in rho.values()), "GALOIS_TRACE_REAL", "conjugate trace is not real")
    require(int_value(character_inner(chi,chi,elements), "GALOIS_REP_IRREDUCIBLE") == 1,
            "GALOIS_REP_IRREDUCIBLE", "conjugate representation is not irreducible")
    return {"representation_homomorphism":True,"dimension":2,"determinant_one":True,
            "faithful":True,"irreducible":True,"sigma_fixes_i":True}


def character_inner(left: dict[Q,Q5], right: dict[Q,Q5], elements: set[Q]) -> Q5:
    total = ZERO
    for g in elements: total = total + left[g]*right[g]
    return total / Q5(len(elements))


def int_value(value: Q5, gate: str) -> int:
    require(value.b == 0 and value.a.denominator == 1, gate, f"expected integer, got {value.text()}")
    return value.a.numerator


def sigma_character(character: dict[Q,Q5]) -> dict[Q,Q5]:
    # For a fixed abstract group element, apply the field automorphism to the value.
    return {g:sigma_q5(v) for g,v in character.items()}


def char_equal(a, b, elements): return all(a[g] == b[g] for g in elements)


def symmetric_power_characters(chi: dict[Q,Q5], elements: set[Q], order: int):
    one = {g:ONE for g in elements}
    out = [one, chi]
    # The SU(2) determinant-one recurrence follows from the characteristic
    # polynomial of each 2x2 matrix: S_(n+1)=chi*S_n-S_(n-1).
    while len(out) <= order:
        out.append({g:chi[g]*out[-1][g]-out[-2][g] for g in elements})
    return out


def derive_irreducibles(elements: set[Q], chi: dict[Q,Q5], classes: list[set[Q]]) -> tuple[list[dict], list[str], dict]:
    irreps: list[dict[Q,Q5]] = []
    labels: list[str] = []
    powers = symmetric_power_characters(chi,elements,len(elements))
    for degree, candidate in enumerate(powers):
        for label, char in ((f"Sym^{degree}(V)",candidate),
                            (f"sigma(Sym^{degree}(V))",sigma_character(candidate))):
            residual=dict(char)
            for known in irreps:
                multiplicity=int_value(character_inner(residual,known,elements), "IRREP_PROJECTION_INTEGER")
                require(multiplicity>=0, "IRREP_NONNEGATIVE_PROJECTION", f"negative known-constituent multiplicity in {label}")
                if multiplicity:
                    residual={g:residual[g]-Q5(multiplicity)*known[g] for g in elements}
            if all(value.is_zero() for value in residual.values()): continue
            norm=int_value(character_inner(residual,residual,elements), "IRREP_RESIDUAL_NORM")
            # Norm-one residuals are new irreducible characters. A larger
            # residual norm is left for a later symmetric-power candidate;
            # no expected dimension list or character table guides this step.
            if norm != 1: continue
            require(residual[IDENTITY].a>0 and residual[IDENTITY].b==0 and residual[IDENTITY].a.denominator==1,
                    "IRREP_POSITIVE_INTEGRAL_DIMENSION", f"{label} residual has invalid degree")
            require(all(int_value(character_inner(residual,known,elements), "IRREP_ORTHOGONALITY_INTEGER") == 0 for known in irreps),
                    "IRREP_ORTHOGONALITY", f"{label} residual has nonzero inner product with a prior irreducible")
            require(all(all(residual[g] == residual[next(iter(cls))] for g in cls) for cls in classes),
                    "IRREP_CLASS_FUNCTION", f"{label} residual is not class constant")
            irreps.append(residual); labels.append(label+" residual")
            dims=[int_value(char0[IDENTITY], "IRREP_DIMENSION_INTEGER") for char0 in irreps]
            if sum(d*d for d in dims)==len(elements):
                require(all(int_value(character_inner(a,b,elements), "IRREP_PAIRING_INTEGER")==(1 if i==j else 0)
                            for i,a in enumerate(irreps) for j,b in enumerate(irreps)),
                        "IRREP_COMPLETE_ORTHONORMAL", "completed system is not orthonormal")
                return irreps,labels,{"count":len(irreps),"dimensions":dims,
                    "sum_squared_dimensions":sum(d*d for d in dims),"orthonormal":True,
                    "derivation":"symmetric powers, exact constituent subtraction, inner products, and field conjugation",
                    "hardcoded_character_table":False,"hardcoded_dimension_list":False}
    raise ValueError("symmetric-power/Galois orbit did not reach the regular-character dimension sum")


def conjugacy_classes(elements: set[Q]) -> list[set[Q]]:
    remaining=set(elements); classes=[]
    while remaining:
        representative=min(remaining,key=qkey)
        cls={qmul(qmul(x,representative),qconj(x)) for x in elements}
        require(cls <= remaining, "CLASS_PARTITION", "conjugacy classes overlap")
        classes.append(cls); remaining-=cls
    require(sum(map(len,classes))==len(elements), "CLASS_COVER", "conjugacy classes do not cover group")
    return classes


def fusion_matrix(elements: set[Q], chars: list[dict[Q,Q5]], faithful: dict[Q,Q5]):
    matrix=[]
    dims=[int_value(char[IDENTITY], "FUSION_DIM") for char in chars]
    product_chars=[]
    for left in chars:
        tensor={g:faithful[g]*left[g] for g in elements}
        product_chars.append(tensor)
        row=[int_value(character_inner(tensor,right,elements), "FUSION_MULTIPLICITY") for right in chars]
        require(all(x>=0 for x in row), "FUSION_NONNEGATIVE", "negative tensor multiplicity")
        require(sum(m*d for m,d in zip(row,dims))==2*dims[len(matrix)], "FUSION_DIMENSION_CONSERVATION", "row dimension conservation failed")
        matrix.append(row)
    require(all(matrix[i][j]==matrix[j][i] for i in range(len(chars)) for j in range(len(chars))),
            "FUSION_SYMMETRY", "fusion matrix is not symmetric")
    return matrix,dims


def graph_data(matrix: list[list[int]], dims: list[int], gate: str):
    n=len(matrix)
    require(all(matrix[i][i]==0 for i in range(n)), gate+"_NO_LOOPS", "derived graph has a loop")
    require(all(matrix[i][j] in (0,1) for i in range(n) for j in range(n)), gate+"_SIMPLE", "derived graph is not simply-laced")
    edges=[[i,j] for i in range(n) for j in range(i+1,n) if matrix[i][j]]
    seen={0}; todo=[0]
    while todo:
        i=todo.pop()
        for j in range(n):
            if matrix[i][j] and j not in seen: seen.add(j);todo.append(j)
    connected=len(seen)==n
    tree=connected and len(edges)==n-1
    eig=all(sum(matrix[i][j]*dims[j] for j in range(n))==2*dims[i] for i in range(n))
    return {"vertices":n,"edges":edges,"edge_count":len(edges),"connected":connected,
            "tree":tree,"simply_laced":True,"dimension_vector":dims,
            "dimension_vector_eigenvalue":2 if eig else None,"dimension_vector_2_eigenvector":eig}


def f5_mul(x,y):
    a,b,c,d=x; e,f,g,h=y
    return ((a*e+b*g)%5,(a*f+b*h)%5,(c*e+d*g)%5,(c*f+d*h)%5)


def enumerate_sl2f5() -> dict:
    matrices={(a,b,c,d) for a,b,c,d in product(range(5),repeat=4) if (a*d-b*c)%5==1}
    ident=(1,0,0,1); minus=(4,0,0,4)
    require(ident in matrices and minus in matrices, "SL2F5_IDENTITIES", "identity or negative identity missing")
    require(all(f5_mul(x,y) in matrices for x in matrices for y in matrices), "SL2F5_CLOSURE", "matrix set is not closed")
    center={x for x in matrices if all(f5_mul(x,y)==f5_mul(y,x) for y in matrices)}
    require(center=={ident,minus}, "SL2F5_CENTER", "center is not exactly +/-I")
    require(all((a*d-b*c)%5==1 for a,b,c,d in matrices), "SL2F5_DETERMINANT", "determinant condition failed")
    require(all((d%5,(-b)%5,(-c)%5,a%5) in matrices for a,b,c,d in matrices),
            "SL2F5_INVERSE", "matrix inverse closure failed")
    return {"definition":"all 2x2 matrices over F5 with determinant one, enumerated from entries",
            "order":len(matrices),"center_size":len(center),"center":["+I","-I"],
            "closure":True,"not_identified_by_order_alone":True}


def main_certificate() -> dict:
    elements=enumerate_2i()
    group=exact_group_checks(elements)
    sl2f5=enumerate_sl2f5()
    classes=conjugacy_classes(elements)
    rep=exact_representation_checks(elements,classes)
    chi=chars_from_doublet(elements)
    irreps,labels,irrep_report=derive_irreducibles(elements,chi,classes)
    fusion,dims=fusion_matrix(elements,irreps,chi)
    graph=graph_data(fusion,dims,"V_FUSION")
    galois_rep=exact_galois_representation_checks(elements)
    chisigma=sigma_character(chi)
    require(not char_equal(chi,chisigma,elements), "GALOIS_DOUBLEt_DISTINCT", "Galois-conjugate trace character equals original")
    require(int_value(character_inner(chisigma,chisigma,elements), "GALOIS_DOUBLEt_NORM") == 1,
            "GALOIS_DOUBLEt_IRREDUCIBLE", "conjugate character is not irreducible")
    sigma_fusion,sigma_dims=fusion_matrix(elements,irreps,chisigma)
    sigma_graph=graph_data(sigma_fusion,sigma_dims,"SIGMA_FUSION")
    require(sigma_dims==dims, "GALOIS_DIMENSIONS", "Galois fusion system dimensions changed")
    traces=list(chi.values())
    require(all(isinstance(v,Q5) for v in traces), "CHARACTER_FIELD_CONTAINMENT", "trace escaped Q(sqrt(5))")
    witness=next((g for g in elements if chi[g].b != 0),None)
    require(witness is not None, "CHARACTER_FIELD_GENERATION", "all traces are rational")
    require((PHI*PHI-PHI-ONE).is_zero(), "PHI_POLYNOMIAL", "phi fails its exact polynomial")
    require(sigma_q5(PHI)==Q5(F(1,2),F(-1,2)), "PHI_GALOIS", "sqrt(5) conjugation did not map phi to psi")
    require(PHI in [chi[g] for g in elements] or -PHI in [chi[g] for g in elements] or
            PHI_INV in [chi[g] for g in elements] or -PHI_INV in [chi[g] for g in elements],
            "GOLDEN_TRACE_WITNESS", "nonrational trace does not exhibit the golden quadratic generator")
    class_sizes=sorted(len(c) for c in classes)
    require(sum(class_sizes)==len(elements), "CLASS_SIZE_SUM", "class sizes fail to sum to group order")
    return {
        "schema":"arithmon.mckay_golden_field.raw.v1",
        "source_group":group,
        "faithful_doublet":rep,
        "conjugacy_classes":{"count":len(classes),"sizes":class_sizes,"partition_verified":True},
        "irreducibles":{**irrep_report,"labels":labels},
        "fusion":{"matrix":fusion,**graph},
        "galois_conjugate_fusion":{"matrix":sigma_fusion,**sigma_graph,
            "automorphism":"sqrt(5)->-sqrt(5); i fixed","doublet_distinct":True,
            **galois_rep,"doublet_irreducible":True,"doublet_faithful":True,"fusion_recomputed":True,
            "same_labeled_graph":fusion==sigma_fusion},
        "golden_character_field":{"all_traces_in_Qsqrt5":True,
            "nonrational_trace_witness":{"quaternion":[v.text() for v in witness],"trace":chi[witness].text()},
            "field_exactly":"Q(sqrt(5))", "witness_minimal_polynomial":quadratic_minimal_polynomial_text(chi[witness])},
        "sl2f5_cross_identification":{**sl2f5,"explicit_isomorphism":False,
            "status":"separately constructed with matching order and center; explicit isomorphism not proved"},
        "golden_embedding":{"phi":PHI.text(),"psi":sigma_q5(PHI).text(),"sigma_phi_equals_psi":sigma_q5(PHI)==Q5(F(1,2),F(-1,2))},
        "primary_verdict":"EXACT_DATA_DERIVED_GRAPH_IDENTIFICATION_PENDING",
        "claim_boundary":["no OPH derivation input","no Koide input","no experimental data",
            "no physical particle identification","no real embedding selector"]
    }


def require(condition: bool, code: str, message: str):
    if not condition: raise ValueError(f"{code}: {message}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument("--no-write",action="store_true",help="print exact receipt without writing the default artifact")
    args=parser.parse_args()
    result=main_certificate()
    data=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if not args.no_write:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(data)
    print(data,end="")

if __name__=="__main__": main()
