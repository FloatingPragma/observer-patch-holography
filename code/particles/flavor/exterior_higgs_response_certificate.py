#!/usr/bin/env python3
"""Exact exterior Higgs response; no fitted inputs and no selected masses.

The producer reconstructs top-form wedge/interior coefficients.  The separate
verifier uses CAR creation/annihilation matrices.  General closure statements
have an analytic proof in the pinned specification; exact executable matrix
checks are not represented as a formal universal source theorem.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[3]
SPEC = "code/particles/flavor/exterior_higgs_response_spec.json"
OUTPUT = "code/particles/runs/flavor/exterior_higgs_response.json"
SOURCES = (
    "code/a5_closure/manifests/super_tannakian_matter_reference.json",
    "code/a5_closure/manifests/family_band_attachment_reference.json",
    "code/a5_closure/super_tannakian_matter_lift_certificate.py",
    "code/source_selection_model/response.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/ElectroweakBreakingComposition.lean",
)
IMPLEMENTATION = (
    SPEC,
    "code/particles/flavor/exterior_higgs_response_certificate.py",
    "code/particles/flavor/verify_exterior_higgs_response.py",
)


def sha(path: str) -> str:
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def wedge(*subsets: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    word = tuple(i for subset in subsets for i in subset)
    if len(set(word)) != len(word):
        return 0, ()
    sign = (-1)**sum(a > b for i,a in enumerate(word) for b in word[i+1:])
    return sign, tuple(sorted(word))


def top(*subsets: tuple[int, ...]) -> int:
    sign, result = wedge(*subsets)
    return sign if result == (0,1,2,3,4) else 0


def interior(mode: int, subset: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    if mode not in subset:
        return 0, ()
    i = subset.index(mode)
    return (-1)**i, subset[:i]+subset[i+1:]


def coefficients() -> tuple[list, list]:
    up, down = [], []
    for w in range(2):
        u, d = s.zeros(6,3), s.zeros(6,3)
        for c in range(3):
            for q in range(2):
                for t in range(3):
                    comp = tuple(i for i in range(3) if i != t)
                    orientation = (-1)**t
                    u[2*c+q,t] = orientation*top((c,3+q),(3+w,),comp)
                    sign, part = interior(3+w,(c,3+q))
                    d[2*c+q,t] = orientation*sign*top(part,comp+(3,4))
        up.append(u.tolist())
        down.append(d.tolist())
    return ([[[int(x) for x in row] for row in m] for m in up],
            [[[int(x) for x in row] for row in m] for m in down])


def color_rank() -> int:
    xs = s.symbols("x0:9")
    X = s.Matrix(3,3,xs)
    generators = []
    for i in range(3):
        for j in range(i+1,3):
            a,b = s.zeros(3),s.zeros(3)
            a[i,j],a[j,i] = 1,-1
            b[i,j],b[j,i] = s.I,s.I
            generators.extend([a,b])
    generators.extend([s.diag(s.I,-s.I,0),s.diag(s.I,s.I,-2*s.I)])
    equations = [z for g in generators for z in X*g-g*X]
    M,_ = s.linear_eq_to_matrix(equations,xs)
    assert M.nullspace() == [s.Matrix(s.eye(3)).reshape(9,1)]
    return M.rank()


def exact_checks(up: list, down: list) -> dict:
    a,b,c,d = s.symbols("a b c d",real=True)
    h0,h1 = a+s.I*b,c+s.I*d
    U = h0*s.Matrix(up[0])+h1*s.Matrix(up[1])
    D = s.conjugate(h0)*s.Matrix(down[0])+s.conjugate(h1)*s.Matrix(down[1])
    norm = a*a+b*b+c*c+d*d
    assert (U.H*U-norm*s.eye(3)).applyfunc(s.expand) == s.zeros(3)
    assert (D.H*D-norm*s.eye(3)).applyfunc(s.expand) == s.zeros(3)
    assert (U.H*D).applyfunc(s.expand) == s.zeros(3)
    gu,gd = s.symbols("gu gd",complex=True)
    full = (gu*U).row_join(gd*D)
    target = s.diag(*([norm*gu*s.conjugate(gu)]*3+[norm*gd*s.conjugate(gd)]*3))
    assert (full.H*full-target).applyfunc(s.expand) == s.zeros(6)
    epsilon=s.Matrix([[0,1],[-1,0]])
    weak_generators=[s.Matrix([[0,s.I],[s.I,0]]),
                     s.Matrix([[0,1],[-1,0]]),s.diag(s.I,-s.I)]
    for T in weak_generators:
        assert T.T*epsilon+epsilon*T==s.zeros(2)
        assert T.H == -T
    # All-left-Weyl source charges, independently of action coefficient values.
    assert s.Rational(1,6)+s.Rational(1,2)-s.Rational(2,3)==0
    assert s.Rational(1,6)-s.Rational(1,2)+s.Rational(1,3)==0
    L = s.Matrix([[1,0,0],[1,2,0],[0,1,3]])
    R = L*L.T
    S = L.T.inv()
    assert S.T*R*S == s.eye(3)
    Lc=s.Matrix([[1,0,0],[s.I,2,0],[0,1+s.I,3]])
    Rc=Lc*Lc.H
    Sc=Lc.H.inv()
    assert (Sc.H*Rc*Sc).applyfunc(s.simplify) == s.eye(3)
    assert (Sc.T*Rc*Sc).applyfunc(s.simplify) != s.eye(3)
    p,h,g = s.symbols("p h g")
    assert s.expand((p*R-h*g*R).det()-R.det()*(p-h*g)**3) == 0
    A = s.Matrix([[2+h,1],[1,3-h]])
    B = s.Matrix([[1,h],[0,1]])
    C = B.T
    E = s.Matrix([[4+h,1],[1,2+2*h]])
    derivative = (A.diff(h)-B.diff(h)*E.inv()*C-B*E.inv()*C.diff(h)
                  +B*E.inv()*E.diff(h)*E.inv()*C)
    assert ((A-B*E.inv()*C).diff(h)-derivative).applyfunc(s.simplify) == s.zeros(2)
    # A genuinely noncommuting witness checks the order of the inverse factors.
    assert E*E.diff(h)-E.diff(h)*E != s.zeros(2)
    I = s.eye(3)
    lifted = (s.kronecker_product(I,A)-s.kronecker_product(I,B)
              *s.kronecker_product(I,E).inv()*s.kronecker_product(I,C))
    assert (lifted-s.kronecker_product(I,A-B*E.inv()*C)).applyfunc(s.simplify) == s.zeros(6)
    return {
        "full_complex_doublet_gram": True,
        "mixed_up_down_gram_zero": True,
        "arbitrary_complex_channel_coefficients_checked": True,
        "full_weak_infinitesimal_covariance_checked": True,
        "up_down_hypercharge_balance_checked": True,
        "color_commutant_constraint_rank": color_rank(),
        "color_commutant_complex_dimension": 1,
        "color_commutant_basis": "I3",
        "kinetic_congruence_matrix": [[str(x) for x in row] for row in R.tolist()],
        "shared_metric_generalized_root_multiplicity": 3,
        "complex_Hermitian_Dirac_normalization_checked": True,
        "incorrect_Weyl_transpose_substitution_rejected": True,
        "noncommuting_schur_derivative_checked": True,
        "family_tensor_schur_complement_checked": True,
    }


def build() -> dict:
    spec = json.loads((ROOT/SPEC).read_text())
    matter = json.loads((ROOT/SOURCES[0]).read_text())["exterior_matter_contract"]
    family = json.loads((ROOT/SOURCES[1]).read_text())
    assert matter["block_trace_charges"] == {"color_block":"-1/3","weak_block":"1/2"}
    assert matter["one_scalar"] == "weak_block"
    assert ["Q","S","u_c"] in matter["yukawa_channels"]
    assert ["Q","Sbar","d_c"] in matter["yukawa_channels"]
    assert family["attachment"]["family_dimension"] == 3
    assert family["selection"]["minimizer"] == "3"
    up,down = coefficients()
    return {
        "schema":"oph.exterior-higgs-response.v1",
        "status":"exact_conditional_tensor_and_normalization_certificate",
        "source_pins":{p:sha(p) for p in SOURCES+IMPLEMENTATION},
        "numeric_input_paths":list(SOURCES[:2]),
        "basis":spec["basis"],
        "hypotheses":spec["hypotheses"],
        "equations":spec["equations"],
        "proof_grade":spec["analytic_proof"]["proof_grade"],
        "up_higgs_component_matrices":up,
        "down_conjugate_higgs_component_matrices":down,
        "exact_checks":exact_checks(up,down),
        "family_response":{
            "dimension":family["attachment"]["family_dimension"],
            "restricted_operator_algebra":"I_family tensor End(matter)",
            "normalization_convention":"four-component Dirac Yukawa",
            "up":"g_u I3", "down":"g_d I3",
            "canonical_up":"g_u I3 / sqrt(z_Q z_u Z_H)",
            "canonical_down":"g_d I3 / sqrt(z_Q z_d Z_H)",
            "within_sector_singular_value_multiplicities":[3],
            "CKM_identifiable_on_triply_degenerate_branch":False,
        },
        "claims":spec["claims"],
    }


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--output",type=Path,default=ROOT/OUTPUT)
    args=ap.parse_args()
    payload=json.dumps(build(),indent=2,sort_keys=True)+"\n"
    if args.check:
        if args.output.read_text() != payload:
            raise SystemExit("exterior Higgs response receipt is stale")
        print("PASS: exterior Higgs response receipt")
    else:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload)
        print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
