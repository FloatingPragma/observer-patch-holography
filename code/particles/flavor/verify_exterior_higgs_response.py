#!/usr/bin/env python3
"""Independent exact CAR verification of the exterior Higgs response.

No producer import, measured quark value, RSCC expression, or mass-scheme
fixture enters this verifier.  Analytic closure and normalization arguments
are pinned separately from the executable finite-dimensional checks.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

ROOT=Path(__file__).resolve().parents[3]
SPEC="code/particles/flavor/exterior_higgs_response_spec.json"
SPEC_HASH="260c99f7c768c6507f10e993152def3af81a6940342dc4ce047c4133e2b82dac"
RECEIPT="code/particles/runs/flavor/exterior_higgs_response.json"
INPUTS=("code/a5_closure/manifests/super_tannakian_matter_reference.json",
        "code/a5_closure/manifests/family_band_attachment_reference.json")
PIN_PATHS=INPUTS+(
    "code/a5_closure/super_tannakian_matter_lift_certificate.py",
    "code/source_selection_model/response.py",
    "Lean/Screen/WeylYukawaConventions.lean",
    "Lean/Screen/ElectroweakBreakingComposition.lean",
    SPEC,
    "code/particles/flavor/exterior_higgs_response_certificate.py",
    "code/particles/flavor/verify_exterior_higgs_response.py",
)


class VerificationError(ValueError):
    pass


def require(condition: bool, code: str) -> None:
    if not condition:
        raise VerificationError(code)


def canonical(value) -> str:
    return json.dumps(value,sort_keys=True,separators=(",",":"))


@lru_cache(maxsize=1)
def car_coefficients() -> tuple[list,list]:
    creators=[]
    for mode in range(5):
        M=s.MutableSparseMatrix(32,32,{})
        for mask in range(32):
            if not mask & (1 << mode):
                preceding=(mask & ((1 << mode)-1)).bit_count()
                M[mask | (1 << mode),mask]=(-1)**preceding
        creators.append(M)
    for i in range(5):
        for j in range(5):
            require(creators[i]*creators[j]+creators[j]*creators[i]==s.zeros(32),"CAR_CREATION")
            require(creators[i].T*creators[j]+creators[j]*creators[i].T
                    ==(s.eye(32) if i==j else s.zeros(32)),"CAR_MIXED")
    vacuum=s.zeros(32,1)
    vacuum[0]=1
    up,down=[],[]
    for higgs in range(2):
        u,d=s.zeros(6,3),s.zeros(6,3)
        for color in range(3):
            for weak in range(2):
                q=creators[color]*creators[3+weak]*vacuum
                contracted=creators[3+higgs].T*q
                for right in range(3):
                    others=[k for k in range(3) if k!=right]
                    uc=(-1)**right*creators[others[0]]*creators[others[1]]*vacuum
                    dc=(-1)**right*creators[others[0]]*creators[others[1]]*creators[3]*creators[4]*vacuum
                    u[2*color+weak,right]=(creators[color]*creators[3+weak]*creators[3+higgs]*uc)[31]
                    d[2*color+weak,right]=sum(contracted[1<<mode]*(creators[mode]*dc)[31]
                                                 for mode in range(5))
        up.append([[int(x) for x in row] for row in u.tolist()])
        down.append([[int(x) for x in row] for row in d.tolist()])
    return up,down


@lru_cache(maxsize=1)
def independent_checks() -> dict:
    up,down=car_coefficients()
    # Check all component Gram coefficients.  These identities imply the
    # general complex-doublet formulas, not merely selected Higgs directions.
    U,D=list(map(s.Matrix,up)),list(map(s.Matrix,down))
    for i in range(2):
        for j in range(2):
            expected=s.eye(3) if i==j else s.zeros(3)
            require(U[i].T*U[j]==expected,"UP_COMPONENT_GRAM")
            require(D[i].T*D[j]==expected,"DOWN_COMPONENT_GRAM")
    # U(H)^dagger D(H) has monomials conjugate(h_i) conjugate(h_j).
    require(U[0].T*D[0]==s.zeros(3) and U[1].T*D[1]==s.zeros(3),"MIXED_DIAGONAL")
    require(U[0].T*D[1]+U[1].T*D[0]==s.zeros(3),"MIXED_CROSS")
    for T in [s.Matrix([[0,s.I],[s.I,0]]),s.Matrix([[0,1],[-1,0]]),s.diag(s.I,-s.I)]:
        Qaction=s.kronecker_product(s.eye(3),T)
        for j in range(2):
            require(sum((U[i]*T[i,j] for i in range(2)),s.zeros(6,3))+Qaction.T*U[j]==s.zeros(6,3),"UP_WEAK_WARD")
            require(sum((D[i]*(-T.T)[i,j] for i in range(2)),s.zeros(6,3))+Qaction.T*D[j]==s.zeros(6,3),"DOWN_WEAK_WARD")
    # Matrix units generate the complexified traceless color algebra.
    basis=[]
    for i in range(3):
        for j in range(3):
            if i!=j:
                M=s.zeros(3); M[i,j]=1; basis.append(M)
    basis.extend([s.diag(1,-1,0),s.diag(0,1,-1)])
    linear=s.Matrix.vstack(*[s.kronecker_product(s.eye(3),M)-s.kronecker_product(M.T,s.eye(3)) for M in basis])
    require(linear.rank()==8,"COLOR_WARD_RANK")
    R=s.Matrix([[1,1,0],[1,5,2],[0,2,10]])
    require([R[:i,:i].det() for i in range(1,4)]==[1,4,36],"KINETIC_POSITIVITY")
    Rc=s.Matrix([[1,-s.I,0],[s.I,5,2-2*s.I],[0,2+2*s.I,11]])
    require(Rc.H==Rc and [Rc[:i,:i].det() for i in range(1,4)]==[1,4,36],"COMPLEX_KINETIC_POSITIVITY")
    Sc=Rc.cholesky().H.inv()
    require((Sc.H*Rc*Sc).applyfunc(s.simplify)==s.eye(3),"DIRAC_CONGRUENCE")
    require((Sc.T*Rc*Sc).applyfunc(s.simplify)!=s.eye(3),"WEYL_TRANSPOSE_CONTROL")
    p,z=s.symbols("p z")
    require(s.expand((p*R-z*R).det()-36*(p-z)**3)==0,"GENERALIZED_PENCIL")
    # Independent elimination using an augmented block solve. The inverse
    # derivative term is checked by differentiating that solution directly.
    h=s.symbols("h")
    A=s.Matrix([[2+h,1],[1,3-h]])
    B=s.Matrix([[1,h],[0,1]])
    E=s.Matrix([[4+h,1],[1,2+2*h]])
    C=B.T
    X=E.inv()*C
    Xprime=E.inv()*(C.diff(h)-E.diff(h)*X)
    require((X.diff(h)-Xprime).applyfunc(s.simplify)==s.zeros(2),"HEAVY_SOLUTION_DERIVATIVE")
    effective=A-B*X
    require((effective.diff(h)-(A.diff(h)-B.diff(h)*X-B*Xprime)).applyfunc(s.simplify)==s.zeros(2),"EFFECTIVE_DERIVATIVE")
    require(E*E.diff(h)!=E.diff(h)*E,"NONCOMMUTING_CONTROL")
    # Direct block-diagonal elimination leaves exactly three repeated copies.
    lifted=s.diag(E,E,E).inv()*s.diag(C,C,C)
    require(lifted==s.diag(X,X,X),"FAMILY_HEAVY_ELIMINATION")
    return {
        "full_complex_doublet_gram":True,
        "mixed_up_down_gram_zero":True,
        "arbitrary_complex_channel_coefficients_checked":True,
        "full_weak_infinitesimal_covariance_checked":True,
        "up_down_hypercharge_balance_checked":True,
        "color_commutant_constraint_rank":8,
        "color_commutant_complex_dimension":1,
        "color_commutant_basis":"I3",
        "kinetic_congruence_matrix":[[str(x) for x in row] for row in R.tolist()],
        "shared_metric_generalized_root_multiplicity":3,
        "complex_Hermitian_Dirac_normalization_checked":True,
        "incorrect_Weyl_transpose_substitution_rejected":True,
        "noncommuting_schur_derivative_checked":True,
        "family_tensor_schur_complement_checked":True,
    }


def verify(receipt: dict,root: Path=ROOT) -> dict:
    expected_keys={"schema","status","source_pins","numeric_input_paths","basis","hypotheses",
                   "equations","proof_grade","up_higgs_component_matrices",
                   "down_conjugate_higgs_component_matrices","exact_checks","family_response","claims"}
    require(set(receipt)==expected_keys,"RECEIPT_KEYS")
    require(receipt["schema"]=="oph.exterior-higgs-response.v1","SCHEMA")
    require(receipt["status"]=="exact_conditional_tensor_and_normalization_certificate","STATUS")
    require(receipt["numeric_input_paths"]==list(INPUTS),"NUMERIC_INPUT_PATHS")
    require(set(receipt["source_pins"])==set(PIN_PATHS),"PIN_PATHS")
    for path in PIN_PATHS:
        require(hashlib.sha256((root/path).read_bytes()).hexdigest()==receipt["source_pins"][path],"SOURCE_PIN:"+path)
    require(receipt["source_pins"][SPEC]==SPEC_HASH,"SPECIFICATION_PIN")
    spec=json.loads((root/SPEC).read_text())
    for key in ("basis","hypotheses","equations","claims"):
        require(canonical(receipt[key])==canonical(spec[key]),"SPEC_PROJECTION:"+key)
    require(receipt["proof_grade"]==spec["analytic_proof"]["proof_grade"],"PROOF_GRADE")
    matter=json.loads((root/INPUTS[0]).read_text())["exterior_matter_contract"]
    charges=matter["block_trace_charges"]
    c,w=Fraction(charges["color_block"]),Fraction(charges["weak_block"])
    require((c,w)==(Fraction(-1,3),Fraction(1,2)) and 3*c+2*w==0,"SOURCE_CHARGES")
    require((c+w)+w+2*c==0 and (c+w)-w+(2*c+2*w)==0,"YUKAWA_CHARGES")
    require(matter["one_scalar"]=="weak_block","SOURCE_SCALAR")
    require(["Q","S","u_c"] in matter["yukawa_channels"] and ["Q","Sbar","d_c"] in matter["yukawa_channels"],"SOURCE_CHANNELS")
    family=json.loads((root/INPUTS[1]).read_text())
    require(family["attachment"]["family_dimension"]==3 and family["selection"]["minimizer"]=="3","SOURCE_FAMILY")
    up,down=car_coefficients()
    require(canonical(receipt["up_higgs_component_matrices"])==canonical(up),"UP_TENSOR")
    require(canonical(receipt["down_conjugate_higgs_component_matrices"])==canonical(down),"DOWN_TENSOR")
    require(canonical(receipt["exact_checks"])==canonical(independent_checks()),"EXACT_CHECKS")
    expected_family={"dimension":3,"restricted_operator_algebra":"I_family tensor End(matter)",
        "normalization_convention":"four-component Dirac Yukawa",
        "up":"g_u I3","down":"g_d I3",
        "canonical_up":"g_u I3 / sqrt(z_Q z_u Z_H)",
        "canonical_down":"g_d I3 / sqrt(z_Q z_d Z_H)",
        "within_sector_singular_value_multiplicities":[3],
        "CKM_identifiable_on_triply_degenerate_branch":False}
    require(canonical(receipt["family_response"])==canonical(expected_family),"FAMILY_RESPONSE")
    return {"verified":True,"scope":"full exterior tensors and conditional normalization; no source-selected masses"}


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("receipt",type=Path,nargs="?",default=ROOT/RECEIPT)
    args=ap.parse_args()
    print(json.dumps(verify(json.loads(args.receipt.read_text())),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
