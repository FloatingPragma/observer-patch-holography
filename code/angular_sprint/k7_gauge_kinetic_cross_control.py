#!/usr/bin/env python3
"""Fail-closed exact cross-control for the pinned K7 gauge claims and OPH plane.

This packet audits whether a K7 gauge object may be typed as OPH's sealed
alpha_inverse column. It reads and pins the vendored source bytes, including
measurement-comparison text in the paper, but consumes no public measurement
values as cross-inputs and never opens the sealed OPH comparison column.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
VENDOR = HERE / "external" / "k7_gauge"
OUT = HERE / "runtime" / "k7_gauge_kinetic_cross_control.json"
SCHEMA = "arithmon.oph.k7_gauge_kinetic_cross.v1"
OPH_SHA = "0660c94573c8f79c1b86955fbd1d7d927a0e4fa2"
K7_SHA = "210480b2ba8a13c98e6a6f146fbfd965e4904c00"
K7L_SHA = "2553ed170d04db8777f84fb2f9ef2b66b386f7fb"
TYPE1_PROVENANCE = "UNKNOWN; pinned honest ledger classifies target-exposed Type-I relations as historical retrodictions/calibration-sensitive unless a target-blind origin is documented"
BTEST_PROVENANCE = "UNKNOWN; pinned sources do not establish a target-blind origin for this historical conditional B-test relation"
OPH_PINS = {
    "Lean/Screen/KineticFormDichotomy.lean": "671578153bb70bf590e380fb7fd035062a341155",
    "Lean/Screen/RGRepresentationFrontier.lean": "c35eb3bc1f8e5f8bd4633bdea07911049bb985d0",
    "code/angular_sprint/kinetic_form_selection_certificate.py": "9cadce4736c666dc48a242ccfb46011d00b24189",
    "code/angular_sprint/runtime/kinetic_form_selection_receipt.json": "7469c729fd54cc7f558830c03fbae1bf510e5f32",
}
SOURCES = [
    ("Arithmon/K7", K7_SHA, "publications/ERRATUM_v3.5.md", "2bcd161a6684f93f83e3541cbd0682738fa21b70", "ERRATUM_v3.5.md"),
    ("Arithmon/K7", K7_SHA, "publications/papers/markdown/k7_framework_3_5_main.md", "3a0895fdaa9624eaed9d083e45566e75b3b1036b", "k7_framework_3_5_main.md"),
    ("Arithmon/K7", K7_SHA, "docs/openwave-candidate/honest_ledger.md", "39efc4eace29269f99a37fc94662fd757574efeb", "honest_ledger.md"),
    ("Arithmon/K7-Lean", K7L_SHA, "GIFT/Relations/GaugeSector.lean", "664e6e620de8f11d628a2c0beb3c5c90d6c190f2", "GaugeSector.lean"),
]


def canon(v):
    return (json.dumps(v, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode()


def rat(x):
    return str(Q(x))


def qdet(a, b, c):
    # columns are triples; calculate determinant without target coefficients.
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])


def quadratic(a=Q(0), b=Q(0)):
    return {"rational_part": str(a), "sqrt2_part": str(b)}


def q2add(x, y):
    return quadratic(Q(x["rational_part"]) + Q(y["rational_part"]), Q(x["sqrt2_part"]) + Q(y["sqrt2_part"]))


def q2mul(x, y):
    a, b = Q(x["rational_part"]), Q(x["sqrt2_part"])
    c, d = Q(y["rational_part"]), Q(y["sqrt2_part"])
    return quadratic(a*c + 2*b*d, a*d + b*c)


def q2scale(x, s):
    return quadratic(Q(x["rational_part"])*s, Q(x["sqrt2_part"])*s)


def q2sub(x, y):
    return q2add(x, q2scale(y, -1))


def q2zero(x):
    return Q(x["rational_part"]) == 0 and Q(x["sqrt2_part"]) == 0


def custody():
    items = []
    for repo, commit, path, blob, name in SOURCES:
        raw = (VENDOR / name).read_bytes()
        if git_blob(raw) != blob:
            raise ValueError(f"vendored K7 source does not match frozen git blob: {path}")
        items.append({"repository": repo, "commit": commit, "path": path, "git_blob": blob,
                      "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw), "vendored_path": f"code/angular_sprint/external/k7_gauge/{name}"})
    return items


def git_blob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def verify_oph_sources():
    for rel, expected in OPH_PINS.items():
        raw = (ROOT / rel).read_bytes()
        if git_blob(raw) != expected:
            raise ValueError(f"OPH frozen source changed: {rel}")


def build():
    verify_oph_sources()
    # Values are independently rederived from the pinned Lean definitions.
    k = (Q(10, 3), Q(2), Q(2))
    nG, nH = 3, 1
    b = (Q(20, 9)*nG + Q(1, 6)*nH,
         -Q(22, 3) + Q(4, 3)*nG + Q(1, 6)*nH,
         -Q(11) + Q(4, 3)*nG)
    cof = (k[1]*b[2]-k[2]*b[1], k[2]*b[0]-k[0]*b[2], k[0]*b[1]-k[1]*b[0])
    scale = Q(1)
    integer = tuple(int(c*9) for c in cof)
    expected = (Q(-23,3), Q(37), Q(-218,9))
    if b != (Q(41,6), Q(-19,6), Q(-7)) or cof != expected or integer != (-69,333,-218):
        raise ValueError("OPH source arithmetic drifted")

    # Low-energy electromagnetic alpha, MZ angle and MZ strong coupling are
    # deliberately assembled only as a forbidden plumbing/control vector.
    A0 = Q(267489, 1952)
    mixed_r = (Q(10,13)*A0, Q(3,13)*A0, Q(0))
    mixed_s = (Q(0), Q(0), Q(6))
    Dmixed = quadratic(sum(c*x for c, x in zip((69,-333,218), mixed_r)), sum(c*x for c, x in zip((69,-333,218), mixed_s)))
    if Dmixed != quadratic(Q(-82654101,25376), Q(1308)):
        raise ValueError(f"mixed-scale exact arithmetic mismatch: {Dmixed}")
    # Derive the B-test scale from alpha_em^-1=alpha_2^-1+alpha_Y^-1,
    # alpha_Y^-1=(5/3)alpha_1^-1, and the pinned GUT ratio 14:7:2.
    # Thus alpha_em^-1=(7+70/3)t=(91/3)t=91 sqrt(2), so t=3 sqrt(2).
    t = q2scale(quadratic(0, 91), Q(3, 91))
    if t != quadratic(0, 3):
        raise ValueError("B-test common-factor reconstruction failed")
    gutvec = [q2scale(t, Q(n)) for n in (14, 7, 2)]
    bvec = [q2scale(t, Q(5,3)*14), gutvec[1], gutvec[2]]
    if bvec != [quadratic(0,n) for n in (70,21,6)]:
        raise ValueError("B-test hypercharge conversion failed")
    Db = quadratic(0, sum(c*Q(x["sqrt2_part"]) for c,x in zip((69,-333,218), bvec)))
    if Db != quadratic(0, -855):
        raise ValueError(f"B-test exact arithmetic mismatch: {Db}")
    xzero=(Q(333),Q(69),Q(0))
    partial_x=qdet((Q(3,5)*xzero[0],xzero[1],xzero[2]),k,b)
    partial_k=qdet(xzero,(Q(3,5)*k[0],k[1],k[2]),b)
    partial_b=qdet(xzero,k,(Q(3,5)*b[0],b[1],b[2]))
    if (partial_x,partial_k,partial_b)!=(Q(5106,5),Q(-644),Q(-1886,5)):
        raise ValueError("partial normalization hostile controls drifted")
    gut_zero=qdet((Q(3,5)*xzero[0],xzero[1],xzero[2]),
                  (Q(3,5)*k[0],k[1],k[2]),
                  (Q(3,5)*b[0],b[1],b[2]))
    if qdet(xzero,k,b)!=0 or gut_zero!=Q(3,5)*qdet(xzero,k,b) or gut_zero!=0:
        raise ValueError("coherent Y/GUT normalization zero equivalence failed")

    src = custody()
    by_name = {Path(x["path"]).name: (VENDOR / x["vendored_path"].split("/")[-1]).read_text(encoding="utf-8") for x in src}
    erratum = by_name["ERRATUM_v3.5.md"]
    paper = by_name["k7_framework_3_5_main.md"]
    lean = by_name["GaugeSector.lean"]
    if "α⁻¹(0)" not in erratum or "267489/1952" not in erratum:
        raise ValueError("pinned erratum no longer types 267489/1952 as the Thomson/low-energy relation")
    if "91\\sqrt{2}" not in paper or "M_GUT" not in paper or "MSSM" not in paper:
        raise ValueError("pinned paper B-test provenance text missing")
    if "alpha_inv_complete_num = 267489" not in lean or "alpha_s_sq_structure" not in lean:
        raise ValueError("pinned Lean gauge identities missing")

    type1 = {
        "alpha_em_inverse": {"exact":"267489/1952", "field_of_definition":"electromagnetic fine-structure inverse", "formal_source":"pinned K7-Lean GaugeSector.lean; source-typed by pinned K7 erratum", "physical_interpretation":"low-energy Thomson-limit structural relation", "scale":"alpha_em^-1(0)", "scheme":"Thomson definition; running-scheme fields UNKNOWN", "u1_normalization":"electromagnetic; not an isolated U(1)_Y coupling", "matter_content_beta_branch":"not an RGE output", "derivation_class":"Type I structural relation with imposed det(g) normalization", "public_data_used_in_construction":TYPE1_PROVENANCE, "status":"historical structural relation"},
        "sin2_theta_W": {"exact":"3/13", "field_of_definition":"weak mixing angle", "formal_source":"pinned K7-Lean GaugeSector.lean; main paper comparison convention", "physical_interpretation":"topological relation compared at M_Z; separate GUT-scale reading discussed", "scale":"M_Z comparison; scale of structural relation otherwise UNKNOWN", "scheme":"UNKNOWN", "u1_normalization":"mixing-angle relation alone does not type alpha_Y", "matter_content_beta_branch":"not derived by an RGE lane", "derivation_class":"Type I structural relation", "public_data_used_in_construction":TYPE1_PROVENANCE, "status":"historical structural relation"},
        "alpha_s": {"exact":"sqrt(2)/12", "field_of_definition":"strong coupling alpha_s", "formal_source":"pinned K7-Lean GaugeSector.lean; paper gauge-structure discussion", "physical_interpretation":"structural relation presented alongside M_Z gauge-coupling comparisons", "scale":"M_Z comparison in paper; structural formula's scale assignment UNKNOWN", "scheme":"UNKNOWN", "u1_normalization":"not applicable", "matter_content_beta_branch":"not an RGE output", "derivation_class":"Type I structural relation", "public_data_used_in_construction":TYPE1_PROVENANCE, "status":"historical structural relation"},
        "assembled_vector": {"values":{"alpha_Y_inverse":"(10/13)*(267489/1952)","alpha_2_inverse":"(3/13)*(267489/1952)","alpha_3_inverse":"6*sqrt(2)"}, "classification":"MIXED_SCALE_VECTOR", "same_object_contract":False, "rejection_reasons":["alpha_em^-1 is explicitly alpha_em^-1(0)","sin2 and alpha_s have M_Z comparison/interpretation while alpha_em^-1(0) is low energy","no pinned common-scale identification","scheme and thresholds are UNKNOWN"]}
    }
    btest = {"inputs":{"alpha_em_inverse":"91*sqrt(2)","sin2_theta_W":"3/13","alpha_s":"sqrt(2)/12"},"normalization_derivation":{"GUT_ratio":"14:7:2","alpha_em_inverse_relation":"alpha_2_inverse + alpha_Y_inverse","alpha_Y_inverse_from_alpha_1_GUT_inverse":"(5/3)*alpha_1_GUT_inverse","common_factor":"alpha_em_inverse/(7+70/3)=3*sqrt(2)","GUT_inverse_vector":["42*sqrt(2)","21*sqrt(2)","6*sqrt(2)"],"Y_inverse_vector":["70*sqrt(2)","21*sqrt(2)","6*sqrt(2)"]}, "values":{"alpha_1_GUT_inverse":"42*sqrt(2)","alpha_2_inverse":"21*sqrt(2)","alpha_3_inverse":"6*sqrt(2)","alpha_Y_inverse":"70*sqrt(2)","alpha_2_inverse_Y":"21*sqrt(2)","alpha_3_inverse_Y":"6*sqrt(2)"},
             "field_of_definition":"conditional inverse-coupling ray reconstructed from B-test identity", "formal_source":"pinned K7 main paper §5.4 and cited holonomy sequence; arithmetic formula is paper-level, not solely Lean theorem", "physical_interpretation":"conditional MSSM B=7/5 relation; main paper calls the ray a GUT-scale quantity", "derivation_class":"conditional B-test / historical holonomy-sequence relation", "public_data_used_in_construction":BTEST_PROVENANCE,
             "ratio_GUT_normalized":"14:7:2", "normalization":"GUT alpha_1, converted coherently to Y for OPH control", "scale":"CONFLICTED: B-test theorem display says M_Z; proof note and holonomy-sequence discussion call it the GUT/M_GUT exact scale", "scheme":"UNKNOWN", "field_content":"MSSM, three generations, one Higgs doublet pair", "rge_branch":"one-loop MSSM B-test; exact holonomy ray is not an OPH SM boundary vector", "physical_role":"conditional boundary ray / holonomy sequence", "same_object_contract":False, "rejection_reasons":["pinned paper conflicts internally on M_Z versus M_GUT scale for this ray","MSSM beta branch differs from frozen OPH SM (nG,nH)=(3,1)","scheme and threshold matching UNKNOWN","boundary scale/direction is not OPH's unresolved alpha_inverse comparison object"]}
    doc = {"schema":SCHEMA,"pins":{"oph_repository":"Arithmon/observer-patch-holography","oph_commit":OPH_SHA,"oph_sources":OPH_PINS,"k7_repository":"Arithmon/K7","k7_commit":K7_SHA,"k7_fetched_main_head":K7_SHA,"k7_lean_repository":"Arithmon/K7-Lean","k7_lean_commit":K7L_SHA,"k7_lean_fetched_main_head":K7L_SHA,"k7_sources":src},
           "oph_contract":{"vector_role":"alpha_inverse=(alpha_Y^-1,alpha_2^-1,alpha_3^-1); issue-639 column sealed and not read", "u1_normalization":"declared hypercharge Y", "kinetic_column":[rat(x) for x in k],"beta_column":[rat(x) for x in b],"nG":nG,"nH":nH,"cofactors":[rat(x) for x in cof],"integer_zero_locus":[69,-333,218],"matter_branch_selected_physical":False,"scheme":"UNKNOWN","thresholds":"single-threshold premise of frozen imported one-loop law; physical threshold realization open"},
           "k7_candidates":{"type1_structural_bundle":type1,"b_test_holonomy_bundle":btest,"rge_bundle":{"found":True,"classification":"CONDITIONAL_NUMERICAL_RGE_OUTPUT__NO_EXACT_SOURCE_CLEAN_REPRODUCER","evidence_scope":["vendored K7 publications/ERRATUM_v3.5.md","vendored K7 publications/papers/markdown/k7_framework_3_5_main.md","vendored K7 docs/openwave-candidate/honest_ledger.md","vendored K7-Lean GIFT/Relations/GaugeSector.lean"],"record":"The pinned K7 source surfaces consumed by this packet establish a common-M_Z numerical output row, but do not establish an exact source-clean triple or an RGE reproducer.","displayed_outputs":{"sin2_theta_W_MZ":"0.2377","alpha_em_inverse_MZ":"131.19","alpha_s_all_MSSM":"0.1038","alpha_s_split_spectrum":"0.1224"},"derived_display_precision_alpha_Y_inverse":"about 100.006137 from (1-0.2377)*131.19; GUT-normalized alpha_1 used by MSSM lane", "derived_display_precision_alpha_2_inverse":"about 31.183863 from 0.2377*131.19", "alpha_3_inverse_split_spectrum":"about 1/0.1224; exact vector unavailable", "exact_vector":"NOT AVAILABLE; displayed decimal outputs are rounded", "field_of_definition":"running inverse gauge-coupling triple reconstructed from the common-M_Z output row", "scale":"M_Z per §5.3 output table", "scheme":"MS-bar per paper-wide M_Z gauge-coupling convention; implementation details UNKNOWN", "u1_normalization":"GUT-normalized alpha_1 in MSSM beta convention", "field_content":"MSSM above M_squark=3165 GeV; SM plus light gauginos below, b3=-5; tan(beta)=2", "rge_branch":"two-loop MSSM with split-spectrum step-function matching; alternate all-MSSM alpha_s output also shown", "boundary_role":"starts from alpha_GUT^-1=25.3 and sin2(theta_W)=3/13 at M_GUT; paper calls these topological boundary values", "direction_of_running":"M_GUT to M_Z", "derivation_class":"Type III conditional numerical RGE output", "public_data_used_in_construction":"UNKNOWN from these pinned source surfaces; output table also contains experimental comparison values, but none are parsed or consumed as cross-inputs", "reproducible_artifact":"NOT ESTABLISHED BY PINNED SOURCE SURFACES", "same_object_contract":False,"rejection_reasons":["MSSM/split-spectrum field-content and two-loop branch differ from frozen OPH SM one-loop branch","displayed values are rounded and have no exact source-clean reproducer","the pinned source surfaces do not resolve whether experimental data entered the numerical derivation"]}},
           "controls":{"mixed_scale_probe":{"xY":["(10/13)*(267489/1952)","(3/13)*(267489/1952)","6*sqrt(2)"],"exact_D":quadratic(Q(-82654101,25376),Q(1308)),"zero_status":"EXACT_NONZERO","classification":"NONPHYSICAL_MIXED_SCALE_CONTROL","eligible_for_framework_cross_verdict":False},"b_test_branch_mismatch_probe":{"xY":["70*sqrt(2)","21*sqrt(2)","6*sqrt(2)"],"exact_D":Db,"zero_status":"EXACT_NONZERO","classification":"RGE_BRANCH_MISMATCH_CONTROL","eligible_for_framework_cross_verdict":False},"normalization_covariance":{"Y_to_GUT_replayed":True,"coherent_first_row_factor":"3/5","all_x_k_b_first_entries_transformed":True,"raw_determinant_row_factor":"3/5","zero_test_input":{"xY":["333","69","0"],"integer_statistic":"0"},"Y_and_GUT_zero_status":"exact zero in both conventions","partial_row_transform_determinants":{"x_only":"5106/5","k_only":"-644","b_only":"-1886/5"},"zero_status_preserved":True,"hostile_partial_row_transform_rejected":True}},
           "comparison":{"admissible_candidate":None,"exact_determinant":None,"exact_zero_status":None,"same_object_contract":False,"reason":"no K7 candidate closes all common-object fields"},
           "primary_verdict":"NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES",
           "trust_boundary":{"public_measurement_consumed_as_cross_input":False,"public_measurement_values_parsed_by_cross_control":False,"vendored_main_paper_contains_measurement_comparison_text":True,"oph_sealed_comparison_opened":False,"k7_modified":False,"oph_matter_branch_selected":False,"physical_cross_framework_identity_claimed":False,"oph_port_response_identified_with_k7":False,"matter_trace_physically_selected":False,"off_plane_does_not_falsify_oph_globally":True}}
    return doc


def main():
    p=argparse.ArgumentParser(); p.add_argument("--check",action="store_true"); a=p.parse_args()
    data=canon(build())
    if a.check:
        if not OUT.exists() or OUT.read_bytes()!=data: raise SystemExit("receipt is stale; regenerate with this producer")
        print("PASS k7 gauge kinetic cross producer check")
    else:
        OUT.parent.mkdir(parents=True,exist_ok=True); OUT.write_bytes(data); print(f"wrote {OUT.relative_to(ROOT)}")

if __name__=="__main__": main()
