#!/usr/bin/env python3
"""Replay the frozen retrospective RSCC law at independently calibrated P.

The output is conditional native mass-chart coordinates. It is neither a
common-scale Yukawa matrix nor a source derivation of the flavor law. The
seven original RSCC hypotheses, representation attachment, mass charts and
absolute unit are stated explicitly. FULL is always the primary law; the
previously defined lower-order ablation is retained only as a control.

Default generation checks a cached gauge root directly against all D10
equations. ``--refresh-source`` also reruns the original gauge-only root solve;
quark and Thomson-continuation methods are disabled. No quark comparison
values enter either path. No legacy RSCC artifact is rewritten.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import sys
from typing import Any

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_source_w5_response_constraints import verify as verify_w5_constraints

ROOT = Path(__file__).resolve().parents[3]
SPEC_PATH = Path(__file__).with_name("conditional_quark_mass_spec.json")
OUT_PATH = ROOT / "code/particles/runs/flavor/conditional_quark_mass_replay.json"
ALPHA_PATH = ROOT / "code/P_derivation/codata_2022_alpha_fixture.json"
GAUGE_PATH = ROOT / "code/P_derivation/paper_math.py"
RESPONSE_PATH = ROOT / "code/particles/runs/flavor/source_w5_response_constraints.json"
RESPONSE_DEPENDENCIES = [
    "code/particles/flavor/source_w5_response_constraints.py",
    "code/particles/flavor/source_w5_response_constraints_spec.json",
    "code/particles/flavor/verify_source_w5_response_constraints.py",
    "code/a5_closure/manifests/echosahedral_federation_reference.json",
    "code/a5_closure/manifests/a3_scheduler_kernel_reference.json",
    "code/a5_closure/manifests/record_counting_mechanism_reference.json",
    "code/a5_closure/source_repair_generator_certificate.py",
]
RESPONSE_CONTRACT = {
    "receipt": "code/particles/runs/flavor/source_w5_response_constraints.json",
    "projection": {"family_dimension": 5, "line_complement_dimension": 4,
                   "scalar_fraction": "1/5", "centered_fraction": "4/5",
                   "common_exposure": "4/15"},
    "physical_boundary": "The dimensions and normalized native projector traces are certified. Choosing the face-axis background, Casimir-line common-exposure readout, scalar budget assignment and its physical quark attachment remains conditional. Composite modules, signed responses and Gaussian/contraction laws are not selected by these coefficients.",
}
HISTORICAL_PATHS = [
    "code/particles/flavor/quark_rscc_completion_candidate.py",
    "code/particles/flavor/verify_quark_rscc_module_arithmetic.py",
    "code/particles/runs/flavor/quark_rscc_completion_candidate_audit.json",
]
SCHEMA = "oph.conditional_quark_mass_replay.v1"
CLAIM = "conditional_retrospective_RSCC_native_mass_chart_replay"
STATUS = "conditional_forward_replay_not_source_selected_or_independent_prediction"
COUNTS = {"N_g": 3, "N_c": 3, "d_std_S3": 2, "order_S3": 6,
          "beta_EW": 4, "oriented_repair_slots": 24, "orientation_dimension": 2}
CHARTS = {
    **{q: {"scheme": "MSbar", "scale": "2 GeV"} for q in ("u", "d", "s")},
    **{q: {"scheme": "MSbar", "scale": "self_scale_mu_equals_running_mass"} for q in ("c", "b")},
    "t": {"scheme": "declared_pole_coordinate", "scale": "not_a_running_scale",
          "boundary": "pole_attachment_is_assumed_not_derived; experimental_extraction_matching_is_separate"},
}
GUARDS = {
    "conditional_native_chart_replay": True,
    "formula_discovery_target_informed": True,
    "module_and_sign_ledger_target_informed": True,
    "quark_reference_values_used_at_runtime": False,
    "old_mixed_branch_scale_used": False,
    "source_selected_flavor_law": False,
    "physical_quark_mass_derivation": False,
    "common_scale_yukawa_matrix": False,
    "top_pole_attachment_assumed": True,
    "independent_prediction": False,
    "new_blind_evidential_weight": False,
    "controls_used_for_model_selection": False,
}
HYPOTHESES = {
    "H1_family_response_carrier": "The physical response is identified with the native A5-irreducible W5 and a supplied face-axis D3 background. Its restriction F=1+2+2 and dimensions5,4 are certified; this physical identification and background are assumed.",
    "H2_module_incidence": "The eight declared composite module incidences are the physical channels.",
    "H3_effects_and_signs": "The remaining composite effect ranks, orientations and sector signs are physical. The common exposure specifically uses the normalized quadratic color-Casimir observable on the native D3 line, not a rank-four invariant color projector; that physical response choice is assumed.",
    "H4_cumulant_law": "A5 symmetry suffices for normalized trace on native W5. Isotropy on the composite modules, Gaussian two-cumulant truncation and the order of color contraction versus nonlinear response remain hypotheses; negative terms require a separate signed-response law.",
    "H5_family_attachment": "The count-only 24-slot register has an additional physical family non-singlet attachment.",
    "H6_readout_laws": "The frozen heat-time, even-response, L/Q log-spectrum and affine-mean equations are the physical readout laws.",
    "H7_residual_selection": "The declared residual functional is physically minimized; this is not inferred from the existence of its mathematical minimum.",
    "H8_chart_attachment": "The six native coordinates have exactly the declared light-MSbar, heavy-self-scale and top-pole charts; a physical top pole and extraction relation are not derived.",
    "H9_normalization": "The D10 transmutation law supplies v/E_star and the external E_star GeV calibration supplies the absolute unit.",
    "H10_gauge_trunk": "The finite SU(2)/SU(3) heat closure, one-loop coefficients (33/5,1,-3), tree weak fixed point and beta_EW=4 define the declared D10 branch.",
    "H11_representation_attachment": "The selected D3 group heat law is represented faithfully by six conjugation channels on native W5, without an added regular submodule. Selecting this background, transposition-only jumps and the physical heat/three-generation readout remains assumed; it is not the native full-seam generator or the A5-triplet restriction.",
    "H12_ordering": "Increasing entries of the up/down output vectors label (u,c,t)/(d,s,b); no CKM or right-family frame is inferred.",
}
MENU = {
    "mean_u_linear": "79/15", "mean_d_linear": "4/15",
    "a_u_linear": "29/5", "a_d_linear": "33/32",
    "mean_u_quadratic": "1/29", "mean_d_quadratic": "-1/432",
    "a_u_quadratic": "1/22", "a_d_quadratic": "1/420",
    "delta_g_P": "1/1008", "delta_g_quadratic_color": "1/432",
    "delta_g_quadratic_relabel": "-1/1584", "q_u": "-rho/10", "q_d": "-1/4",
    "module_dimensions": [29, 432, 22, 32, 840, 1008, 432, 1584],
}
EQUATIONS = {
    "root": "P=phi+sqrt(pi)/alpha_inverse; w=pi*alpha_U(P)",
    "heat": "tau=P/4-w/5; r=exp(-3*tau); rho=3/(2+r); x=(r-1)/(r+1)",
    "mean_up": "3*P+79*w/15+w^2/29", "mean_down": "2*P+4*w/15-w^2/432",
    "a_up": "3*P+29*w/5+w^2/22", "a_down": "2*P+33*w/32+w^2/420",
    "delta_g": "P/1008+w^2/432-w^2/1584", "common": "g_ch/v=2*exp(-2*pi+delta_g)",
    "b_up": "a_u*(-rho*x+rho-x-1)/((1+rho)*(x^2-1))-w*rho/10",
    "b_down": "a_d*(-rho*x-rho-x+1)/((1+rho)*(x^2-1))-w/4",
    "basis": "L=ctr(-1,x,1); Q=ctr(1,x^2,1)",
    "affine": "A=1/(2*(1+rho-x^2)); B=1/(2*(1-x^2-x^2/(1+rho))); sigma=(mean_u+mean_d)/2; eta=(mean_u-mean_d)/2",
    "up": "m_up/v=(g_ch/v)*exp(-A*sigma+B*eta+a_u*L+b_u*Q)",
    "down": "m_down/v=(g_ch/v)*exp(-A*sigma-B*eta+a_d*L+b_d*Q)",
    "unit": "v/E_star=P^(-1/2)*exp(-2*pi/(4*alpha_U)); m_GeV=(m/v)*(v/E_star)*E_star_GeV",
    "control": "Remove only the four quadratic corrections in mean/a and all delta_g; keep the same remaining law and inputs.",
}


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_spec(path: Path = SPEC_PATH) -> dict[str, Any]:
    spec = json.loads(path.read_text())
    expected = {"schema", "claim_class", "status", "calibrations", "structural_counts",
                "charts", "guards", "hypotheses", "coefficient_menu", "equations",
                "gauge_solver", "source_cache", "historical_law_sha256", "source_response"}
    require(set(spec) == expected, "unexpected spec field or hidden input")
    for key, value in {"schema": SCHEMA, "claim_class": CLAIM, "status": STATUS,
                       "structural_counts": COUNTS, "charts": CHARTS, "guards": GUARDS,
                       "hypotheses": HYPOTHESES, "coefficient_menu": MENU, "equations": EQUATIONS,
                       "source_response": RESPONSE_CONTRACT}.items():
        require(spec[key] == value, f"frozen {key} changed")
    require(spec["calibrations"] == {
        "alpha_inverse": {"value": "137.035999177", "role": "external_measured_calibration"},
        "E_star_GeV": {"value": "1.220890e19", "role": "external_dimensionful_unit_calibration"}}, "calibration or unit changed")
    fixture = json.loads(ALPHA_PATH.read_text())
    require(fixture["inverse_fine_structure_constant"]["value"] == "137.035999177", "alpha fixture drift")
    require(spec["gauge_solver"] == {"precision": 50, "su2_cutoff": 80, "su3_cutoff": 60,
                                     "alpha_U_domain": ["0.02", "0.08"]}, "gauge recipe changed")
    require(set(spec["historical_law_sha256"]) == set(HISTORICAL_PATHS), "historical ancestry changed")
    for rel, digest in spec["historical_law_sha256"].items():
        require(sha(ROOT / rel) == digest, "frozen historical law or audit changed")
    return spec


def _number(value: str) -> mp.mpf:
    require(isinstance(value, str), "numeric inputs must be decimal strings")
    number = mp.mpf(value)
    require(bool(mp.isfinite(number)), "nonfinite numeric value")
    return number


def verify_source_packet(packet: dict[str, Any], spec: dict[str, Any]) -> None:
    """Direct D10 equation checks; does not trust a cached root label."""
    keys = {"P", "alpha_U", "v_over_E_star", "mz_over_E_star", "mu_U_over_E_star",
            "alpha1_mz", "alpha2_mz", "alpha3_mz"}
    require(set(packet) == keys, "source packet field or target injection")
    with mp.workdps(70):
        z = {key: _number(value) for key, value in packet.items()}
        p = (1+mp.sqrt(5))/2 + mp.sqrt(mp.pi)/_number(spec["calibrations"]["alpha_inverse"]["value"])
        def close(actual: mp.mpf, expected: mp.mpf, label: str) -> None:
            require(abs(actual-expected) <= mp.mpf("2e-16")*max(abs(expected), mp.mpf("1e-40")), f"inconsistent gauge {label}")
        close(z["P"], p, "P")
        require(mp.mpf(".02") < z["alpha_U"] < mp.mpf(".08"), "root outside declared domain")
        v = p**(-mp.mpf(".5"))*mp.exp(-2*mp.pi/(4*z["alpha_U"]))
        mu = mp.exp(-2*mp.pi)*p**(mp.mpf(1)/6)
        close(z["v_over_E_star"], v, "transmutation scale")
        close(z["mu_U_over_E_star"], mu, "unification scale")
        for key, beta in zip(("alpha1_mz", "alpha2_mz", "alpha3_mz"), (mp.mpf(33)/5, mp.mpf(1), mp.mpf(-3))):
            value = 1/(1/z["alpha_U"] + beta/(2*mp.pi)*mp.log(mu/z["mz_over_E_star"]))
            close(z[key], value, key)
        mz = v*mp.sqrt(4*mp.pi*(mp.mpf(3)/5*z["alpha1_mz"]+z["alpha2_mz"]))/2
        close(z["mz_over_E_star"], mz, "weak fixed point")
        t2, t3 = 4*mp.pi**2*z["alpha2_mz"], 4*mp.pi**2*z["alpha3_mz"]
        su2 = [(mp.mpf(n+1),mp.mpf(n)*(n+2)/4) for n in range(81)]
        su3 = [(mp.mpf((a+1)*(b+1)*(a+b+2))/2, mp.mpf(a*a+b*b+a*b+3*a+3*b)/3)
               for a in range(61) for b in range(61)]
        def ell(terms: list, time: mp.mpf) -> mp.mpf:
            weights = [dim*mp.exp(-time*c) for dim,c in terms]
            return mp.fsum(weight*mp.log(term[0]) for weight,term in zip(weights,terms))/mp.fsum(weights)
        require(abs(ell(su2,t2)+ell(su3,t3)-p/4) < mp.mpf("2e-16"), "D10 heat closure failed")


def refresh_source_packet(spec: dict[str, Any]) -> dict[str, str]:
    module_spec = importlib.util.spec_from_file_location("_oph_conditional_gauge_math", GAUGE_PATH)
    require(module_spec is not None and module_spec.loader is not None, "gauge source unavailable")
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[module_spec.name] = module
    module_spec.loader.exec_module(module)
    class GaugeOnlyContext(module.PaperMathContext):
        def _derive_stage5_integer_vectors(self):
            return {}
        def diagonal_quark_masses(self, *args, **kwargs):
            raise RuntimeError("quark continuation forbidden in conditional gauge input")
        def structured_thomson_running(self, *args, **kwargs):
            raise RuntimeError("Thomson closure forbidden in independently calibrated P input")
        def structured_thomson_running_asymptotic(self, *args, **kwargs):
            raise RuntimeError("Thomson closure forbidden in independently calibrated P input")
        def solve_closure(self, *args, **kwargs):
            raise RuntimeError("P closure forbidden: P is independently calibrated")
    ctx = GaugeOnlyContext(precision=50, su2_cutoff=80, su3_cutoff=60)
    p = ctx.p_from_inverse_alpha(spec["calibrations"]["alpha_inverse"]["value"])
    point = ctx.build_d10_from_p(p)
    mapping = {"P":"p", "alpha_U":"alpha_u", "v_over_E_star":"v", "mz_over_E_star":"mz_run",
               "mu_U_over_E_star":"mu_u", "alpha1_mz":"alpha1_mz", "alpha2_mz":"alpha2_mz", "alpha3_mz":"alpha3_mz"}
    packet = {key:str(getattr(point,attr)) for key,attr in mapping.items()}
    verify_source_packet(packet,spec)
    return packet


def load_source_response(path: Path = RESPONSE_PATH) -> dict[str, Any]:
    certificate = json.loads(path.read_text())
    verify_w5_constraints(certificate)
    carrier, color = certificate["carrier"], certificate["color"]
    projection = {"family_dimension": carrier["dimension"],
                  "line_complement_dimension": carrier["D3_line_complement_dimension"],
                  "scalar_fraction": carrier["normalized_native_projector_weights"][0],
                  "centered_fraction": carrier["normalized_native_projector_weights"][1],
                  "common_exposure": color["normalized_trace"]}
    require(projection == RESPONSE_CONTRACT["projection"], "source-response coefficient projection changed")
    return projection


def _fmt(value: mp.mpf) -> str:
    return mp.nstr(value, 40)


def evaluate_law(packet: dict[str, str], *, full: bool, energy: str, response: dict[str, Any]) -> dict[str, Any]:
    with mp.workdps(80):
        p, alpha, v, unit = map(_number, (packet["P"],packet["alpha_U"],packet["v_over_E_star"],energy))
        def rational(value):
            q = Fraction(value)
            return mp.mpf(q.numerator)/q.denominator
        nf = response["family_dimension"]; f0 = response["line_complement_dimension"]
        scalar = rational(response["scalar_fraction"])
        centered = rational(response["centered_fraction"])
        exposure = rational(response["common_exposure"])
        w=mp.pi*alpha; tau=p/f0-w*scalar; r=mp.exp(-3*tau); rho=3/(2+r); x=(r-1)/(r+1)
        mean_u=3*p+(nf+exposure)*w; mean_d=2*p+exposure*w
        a_u=3*p+(nf+centered)*w; a_d=2*p+mp.mpf(33)*w/32; delta_g=mp.mpf(0)
        if full:
            mean_u+=w*w/29; mean_d-=w*w/432; a_u+=w*w/22; a_d+=w*w/420
            delta_g=p/1008+w*w/432-w*w/1584
        b_u_ray=a_u*(-rho*x+rho-x-1)/((1+rho)*(x*x-1))
        b_d_ray=a_d*(-rho*x-rho-x+1)/((1+rho)*(x*x-1))
        q_u=-w*rho/10; q_d=-w/4; b_u=b_u_ray+q_u; b_d=b_d_ray+q_d
        A=1/(2*(1+rho-x*x)); B=1/(2*(1-x*x-x*x/(1+rho)))
        sigma=(mean_u+mean_d)/2; eta=(mean_u-mean_d)/2
        linear=[z-x/3 for z in (-mp.mpf(1),x,mp.mpf(1))]
        quadratic=[z-(2+x*x)/3 for z in (mp.mpf(1),x*x,mp.mpf(1))]
        common=2*mp.exp(-2*mp.pi+delta_g)
        up=[common*mp.exp(-A*sigma+B*eta+a_u*linear[i]+b_u*quadratic[i]) for i in range(3)]
        down=[common*mp.exp(-A*sigma-B*eta+a_d*linear[i]+b_d*quadratic[i]) for i in range(3)]
        require(0<up[0]<up[1]<up[2] and 0<down[0]<down[1]<down[2], "species ordering failed")
        masses=dict(zip(("u","c","t","d","s","b"),up+down))
        coefficient_names=("w","tau","r","rho","x","mean_u","mean_d","a_u","a_d","b_u_ray","b_d_ray","q_u","q_d","b_u","b_d","A","B","delta_g","g_ch_over_v")
        values=(w,tau,r,rho,x,mean_u,mean_d,a_u,a_d,b_u_ray,b_d_ray,q_u,q_d,b_u,b_d,A,B,delta_g,common)
        return {"coefficients":dict(zip(coefficient_names,map(_fmt,values))),
                "basis":{"linear":list(map(_fmt,linear)),"quadratic":list(map(_fmt,quadratic))},
                "mass_over_v":{k:_fmt(z) for k,z in masses.items()},
                "mass_over_E_star":{k:_fmt(z*v) for k,z in masses.items()},
                "native_mass_chart_GeV":{k:_fmt(z*v*unit) for k,z in masses.items()},
                "light_ratios":{"u_over_d":_fmt(masses["u"]/masses["d"]),
                                "s_over_mean_ud":_fmt(2*masses["s"]/(masses["u"]+masses["d"]))}}


def build_payload(*, spec_path: Path = SPEC_PATH, refresh_source: bool = False) -> dict[str, Any]:
    spec=load_spec(spec_path)
    response=load_source_response()
    packet=copy.deepcopy(spec["source_cache"])
    verify_source_packet(packet,spec)
    if refresh_source:
        fresh=refresh_source_packet(spec)
        require(fresh==packet,"fresh gauge solve differs from frozen source cache")
    energy=spec["calibrations"]["E_star_GeV"]["value"]
    source_paths=[Path(__file__).resolve(),spec_path.resolve(),ALPHA_PATH,GAUGE_PATH]
    source_paths += [ROOT/rel for rel in HISTORICAL_PATHS+RESPONSE_DEPENDENCIES]+[RESPONSE_PATH]
    return {"schema":SCHEMA,"claim_class":CLAIM,"status":STATUS,
            "inputs":{**copy.deepcopy(spec["calibrations"]),
                      "P":{"value":packet["P"],"definition":"phi+sqrt(pi)/alpha_inverse"},
                      "structural_counts":COUNTS},
            "source_packet":packet,
            "model":{"primary":"full_rscc","control":"lower_order_ablation_not_selected",
                     "charts":CHARTS,"coefficient_menu":MENU,"equations":EQUATIONS,
                     "gauge_solver":spec["gauge_solver"],
                     "source_response":RESPONSE_CONTRACT,
                     "numeric_input_paths":[ALPHA_PATH.relative_to(ROOT).as_posix(),SPEC_PATH.relative_to(ROOT).as_posix(),
                                            RESPONSE_PATH.relative_to(ROOT).as_posix()]+[p for p in RESPONSE_DEPENDENCIES if p.endswith(".json")],
                     "historical_law_paths":HISTORICAL_PATHS,
                     "precision_boundary":"Printed digits are deterministic numerical output, not certified accuracy, model error or calibration uncertainty."},
            "full_rscc":evaluate_law(packet,full=True,energy=energy,response=response),
            "lower_order_control":evaluate_law(packet,full=False,energy=energy,response=response),
            "hypotheses":HYPOTHESES,"guards":GUARDS,
            "source_pins":{p.relative_to(ROOT).as_posix():sha(p) for p in source_paths}}


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--refresh-source",action="store_true")
    parser.add_argument("--out",type=Path,default=OUT_PATH)
    args=parser.parse_args()
    text=json.dumps(build_payload(refresh_source=args.refresh_source),indent=2,sort_keys=True)+"\n"
    if args.check:
        require(args.out.read_text()==text,"conditional replay receipt is stale or altered")
        print("conditional quark mass replay: PASS")
    else:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(text)
        print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
