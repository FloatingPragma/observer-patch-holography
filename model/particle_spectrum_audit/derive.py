"""Recompute conditional coefficients and arithmetic comparison boundaries.

The input packet fixes a stipulated model branch. The multiloop outputs are
imported numerical inputs; replay_smdr.py checks their heavy-particle part.
No finite theorem identifies the lepton coefficients with physical poles.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
ORDER = ("e", "mu", "tau", "u", "d", "s", "c", "b", "t", "W", "Z", "H")
EXPONENTS = {"u": 6, "d": 6, "s": 4, "c": 3, "b": 2}
GEV = {"GeV": F(1), "MeV": F(1, 1000), "eV": F(1, 10**9)}


def read(name):
    return json.loads((HERE / name).read_text())


def forward_coefficients(packet):
    """Independent closed-form evaluation of the declared three-register law."""
    with mp.workdps(80):
        p, u, v = (mp.mpf(str(packet[k])) for k in
                   ("P_decimal", "alpha_U_decimal", "v_GeV"))
        phi = (1 + mp.sqrt(5))/2
        kappa = (u/50)/(1 + u/31 + u*u/310)
        chi = (u*u/512)/(1 - u/77)
        zeta = ((p-phi)**2/(21*mp.pi))/(1-u/27-u*u/135)
        rho, delta = mp.sqrt(2)*mp.exp(-chi), mp.mpf(2)/9+zeta
        roots = [1+rho*mp.cos(delta+2*mp.pi*j/3) for j in (1, 2, 0)]
        if min(roots) <= 0:
            raise ValueError("The positive-root circulant identity needs positive roots")
        squares = [x*x for x in roots]
        determinant = mp.exp(3*kappa)/(2*6**14)
        scale = (determinant/mp.fprod(squares))**(mp.mpf(1)/3)
        masses = [v*scale*x for x in squares]
        q = sum(masses)/sum(mp.sqrt(m) for m in masses)**2
        if abs(q-(1+mp.exp(-2*chi))/3) > mp.mpf("1e-70"):
            raise ArithmeticError("Circulant invariant identity failed")
        if abs(mp.fprod(masses)/v**3/determinant-1) > mp.mpf("1e-70"):
            raise ArithmeticError("Determinant identity failed")
        return {"masses_GeV": dict(zip(("e", "mu", "tau"), map(float, masses))),
                "kappa": float(kappa), "chi": float(chi), "zeta": float(zeta),
                "rho": float(rho), "delta": float(delta), "koide_Q": float(q),
                "determinant_over_v_cubed": float(determinant),
                "physical_pole_identification": "adopted comparison premise",
                "target_informed_architecture": True}


def shape(masses):
    with mp.workdps(60):
        roots = [mp.sqrt(str(x)) for x in masses]
        avg = sum(roots)/3
        q = sum(x*x for x in roots)/sum(roots)**2
        z = sum((x/avg-1)*mp.exp(-2j*mp.pi*j/3)
                for x, j in zip(roots, (1, 2, 0)))
        return {"koide_Q": float(q), "delta": float(mp.arg(z)),
                "rho": float(mp.sqrt(2*(3*q-1)))}


def reference_interval(row):
    scale = GEV[row["unit"]]
    center = F(row["central"])*scale
    return center, ((F(row["central"])-F(row["error_minus"]))*scale,
                    (F(row["central"])+F(row["error_plus"]))*scale)


def comparison(value, reference):
    value = F(str(value))
    center, (lo, hi) = reference_interval(reference)
    return {"value_GeV": float(value), "reference_GeV": float(center),
            "reference_interval_GeV": [float(lo), float(hi)],
            "residual_percent": float(100*(value/center-1)),
            "residual_ppm": float(10**6*(value/center-1)),
            "within_printed_interval": lo <= value <= hi}


def factor_intersection(values, references):
    """Exact rectangle feasibility for a common positive output multiplier."""
    bounds = {}
    for key, value in values.items():
        if value <= 0:
            raise ValueError("A common positive-factor check needs positive values")
        _, (lo, hi) = reference_interval(references[key])
        bounds[key] = lo/F(str(value)), hi/F(str(value))
    low_key = max(bounds, key=lambda k: bounds[k][0])
    high_key = min(bounds, key=lambda k: bounds[k][1])
    lo, hi = bounds[low_key][0], bounds[high_key][1]
    return {"lower_bound": float(lo), "upper_bound": float(hi),
            "lower_bound_species": low_key, "upper_bound_species": high_key,
            "intersection_nonempty": lo <= hi,
            "scope": "Printed marginal interval rectangle; no joint confidence level or theory-error bound inferred."}


def frozen_transport_exponent(n, output, target, base=6):
    """Invert output=A*base**(-n) with A held fixed, including transport."""
    if base <= 1 or output <= 0 or target <= 0:
        raise ValueError("Positive masses and base > 1 required")
    return n+math.log(output/target, base)


def ratio_interval(numerator, denominator):
    _, (nlo, nhi) = reference_interval(numerator)
    _, (dlo, dhi) = reference_interval(denominator)
    if dlo <= 0:
        raise ValueError("A positive denominator interval is required")
    return [float(nlo/dhi), float(nhi/dlo)]


def build(inputs=None, reference=None):
    inputs = read("inputs.json") if inputs is None else inputs
    reference = read("pdg2026_reference.json") if reference is None else reference
    refs = reference["masses"]
    forward = forward_coefficients(inputs["boundary_packet"])
    values = {**forward["masses_GeV"], **{k: row["value_GeV"]
              for k, row in inputs["imported_readouts"].items()}}
    replay_error = max(abs(values[k]/inputs["lepton_replay_reference_GeV"][k]-1)
                       for k in ("e", "mu", "tau"))
    if replay_error > 1e-13:
        raise ArithmeticError("Independent lepton law does not reproduce input packet")
    rows = []
    for key in ORDER:
        row = {"parameter": key, **comparison(values[key], refs[key]),
               "reference_convention": refs[key]["convention"],
               "reference_confidence_level": refs[key]["confidence_level"],
               "reference_source": reference["sources"][refs[key]["source"]]["URL"]}
        row["convention"] = ("Pole identification adopted for the target-informed circulant coefficient"
                             if key in forward["masses_GeV"] else inputs["imported_readouts"][key]["convention"])
        rows.append(row)
    pole = [float(reference_interval(refs[k])[0]) for k in ("e", "mu", "tau")]
    native = [values[k] for k in ("e", "mu", "tau")]
    relation = inputs["one_loop_lepton_relation"]
    controls = {}
    for k in ("e", "mu", "tau"):
        m = values[k]
        converted = m*(1+relation["delta_EW"]+relation["alpha_full"]/math.pi*
                       (1+1.5*math.log(relation["Q0_GeV"]/m)))
        if not math.isclose(converted, inputs["lepton_msbar_attachment_control_GeV"][k], rel_tol=1e-13):
            raise ArithmeticError("Retained one-loop control disagrees")
        controls[k] = comparison(converted, refs[k])
    with mp.workdps(60):
        exact_shape = [float((1+mp.sqrt(2)*mp.cos(mp.mpf(2)/9+2*mp.pi*j/3))**2)
                       for j in (1, 2, 0)]
    ratios = {}
    for label, i, j in (("mu_over_e", 1, 0), ("tau_over_e", 2, 0)):
        ratio, observed = native[i]/native[j], pole[i]/pole[j]
        keys = ("e", "mu", "tau")
        interval = ratio_interval(refs[keys[i]], refs[keys[j]])
        ideal = exact_shape[i]/exact_shape[j]
        ratios[label] = {"value": ratio, "reference": observed,
                         "residual_ppm": 1e6*(ratio/observed-1),
                         "printed_rectangle_interval": interval,
                         "within_printed_rectangle": interval[0] <= ratio <= interval[1],
                         "exact_Koide_Brannen_ratio": ideal,
                         "exact_Koide_Brannen_residual_ppm": 1e6*(ideal/observed-1),
                         "exact_point_inside_rectangle": interval[0] <= ideal <= interval[1]}
    quarks = {}
    for key, n in EXPONENTS.items():
        target = float(reference_interval(refs[key])[0])
        needed = frozen_transport_exponent(n, values[key], target)
        if not math.isclose(values[key]*6**(-(needed-n)), target, rel_tol=1e-14):
            raise ArithmeticError("Exponent inversion round trip failed")
        quarks[key] = {"declared_exponent": n, "frozen_transport_exponent": needed,
                       "shift": needed-n}
    sd_interval = ratio_interval(refs["s"], refs["d"])
    ud_ref = reference["additional_observables"]["u_over_d"]
    ud_interval = [float(F(ud_ref["central"])-F(ud_ref["error_minus"])),
                   float(F(ud_ref["central"])+F(ud_ref["error_plus"]))]
    additional = reference["additional_observables"]
    dm21 = additional["Delta_m21_squared"]
    dm32 = additional["Delta_m32_squared_NO"]
    a = (1-(95/96)*math.sqrt(653/848))/2
    gap_rectangle = [(float(dm21["central"])-float(dm21["error_minus"])) /
                     (float(dm32["central"])+float(dm32["error_plus"])),
                     (float(dm21["central"])+float(dm21["error_plus"])) /
                     (float(dm32["central"])-float(dm32["error_minus"]))]
    by_output = {"W_BW_GeV": "W", "Z_BW_GeV": "Z", "Higgs_pole_GeV": "H",
                 "top_QCD4_EW2_pole_GeV": "t"}
    vev_controls = []
    for control in inputs["substitution_controls"]:
        vev_controls.append({"substituted_reference_inputs": control["substituted_from_SM_reference"],
                             "outputs": {key: comparison(control["outputs"][out], refs[key])
                                         for out, key in by_output.items()},
                             "prospective_prediction": False})
    return {"schema": "oph.particle_spectrum_audit.result.v1", "forward_leptons": forward,
            "lepton_replay_max_relative_error": replay_error, "rows": rows,
            "outside_printed_intervals": sum(not r["within_printed_interval"] for r in rows),
            "lepton_shapes": {"declared": shape(native), "PDG_poles": shape(pole)},
            "lepton_ratios": ratios, "msbar_attachment_one_loop_controls": controls,
            "common_factor": {"all": factor_intersection(values, refs),
                              "leptons": factor_intersection({k: values[k] for k in ("e", "mu", "tau")}, refs),
                              "WZ": factor_intersection({k: values[k] for k in ("W", "Z")}, refs)},
            "quark_exponent_diagnostics": quarks,
            "quark_exponent_scope": "Multiplicative transport held fixed. Changing a Yukawa changes threshold feedback and c/b self-scale roots; these numbers are not exact source exponents.",
            "quark_same_scale_ratios": {
                "s_over_d": {"source": values["s"]/values["d"], "reference_rectangle": sd_interval,
                             "source_in_rectangle": sd_interval[0] <= values["s"]/values["d"] <= sd_interval[1]},
                "u_over_d": {"source": values["u"]/values["d"], "direct_reference_interval": ud_interval,
                             "source_in_interval": ud_interval[0] <= values["u"]/values["d"] <= ud_interval[1]}},
            "neutral": {"minimal_action_masses_GeV": [0, 0, 0], "minimal_action_squared_gaps_eV2": [0, 0],
                        "individual_neutrino_masses_measured": False,
                        "normal_gap_ratio_reference": float(dm21["central"])/float(dm32["central"]),
                        "normal_gap_ratio_reference_rectangle": gap_rectangle,
                        "conditional_unshifted_phase_bound": a*a/(1-2*a),
                        "bound_hypotheses": "p in [0,1/12], 0<h<=1/2, same-scale nonnegative unshifted eigenphases identified with masses; see particle_finite_response.",
                        "arbitrary_neutral_models_excluded": False},
            "massless_gauge": {
                "photon": {"conditional_mass_GeV": 0,
                           "reference_kind": refs["photon"]["kind"],
                           "reference_upper_limit_eV": refs["photon"]["upper_limit"],
                           "zero_is_experimentally_resolved": False},
                "gluon": {"conditional_mass_GeV": 0,
                          "reference_kind": refs["gluon"]["kind"],
                          "separate_free_color_pole_masses_measured": False}},
            "vev_substitution_controls": vev_controls,
            "vev_controls_scope": "Reference substitutions diagnose sensitivity; they do not select an attachment scale. Controls omit the separately stored tiny finite-light W/Z shifts.",
            "strong_coupling": inputs["strong_coupling"],
            "scope": {"physical_spectrum_identity_established": False,
                      "prospective_prediction": False, "theory_error_bound": None,
                      "joint_statistical_significance_computed": False,
                      "particle_abundance_distribution_computed": False,
                      "exact_equality_to_experimental_central_values_required": False},
            "source_sha256": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in ("derive.py", "inputs.json", "pdg2026_reference.json")}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build()
    path = HERE/"results.json"
    if args.check:
        if result != json.loads(path.read_text()):
            raise SystemExit("Spectrum comparison receipt differs; regenerate and review")
        print("Particle spectrum arithmetic and source bindings verified")
    else:
        path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
        print(path)


if __name__ == "__main__":
    main()
