"""Replay the frozen branch with a supplied, built SMDR source tree.

The boundary packet is stipulated. This replay does not solve the pixel
selection problem or identify finite coefficients with physical particles.
"""
from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

from scipy.optimize import brentq

from .derive import HERE, read


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def invoke(binary, *args, cwd=None):
    argv = [format(x, ".18g") if isinstance(x, float) else str(x) for x in args]
    done = subprocess.run([str(binary), *argv], cwd=cwd, check=True,
                          capture_output=True, text=True, timeout=120)
    row = json.loads(done.stdout.strip().splitlines()[-1])
    def finite(x):
        if isinstance(x, dict):
            return all(finite(v) for v in x.values())
        if isinstance(x, list):
            return all(map(finite, x))
        return not isinstance(x, float) or math.isfinite(x)
    if not finite(row):
        raise ArithmeticError("SMDR returned a non-finite field")
    return row


def compile_driver(smdr, work, name, sources):
    binary = work/name
    command = ["cc", "-O2", "-fcommon", "-I", str(smdr), "-I", str(smdr/"src"),
               "-I", str(HERE), *[str(HERE/p) for p in sources], "-L", str(smdr),
               "-lsmdr", "-l3vil", "-ltsil", "-lm", "-o", str(binary)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    return binary


def finite_mass_threshold(r, nq, alpha_s):
    """Retained small-mass terms of Martin, PRD 99, 033007, eq. 3.15."""
    if r == 0:
        return 0.0
    log_r = math.log(r)
    d2 = r*(8*log_r/15-76/75)+r*r*(9*log_r/70-1389/9800)
    d3 = r*((64*nq/135-451/81)*log_r**2+(84887/7290-128*nq/135)*log_r+2.77670-.22452*nq)
    d3 += r*r*((4*nq/35-239/270)*log_r**2+(580157/396900-6*nq/35)*log_r+.52092+.03556*nq)
    return alpha_s**2*d2/(12*math.pi**2)+alpha_s**3*d3/(64*math.pi**3)


def low_energy(light, state):
    """Five-loop running, four-loop matching, and source self-scale roots."""
    def run(current, target, nu=2, nd=3):
        result = invoke(light, "run", current["Q"], target, nu, nd, 3, 5,
                        current["alphaS"], current["alpha"])
        if result["code"] != 0 or not 0 < result["alphaS"] < 1:
            raise ArithmeticError("QCD+QED evolution failed")
        group = {"u": "u", "c": "u", "d": "d", "s": "d", "b": "d",
                 "e": "e", "mu": "e", "tau": "e"}
        return {"Q": target, "alphaS": result["alphaS"], "alpha": result["alpha"],
                "masses": {k: m*result["c"+group[k]] for k, m in current["masses"].items()}}

    def self_scale(current, species, nd):
        high = current["Q"]
        low = high
        for _ in range(80):
            low *= .8
            if run(current, low, nd=nd)["masses"][species] >= low:
                break
        else:
            raise ArithmeticError("No self-scale bracket")
        root = brentq(lambda q: run(current, q, nd=nd)["masses"][species]-q,
                      low, high, xtol=2e-13, rtol=1e-13)
        result = run(current, root, nd=nd)
        if abs(result["masses"][species]-root) > 1e-10:
            raise ArithmeticError("Self-scale root residual is too large")
        return result

    bottom = self_scale(state, "b", 3)
    match = invoke(light, "threshold", "d", bottom["Q"], bottom["masses"]["b"],
                   bottom["alphaS"], bottom["alpha"], 4, 4)
    masses = {k: v for k, v in bottom["masses"].items() if k != "b"}
    factors = {k: match["zu" if k in ("u", "c") else "zd" if k in ("d", "s") else "ze"]
               for k in masses}
    for key in ("c", "s", "u", "d"):
        factors[key] += finite_mass_threshold((masses[key]/bottom["masses"]["b"])**2,
                                              4, bottom["alphaS"])
    below = {"Q": bottom["Q"], "alphaS": match["alphaS"], "alpha": match["alpha"],
             "masses": {k: m*factors[k] for k, m in masses.items()}}
    at_two = run(below, 2.0, nd=2)
    charm = self_scale(below, "c", 2)
    return {**{k: at_two["masses"][k] for k in ("u", "d", "s")},
            "c": charm["masses"]["c"], "b": bottom["masses"]["b"]}


def replay(smdr):
    smdr = Path(smdr).resolve()
    pin = read("smdr_source_pin.json")
    for name, expected in pin["files_sha256"].items():
        if digest(smdr/name) != expected:
            raise ValueError(f"SMDR source pin mismatch: {name}")
    inputs = read("inputs.json")
    packet, g = inputs["boundary_packet"], inputs["couplings_Q0"]
    q0, v = packet["Q0_GeV"], packet["v_GeV"]
    order = ("yu", "yc", "yt", "yd", "ys", "yb", "ye", "ymu", "ytau", "lambda")
    native = [q0, v, g["g1_GUT"]*math.sqrt(3/5), g["g2"], g["g3"], *[g[k] for k in order]]
    with tempfile.TemporaryDirectory(prefix="oph-spectrum-") as directory:
        work = Path(directory)
        heavy = compile_driver(smdr, work, "heavy", ["smdr_flavor_driver.c"])
        light = compile_driver(smdr, work, "light", ["light_driver.c", "corrected_decouple_thZW.c"])
        wz = compile_driver(smdr, work, "wz", ["light_fermion_wz.c"])
        precise = compile_driver(smdr, work, "rg", ["rg_precise.c"])
        matched = invoke(heavy, "match", *native)["rows"][2]
        shifts = invoke(wz, "light", *native)["delta_complex_pole_squared_GeV2"]
        outputs = {"H": matched["higgs"]["mass_GeV"], "t": matched["top_QCD4_EW2"]["mass_GeV"]}
        for key in ("W", "Z"):
            row = matched[key]
            pole_squared = complex(row["mass_GeV"], -row["width_GeV"]/2)**2+complex(*shifts[key])
            outputs[key] = math.sqrt(pole_squared.real+pole_squared.imag**2/pole_squared.real)
        decoupled = invoke(light, "fixed", *native)
        outputs.update(low_energy(light, decoupled))
        errors = {k: outputs[k]/row["value_GeV"]-1 for k, row in inputs["imported_readouts"].items()}
        if max(map(abs, errors.values())) > 1e-9:
            raise ArithmeticError(f"Frozen mass replay mismatch: {errors}")
        strong = inputs["strong_coupling"]
        alpha = invoke(light, "run", q0, strong["MZ_comparison_scale_GeV"], 2, 3, 3, 5,
                       decoupled["alphaS"], decoupled["alpha"])
        if not math.isclose(alpha["alphaS"], strong["alpha_s_5_at_PDG_MZ"], rel_tol=1e-11):
            raise ArithmeticError("Strong-coupling scheme conversion failed")
        names = {"Q": q0, "v": v, "gp": native[2], "g": g["g2"], "g3": g["g3"],
                 **{k: g[k] for k in order}, "m2": -g["lambda"]*v*v, "Lambda": 0}
        dat = work/"input.dat"
        dat.write_text("\n".join(f"SMDR_{k}_in = {x!r};" for k, x in names.items())+"\n")
        critical = {str(n): invoke(precise, packet["UV_GeV"], n, dat.name, cwd=work) for n in (2, 3)}
        if abs(critical["2"]["lambda"]) > 3e-8 or abs(critical["2"]["beta_lambda"]) > 1e-9:
            raise ArithmeticError("The two criticality conditions fail the frozen replay tolerance")
        controls = []
        for control in inputs["substitution_controls"]:
            c = control["inputs"]
            args = [q0, c["v"], c["gY"], c["g2"], c["g3"],
                    *[c[k] for k in order[:-1]], c["lam"]]
            row = invoke(heavy, "match", *args)["rows"][2]
            got = {"W_BW_GeV": row["W"]["BW_mass_GeV"], "Z_BW_GeV": row["Z"]["BW_mass_GeV"],
                   "Higgs_pole_GeV": row["higgs"]["mass_GeV"], "top_QCD4_EW2_pole_GeV": row["top_QCD4_EW2"]["mass_GeV"]}
            error = max(abs(val/control["outputs"][key]-1) for key, val in got.items())
            if error > 1e-9:
                raise ArithmeticError("Reference-substitution replay failed")
            controls.append({"substituted_reference_inputs": control["substituted_from_SM_reference"],
                             "outputs": got, "max_relative_replay_error": error})
    return {"schema": "oph.particle_spectrum_audit.smdr_replay.v1", "all_checks_pass": True,
            "outputs_GeV": outputs, "relative_mass_replay_errors": errors,
            "alpha_s_5_at_comparison_MZ": alpha["alphaS"], "criticality_boundary": critical,
            "reference_substitution_controls": controls,
            "inputs_sha256": digest(HERE/"inputs.json"),
            "sources_sha256": {p.name: digest(p) for p in sorted(HERE.iterdir())
                               if p.suffix in (".py", ".c", ".h")},
            "library_sha256": {name: digest(smdr/name) for name in ("libsmdr.a", "lib3vil.a", "libtsil.a")},
            "scope": "Replay of a stipulated frozen boundary and retained perturbative approximations; no physical identity or complete theory-error bound."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smdr-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE/"smdr_replay.json")
    args = parser.parse_args()
    result = replay(args.smdr_dir)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print("Frozen nine-particle multiloop outputs, criticality and substitution controls verified")


if __name__ == "__main__":
    main()
