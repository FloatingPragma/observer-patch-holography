"""Audit an explicitly calibrated SM parameter point with external SMDR.

The ordinary --check verifies custody and recomputes comparisons from retained
external outputs. --smdr-root additionally compiles the C adapter and repeats
the forward calculations, without network access. The external perturbative
SM calculation is neither native observer evolution nor a mass derivation.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from tools import strict_json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "calibrated_matching_results.json"
SMDR_COMMIT = "5c4dd92c918cadb7f91ad3795694f8398b5c287a"
PACKET = HERE / "calibrated_matter_parameters.json"
ADAPTER = HERE / "calibrated_matching.c"

# Explicit calibration targets from the pinned SMDR reference, with the
# measurement conventions and uncertainties in its src/smdr_pdg.h. They are
# not an all-data likelihood. W is excluded from the inverse calibration.
TARGETS = {
    "SMDR_Mt_pole": ("172.4", "0.7", "top pole extraction"),
    "SMDR_Mh_pole": ("125.20", "0.11", "Higgs pole"),
    "SMDR_MZ_PDG": ("91.1876", "0.0021", "LEP-only Breit-Wigner Z"),
    "SMDR_mbmb": ("4.183", "0.004", "MSbar mb(mb), five-flavour QCD+QED"),
    "SMDR_mcmc": ("1.2730", "0.0028", "MSbar mc(mc), four-quark/two-lepton QCD+QED"),
    "SMDR_ms_2GeV": ("0.0935", "0.0005", "MSbar strange, 2 GeV"),
    "SMDR_mu_2GeV": ("0.00216", "0.00004", "MSbar up, 2 GeV"),
    "SMDR_md_2GeV": ("0.00470", "0.00004", "MSbar down, 2 GeV"),
    "SMDR_Mtau_pole": ("1.77693", "0.00009", "tau pole"),
    "SMDR_Mmuon_pole": ("0.1056583755", "0.0000000023", "muon pole"),
    "SMDR_Melectron_pole": ("0.00051099895000", "0.00000000000015", "electron pole"),
    "SMDR_alpha": ("1/137.035999178", "0.00000000000043", "Sommerfeld alpha"),
    "SMDR_alphaS_5_MZ": ("0.1180", "0.0009", "five-flavour alphaS at MZ"),
    "SMDR_GFermi": ("0.000011663782", "0.000000000006", "SMDR Fermi-constant convention"),
    "SMDR_MW_PDG": ("80.3692", "0.0133", "Breit-Wigner W reference excluding CDF; not fitted"),
}
CHARGED = ("SMDR_Mt_pole", "SMDR_mbmb", "SMDR_mcmc", "SMDR_ms_2GeV",
           "SMDR_mu_2GeV", "SMDR_md_2GeV", "SMDR_Mtau_pole",
           "SMDR_Mmuon_pole", "SMDR_Melectron_pole")
OUTPUT_NAMES = set(TARGETS) | {"SMDR_MZ_pole", "SMDR_MW_pole", "SMDR_Gammat_pole",
    "SMDR_Gammah_pole", "SMDR_GammaZ_pole", "SMDR_GammaW_pole", "sizeof_long_double"}
# Independent shipped SMDR reference at Q=172.4 GeV; not the calibrated
# common-MZ point. Comparing both detects transport and adapter mistakes.
REFERENCE = dict(zip(
    ("Q", "g3", "g", "gp", "v", "lambda", "yt", "yb", "yc", "ys", "yd", "yu", "ytau", "ymu", "ye", "Delta_alpha_had_5_MZ"),
    ("172.4", "1.1635013021305831", "0.6475934615902679", "0.3585947982153167",
     "246.60029586365829", "0.1262969232875777", "0.9309715465595599",
     "0.0155080080991419", "0.003400275709585654", "0.000292915817833262",
     "0.00001472409137072013", "0.00000673214356538510", "0.0099949323705356",
     "0.00058838465029327", "0.00000279299434070567", "0.02783")))


def _fraction(value):
    if type(value) is not str:
        raise ValueError("finite decimal or decimal-ratio string required")
    pieces = value.split("/")
    if len(pieces) > 2:
        raise ValueError("one optional rational division required")
    result = F(pieces[0])
    if len(pieces) == 2:
        result /= F(pieces[1])
    return result


def repository_key(path, root):
    """Repository custody keys use forward slashes on every host."""
    return path.relative_to(root).as_posix()


def source_bindings():
    return {repository_key(p, ROOT): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__), ADAPTER, PACKET)}


def smdr_inputs(packet=None, *, rounded=False):
    packet = strict_json.load(PACKET) if packet is None else packet
    base = packet["rounded_table_control"] if rounded else packet
    values, scale = base["matched_SM"], base["common_scale"]
    y = values["yukawa"]
    with localcontext() as ctx:
        ctx.prec = 55
        result = {"Q": scale["reference_MZ_GeV"], "v": values["v_SM_GeV"],
            "lambda": str(Decimal(values["lambda_paper"])/4),
            "g": values["g2"], "g3": values["g3"],
            "gp": str(Decimal(values["g1_GUT"])*(Decimal(3)/5).sqrt()),
            "Delta_alpha_had_5_MZ": "0.02783"}
    for channel, names in (("up", ("yu", "yc", "yt")), ("down", ("yd", "ys", "yb")),
                           ("charged_lepton", ("ye", "ymu", "ytau"))):
        if len(y[channel]) != 3:
            raise ValueError("three family coordinates required")
        result.update(zip(names, y[channel]))
    if set(result) != set(REFERENCE) or any(_fraction(v) <= 0 for v in result.values()):
        raise ValueError("complete positive SMDR boundary required")
    return result


def input_text(values):
    if set(values) != set(REFERENCE):
        raise ValueError("complete SMDR input coordinates required")
    for value in values.values():
        if _fraction(value) <= 0:
            raise ValueError("positive SMDR input required")
    return "".join(f"SMDR_{key}_in = {values[key]};\n" for key in sorted(values))


def parse_output(text):
    rows = [x.removeprefix("OPH_MATCHING_JSON=") for x in text.splitlines()
            if x.startswith("OPH_MATCHING_JSON=")]
    if len(rows) != 1:
        raise ValueError("one complete adapter output required")
    result = strict_json.loads(rows[0])
    if set(result) != OUTPUT_NAMES:
        raise ValueError("incomplete or extra forward output")
    if type(result["sizeof_long_double"]) is not int or result["sizeof_long_double"] not in (8, 12, 16):
        raise ValueError("explicit floating-point representation required")
    if any(_fraction(v) <= 0 for k, v in result.items() if k != "sizeof_long_double"):
        raise ValueError("positive finite matched outputs required")
    return result


def comparison(outputs):
    # Validate stored data with exactly the same strict parser as fresh output.
    outputs = parse_output("OPH_MATCHING_JSON="+json.dumps(outputs, allow_nan=False))
    rows = {}
    for name, (target, sigma, convention) in TARGETS.items():
        residual = (_fraction(outputs[name])-_fraction(target))/_fraction(sigma)
        rows[name] = {"target": target, "quoted_standard_uncertainty": sigma,
            "convention": convention, "used_for_calibration": name != "SMDR_MW_PDG",
            "matched_value": outputs[name], "residual_over_quoted_uncertainty": float(residual),
            "absolute_residual_le_one_quoted_uncertainty": abs(residual) <= 1}
    return rows


def build_result(outputs, external_source):
    if set(outputs) != {"calibrated", "rounded_table", "shipped_reference"}:
        raise ValueError("full-precision, rounded and reference controls required")
    rows = {name: comparison(value) for name, value in outputs.items()}
    inputs = {"calibrated": smdr_inputs(), "rounded_table": smdr_inputs(rounded=True),
              "shipped_reference": REFERENCE}
    return {
        "schema": "oph.calibrated_matching.v1",
        "role": "external_SM_calibration_sanity_check_not_mass_derivation",
        "source_bindings": source_bindings(), "external_source": external_source,
        "inputs": inputs, "external_outputs": outputs, "comparisons": rows,
        "summary": {
            "charged_mass_inputs_reproduced_within_quoted_one_sigma": all(
                rows["calibrated"][k]["absolute_residual_le_one_quoted_uncertainty"] for k in CHARGED),
            "all_fitted_inputs_reproduced_within_quoted_one_sigma": all(
                v["absolute_residual_le_one_quoted_uncertainty"] for v in rows["calibrated"].values()
                if v["used_for_calibration"]),
            "rounded_muon_and_electron_miss_their_precision": all(
                not rows["rounded_table"][k]["absolute_residual_le_one_quoted_uncertainty"]
                for k in ("SMDR_Mmuon_pole", "SMDR_Melectron_pole")),
            "W_reference_was_not_fitted": True,
            "W_within_one_quoted_sigma": rows["calibrated"]["SMDR_MW_PDG"]["absolute_residual_le_one_quoted_uncertainty"],
        },
        "scope": {
            "public_SMDR_source_commit": SMDR_COMMIT,
            "fixed_point": "All charged Yukawas, three gauge couplings, Higgs quartic and v at 91.1876 GeV; one common action point",
            "fit": "SMDR calc_fit -i ReferenceModel.dat -Q 91.1876 -e 1e-14 -o CalibratedFull.dat",
            "fit_stopping": "12-iteration cap reached; requested 1e-14 fractional tolerance not attained. Last printed aggregate fractional residual 1.51230e-10; forward observables checked separately.",
            "cache_control": "Adapter transports the supplied point to 200 GeV before forward matching to avoid SMDR1.3's zero-length RG return with stale derived coupling powers; no refit in this step",
            "additional_empirical_input": "Delta_hadronic^(5) alpha(MZ)=0.02783; pinned SMDR threshold scales and PDG reference conventions",
            "error_meaning": "Observed floating-point calibration residuals only. Quoted marginal uncertainties are not hard error bounds, a joint likelihood or a perturbative-theory uncertainty budget.",
            "independent_data": "None of the fitted masses is held out. The W reference comparison is retained and is not within one quoted standard uncertainty.",
            "native_attachment": "No native interacting self-energy, confinement, decay widths, quantum gauge dynamics or full Standard Model observer history is executed by this external calculation.",
            "neutrinos": "This SMDR calculation contains no neutrino mass or Weinberg matching.",
            "flavour_approximation": "SMDR neglects CKM mixing in this precision matching and diagonal-Yukawa running (arXiv:1907.02500, footnote following Eq.1.4); the separate calibrated action retains its supplied CKM matrix. The omitted contribution is not assigned a certified error here.",
            "offline_check": "Recomputes comparisons and checks custody of the retained external output; does not rerun SMDR.",
            "external_replay": "python -m model.common_source.calibrated_matching --check --smdr-root /path/to/built/SMDR",
        },
    }


def external_identity(directory):
    directory = Path(directory).resolve()
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=directory, text=True).strip()
    if commit != SMDR_COMMIT:
        raise ValueError("the declared SMDR source revision is required")
    subprocess.run(["git", "diff", "--exit-code", "HEAD", "--"], cwd=directory,
                   check=True, stdout=subprocess.DEVNULL)
    names = subprocess.check_output(["git", "ls-files", "-z"], cwd=directory).split(b"\0")
    digest = hashlib.sha256()
    for raw in sorted(n for n in names if n):
        data = (directory/raw.decode()).read_bytes()
        digest.update(raw+b"\0"+hashlib.sha256(data).digest())
    return {"repository": "https://github.com/davidgrobertson/SMDR", "commit": commit,
            "tracked_tree_sha256": digest.hexdigest(), "version": "1.3"}


def execute(directory, compiler="cc"):
    directory = Path(directory).resolve()
    identity = external_identity(directory)
    values = {"calibrated": smdr_inputs(), "rounded_table": smdr_inputs(rounded=True),
              "shipped_reference": REFERENCE}
    with tempfile.TemporaryDirectory(prefix="oph-sm-matching-") as temporary:
        temporary = Path(temporary)
        executable = temporary/"matching"
        subprocess.run([compiler, "-O2", "-I"+str(directory), str(ADAPTER), "-L"+str(directory),
                        "-lsmdr", "-l3vil", "-ltsil", "-lm", "-o", str(executable)],
                       check=True, timeout=120)
        outputs = {}
        for name, data in values.items():
            path = temporary/(name+".dat")
            path.write_text(input_text(data))
            run = subprocess.run([str(executable), str(path)], cwd=directory, text=True,
                                 capture_output=True, check=True, timeout=120)
            outputs[name] = parse_output(run.stdout)
    return build_result(outputs, identity)


def check(result):
    if result["external_source"]["commit"] != SMDR_COMMIT:
        raise ValueError("wrong retained external revision")
    rebuilt = build_result(result["external_outputs"], result["external_source"])
    if json.dumps(result, sort_keys=True) != json.dumps(rebuilt, sort_keys=True):
        raise ValueError("matching receipt, inputs, comparisons or source custody changed")
    return rebuilt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path, nargs="?", const=OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--smdr-root", type=Path)
    args = parser.parse_args()
    if args.check:
        retained = check(strict_json.load(args.check))
        if args.smdr_root:
            fresh = execute(args.smdr_root)
            if fresh["external_source"] != retained["external_source"]:
                raise ValueError("external source tree changed")
            # Numerical replay tolerance is explicit and distinct from each
            # measured uncertainty. Fresh classifications must also agree.
            for name, outputs in fresh["external_outputs"].items():
                for key, value in outputs.items():
                    if key != "sizeof_long_double":
                        old = _fraction(retained["external_outputs"][name][key])
                        if abs(_fraction(value)-old) > abs(old)*F(2, 10**8):
                            raise ArithmeticError("external forward replay changed: "+name+"/"+key)
            if fresh["summary"] != retained["summary"]:
                raise ArithmeticError("external forward comparison classification changed")
            print("SMDR forward replay and calibration comparisons passed")
        else:
            print("retained SMDR output custody and exact comparison arithmetic passed (offline)")
    else:
        if args.smdr_root is None:
            parser.error("writing results requires --smdr-root; retained data is never fabricated")
        result = execute(args.smdr_root)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
        print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
