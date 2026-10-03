#!/usr/bin/env python3
"""Independent checks of the conditional, retrospective six-mass replay.

No producer or D10 implementation is imported. The gauge packet is checked
against heat sums and defining equations. The spectrum is checked through
its geometric means and two log-difference invariants in each sector,
rather than by reusing the producer's centered exponential evaluation.
These are numerical consistency checks, not a physical attachment theorem.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import verify_source_w5_response_constraints as source_w5

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_RECEIPT = HERE.parent / "runs/flavor/conditional_quark_mass_replay.json"
SPEC = HERE / "conditional_quark_mass_spec.json"
PARTICLES = ("u", "c", "t", "d", "s", "b")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def number(value: Any, label: str) -> float:
    require(isinstance(value, str) and len(value) < 200,
            f"{label}: expected a finite decimal string")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{label}: invalid decimal") from exc
    require(math.isfinite(result), f"{label}: nonfinite decimal")
    return result


def close(actual: float, expected: float, label: str,
          *, rtol: float = 2e-11, atol: float = 2e-13) -> None:
    require(math.isfinite(actual) and math.isfinite(expected)
            and abs(actual - expected) <= atol + rtol * abs(expected),
            f"{label}: defining equation failed")


def heat_log_dimension(coupling: float, rank: int) -> float:
    """Separate complete-enough finite sums; no producer heat routine."""
    require(coupling > 0, "heat coupling must be positive")
    t = 4 * math.pi**2 * coupling
    pairs = []
    if rank == 2:
        for dimension in range(1, 82):
            casimir = (dimension * dimension - 1) / 4
            pairs.append((dimension, casimir))
    else:
        require(rank == 3, "unknown heat rank")
        for p in range(61):
            for q in range(61):
                dimension = (p + 1) * (q + 1) * (p + q + 2) // 2
                casimir = ((p + q) ** 2 - p * q + 3 * (p + q)) / 3
                pairs.append((dimension, casimir))
    weights = [dimension * math.exp(-t * casimir)
               for dimension, casimir in pairs]
    return math.fsum(w * math.log(pair[0]) for w, pair in zip(weights, pairs)) / math.fsum(weights)


def check_source_packet(packet: dict[str, Any], p: float) -> dict[str, float]:
    keys = ("P", "alpha_U", "v_over_E_star", "mz_over_E_star",
            "mu_U_over_E_star", "alpha1_mz", "alpha2_mz", "alpha3_mz")
    values = {key: number(packet[key], f"source_packet.{key}") for key in keys}
    close(values["P"], p, "source P", rtol=2e-14, atol=0)
    alpha = values["alpha_U"]
    require(0.02 < alpha < 0.08, "D10 alpha_U outside its declared search domain")
    scale = math.exp(-2 * math.pi) * p ** (1 / 6)
    v = math.exp(-math.pi / (2 * alpha)) / math.sqrt(p)
    close(values["mu_U_over_E_star"], scale, "D10 unification scale", atol=0)
    close(values["v_over_E_star"], v, "D10 transmutation", atol=0)
    mz = values["mz_over_E_star"]
    require(mz > 0 and v > 0 and mz < scale, "D10 scale ordering")
    for index, slope in enumerate((33 / 5, 1, -3), 1):
        expected = 1 / (1 / alpha + slope * math.log(scale / mz) / (2 * math.pi))
        close(values[f"alpha{index}_mz"], expected, f"D10 gauge {index}")
    ay = 3 * values["alpha1_mz"] / 5
    close(mz, v * math.sqrt(math.pi * (ay + values["alpha2_mz"])),
          "D10 Z fixed point", atol=0)
    heat = heat_log_dimension(values["alpha2_mz"], 2) + heat_log_dimension(values["alpha3_mz"], 3)
    close(heat, p / 4, "D10 heat closure", rtol=2e-11, atol=2e-13)
    return values


def expected_coefficients(p: float, alpha: float, full: bool) -> dict[str, float]:
    """The previously frozen rational law, independently transcribed."""
    w = math.pi * alpha
    tau = p / 4 - w / 5
    r = math.exp(-3 * tau)
    rho = 3 / (2 + r)
    x = -math.tanh(3 * tau / 2)
    mu_u = 3 * p + 79 * w / 15
    mu_d = 2 * p + 4 * w / 15
    a_u = 3 * p + 29 * w / 5
    a_d = 2 * p + 33 * w / 32
    delta = 0.0
    if full:
        mu_u += w * w / 29
        mu_d -= w * w / 432
        a_u += w * w / 22
        a_d += w * w / 420
        delta = p / 1008 + w * w * (1 / 432 - 1 / 1584)
    # The ray slopes are simplified before evaluation; cancellation differs
    # from the producer's unsimplified rational form.
    bu_ray = a_u * ((rho - 1) - x * (rho + 1)) / ((rho + 1) * (x * x - 1))
    bd_ray = a_d * ((1 - rho) - x * (rho + 1)) / ((rho + 1) * (x * x - 1))
    qu, qd = -w * rho / 10, -w / 4
    aa = 1 / (2 * (1 + rho - x * x))
    bb = (1 + rho) / (2 * ((1 + rho) - (2 + rho) * x * x))
    return dict(w=w, tau=tau, r=r, rho=rho, x=x, mean_u=mu_u,
                mean_d=mu_d, a_u=a_u, a_d=a_d, b_u_ray=bu_ray,
                b_d_ray=bd_ray, q_u=qu, q_d=qd, b_u=bu_ray + qu,
                b_d=bd_ray + qd, A=aa, B=bb, delta_g=delta,
                g_ch_over_v=2 * math.exp(-2 * math.pi + delta))


def check_spectrum(row: dict[str, Any], p: float, alpha: float,
                   v: float, energy: float, full: bool) -> None:
    tag = "full_rscc" if full else "lower_order_control"
    require(set(row) == {"coefficients", "basis", "mass_over_v", "mass_over_E_star",
                         "native_mass_chart_GeV", "light_ratios"},
            f"{tag}: unknown scientific field or hidden input")
    require(set(row["basis"]) == {"linear", "quadratic"}, f"{tag}: basis keys")
    require(set(row["light_ratios"]) == {"u_over_d", "s_over_mean_ud"},
            f"{tag}: light ratio keys")
    expected = expected_coefficients(p, alpha, full)
    require(set(row["coefficients"]) == set(expected), f"{tag}: coefficient keys")
    for key, value in expected.items():
        close(number(row["coefficients"][key], key), value, f"{tag}.{key}")
    x = expected["x"]
    for key, values in {"linear": [-1 - x / 3, 2 * x / 3, 1 - x / 3],
                        "quadratic": [(1 - x * x) / 3,
                                      2 * (x * x - 1) / 3,
                                      (1 - x * x) / 3]}.items():
        require(len(row["basis"][key]) == 3, f"{tag}: basis length")
        for actual, value in zip(row["basis"][key], values):
            close(number(actual, key), value, f"{tag}.{key}")
    maps = ("mass_over_v", "mass_over_E_star", "native_mass_chart_GeV")
    for name in maps:
        require(set(row[name]) == set(PARTICLES), f"{tag}: spectrum particle keys")
    masses = {name: number(value, f"{tag}.{name}")
              for name, value in row["mass_over_v"].items()}
    require(all(m > 0 for m in masses.values()), f"{tag}: nonpositive mass")
    for name in PARTICLES:
        over_e = number(row["mass_over_E_star"][name], name)
        gev = number(row["native_mass_chart_GeV"][name], name)
        close(over_e, masses[name] * v, f"{tag}.{name} normalization", atol=0)
        close(gev, over_e * energy, f"{tag}.{name} unit conversion", atol=0)
    sigma = (expected["mean_u"] + expected["mean_d"]) / 2
    eta = (expected["mean_u"] - expected["mean_d"]) / 2
    for sector, names, sign in (("u", ("u", "c", "t"), 1),
                                ("d", ("d", "s", "b"), -1)):
        logs = [math.log(masses[name]) for name in names]
        a, b = expected[f"a_{sector}"], expected[f"b_{sector}"]
        center = math.log(expected["g_ch_over_v"]) - expected["A"] * sigma + sign * expected["B"] * eta
        close(math.fsum(logs) / 3, center, f"{tag}.{sector} geometric mean")
        close(logs[2] - logs[0], 2 * a, f"{tag}.{sector} endpoint log ratio")
        close(logs[1] - (logs[0] + logs[2]) / 2, a * x + b * (x * x - 1),
              f"{tag}.{sector} centered log curvature")
    light = row["light_ratios"]
    close(number(light["u_over_d"], "u_over_d"), masses["u"] / masses["d"], f"{tag} u/d")
    close(number(light["s_over_mean_ud"], "s_over_mean_ud"),
          2 * masses["s"] / (masses["u"] + masses["d"]), f"{tag} s/mean")


GAUGE_METHODS = {
    "__init__", "_build_su2_terms", "_build_su3_terms",
    "_derive_stage5_integer_vectors", "alpha_em_from_alpha1_alpha2",
    "alpha_run_1loop", "build_d10_from_p", "e_cell", "ellbar_su2",
    "ellbar_su3", "mz_tree_from_v_and_couplings", "outer_p_from_alpha",
    "p_from_inverse_alpha", "pixel_residual", "run_alphas_from_unification",
    "solve_alpha_u_from_p", "solve_mz_fixed_point_tree",
    "unification_scale_gev", "v_from_transmutation",
}


def check_gauge_call_graph(producer_text: str, gauge_text: str) -> None:
    """Audit the optional source refresh's actual inherited method reachability.

    The historical module contains a flavor closure elsewhere. Presence of
    that code is distinct from consuming it: this verifies the gauge-only
    subclass and the reachable self-method graph instead of treating the
    whole historical module as a numerical quark ancestor.
    """
    tree = ast.parse(producer_text)
    subclasses = [n for n in ast.walk(tree) if isinstance(n, ast.ClassDef)
                  and any((isinstance(b, ast.Name) and b.id == "PaperMathContext")
                          or (isinstance(b, ast.Attribute) and b.attr == "PaperMathContext")
                          for b in n.bases)]
    require(len(subclasses) == 1, "exactly one gauge-only context required")
    override = {n.name: n for n in subclasses[0].body if isinstance(n, ast.FunctionDef)}
    require("_derive_stage5_integer_vectors" in override,
            "Stage-5 constructor override missing")
    body = override["_derive_stage5_integer_vectors"].body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    require(len(body) == 1 and isinstance(body[0], ast.Return)
            and isinstance(body[0].value, ast.Dict) and not body[0].value.keys,
            "Stage-5 constructor must be disabled, not evaluated")
    require("diagonal_quark_masses" in override,
            "quark mass guard missing")
    require(any(isinstance(n, ast.Raise) for n in ast.walk(override["diagonal_quark_masses"])),
            "quark mass guard must raise")
    gauge_tree = ast.parse(gauge_text)
    classes = [n for n in gauge_tree.body if isinstance(n, ast.ClassDef)
               and n.name == "PaperMathContext"]
    require(len(classes) == 1, "gauge context missing")
    methods = {n.name: n for n in classes[0].body if isinstance(n, ast.FunctionDef)}
    seen: set[str] = set()
    pending = ["__init__", "p_from_inverse_alpha", "build_d10_from_p"]
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        require(name in GAUGE_METHODS, f"nongauge method in refresh ancestry: {name}")
        if name == "_derive_stage5_integer_vectors":
            continue
        require(name in methods, f"missing gauge method: {name}")
        for node in ast.walk(methods[name]):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
                require(node.func.attr not in {"read_text", "read_bytes", "open"},
                        "external input in mathematical gauge ancestry")
                if (isinstance(node.func.value, ast.Name)
                        and node.func.value.id == "self" and node.func.attr in methods):
                    pending.append(node.func.attr)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                require(node.func.id not in {"open", "eval", "exec", "__import__"},
                        "dynamic input in mathematical gauge ancestry")
    require(seen == GAUGE_METHODS, "reviewed gauge ancestry changed")


COUNTS = {"N_g": 3, "N_c": 3, "d_std_S3": 2, "order_S3": 6,
          "beta_EW": 4, "oriented_repair_slots": 24, "orientation_dimension": 2}
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
MENU = {
    "mean_u_linear": "79/15", "mean_d_linear": "4/15",
    "a_u_linear": "29/5", "a_d_linear": "33/32",
    "mean_u_quadratic": "1/29", "mean_d_quadratic": "-1/432",
    "a_u_quadratic": "1/22", "a_d_quadratic": "1/420",
    "delta_g_P": "1/1008", "delta_g_quadratic_color": "1/432",
    "delta_g_quadratic_relabel": "-1/1584", "q_u": "-rho/10", "q_d": "-1/4",
    "module_dimensions": [29, 432, 22, 32, 840, 1008, 432, 1584],
}
CHARTS = {
    **{q: {"scheme": "MSbar", "scale": "2 GeV"} for q in ("u", "d", "s")},
    **{q: {"scheme": "MSbar", "scale": "self_scale_mu_equals_running_mass"}
       for q in ("c", "b")},
    "t": {"scheme": "declared_pole_coordinate", "scale": "not_a_running_scale",
          "boundary": "pole_attachment_is_assumed_not_derived; experimental_extraction_matching_is_separate"},
}
CALIBRATIONS = {
    "alpha_inverse": {"value": "137.035999177", "role": "external_measured_calibration"},
    "E_star_GeV": {"value": "1.220890e19", "role": "external_dimensionful_unit_calibration"},
}
GAUGE_SOLVER = {"precision": 50, "su2_cutoff": 80, "su3_cutoff": 60,
                "alpha_U_domain": ["0.02", "0.08"]}
SCHEMA = "oph.conditional_quark_mass_replay.v1"
CLAIM = "conditional_retrospective_RSCC_native_mass_chart_replay"
STATUS = "conditional_forward_replay_not_source_selected_or_independent_prediction"
HISTORY = [
    "code/particles/flavor/quark_rscc_completion_candidate.py",
    "code/particles/flavor/verify_quark_rscc_module_arithmetic.py",
    "code/particles/runs/flavor/quark_rscc_completion_candidate_audit.json",
]
NUMERIC_INPUTS = ["code/P_derivation/codata_2022_alpha_fixture.json",
                  "code/particles/flavor/conditional_quark_mass_spec.json",
                  "code/particles/runs/flavor/source_w5_response_constraints.json",
                  "code/particles/flavor/source_w5_response_constraints_spec.json",
                  "code/a5_closure/manifests/echosahedral_federation_reference.json",
                  "code/a5_closure/manifests/a3_scheduler_kernel_reference.json",
                  "code/a5_closure/manifests/record_counting_mechanism_reference.json"]
RESPONSE_CODE_PATHS = ["code/particles/flavor/source_w5_response_constraints.py",
                       "code/particles/flavor/verify_source_w5_response_constraints.py",
                       "code/a5_closure/source_repair_generator_certificate.py"]
RESPONSE_CONTRACT = {'physical_boundary': 'The dimensions and normalized native projector traces are certified. Choosing the face-axis background, Casimir-line common-exposure readout, scalar budget assignment and its physical quark attachment remains conditional. Composite modules, signed responses and Gaussian/contraction laws are not selected by these coefficients.', 'projection': {'centered_fraction': '4/5', 'common_exposure': '4/15', 'family_dimension': 5, 'line_complement_dimension': 4, 'scalar_fraction': '1/5'}, 'receipt': 'code/particles/runs/flavor/source_w5_response_constraints.json'}

PRODUCER = "code/particles/flavor/conditional_quark_mass_replay.py"
GAUGE = "code/P_derivation/paper_math.py"
STATEMENT_HASHES = {
    "hypotheses": "51a8a0dc8c959237c76296f2113b04ecbca50beef7229eaa88f31b3bdb56abac",
    "equations": "2df54fbe8650d0839023afbdddecf5172cc6134820c819ba6e20608dbf04543f",
}


def content_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def check_producer_inputs(source: str) -> None:
    """Constrain executable imports and reads to the reviewed numerical path.

    Historical source files are hashed, never imported or parsed as inputs.
    Together with exact pins, fixed inputs and the independent numerical
    equations this prevents a comparison fixture becoming a hidden input.
    """
    tree = ast.parse(source)
    parents = {child: node for node in ast.walk(tree) for child in ast.iter_child_nodes(node)}

    def owner(node: ast.AST) -> str:
        while node in parents:
            node = parents[node]
            if isinstance(node, ast.FunctionDef):
                return node.name
        return "<module>"

    allowed_imports = {"__future__", "argparse", "copy", "hashlib", "importlib.util",
                       "json", "pathlib", "sys", "typing", "mpmath", "fractions",
                       "verify_source_w5_response_constraints"}
    read_text_receivers = set()
    read_count = byte_count = parse_count = 0
    dynamic_loader_calls = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            require(all(alias.name in allowed_imports for alias in node.names),
                    "unreviewed producer import")
        if isinstance(node, ast.ImportFrom):
            require(node.module in allowed_imports, "unreviewed producer import")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            require(node.func.id not in {"open", "eval", "exec", "__import__"},
                    "dynamic producer input")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            name = node.func.attr
            if name == "read_text":
                receiver = ast.unparse(node.func.value)
                require(receiver in {"path", "ALPHA_PATH", "args.out"},
                        "unreviewed numerical input read")
                allowed_owners = {"main"} if receiver == "args.out" else ({"load_spec", "load_source_response"} if receiver == "path" else {"load_spec"})
                require(owner(node) in allowed_owners,
                        "historical or hidden data parsed outside input loader")
                read_text_receivers.add(receiver)
                read_count += 1
            if name == "read_bytes":
                require(ast.unparse(node.func.value) == "path" and owner(node) == "sha",
                        "historical source bytes used outside hash function")
                byte_count += 1
            if name == "loads":
                require(isinstance(node.func.value, ast.Name) and node.func.value.id == "json"
                        and owner(node) in {"load_spec", "load_source_response"} and len(node.args) == 1
                        and ast.unparse(node.args[0]) in {"path.read_text()", "ALPHA_PATH.read_text()"},
                        "historical or observed target parsed as numerical input")
                parse_count += 1
            if name == "spec_from_file_location":
                require(len(node.args) == 2 and isinstance(node.args[1], ast.Name)
                        and node.args[1].id == "GAUGE_PATH", "unreviewed dynamic source import")
                dynamic_loader_calls.append(node)
    require(read_text_receivers == {"path", "ALPHA_PATH", "args.out"},
            "reviewed producer read surface changed")
    require((read_count, byte_count, parse_count) == (4, 1, 3),
            "reviewed producer input surface changed")
    require(len(dynamic_loader_calls) == 1, "gauge loader ancestry changed")


def verify_payload(payload: dict[str, Any], *, repo: Path = REPO) -> dict[str, Any]:
    require(isinstance(payload, dict), "receipt must be an object")
    require(set(payload) == {"schema", "claim_class", "status", "inputs", "source_packet",
                             "model", "full_rscc", "lower_order_control", "hypotheses",
                             "guards", "source_pins"},
            "unknown receipt field or hidden input")
    for key, value in {"schema": SCHEMA, "claim_class": CLAIM, "status": STATUS}.items():
        require(payload[key] == value, f"noncanonical {key}")
    require(payload["guards"] == GUARDS
            and all(type(value) is bool for value in payload["guards"].values()),
            "scientific claim or provenance guard changed")
    require(content_hash(payload["hypotheses"]) == STATEMENT_HASHES["hypotheses"],
            "physical hypotheses removed or promoted")
    inputs = payload["inputs"]
    require(set(inputs) == {"alpha_inverse", "E_star_GeV", "P", "structural_counts"},
            "hidden numeric input")
    for key, calibration in CALIBRATIONS.items():
        require(inputs[key] == calibration, f"external calibration changed: {key}")
    require(inputs["structural_counts"] == COUNTS, "structural incidence changed")
    require(set(inputs["P"]) == {"value", "definition"}
            and inputs["P"]["definition"] == "phi+sqrt(pi)/alpha_inverse",
            "independent P definition changed")
    p = (1 + math.sqrt(5)) / 2 + math.sqrt(math.pi) / number(inputs["alpha_inverse"]["value"], "alpha inverse")
    close(number(inputs["P"]["value"], "P"), p, "outer P identity", rtol=2e-14, atol=0)
    packet = payload["source_packet"]
    require(set(packet) == {"P", "alpha_U", "v_over_E_star", "mz_over_E_star",
                            "mu_U_over_E_star", "alpha1_mz", "alpha2_mz", "alpha3_mz"},
            "source packet hidden input")
    require(inputs["P"]["value"] == packet["P"], "mixed source P branches")
    source_values = check_source_packet(packet, p)
    model = payload["model"]
    require(set(model) == {"primary", "control", "charts", "coefficient_menu", "equations",
                           "gauge_solver", "numeric_input_paths", "historical_law_paths",
                           "precision_boundary", "source_response"}, "unknown model field or hidden selector")
    expected_model = {"primary": "full_rscc", "control": "lower_order_ablation_not_selected",
                      "charts": CHARTS, "coefficient_menu": MENU, "gauge_solver": GAUGE_SOLVER,
                      "numeric_input_paths": NUMERIC_INPUTS, "historical_law_paths": HISTORY,
                      "source_response": RESPONSE_CONTRACT,
                      "precision_boundary": "Printed digits are deterministic numerical output, not certified accuracy, model error or calibration uncertainty."}
    for key, value in expected_model.items():
        require(model[key] == value, f"model selection, chart or ancestry changed: {key}")
    require(content_hash(model["equations"]) == STATEMENT_HASHES["equations"],
            "frozen equations changed")
    energy = number(inputs["E_star_GeV"]["value"], "external GeV unit")
    for name, full in (("full_rscc", True), ("lower_order_control", False)):
        check_spectrum(payload[name], p, source_values["alpha_U"],
                       source_values["v_over_E_star"], energy, full)
    pins = payload["source_pins"]
    require(set(pins) == {PRODUCER, GAUGE, *NUMERIC_INPUTS, *HISTORY, *RESPONSE_CODE_PATHS},
            "unexpected source ancestry, target input, or missing source pin")
    for relative, digest in pins.items():
        require(isinstance(digest, str) and len(digest) == 64,
                "invalid source digest")
        require(hashlib.sha256((repo / relative).read_bytes()).hexdigest() == digest,
                f"source pin drift: {relative}")
    spec = json.loads((repo / NUMERIC_INPUTS[1]).read_text())
    require(set(spec) == {"schema", "claim_class", "status", "calibrations", "structural_counts",
                          "charts", "guards", "hypotheses", "coefficient_menu", "equations",
                          "gauge_solver", "source_cache", "historical_law_sha256", "source_response"},
            "spec hidden input")
    projections = {"schema": SCHEMA, "claim_class": CLAIM, "status": STATUS,
                   "calibrations": CALIBRATIONS, "structural_counts": COUNTS,
                   "charts": CHARTS, "guards": GUARDS, "hypotheses": payload["hypotheses"],
                   "coefficient_menu": MENU, "equations": model["equations"],
                   "gauge_solver": GAUGE_SOLVER, "source_cache": packet,
                   "source_response": RESPONSE_CONTRACT,
                   "historical_law_sha256": {path: pins[path] for path in HISTORY}}
    require(spec == projections, "spec and independently checked receipt disagree")
    fixture = json.loads((repo / NUMERIC_INPUTS[0]).read_text())
    require(fixture["inverse_fine_structure_constant"]["value"] == "137.035999177",
            "independent alpha calibration drift")
    source_certificate = json.loads((repo / NUMERIC_INPUTS[2]).read_text())
    source_w5.verify(source_certificate, repo=repo)
    carrier, color = source_certificate["carrier"], source_certificate["color"]
    dimension = carrier["dimension"]
    projection = {"family_dimension": dimension,
                  "line_complement_dimension": dimension-1,
                  "scalar_fraction": str(Fraction(1,dimension)),
                  "centered_fraction": str(Fraction(dimension-1,dimension)),
                  "common_exposure": str(Fraction(3*3-1,2*3)/dimension)}
    require(projection == RESPONSE_CONTRACT["projection"], "independent source coefficient projection")
    require(color["normalized_trace"] == projection["common_exposure"]
            and carrier["D3_line_complement_dimension"] == projection["line_complement_dimension"],
            "source certificate projection mismatch")
    require(str(Fraction(dimension)+Fraction(projection["common_exposure"])) == MENU["mean_u_linear"]
            and projection["common_exposure"] == MENU["mean_d_linear"]
            and str(Fraction(dimension)+Fraction(projection["centered_fraction"])) == MENU["a_u_linear"],
            "certified coefficients do not match the unchanged law")
    producer_text = (repo / PRODUCER).read_text()
    check_producer_inputs(producer_text)
    check_gauge_call_graph(producer_text, (repo / GAUGE).read_text())
    return {"status": "PASS", "rows_verified": 12,
            "scope": "conditional retrospective arithmetic and stated ancestry; no physical mass theorem"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    result = verify_payload(json.loads(args.receipt.read_text()))
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
