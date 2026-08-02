#!/usr/bin/env python3
"""Independent verifier for the issue-646 strict-pole scale quotient.

This verifier imports neither the quotient producer nor the strict-pole
consumer.  It checks the pinned consumer files, reconstructs the exact
rational-complex quotient and mutation family, and checks every frozen
formula and boundary in the receipt.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
import re
import subprocess
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
DEFAULT_RECEIPT = HERE / "runtime" / "electroweak_pole_quotient_receipt.json"
SCHEMA = "oph.electroweak_pole_scale_quotient.v1"
EXPECTED_STATUS = (
    "EXACT_CONSUMER_SCALE_QUOTIENT__NORMALIZED_SELF_ENERGIES_OPEN__"
    "NO_NUMERICAL_POSTDICTION"
)
EXPECTED_COMMAND = (
    "python3 code/angular_sprint/verify_electroweak_pole_quotient.py "
    "--receipt code/angular_sprint/runtime/"
    "electroweak_pole_quotient_receipt.json"
)
EXPECTED_TOP_LEVEL_KEYS = {
    "schema",
    "issue",
    "status",
    "consumer_convention",
    "global_quotient",
    "symbolic_log_jacobian",
    "frozen_output_vocabulary",
    "upstream_dependency_factorization",
    "minimal_surviving_source_contract",
    "counterfamilies",
    "comparison_boundary",
    "lean_check",
    "verification",
    "parent_pins",
    "receipt_sha256",
}
EXPECTED_LEAN_THEOREMS = [
    "common_scale_cancels",
    "common_rescaling_invariant",
    "born_ratio",
    "normalized_w_correction_counterfamily",
    "common_scale_cancels_massWidthVector",
]
EXPECTED_PARENT_PINS = [
    {
        "path": "code/particles/calibration/strict_one_loop_pole_map/src/wz_pole_map.py",
        "bytes": 13807,
        "sha256": "sha256:30679960ecd2baa33458773ff64579d6402e9f7d7fd52123ea7c2051181bdc69",
    },
    {
        "path": "code/particles/calibration/strict_one_loop_pole_map/README.md",
        "bytes": 5546,
        "sha256": "sha256:8b0a4fabb687d86338b33c112bdce8dee391f858ea6226f7dc4c6ed3273bef38",
    },
    {
        "path": (
            "code/particles/calibration/strict_one_loop_pole_map/data/"
            "oph_fj_input_template.json"
        ),
        "bytes": 2418,
        "sha256": "sha256:39daa89a581a0991fe74a200d8d0f3dfcbb435c7c17577dfd3c05273fe5ce1ec",
    },
    {
        "path": "Lean/Screen/ElectroweakPoleScaleQuotient.lean",
        "role": "kernel-checked common-scale cancellation algebra",
        "bytes": 4590,
        "sha256": "sha256:9305d115cfc8b35aa8c1179330a7ea3ec60cb1de9c890651cee59bbd9fff76ae",
        "sorry_free": True,
        "theorems": EXPECTED_LEAN_THEOREMS,
    },
]
EXPECTED_CONSUMER_CONVENTION = {
    "inverse_propagator": "Gamma=s-m2-Delta; s_pole=m2+Delta_at_tree",
    "tree_coordinates": {
        "S": "g*v_F/2",
        "t": "g_prime/g",
        "w": "S^2",
        "z": "S^2*(1+t^2)",
    },
    "normalized_corrections": {
        "dW": "Delta_WW(w)/w",
        "dZ": "Delta_ZZ(z)/z",
    },
    "strict_poles": {
        "sW": "S^2*(1+dW)",
        "sZ": "S^2*(1+t^2)*(1+dZ)",
    },
}
EXPECTED_GLOBAL_QUOTIENT = {
    "identity": "sW/sZ=(1+dW)/((1+t^2)*(1+dZ))",
    "domain": "S!=0, 1+t^2!=0, 1+dZ!=0",
    "common_scale_cancels_exactly": True,
    "born_specialization": "dW=dZ=0 => sW/sZ=1/(1+t^2)",
    "truncation_typing": (
        "exact quotient of the two one-loop-truncated pole coordinates; "
        "it is not itself the quotient re-expanded to strict one-loop order"
    ),
    "strict_one_loop_expansion": (
        "sW/sZ=(1/(1+t^2))*(1+dW-dZ)+O(loop^2)"
    ),
    "physical_boundary": (
        "the identity is exact consumer algebra at fixed normalized "
        "corrections; it covers passive common-unit rescaling but does "
        "not prove that active vacuum-scale changes leave dW and dZ "
        "fixed, and it does not select t, dW, dZ, a pole sheet, or a "
        "physical source action"
    ),
}
EXPECTED_JACOBIAN = {
    "input_differentials": ["d_log_S", "d_log_t", "d_dW", "d_dZ"],
    "coordinate_field": (
        "four complexified factored coordinates; the displayed ratio row "
        "has rank one and kernel dimension three over C"
    ),
    "physical_real_slice": (
        "S and t real positive with dW and dZ complex gives six real input "
        "directions; the complex ratio has generic real rank two and a "
        "four-real-dimensional kernel: common scale plus three correlated "
        "directions"
    ),
    "domain": "S!=0, t!=0, 1+t^2!=0, 1+dW!=0, 1+dZ!=0",
    "rows": {
        "d_log_sW": ["2", "0", "1/(1+dW)", "0"],
        "d_log_sZ": ["2", "2*t^2/(1+t^2)", "0", "1/(1+dZ)"],
        "d_log_sW_over_sZ": [
            "0",
            "-2*t^2/(1+t^2)",
            "1/(1+dW)",
            "-1/(1+dZ)",
        ],
    },
    "exact_null_direction": {
        "direction": ["1", "0", "0", "0"],
        "absolute_pole_image": ["2", "2"],
        "ratio_image": "0",
        "typing": (
            "partial consumer-coordinate derivative at fixed t, dW, "
            "and dZ; not a total source derivative under an active "
            "vacuum-scale change"
        ),
    },
    "noncancellation": (
        "on the declared nonzero domain the t, dW, and dZ columns of "
        "the ratio row are individually nonzero, so none is a "
        "coordinate-axis null direction; over the complexified coordinate "
        "algebra the one-row map has a three-complex-dimensional local "
        "kernel consisting of the common-scale direction and two correlated "
        "t-dW-dZ directions"
    ),
}
EXPECTED_OUTPUT_VOCABULARY = [
    {
        "output": "M_W/M_Z",
        "convention": "exact energy coordinates of the two truncated complex poles",
        "consumer_inputs": ["t", "dW", "dZ", "decaying square-root sheet"],
        "common_scale_dependency": "explicit_S_cancels_at_fixed_normalized_corrections",
    },
    {
        "output": "Gamma_W/M_W",
        "convention": "exact energy coordinates of the truncated W pole",
        "consumer_inputs": ["dW", "decaying square-root sheet"],
        "common_scale_dependency": "explicit_S_cancels_at_fixed_normalized_correction",
    },
    {
        "output": "Gamma_Z/M_Z",
        "convention": "exact energy coordinates of the truncated Z pole",
        "consumer_inputs": ["dZ", "decaying square-root sheet"],
        "common_scale_dependency": "explicit_S_cancels_at_fixed_normalized_correction",
    },
    {
        "output": "pole-residue ratios",
        "convention": "matrix Laurent data in the same convention-separated inverse propagator",
        "consumer_inputs": [
            "charged inverse-propagator derivative at the W pole",
            "full neutral inverse/self-energy matrix and derivative at the Z pole",
            "neutral adjugate or normalized left/right null-vector data",
            "current and field normalization",
        ],
        "status": "not_emitted_by_strict_pole_consumer",
    },
    {
        "output": "asymmetry combinations",
        "convention": "effective-current form factors kept separate from pole coordinates",
        "consumer_inputs": [
            "complete neutral residue packet",
            "dressed renormalized neutral-current vertices",
            "external-state and effective-angle convention",
        ],
        "status": "not_emitted_by_strict_pole_consumer",
    },
]
EXPECTED_DEPENDENCY_FACTORIZATION = {
    "passive_common_unit_or_clock_scale": {
        "interface_coordinates": ["S"],
        "dimensionless_output_status": "proved_canceled",
    },
    "active_vacuum_scale": {
        "interface_coordinates": ["S", "dW", "dZ", "derivatives", "vertices"],
        "dimensionless_output_status": (
            "explicit_S_factor_cancels; total derivative unresolved "
            "because thresholds can move normalized corrections"
        ),
    },
    "g_prime_over_g": {
        "interface_coordinates": ["t"],
        "dimensionless_output_status": "survives",
    },
    "scalar_and_yukawa": {
        "interface_coordinates": ["dW", "dZ", "derivatives", "vertices"],
        "dimensionless_output_status": "unresolved_upstream_map",
    },
    "field_normalization": {
        "interface_coordinates": ["residues", "current vertices"],
        "dimensionless_output_status": "unresolved_outside_pole_locations",
    },
    "scheme_and_threshold": {
        "interface_coordinates": ["t", "dW", "dZ", "derivatives", "vertices"],
        "dimensionless_output_status": "unresolved_upstream_map",
    },
    "family_and_scalar_census": {
        "interface_coordinates": ["dW", "dZ", "derivatives", "vertices"],
        "dimensionless_output_status": "unresolved_upstream_map",
    },
    "repair_selector_and_continuum": {
        "interface_coordinates": ["t", "dW", "dZ", "derivatives", "vertices"],
        "dimensionless_output_status": "unresolved_upstream_map",
    },
}
EXPECTED_SOURCE_CONTRACT = {
    "mass_width_triple": [
        "source-selected t^2=(g_prime/g)^2 in one declared scheme",
        "source-selected normalized complex correction dW",
        "source-selected normalized complex correction dZ",
        "simple decaying-pole and shared convention certificate",
    ],
    "residue_extension": [
        "charged inverse-propagator derivative at the W pole",
        "full neutral inverse/self-energy matrix at the Z pole and its derivative",
        "neutral adjugate data or normalized left/right null vectors",
        "source-selected current and field normalization",
    ],
    "asymmetry_extension": [
        "complete residue-extension packet",
        "dressed renormalized neutral-current vertex form factors",
        "external-state and effective-angle convention",
    ],
    "absolute_scale_required": False,
    "absolute_scale_boundary": (
        "false at the consumer interface once dimensionless dW and dZ "
        "are supplied; their production can require source-selected "
        "threshold and renormalization-scale ratios"
    ),
    "result_if_supplied": (
        "one correlated dimensionless vector (M_W/M_Z, Gamma_W/M_W, "
        "Gamma_Z/M_Z) becomes consumer-evaluable without a physical clock"
    ),
}
EXPECTED_COMPARISON_BOUNDARY = {
    "public_measurement_read": False,
    "comparison_permitted": False,
    "candidate_selected_by_target_proximity": False,
}
EXPECTED_VERIFICATION = {
    "command": EXPECTED_COMMAND,
    "classification": (
        "independent exact rational-complex replay and symbolic differentiation "
        "of the quotient, parent pins, and factored-coordinate examples"
    ),
}
EXPECTED_COUNTERFAMILY_TYPING = (
    "exact counterfamilies on unconstrained factored coordinates; they are not "
    "accepted strict-consumer packets, source-admissible OPH actions, or "
    "physical EFT completions"
)
EXPECTED_BASE_INPUTS = {
    "common_scale_squared": "4",
    "t_squared": "1/4",
    "normalized_pole_factor_w": ["1599/1600", "-1/20"],
    "normalized_pole_factor_z": ["2499/2500", "-1/25"],
}
EXPECTED_MUTATION_UPDATES = {
    "common_scale": {"common_scale_squared": "9"},
    "weak_coupling_ratio": {"t_squared": "1/9"},
    "w_real_pole_factor": {"w_mass_factor": "6/5"},
    "w_absorptive_pole_factor": {"w_width_over_mass": "1/10"},
    "z_real_pole_factor": {"z_mass_factor": "7/6"},
    "z_absorptive_pole_factor": {"z_width_over_mass": "2/25"},
}
EXPECTED_MUTATION_INPUTS = {
    "common_scale": {
        **EXPECTED_BASE_INPUTS,
        "common_scale_squared": "9",
    },
    "weak_coupling_ratio": {
        **EXPECTED_BASE_INPUTS,
        "t_squared": "1/9",
    },
    "w_real_pole_factor": {
        **EXPECTED_BASE_INPUTS,
        "normalized_pole_factor_w": ["14391/10000", "-9/125"],
    },
    "w_absorptive_pole_factor": {
        **EXPECTED_BASE_INPUTS,
        "normalized_pole_factor_w": ["399/400", "-1/10"],
    },
    "z_real_pole_factor": {
        **EXPECTED_BASE_INPUTS,
        "normalized_pole_factor_z": ["40817/30000", "-49/900"],
    },
    "z_absorptive_pole_factor": {
        **EXPECTED_BASE_INPUTS,
        "normalized_pole_factor_z": ["624/625", "-2/25"],
    },
}


class VerificationError(ValueError):
    """The independent verifier rejected an input."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def strict_json_bytes(payload: bytes, *, label: str) -> dict[str, Any]:
    def object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for key, value in pairs:
            if key in out:
                raise VerificationError(f"{label}: duplicate JSON key {key!r}")
            out[key] = value
        return out

    def reject_constant(value: str) -> Any:
        raise VerificationError(f"{label}: non-finite JSON constant {value!r}")

    try:
        value = json.loads(
            payload.decode("ascii"),
            object_pairs_hook=object_pairs,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f"{label}: invalid JSON: {exc}") from exc
    require(isinstance(value, dict), f"{label}: top-level value is not an object")
    return value


def strict_load(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    value = strict_json_bytes(payload, label=path.as_posix())
    require(
        payload == canonical_json_bytes(value),
        f"{path.as_posix()}: noncanonical JSON",
    )
    return value


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(payload: bytes) -> str:
    return "sha256:" + hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class QI:
    real: Fraction
    imaginary: Fraction = Fraction(0)

    @staticmethod
    def of(real: int | str | Fraction, imaginary: int | str | Fraction = 0) -> "QI":
        return QI(Fraction(real), Fraction(imaginary))

    @staticmethod
    def pair(value: Any) -> "QI":
        require(isinstance(value, list) and len(value) == 2, "Q(i) pair malformed")
        return QI.of(value[0], value[1])

    def __add__(self, other: "QI") -> "QI":
        return QI(self.real + other.real, self.imaginary + other.imaginary)

    def __mul__(self, other: "QI") -> "QI":
        return QI(
            self.real * other.real - self.imaginary * other.imaginary,
            self.real * other.imaginary + self.imaginary * other.real,
        )

    def __truediv__(self, other: "QI") -> "QI":
        denominator = other.real * other.real + other.imaginary * other.imaginary
        require(denominator != 0, "division by zero in Q(i)")
        return QI(
            (self.real * other.real + self.imaginary * other.imaginary)
            / denominator,
            (self.imaginary * other.real - self.real * other.imaginary)
            / denominator,
        )


@lru_cache(maxsize=1)
def require_lean_kernel_build(repo_root: Path) -> None:
    """Fail unless Lake compiles the independently pinned Lean module."""

    completed = subprocess.run(
        ["lake", "build", "ElectroweakPoleScaleQuotient"],
        cwd=repo_root / "Lean",
        check=False,
        capture_output=True,
        text=True,
    )
    require(
        completed.returncode == 0,
        "electroweak quotient Lean build failed:\n"
        + completed.stdout[-4000:]
        + completed.stderr[-4000:],
    )


def check_parent_pins(receipt: dict[str, Any], repo_root: Path) -> None:
    pins = receipt.get("parent_pins")
    require(pins == EXPECTED_PARENT_PINS, "parent pin inventory drift")
    for row in EXPECTED_PARENT_PINS:
        path_text = row["path"]
        path = repo_root / path_text
        payload = path.read_bytes()
        require(len(payload) == row["bytes"], f"parent byte drift: {path_text}")
        require(
            tagged_sha256(payload) == row["sha256"],
            f"parent hash drift: {path_text}",
        )
    require_lean_kernel_build(repo_root)
    lean = EXPECTED_PARENT_PINS[-1]
    text = (repo_root / lean["path"]).read_text(encoding="utf-8")
    require("sorry" not in text and "admit" not in text, "Lean placeholder present")
    for theorem in EXPECTED_LEAN_THEOREMS:
        require(
            re.search(rf"\btheorem\s+{re.escape(theorem)}\b", text) is not None,
            f"missing Lean theorem: {theorem}",
        )


def check_exact_readout(readout: dict[str, Any]) -> dict[str, str | list[str]]:
    require(
        set(readout) == {"inputs", "absolute_poles", "dimensionless_outputs"},
        "readout key inventory drift",
    )
    inputs = readout["inputs"]
    outputs = readout["dimensionless_outputs"]
    absolute = readout["absolute_poles"]
    require(
        set(inputs)
        == {
            "common_scale_squared",
            "t_squared",
            "normalized_pole_factor_w",
            "normalized_pole_factor_z",
        },
        "readout input inventory drift",
    )
    require(set(absolute) == {"s_w", "s_z"}, "absolute-pole inventory drift")
    require(
        set(outputs)
        == {
            "complex_pole_ratio_s_w_over_s_z",
            "mass_ratio_squared",
            "gamma_w_over_m_w",
            "gamma_z_over_m_z",
        },
        "dimensionless-output inventory drift",
    )
    scale = Fraction(inputs["common_scale_squared"])
    t2 = Fraction(inputs["t_squared"])
    require(scale > 0 and t2 > 0, "nonpositive readout coordinate")
    p_w = QI.pair(inputs["normalized_pole_factor_w"])
    p_z = QI.pair(inputs["normalized_pole_factor_z"])
    s_w = QI.pair(absolute["s_w"])
    s_z = QI.pair(absolute["s_z"])
    require(s_w == QI.of(scale) * p_w, "absolute W pole mismatch")
    require(s_z == QI.of(scale * (1 + t2)) * p_z, "absolute Z pole mismatch")
    ratio = s_w / s_z
    require(
        ratio == p_w / (QI.of(1 + t2) * p_z),
        "global quotient mismatch",
    )
    require(
        ratio == QI.pair(outputs["complex_pole_ratio_s_w_over_s_z"]),
        "stored complex-pole ratio mismatch",
    )
    q_w = Fraction(outputs["gamma_w_over_m_w"])
    q_z = Fraction(outputs["gamma_z_over_m_z"])
    require(q_w >= 0 and q_z >= 0, "negative normalized width")

    def root_mass_squared(pole: QI, width_ratio: Fraction, label: str) -> Fraction:
        mass_squared = (
            -pole.imaginary / width_ratio
            if width_ratio != 0
            else pole.real
        )
        expected = QI(
            mass_squared * (1 - width_ratio * width_ratio / 4),
            -mass_squared * width_ratio,
        )
        require(mass_squared > 0 and pole == expected, f"{label} pole-factor mismatch")
        return mass_squared

    mass_w_squared = root_mass_squared(p_w, q_w, "W")
    mass_z_squared = root_mass_squared(p_z, q_z, "Z")
    require(
        Fraction(outputs["mass_ratio_squared"])
        == mass_w_squared / ((1 + t2) * mass_z_squared),
        "stored mass-ratio mismatch",
    )
    return outputs


def check_symbolic_jacobian(jacobian: dict[str, Any]) -> None:
    """Derive the logarithmic rows without trusting stored formula strings."""

    scale, coupling_ratio, d_w, d_z = sp.symbols("S t dW dZ")
    s_w = scale**2 * (1 + d_w)
    s_z = scale**2 * (1 + coupling_ratio**2) * (1 + d_z)
    variables = (scale, coupling_ratio, d_w, d_z)
    logarithmic_weights = (scale, coupling_ratio, sp.Integer(1), sp.Integer(1))

    def row(expression: sp.Expr) -> list[sp.Expr]:
        return [
            sp.simplify(weight * sp.diff(expression, variable) / expression)
            for variable, weight in zip(variables, logarithmic_weights, strict=True)
        ]

    derived = {
        "d_log_sW": row(s_w),
        "d_log_sZ": row(s_z),
        "d_log_sW_over_sZ": row(sp.cancel(s_w / s_z)),
    }
    parse_locals = {"S": scale, "t": coupling_ratio, "dW": d_w, "dZ": d_z}
    for name, exact_row in derived.items():
        stored_row = [sp.sympify(value, locals=parse_locals) for value in jacobian["rows"][name]]
        require(
            all(sp.simplify(stored - exact) == 0 for stored, exact in zip(stored_row, exact_row, strict=True)),
            f"independent symbolic Jacobian mismatch: {name}",
        )

    ratio_row = sp.Matrix([derived["d_log_sW_over_sZ"]])
    require(ratio_row.rank() == 1, "independent ratio-row rank mismatch")
    require(
        len(ratio_row.nullspace()) == 3,
        "independent complexified ratio-kernel dimension mismatch",
    )
    require(
        all(sp.simplify(value) != 0 for value in derived["d_log_sW_over_sZ"][1:]),
        "independent surviving-coordinate test failed",
    )
    common_scale = sp.Matrix([1, 0, 0, 0])
    require(
        (ratio_row * common_scale)[0] == 0,
        "independent common-scale null-direction mismatch",
    )

    # On the physical slice, S and t are real while dW and dZ are complex.
    # The real-linear block contributed by dW alone is multiplication by
    # 1/(1+dW).  Its determinant is |1/(1+dW)|^2, so the complex ratio has
    # generic real rank two.  Six real inputs therefore leave a four-real-
    # dimensional kernel.
    d_w_real, d_w_imaginary = sp.symbols("dWr dWi", real=True)
    denominator = (1 + d_w_real) ** 2 + d_w_imaginary**2
    coefficient_real = (1 + d_w_real) / denominator
    coefficient_imaginary = -d_w_imaginary / denominator
    d_w_real_block = sp.Matrix(
        [
            [coefficient_real, -coefficient_imaginary],
            [coefficient_imaginary, coefficient_real],
        ]
    )
    require(
        sp.simplify(d_w_real_block.det() - 1 / denominator) == 0,
        "independent physical-slice rank certificate failed",
    )


def verify_receipt(
    receipt: dict[str, Any], *, repo_root: Path = REPO_ROOT
) -> dict[str, Any]:
    reasons: list[str] = []
    try:
        require(set(receipt) == EXPECTED_TOP_LEVEL_KEYS, "top-level key inventory drift")
        require(receipt.get("schema") == SCHEMA, "schema drift")
        require(receipt.get("issue") == 646, "issue drift")
        require(receipt.get("status") == EXPECTED_STATUS, "status drift")
        claimed_hash = receipt.get("receipt_sha256")
        payload = dict(receipt)
        payload.pop("receipt_sha256", None)
        require(
            claimed_hash == tagged_sha256(canonical_json_bytes(payload)),
            "receipt self-hash mismatch",
        )
        check_parent_pins(receipt, repo_root)
        require(
            receipt.get("lean_check")
            == "Lean/Screen/ElectroweakPoleScaleQuotient.lean",
            "Lean check path drift",
        )
        require(
            receipt.get("consumer_convention") == EXPECTED_CONSUMER_CONVENTION,
            "consumer convention drift",
        )

        quotient = receipt["global_quotient"]
        require(quotient == EXPECTED_GLOBAL_QUOTIENT, "global quotient statement drift")

        jacobian = receipt["symbolic_log_jacobian"]
        require(jacobian == EXPECTED_JACOBIAN, "symbolic Jacobian drift")
        check_symbolic_jacobian(jacobian)
        require(
            receipt["frozen_output_vocabulary"] == EXPECTED_OUTPUT_VOCABULARY,
            "output vocabulary drift",
        )
        require(
            receipt["upstream_dependency_factorization"]
            == EXPECTED_DEPENDENCY_FACTORIZATION,
            "dependency factorization drift",
        )

        counter = receipt["counterfamilies"]
        require(
            set(counter) == {"typing", "base", "mutations"},
            "counterfamily key inventory drift",
        )
        require(
            counter["typing"] == EXPECTED_COUNTERFAMILY_TYPING,
            "counterfamily typing drift",
        )
        require(
            counter["base"].get("inputs") == EXPECTED_BASE_INPUTS,
            "counterfamily base drift",
        )
        base_outputs = check_exact_readout(counter["base"])
        expected = {
            "common_scale": [],
            "weak_coupling_ratio": [
                "complex_pole_ratio_s_w_over_s_z",
                "mass_ratio_squared",
            ],
            "w_real_pole_factor": [
                "complex_pole_ratio_s_w_over_s_z",
                "mass_ratio_squared",
            ],
            "w_absorptive_pole_factor": [
                "complex_pole_ratio_s_w_over_s_z",
                "gamma_w_over_m_w",
            ],
            "z_real_pole_factor": [
                "complex_pole_ratio_s_w_over_s_z",
                "mass_ratio_squared",
            ],
            "z_absorptive_pole_factor": [
                "complex_pole_ratio_s_w_over_s_z",
                "gamma_z_over_m_z",
            ],
        }
        names: list[str] = []
        for row in counter["mutations"]:
            require(
                set(row)
                == {
                    "name",
                    "updated_coordinates",
                    "changed_dimensionless_outputs",
                    "unchanged_dimensionless_outputs",
                    "readout",
                },
                "counterfamily row key inventory drift",
            )
            outputs = check_exact_readout(row["readout"])
            changed = sorted(
                key for key in base_outputs if base_outputs[key] != outputs[key]
            )
            name = row["name"]
            names.append(name)
            require(name in expected, f"unknown mutation row: {name}")
            require(
                row["readout"]["inputs"] == EXPECTED_MUTATION_INPUTS[name],
                f"mutation readout-input drift: {name}",
            )
            require(changed == expected[name], f"mutation semantics drift: {name}")
            require(
                changed == row["changed_dimensionless_outputs"],
                f"mutation claim drift: {name}",
            )
            require(
                row["unchanged_dimensionless_outputs"]
                == sorted(key for key in base_outputs if key not in changed),
                f"mutation unchanged-output drift: {name}",
            )
            require(
                row["updated_coordinates"] == EXPECTED_MUTATION_UPDATES[name],
                f"mutation coordinate drift: {name}",
            )
        require(names == list(expected), "mutation inventory or order drift")

        contract = receipt["minimal_surviving_source_contract"]
        require(contract == EXPECTED_SOURCE_CONTRACT, "minimal source contract drift")
        comparison = receipt["comparison_boundary"]
        require(comparison == EXPECTED_COMPARISON_BOUNDARY, "comparison boundary drift")
        require(
            receipt["verification"] == EXPECTED_VERIFICATION,
            "independent verifier metadata drift",
        )
    except (
        ArithmeticError,
        KeyError,
        OSError,
        TypeError,
        UnicodeError,
        ValueError,
        VerificationError,
    ) as exc:
        reasons.append(str(exc))
    return {
        "schema": "oph.electroweak_pole_scale_quotient.verification.v1",
        "receipt": not reasons,
        "reasons": reasons,
        "exact_global_scale_quotient": not reasons,
        "counterfamily_mutations": 6 if not reasons else 0,
        "comparison_data_read": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args()
    result = verify_receipt(strict_load(args.receipt))
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["receipt"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
