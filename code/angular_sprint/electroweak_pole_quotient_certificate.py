#!/usr/bin/env python3
"""Issue #646: exact scale quotient of the strict W/Z pole consumer.

The validated strict-one-loop consumer uses

    w = g^2 v_F^2 / 4,
    z = (g^2 + g'^2) v_F^2 / 4,
    s_W = w + Delta_WW(w),
    s_Z = z + Delta_ZZ(z).

Put ``S = g v_F / 2``, ``t = g'/g``,
``dW = Delta_WW(w)/w`` and ``dZ = Delta_ZZ(z)/z``.  Then

    s_W = S^2 (1+dW),
    s_Z = S^2 (1+t^2) (1+dZ),

so the complex-pole ratio is exactly

    s_W/s_Z = (1+dW) / ((1+t^2)(1+dZ)).

This packet proves the explicit common pole factor cancels globally at fixed
normalized corrections, records the exact local log-Jacobian, and constructs
exact factored-coordinate examples showing that the mass and width ratios vary
algebraically while ``t``, ``dW`` and ``dZ`` are unresolved. A
passive rescaling of all dimensionful units is covered. An active change of
the vacuum scale need not be a null direction because thresholds can change
the normalized corrections. The six mutations are unconstrained algebraic
coordinates, not accepted strict-consumer packets, a global physical
surjectivity theorem, or claimed source-admissible OPH actions. They prevent the
consumer quotient from being mistaken for a numerical postdiction.

Pole residues and asymmetries do not occur in the strict-pole consumer.  The
minimal extension for them is stated explicitly.  In particular, the neutral
residue needs the full neutral inverse-propagator matrix and its Laurent data,
not only a scalar derivative.  Current and vertex form factors must be supplied
by the same source packet before an asymmetry is evaluable.

No measured coupling, mass, width, residue, or asymmetry is read.
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

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "electroweak_pole_quotient_receipt.json"
CONSUMER_SOURCE = (
    REPO_ROOT
    / "code"
    / "particles"
    / "calibration"
    / "strict_one_loop_pole_map"
    / "src"
    / "wz_pole_map.py"
)
CONSUMER_README = (
    REPO_ROOT
    / "code"
    / "particles"
    / "calibration"
    / "strict_one_loop_pole_map"
    / "README.md"
)
CONSUMER_TEMPLATE = (
    REPO_ROOT
    / "code"
    / "particles"
    / "calibration"
    / "strict_one_loop_pole_map"
    / "data"
    / "oph_fj_input_template.json"
)
LEAN_SOURCE = REPO_ROOT / "Lean" / "Screen" / "ElectroweakPoleScaleQuotient.lean"

LEAN_THEOREMS = (
    "common_scale_cancels",
    "common_rescaling_invariant",
    "born_ratio",
    "normalized_w_correction_counterfamily",
    "common_scale_cancels_massWidthVector",
)

SCHEMA = "oph.electroweak_pole_scale_quotient.v1"
STATUS = (
    "EXACT_CONSUMER_SCALE_QUOTIENT__NORMALIZED_SELF_ENERGIES_OPEN__"
    "NO_NUMERICAL_POSTDICTION"
)
INDEPENDENT_VERIFIER = (
    "python3 code/angular_sprint/verify_electroweak_pole_quotient.py "
    "--receipt code/angular_sprint/runtime/"
    "electroweak_pole_quotient_receipt.json"
)


class PoleQuotientError(ValueError):
    """The exact quotient packet refused malformed input."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise PoleQuotientError(message)


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
class GaussianRational:
    """Exact number ``real + imaginary*i`` with rational components."""

    real: Fraction
    imaginary: Fraction = Fraction(0)

    @staticmethod
    def of(real: int | str | Fraction, imaginary: int | str | Fraction = 0) -> "GaussianRational":
        return GaussianRational(Fraction(real), Fraction(imaginary))

    def __add__(self, other: "GaussianRational") -> "GaussianRational":
        return GaussianRational(
            self.real + other.real,
            self.imaginary + other.imaginary,
        )

    def __sub__(self, other: "GaussianRational") -> "GaussianRational":
        return GaussianRational(
            self.real - other.real,
            self.imaginary - other.imaginary,
        )

    def __mul__(self, other: "GaussianRational") -> "GaussianRational":
        return GaussianRational(
            self.real * other.real - self.imaginary * other.imaginary,
            self.real * other.imaginary + self.imaginary * other.real,
        )

    def __truediv__(self, other: "GaussianRational") -> "GaussianRational":
        denominator = other.real * other.real + other.imaginary * other.imaginary
        require(denominator != 0, "division by zero in Q(i)")
        return GaussianRational(
            (self.real * other.real + self.imaginary * other.imaginary)
            / denominator,
            (self.imaginary * other.real - self.real * other.imaginary)
            / denominator,
        )

    def square(self) -> "GaussianRational":
        return self * self

    def pair(self) -> list[str]:
        return [str(self.real), str(self.imaginary)]


ONE = GaussianRational.of(1)


def pole_factor(mass_factor: Fraction, width_over_mass: Fraction) -> GaussianRational:
    """Return ``(mass_factor * (1 - i q/2))^2`` exactly."""

    require(mass_factor > 0, "mass factor must be positive")
    require(width_over_mass >= 0, "width ratio must be nonnegative")
    root = GaussianRational(mass_factor, -mass_factor * width_over_mass / 2)
    return root.square()


def pole_ratio(
    t_squared: Fraction,
    p_w: GaussianRational,
    p_z: GaussianRational,
) -> GaussianRational:
    require(t_squared > 0, "t squared must be positive")
    return p_w / (GaussianRational.of(1 + t_squared) * p_z)


def exact_readout(
    *,
    common_scale_squared: Fraction,
    t_squared: Fraction,
    w_mass_factor: Fraction,
    w_width_over_mass: Fraction,
    z_mass_factor: Fraction,
    z_width_over_mass: Fraction,
) -> dict[str, Any]:
    """Construct one exact decaying-pole coordinate and its scale-free readout."""

    require(common_scale_squared > 0, "common scale squared must be positive")
    p_w = pole_factor(w_mass_factor, w_width_over_mass)
    p_z = pole_factor(z_mass_factor, z_width_over_mass)
    s_w = GaussianRational.of(common_scale_squared) * p_w
    s_z = (
        GaussianRational.of(common_scale_squared * (1 + t_squared)) * p_z
    )
    direct_ratio = s_w / s_z
    factored_ratio = pole_ratio(t_squared, p_w, p_z)
    require(direct_ratio == factored_ratio, "global scale quotient failed")
    mass_ratio_squared = (
        w_mass_factor * w_mass_factor
        / ((1 + t_squared) * z_mass_factor * z_mass_factor)
    )
    return {
        "inputs": {
            "common_scale_squared": str(common_scale_squared),
            "t_squared": str(t_squared),
            "normalized_pole_factor_w": p_w.pair(),
            "normalized_pole_factor_z": p_z.pair(),
        },
        "absolute_poles": {"s_w": s_w.pair(), "s_z": s_z.pair()},
        "dimensionless_outputs": {
            "complex_pole_ratio_s_w_over_s_z": direct_ratio.pair(),
            "mass_ratio_squared": str(mass_ratio_squared),
            "gamma_w_over_m_w": str(w_width_over_mass),
            "gamma_z_over_m_z": str(z_width_over_mass),
        },
    }


def counterfamilies() -> dict[str, Any]:
    """Exact mutations of the unconstrained factored-coordinate directions."""

    base_args = {
        "common_scale_squared": Fraction(4),
        "t_squared": Fraction(1, 4),
        "w_mass_factor": Fraction(1),
        "w_width_over_mass": Fraction(1, 20),
        "z_mass_factor": Fraction(1),
        "z_width_over_mass": Fraction(1, 25),
    }
    base = exact_readout(**base_args)

    def mutation(name: str, **updates: Fraction) -> dict[str, Any]:
        args = dict(base_args)
        args.update(updates)
        mutated = exact_readout(**args)
        before = base["dimensionless_outputs"]
        after = mutated["dimensionless_outputs"]
        changed = sorted(key for key in before if before[key] != after[key])
        unchanged = sorted(key for key in before if before[key] == after[key])
        return {
            "name": name,
            "updated_coordinates": {key: str(value) for key, value in updates.items()},
            "changed_dimensionless_outputs": changed,
            "unchanged_dimensionless_outputs": unchanged,
            "readout": mutated,
        }

    rows = [
        mutation("common_scale", common_scale_squared=Fraction(9)),
        mutation("weak_coupling_ratio", t_squared=Fraction(1, 9)),
        mutation("w_real_pole_factor", w_mass_factor=Fraction(6, 5)),
        mutation("w_absorptive_pole_factor", w_width_over_mass=Fraction(1, 10)),
        mutation("z_real_pole_factor", z_mass_factor=Fraction(7, 6)),
        mutation("z_absorptive_pole_factor", z_width_over_mass=Fraction(2, 25)),
    ]
    expected = {
        "common_scale": [],
        "weak_coupling_ratio": ["complex_pole_ratio_s_w_over_s_z", "mass_ratio_squared"],
        "w_real_pole_factor": ["complex_pole_ratio_s_w_over_s_z", "mass_ratio_squared"],
        "w_absorptive_pole_factor": [
            "complex_pole_ratio_s_w_over_s_z",
            "gamma_w_over_m_w",
        ],
        "z_real_pole_factor": ["complex_pole_ratio_s_w_over_s_z", "mass_ratio_squared"],
        "z_absorptive_pole_factor": [
            "complex_pole_ratio_s_w_over_s_z",
            "gamma_z_over_m_z",
        ],
    }
    for row in rows:
        require(
            row["changed_dimensionless_outputs"] == expected[row["name"]],
            f"counterfamily mutation drift for {row['name']}",
        )
    return {
        "typing": (
            "exact counterfamilies on unconstrained factored coordinates; they "
            "are not accepted strict-consumer packets, source-admissible OPH "
            "actions, or physical EFT completions"
        ),
        "base": base,
        "mutations": rows,
    }


def parent_pin(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
    }


@lru_cache(maxsize=1)
def require_lean_kernel_build() -> None:
    """Require Lake to compile the pinned theorem module successfully."""

    completed = subprocess.run(
        ["lake", "build", "ElectroweakPoleScaleQuotient"],
        cwd=REPO_ROOT / "Lean",
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


def lean_parent_pin() -> dict[str, Any]:
    require_lean_kernel_build()
    payload = LEAN_SOURCE.read_bytes()
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PoleQuotientError("electroweak quotient Lean source is not UTF-8") from error
    require("sorry" not in text and "admit" not in text, "Lean placeholder present")
    for theorem in LEAN_THEOREMS:
        require(
            re.search(rf"\btheorem\s+{re.escape(theorem)}\b", text) is not None,
            f"missing Lean theorem: {theorem}",
        )
    return {
        "path": LEAN_SOURCE.relative_to(REPO_ROOT).as_posix(),
        "role": "kernel-checked common-scale cancellation algebra",
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
        "sorry_free": True,
        "theorems": list(LEAN_THEOREMS),
    }


def build_receipt() -> dict[str, Any]:
    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": 646,
        "status": STATUS,
        "consumer_convention": {
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
        },
        "global_quotient": {
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
        },
        "symbolic_log_jacobian": {
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
                "d_log_sZ": [
                    "2",
                    "2*t^2/(1+t^2)",
                    "0",
                    "1/(1+dZ)",
                ],
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
        },
        "frozen_output_vocabulary": [
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
        ],
        "upstream_dependency_factorization": {
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
        },
        "minimal_surviving_source_contract": {
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
                "one correlated dimensionless vector (M_W/M_Z, "
                "Gamma_W/M_W, Gamma_Z/M_Z) becomes consumer-evaluable "
                "without a physical clock"
            ),
        },
        "counterfamilies": counterfamilies(),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
            "candidate_selected_by_target_proximity": False,
        },
        "lean_check": "Lean/Screen/ElectroweakPoleScaleQuotient.lean",
        "verification": {
            "command": INDEPENDENT_VERIFIER,
            "classification": (
                "independent exact rational-complex replay and symbolic differentiation "
                "of the quotient, parent pins, and factored-coordinate examples"
            ),
        },
        "parent_pins": [
            parent_pin(CONSUMER_SOURCE),
            parent_pin(CONSUMER_README),
            parent_pin(CONSUMER_TEMPLATE),
            lean_parent_pin(),
        ],
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return RECEIPT_PATH


def verify_runtime() -> None:
    require(RECEIPT_PATH.exists(), "runtime receipt is missing")
    require(
        RECEIPT_PATH.read_bytes() == canonical_json_bytes(build_receipt()),
        "runtime receipt is stale",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("ELECTROWEAK_POLE_QUOTIENT_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        print(
            json.dumps(
                {
                    "status": receipt["status"],
                    "identity": receipt["global_quotient"]["identity"],
                    "common_scale_cancels": receipt["global_quotient"][
                        "common_scale_cancels_exactly"
                    ],
                    "minimum_source_inputs": receipt[
                        "minimal_surviving_source_contract"
                    ]["mass_width_triple"],
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
