#!/usr/bin/env python3
"""Exact certificate for GitHub issue #314: the source-bound super-Tannakian matter lift.

The input is a matter-lift manifest.  It declares only:

* the port-current response manifest of issue #566 (by path and hash)
  together with its stored receipt hash - the exact carrier-dynamics
  response source of issue #599 is strictly upstream;
* the measured spin statistics artifact of the simulator (by path and
  hash) - the target-blind transport measurement of the deck lift group,
  its centre, its section obstruction, and the unique spin structure;
* the trace-balanced exterior matter contract: the exact charge pair on
  the color and weak blocks, the one-scalar choice (the weak block
  itself), and the declared invariant Yukawa channel list - the charges
  are validated against the balance derivation and the channel list
  against the compatibility scan;
* the statistics and category contracts, which are validated against the
  scan and forcing derivations below rather than accepted as premises,
  the operator-projector realization contract, the kernel emission
  contract, and the declared candidate matter class with uniqueness promotion
  switched off.

From that packet the verifier derives, rather than assumes:

* the exact algebraic PORT-SPIN-LIFT target: special-unitary lifts of all
  sixty proper implementers whose 120-element lift group has a unique
  involution, so the double cover is genuinely non-split (binary
  icosahedral);
* a faithful skew-adjoint Lie-algebra homomorphism from the twelve
  dimensional current algebra to the matter carrier V = C (+) W;
* the auxiliary CAR/Fock space Lambda^* V with derived fermionic parity,
  super tensor structure (wedge with Koszul signs), and conjugation
  through the pairing into the top line, which is exactly trivial because
  of the declared trace balance;
* the unordered conjugate pair of equivariant projectors (a parity
  projector minus its complete invariant line), each of exact rank
  fifteen; neither Weyl representative is selected by the source;
* exact chirality (disjoint charge spectra force a zero intertwiner
  space), realized perturbative anomaly traces, the even Witten parity,
  and exactly one invariant line per declared Yukawa channel;
* the common action kernel on the simply connected cover R x SU(3) x SU(2),
  emitted as data: infinite cyclic with generator (zeta_6-turn, omega, -1),
  whose sixth power is the unit deck translation (one full central turn) -
  not the identity on the cover - with residual order six modulo the pure
  deck translations; neither the compactification of the central R nor any
  global quotient is formed;
* descent along the declared algebraic carrier tower maps.

Beyond the conditional algebra, the verifier derives the physical typing
at finite source-model scope: the exhaustive 1024-subset anomaly scan
selects the unordered conjugate rank-fifteen pair with the fermionic-
parity grading as an output; the measured lift centre {+1, -1} and the
measured section obstruction over every Klein four-subgroup make the
central -1 of the non-split double cover the unique source
implementation of that grading, with the gauge-centre branch excluded by
the Lean fermion-parity no-go; and the Spin/odd-Weyl super typing is
thereby forced rather than declared.  Scalar existence and economy stay
typed branch premises owned by issue #609, listed as deferred rows that
never enter the passing gate.

Vec-typed, split-sVec, opposite-Weyl, bosonic-statistics, truncated
selection, full-even-module, empty-Gauss, assumed-quotient,
kernel-killing, representation-arithmetic-only, charge-dead, unbalanced,
uniqueness-promoting, and firewall countermodels fail closed with typed
error codes.  No family attachment, scalar-potential, pole-mass, measured
coupling, or global-form choice is accepted in a source manifest.  Every
arithmetic decision is exact in Q(sqrt5); no floating point appears in a
proof step.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as p566  # noqa: E402

SCHEMA = "oph.super_tannakian_matter_manifest.v5"
RECEIPT_SCHEMA = "oph.super_tannakian_matter_receipt.v5"
NEGATIVE_SCHEMA = "oph.super_tannakian_matter_negative_controls.v5"

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

F5 = p566.F5
C5 = p566.C5
ZERO = p566.ZERO
ONE = p566.ONE
czeros = p566.czeros
cmul = p566.cmul
cadd = p566.cadd
csub = p566.csub
cdagger = p566.cdagger
commutator = p566.commutator
ctrace = p566.ctrace
c_is_zero = p566.c_is_zero
rref = p566.rref
nullspace = p566.nullspace
rank = p566.rank

IZERO = C5(ZERO, ZERO)
IONE = C5(ONE, ZERO)

FORBIDDEN_MATTER_TOKENS = (
    "threefamil",
    "familyattachment",
    "familyrank",
    "scalarpotential",
    "polemass",
    "yukawacouplingvalue",
    "ckm",
    "pmns",
    "neutrinomass",
    "globalquotientchoice",
)


# ---------------------------------------------------------------------------
# Exact square roots in Q(sqrt5)
# ---------------------------------------------------------------------------


def sqrt_fraction(value: Fraction) -> Fraction | None:
    """Exact square root of a non-negative rational, or None."""

    if value < 0:
        return None
    num = math.isqrt(value.numerator)
    den = math.isqrt(value.denominator)
    if num * num != value.numerator or den * den != value.denominator:
        return None
    return Fraction(num, den)


def sqrt_f5(value: F5) -> F5 | None:
    """An exact square root of value in Q(sqrt5), or None if none exists."""

    if value.is_zero():
        return F5(0)
    if value.b == 0:
        root = sqrt_fraction(value.a)
        if root is not None:
            return F5(root)
        root = sqrt_fraction(value.a / 5)
        if root is not None:
            return F5(0, root)
        return None
    # (x + y*sqrt5)^2 = value: x^2 + 5y^2 = value.a and 2xy = value.b.
    # x^2 solves t^2 - value.a * t + 5 (value.b/2)^2 = 0.
    disc = value.a * value.a - 5 * value.b * value.b
    disc_root = sqrt_fraction(disc)
    if disc_root is None:
        return None
    for branch in ((value.a + disc_root) / 2, (value.a - disc_root) / 2):
        x = sqrt_fraction(branch)
        if x is None or x == 0:
            continue
        y = value.b / (2 * x)
        candidate = F5(x, y)
        if (candidate * candidate - value).is_zero():
            return candidate
    return None


# ---------------------------------------------------------------------------
# Small exact helpers on C5 matrices
# ---------------------------------------------------------------------------


def cidentity(n: int) -> list[list[C5]]:
    return [[IONE if i == j else IZERO for j in range(n)] for i in range(n)]


def cscale(x: Sequence[Sequence[C5]], s: C5) -> list[list[C5]]:
    return [[s * entry for entry in row] for row in x]


def cdet(m: Sequence[Sequence[C5]]) -> C5:
    n = len(m)
    if n == 0:
        return IONE
    if n == 1:
        return m[0][0]
    total = IZERO
    sign_positive = True
    for j in range(n):
        entry = m[0][j]
        if not entry.is_zero():
            minor = [[m[i][k] for k in range(n) if k != j] for i in range(1, n)]
            term = entry * cdet(minor)
            total = total + (term if sign_positive else C5(-term.re, -term.im))
        sign_positive = not sign_positive
    return total


def flatten_one(matrix: Sequence[Sequence[C5]]) -> list[F5]:
    out: list[F5] = []
    for row in matrix:
        for entry in row:
            out.append(entry.re)
            out.append(entry.im)
    return out


def span_dimension(matrices: Sequence[Sequence[Sequence[C5]]]) -> int:
    flats = [flatten_one(m) for m in matrices]
    return rank([row[:] for row in flats])


def basis_by_pivots(matrices: Sequence[Sequence[Sequence[C5]]]) -> list[Sequence[Sequence[C5]]]:
    """Matrices that introduce new pivots, in order: an exact basis of the span."""

    basis: list[Sequence[Sequence[C5]]] = []
    kept: list[list[F5]] = []
    for m in matrices:
        candidate = [row[:] for row in kept] + [flatten_one(m)]
        if rank(candidate) > len(kept):
            basis.append(m)
            kept.append(flatten_one(m))
    return basis


def complex_nullspace_dimension(rows: list[list[F5]]) -> int:
    """Complex dimension of the joint kernel encoded as interleaved real rows."""

    filtered = [row for row in rows if any(not entry.is_zero() for entry in row)]
    if not filtered:
        # No constraints: the kernel is everything.
        width = len(rows[0]) if rows else 0
        return width // 2
    real_dimension = len(nullspace(filtered))
    require(real_dimension % 2 == 0, "COMPLEX_STRUCTURE", "complex-linear kernel has odd real dimension")
    return real_dimension // 2


def complex_constraint_rows(matrix: Sequence[Sequence[C5]]) -> list[list[F5]]:
    """Real rows encoding matrix @ v = 0 for v with interleaved (re, im) coords."""

    rows: list[list[F5]] = []
    n = len(matrix)
    m = len(matrix[0]) if n else 0
    for i in range(n):
        real_row: list[F5] = []
        imag_row: list[F5] = []
        for j in range(m):
            entry = matrix[i][j]
            real_row.extend([entry.re, -entry.im])
            imag_row.extend([entry.im, entry.re])
        rows.append(real_row)
        rows.append(imag_row)
    return rows


# ---------------------------------------------------------------------------
# Manifest validation
# ---------------------------------------------------------------------------


def enforce_matter_firewall(manifest: Mapping[str, Any]) -> None:
    e565.enforce_source_firewall(manifest)
    hits: list[str] = []
    for path, text in e565.walk_strings(manifest):
        token = e565.normalized_token(text)
        for forbidden in FORBIDDEN_MATTER_TOKENS:
            if forbidden in token:
                hits.append(f"{path}:{text}")
    if hits:
        raise CertificateError(
            "FORBIDDEN_DEPENDENCY",
            "source manifest contains downstream matter target data: " + "; ".join(hits[:4]),
        )


def parse_rational(value: Any, code: str) -> Fraction:
    try:
        return Fraction(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError(code, f"cannot parse exact rational {value!r}") from exc


SELECTION_RULES = (
    "parity_even_minus_derived_invariants",
    "parity_odd_minus_derived_invariants",
    "lambda2_only",
    "even_including_vacuum",
)

FIELD_LABELS = ("Q", "u_c", "e_c", "d_c", "L")
SCALAR_LABELS = ("S", "Sbar")


def validate_manifest(
    manifest: Mapping[str, Any],
    *,
    allow_control_contracts: bool = False,
) -> dict[str, Any]:
    enforce_matter_firewall(manifest)
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")

    exterior = manifest.get("exterior_matter_contract")
    require(isinstance(exterior, Mapping), "EXTERIOR_CONTRACT", "exterior_matter_contract is missing")
    charges = exterior.get("block_trace_charges")
    require(isinstance(charges, Mapping), "EXTERIOR_CONTRACT", "block_trace_charges is missing")
    y_color = parse_rational(charges.get("color_block"), "EXTERIOR_CONTRACT")
    y_weak = parse_rational(charges.get("weak_block"), "EXTERIOR_CONTRACT")
    # Trace balance is checked arithmetically from the declared charge pair;
    # a redundant declared balance flag is not accepted.
    require(
        3 * y_color + 2 * y_weak == 0,
        "TRACE_BALANCE",
        f"declared block charges are not trace balanced: 3*({y_color}) + 2*({y_weak}) != 0",
    )
    require(
        exterior.get("one_scalar") == "weak_block",
        "EXTERIOR_CONTRACT",
        "the declared one-scalar choice must be the weak block itself",
    )
    channels_raw = exterior.get("yukawa_channels")
    require(isinstance(channels_raw, list), "EXTERIOR_CONTRACT", "yukawa_channels must be a list")
    require(
        len(channels_raw) > 0,
        "GAUSS_DATA_EMPTY",
        "the declared Gauss data is empty: at least one invariant Yukawa channel is required",
    )
    channels: list[tuple[str, str, str]] = []
    for row in channels_raw:
        require(
            isinstance(row, list) and len(row) == 3,
            "EXTERIOR_CONTRACT",
            "each Yukawa channel must list a matter field, a scalar, and a matter field",
        )
        left, scalar, right = (str(x) for x in row)
        require(left in FIELD_LABELS and right in FIELD_LABELS, "EXTERIOR_CONTRACT", f"unknown matter field in channel {row}")
        require(scalar in SCALAR_LABELS, "EXTERIOR_CONTRACT", f"unknown scalar in channel {row}")
        channels.append((left, scalar, right))
    extra_scalars_raw = exterior.get("extra_scalars", [])
    require(isinstance(extra_scalars_raw, list), "EXTERIOR_CONTRACT", "extra_scalars must be a list")
    extra_scalars: list[dict[str, Any]] = []
    for row in extra_scalars_raw:
        require(isinstance(row, Mapping), "EXTERIOR_CONTRACT", "each extra scalar must be an object")
        extra_scalars.append(
            {
                "label": str(row.get("label")),
                "charge": parse_rational(row.get("charge"), "EXTERIOR_CONTRACT"),
                "triality": int(row.get("triality", 0)) % 3,
                "duality": int(row.get("duality", 0)) % 2,
            }
        )

    statistics = manifest.get("statistics_contract")
    require(isinstance(statistics, Mapping), "STATISTICS_TYPING", "statistics_contract is missing")
    matter_statistics = statistics.get("matter_statistics")
    allowed_statistics = (
        ("fermionic_odd", "bosonic_even") if allow_control_contracts else ("fermionic_odd",)
    )
    require(
        matter_statistics in allowed_statistics,
        "STATISTICS_TYPING",
        "the production packet must declare fermionic_odd matter statistics",
    )
    require(
        statistics.get("distinct_from_bosonic_records") is True,
        "STATISTICS_TYPING",
        "the statistics contract must keep matter distinct from bosonic records",
    )

    category = manifest.get("category_contract")
    require(isinstance(category, Mapping), "CATEGORY_TYPING", "category_contract is missing")
    typing = category.get("typing")
    allowed_typings = (
        ("spin_odd_weyl_super", "svec", "vec") if allow_control_contracts else ("spin_odd_weyl_super",)
    )
    require(
        typing in allowed_typings,
        "CATEGORY_TYPING",
        "the production packet must declare the spin_odd_weyl_super category typing",
    )
    spin = category.get("spin_lift")
    require(isinstance(spin, Mapping), "CATEGORY_TYPING", "spin_lift contract is missing")
    double_cover = spin.get("double_cover")
    require(double_cover in (True, False), "CATEGORY_TYPING", "spin_lift.double_cover must be boolean")
    require(
        allow_control_contracts or double_cover is True,
        "CATEGORY_TYPING",
        "the production packet must declare the genuine double-cover spin lift",
    )
    declared_selection_rule = category.get("selection_rule")
    require(
        "projector_representative_convention" not in category,
        "SELECTION_RULE",
        "the source packet may not prefer either member of the derived "
        "charge-conjugate projector pair",
    )
    if allow_control_contracts:
        selection_rule = (
            str(declared_selection_rule)
            if declared_selection_rule is not None
            else "parity_even_minus_derived_invariants"
        )
        require(
            selection_rule in SELECTION_RULES,
            "SELECTION_RULE",
            "unknown control-lane selection rule",
        )
    else:
        require(
            declared_selection_rule is None,
            "SELECTION_RULE",
            "the production packet may not choose a matter projector; the "
            "conjugate rank-fifteen projector pair is derived by the verifier",
        )
        selection_rule = "parity_even_minus_derived_invariants"
    realization = category.get("realization")
    require(
        realization in ("operator_projector", "representation_arithmetic"),
        "CATEGORY_TYPING",
        "realization must be 'operator_projector' or 'representation_arithmetic'",
    )
    require(
        realization == "operator_projector",
        "REALIZATION_NOT_OPERATOR",
        "representation arithmetic alone is not physical realization: the packet requires the operator projector receipt",
    )

    kernel = manifest.get("kernel_emission_contract")
    require(isinstance(kernel, Mapping), "KERNEL_EMISSION_CONTRACT", "kernel_emission_contract is missing")
    require(
        kernel.get("emit") is True,
        "KERNEL_EMISSION_CONTRACT",
        "the common action kernel must be emitted as data",
    )
    require(
        kernel.get("assume_global_quotient") is False,
        "KERNEL_EMISSION_CONTRACT",
        "the packet must not assume the final global quotient",
    )

    mar = manifest.get("candidate_matter_class")
    require(isinstance(mar, Mapping), "MATTER_CLASS_CONTRACT", "candidate_matter_class is missing")
    require(
        mar.get("declared") == "one_generation_one_scalar_chiral_anomaly_free",
        "MATTER_CLASS_CONTRACT",
        "the declared candidate class must be the one-generation one-scalar chiral anomaly-free class",
    )
    require(
        mar.get("promote_uniqueness") is False,
        "CLASS_UNIQUENESS_PROMOTION",
        "class uniqueness may not be promoted before the class is proved nonempty and the uniqueness lane runs",
    )

    return {
        "y_color": y_color,
        "y_weak": y_weak,
        "channels": channels,
        "extra_scalars": extra_scalars,
        "matter_statistics": matter_statistics,
        "category_typing": typing,
        "double_cover": bool(double_cover),
        "selection_rule": selection_rule,
    }


def load_upstream(manifest: Mapping[str, Any], base_dir: Path) -> dict[str, Any]:
    current_path_raw = manifest.get("current_manifest_path")
    require(isinstance(current_path_raw, str), "UPSTREAM_REFERENCE", "current_manifest_path is missing")
    current_path = Path(current_path_raw)
    if not current_path.is_absolute():
        current_path = base_dir / current_path
    current_manifest = load_json(current_path)
    digest = sha256_json(current_manifest)
    require(
        manifest.get("current_manifest_sha256") == digest,
        "UPSTREAM_HASH",
        "the #566 current manifest hash does not match the declared pin",
    )

    receipt_path_raw = manifest.get("current_receipt_path")
    require(isinstance(receipt_path_raw, str), "UPSTREAM_REFERENCE", "current_receipt_path is missing")
    receipt_path = Path(receipt_path_raw)
    if not receipt_path.is_absolute():
        receipt_path = base_dir / receipt_path
    current_receipt = load_json(receipt_path)
    require(
        manifest.get("current_receipt_sha256") == sha256_json(current_receipt),
        "UPSTREAM_HASH",
        "the #566 current receipt hash does not match the declared pin",
    )
    require(
        current_receipt.get("schema") == p566.RECEIPT_SCHEMA and current_receipt.get("issue") == 566,
        "UPSTREAM_RECEIPT",
        "the pinned upstream receipt is not a #566 port-current receipt",
    )
    require(
        current_receipt.get("manifest_sha256") == digest,
        "UPSTREAM_RECEIPT",
        "the pinned upstream receipt does not certify the pinned current manifest",
    )
    gate = current_receipt.get("conditional_algebraic_gate")
    require(
        isinstance(gate, Mapping) and gate.get("passed") is True,
        "UPSTREAM_RECEIPT",
        "the pinned upstream receipt did not pass the conditional algebraic gate",
    )
    physical_gate = current_receipt.get("physical_source_gate")
    require(
        isinstance(physical_gate, Mapping)
        and physical_gate.get("passed") is True
        and physical_gate.get("target_blind_impulse_readback_recomputed") is True,
        "UPSTREAM_RECEIPT",
        "the upstream receipt must pass the target-blind impulse/readback "
        "source gate",
    )
    binding = current_receipt.get("semantic_response_binding")
    require(
        isinstance(binding, Mapping)
        and binding.get("sector_structure_recomputed") is True,
        "UPSTREAM_RECEIPT",
        "the upstream physical gate must be justified by a recomputed semantic "
        "response binding",
    )
    return {
        "current_manifest": current_manifest,
        "current_manifest_sha256": digest,
        "current_receipt_sha256": sha256_json(current_receipt),
        "carrier_manifest_sha256": current_receipt.get("carrier_manifest_sha256"),
        "semantic_response_artifact_sha256": binding.get("artifact_sha256"),
        "upstream_physical_source_gate_passed": physical_gate.get("passed") is True,
    }


SPIN_ARTIFACT_SCHEMA = "oph.spin_statistics_semantic_artifact.v1"
BINARY_ICOSAHEDRAL_PROFILE = {"1": 1, "2": 1, "3": 20, "4": 30, "5": 24, "6": 20, "10": 24}


def load_spin_statistics_artifact(
    manifest: Mapping[str, Any], base_dir: Path, upstream: Mapping[str, Any]
) -> dict[str, Any]:
    """Load and verify the hash-pinned measured spin-transport artifact.

    The artifact is produced target-blind by the simulator from the certified
    carrier: it measures the deck group, the exact quaternion lift closure
    with its order profile, unique involution, and two-element centre, the
    section obstruction over every Klein four-subgroup, the unique spin
    structure on the oriented support, and the orientation convention. This
    loader verifies the pin, the self-hash, the schema, and that the artifact
    binds the same certified carrier as the upstream response packet.
    """

    path_raw = manifest.get("spin_statistics_artifact_path")
    require(isinstance(path_raw, str), "SPIN_ARTIFACT", "spin_statistics_artifact_path is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    artifact = load_json(path)
    declared = manifest.get("spin_statistics_artifact_sha256")
    require(
        isinstance(declared, str) and declared == artifact.get("artifact_sha256"),
        "UPSTREAM_HASH",
        "the spin statistics artifact hash does not match the declared pin",
    )
    body = {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    require(
        artifact.get("artifact_sha256") == "sha256:" + sha256_json(body),
        "SPIN_ARTIFACT",
        "the spin statistics artifact self-hash does not recompute",
    )
    require(
        artifact.get("schema") == SPIN_ARTIFACT_SCHEMA and artifact.get("issue") == 314,
        "SPIN_ARTIFACT",
        "the pinned artifact is not a #314 spin statistics artifact",
    )
    require(
        artifact.get("carrier_binding", {}).get("carrier_manifest_sha256")
        == upstream["carrier_manifest_sha256"],
        "SPIN_ARTIFACT",
        "the spin artifact does not bind the same certified carrier as the response packet",
    )
    lift = artifact.get("lift_measurement", {})
    require(
        lift.get("lift_group_order") == 120
        and lift.get("order_profile") == BINARY_ICOSAHEDRAL_PROFILE
        and lift.get("unique_nontrivial_involution") == "-1"
        and lift.get("centre_order") == 2,
        "SPIN_ARTIFACT",
        "the measured lift group is not the binary icosahedral transport class",
    )
    obstruction = artifact.get("section_obstruction", {})
    require(
        obstruction.get("deck_involutions") == 15
        and obstruction.get("klein_four_subgroups") == 5
        and obstruction.get("no_section_over_any_klein_four_subgroup") is True,
        "SPIN_ARTIFACT",
        "the measured section obstruction is not total over the Klein four-subgroups",
    )
    homology = artifact.get("support_homology", {})
    require(
        homology.get("betti_numbers") == [1, 0, 1]
        and homology.get("spin_structure_count") == 1,
        "SPIN_ARTIFACT",
        "the oriented support does not carry a unique spin structure",
    )
    orientation = artifact.get("orientation_convention", {})
    require(
        orientation.get("rotations_preserve_oriented_faces") is True
        and orientation.get("improper_coset_reverses_oriented_faces") is True,
        "SPIN_ARTIFACT",
        "the measured orientation convention is incomplete",
    )
    gate = artifact.get("physical_source_gate", {})
    require(
        gate.get("passed") is True,
        "SPIN_ARTIFACT",
        "the spin artifact's own measured source gate does not pass",
    )
    return {
        "artifact_sha256": artifact["artifact_sha256"],
        "order_profile": lift["order_profile"],
        "centre_order": lift["centre_order"],
        "klein_four_subgroups": obstruction["klein_four_subgroups"],
        "no_section": obstruction["no_section_over_any_klein_four_subgroup"],
        "spin_structure_count": homology["spin_structure_count"],
        "orientation_convention": orientation,
        "laboratory_exchange_measurement": gate.get("laboratory_exchange_measurement"),
    }


def _scan_components(a: int, b: int) -> list[dict[str, Any]]:
    """The ten nontrivial isotypic components of the exterior module over the
    derived integer block charges (a, b), with conjugate pairing data."""

    return [
        {"color": "3", "weak": 1, "q": a, "parity": 1, "pair": 0, "side": 0},
        {"color": "1", "weak": 2, "q": b, "parity": 1, "pair": 1, "side": 0},
        {"color": "3bar", "weak": 1, "q": 2 * a, "parity": 0, "pair": 2, "side": 0},
        {"color": "3", "weak": 2, "q": a + b, "parity": 0, "pair": 3, "side": 0},
        {"color": "1", "weak": 1, "q": 2 * b, "parity": 0, "pair": 4, "side": 0},
        {"color": "3", "weak": 1, "q": -2 * a, "parity": 1, "pair": 2, "side": 1},
        {"color": "3bar", "weak": 2, "q": -a - b, "parity": 1, "pair": 3, "side": 1},
        {"color": "1", "weak": 1, "q": -2 * b, "parity": 1, "pair": 4, "side": 1},
        {"color": "3bar", "weak": 1, "q": -a, "parity": 0, "pair": 0, "side": 1},
        {"color": "1", "weak": 2, "q": -b, "parity": 0, "pair": 1, "side": 1},
    ]


def _scan_anomalies(selection: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    def color_dim(component: Mapping[str, Any]) -> int:
        return 3 if component["color"] in ("3", "3bar") else 1

    grav = sum(color_dim(c) * c["weak"] * c["q"] for c in selection)
    su3 = sum(c["weak"] * c["q"] for c in selection if c["color"] in ("3", "3bar"))
    su2 = sum(color_dim(c) * c["q"] for c in selection if c["weak"] == 2)
    u1_cubed = sum(color_dim(c) * c["weak"] * c["q"] ** 3 for c in selection)
    witten = sum(color_dim(c) for c in selection if c["weak"] == 2) % 2
    return {"grav": grav, "su3": su3, "su2": su2, "u1_cubed": u1_cubed, "witten": witten}


def exterior_selection_scan(y_color: Fraction, y_weak: Fraction) -> dict[str, Any]:
    """Exhaustive matter-selection theorem on the exterior module.

    Over the derived primitive integer block charges, every one of the 1024
    subsets of the ten nontrivial isotypic components is classified. Exactly
    two nonempty chiral subsets are anomaly free (gravitational, SU(3)^2 U(1),
    SU(2)^2 U(1), and U(1)^3), they are exchanged by charge conjugation, and
    each is exactly one fermionic-parity sector minus its invariant line. The
    parity grading of the matter object is therefore an output of the scan,
    not a declared statistics contract, and the unordered conjugate pair is
    selected without any rank, complementarity, or faithfulness assumption;
    Witten parity and faithfulness hold automatically on both survivors.
    """

    normalization = 6
    a = int(y_color * normalization)
    b = int(y_weak * normalization)
    require(
        Fraction(a, normalization) == y_color and Fraction(b, normalization) == y_weak,
        "SELECTION_SCAN",
        "the derived block charges are not sixths-integral",
    )
    require(3 * a + 2 * b == 0, "SELECTION_SCAN", "the derived block charges are not balanced")
    components = _scan_components(a, b)

    def survivors(table: Sequence[Mapping[str, Any]]) -> list[tuple[int, ...]]:
        passing: list[tuple[int, ...]] = []
        for mask in range(1, 1 << len(table)):
            selection = [table[i] for i in range(len(table)) if mask >> i & 1]
            picked_pairs = {(c["pair"], c["side"]) for c in selection}
            if any((p, 0) in picked_pairs and (p, 1) in picked_pairs for p in range(5)):
                continue
            anomalies = _scan_anomalies(selection)
            if (
                anomalies["grav"] == 0
                and anomalies["su3"] == 0
                and anomalies["su2"] == 0
                and anomalies["u1_cubed"] == 0
            ):
                passing.append(tuple(sorted(i for i in range(len(table)) if mask >> i & 1)))
        return passing

    passing = survivors(components)
    require(
        len(passing) == 2,
        "SELECTION_SCAN",
        f"expected exactly two anomaly-free chiral subsets, got {len(passing)}",
    )
    parity_sectors = {
        0: tuple(sorted(i for i, c in enumerate(components) if c["parity"] == 0)),
        1: tuple(sorted(i for i, c in enumerate(components) if c["parity"] == 1)),
    }
    require(
        set(passing) == set(parity_sectors.values()),
        "SELECTION_SCAN",
        "the anomaly-free pair is not the fermionic-parity pair",
    )
    conjugate_of = {
        i: next(
            j
            for j, other in enumerate(components)
            if other["pair"] == c["pair"] and other["side"] == 1 - c["side"]
        )
        for i, c in enumerate(components)
    }
    require(
        tuple(sorted(conjugate_of[i] for i in passing[0])) == passing[1],
        "SELECTION_SCAN",
        "the two survivors are not exchanged by charge conjugation",
    )
    for survivor in passing:
        selection = [components[i] for i in survivor]
        anomalies = _scan_anomalies(selection)
        require(anomalies["witten"] == 0, "SELECTION_SCAN", "a survivor has odd Witten parity")
        require(
            any(c["color"] in ("3", "3bar") for c in selection)
            and any(c["weak"] == 2 for c in selection)
            and any(c["q"] != 0 for c in selection),
            "SELECTION_SCAN",
            "a survivor does not carry a faithful current action",
        )
    # Charge conjugation of the input charges selects the same unordered pair.
    conjugated = survivors(_scan_components(-a, -b))
    require(
        set(conjugated) == set(passing),
        "SELECTION_SCAN",
        "the scan is not invariant under charge conjugation of the derived charges",
    )
    return {
        "subsets_enumerated": 1024,
        "constraints": [
            "chirality (no component together with its conjugate)",
            "gravitational anomaly",
            "SU(3)^2 U(1)",
            "SU(2)^2 U(1)",
            "U(1)^3",
        ],
        "survivor_count": 2,
        "survivors_are_conjugate_pair": True,
        "survivors_equal_parity_sectors": True,
        "witten_and_faithfulness_automatic": True,
        "charge_conjugation_invariant": True,
        "derived_block_charges": {"a": a, "b": b, "normalization": normalization},
        "conclusion": (
            "the unordered conjugate pair {even parity minus vacuum, odd parity "
            "minus top} is the unique nonempty chiral anomaly-free selection of "
            "the exterior module; the fermionic-parity grading is derived by the "
            "scan rather than declared"
        ),
    }


def statistics_forcing_certificate(
    scan: Mapping[str, Any],
    spin: Mapping[str, Any],
    spin_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    """Force the physical statistics and category typing from measured data.

    The derived grading (from the exhaustive selection scan) must be
    implemented by measured central structure to count as source-produced.
    The implementation menu is enumerated exhaustively:

    - gauge centre: no central parameter among the 36 candidates acts by -1
      on every component of the selected matter pair (recomputed here and
      machine-checked in Lean/Screen/Z6Descent.lean as
      no_universal_fermion_minus_one), so the gauge centre cannot implement
      the grading;
    - measured frame-transport centre: the artifact measures the lift centre
      to be exactly {+1, -1}, so the unique nontrivial implementation is the
      central -1 of the measured double cover;
    - a split abstract grading element would need a deck-side involution
      lifting with square +1, which the measured section obstruction refutes
      over every Klein four-subgroup.

    A Vec typing implements the grading trivially and contradicts its
    derived nontriviality; a split sVec typing is not source-realized; the
    opposite-Weyl relabeling reverses the measured orientation convention
    and is a same-reduct symmetry only when composed with charge
    conjugation, which is already recorded as the unordered pair symmetry.
    The Spin/odd-Weyl super typing is therefore forced on the realized
    module at finite source-model scope.
    """

    a = scan["derived_block_charges"]["a"]
    b = scan["derived_block_charges"]["b"]
    components = _scan_components(a, b)
    parity_one = [c for c in components if c["parity"] == 1]
    parity_zero = [c for c in components if c["parity"] == 0]

    def triality(component: Mapping[str, Any]) -> int:
        if component["color"] == "3":
            return 1
        if component["color"] == "3bar":
            return 2
        return 0

    def duality(component: Mapping[str, Any]) -> int:
        return 1 if component["weak"] == 2 else 0

    gauge_minus_one_candidates = []
    for k in range(3):
        for l in range(2):
            for r in range(6):
                for sector in (parity_one, parity_zero):
                    if all(
                        (2 * k * triality(c) + 3 * l * duality(c) + r * (c["q"] % 6)) % 6 == 3
                        for c in sector
                    ):
                        gauge_minus_one_candidates.append((k, l, r))
    require(
        not gauge_minus_one_candidates,
        "GRADING_IMPLEMENTATION",
        "a gauge central element implements fermion parity, contradicting the Lean no-go",
    )
    require(
        spin_artifact["centre_order"] == 2 and spin.get("centre_order") == 2,
        "GRADING_IMPLEMENTATION",
        "the measured and recomputed lift centres must both be {+1, -1}",
    )
    require(
        spin_artifact["no_section"] is True and spin.get("involution_lift_order") == 4,
        "GRADING_IMPLEMENTATION",
        "the measured double cover must be non-split with involution lifts of order four",
    )
    return {
        "derived_grading": (
            "the exhaustive selection scan outputs the fermionic-parity pair, so "
            "the matter grading is derived, nontrivial, and conjugation-symmetric"
        ),
        "implementation_menu": {
            "gauge_centre": {
                "candidates_enumerated": 36,
                "acts_as_minus_one_on_selected_matter": 0,
                "excluded": True,
                "lean": "Z6Descent.no_universal_fermion_minus_one",
            },
            "measured_frame_transport_centre": {
                "centre": ["+1", "-1"],
                "unique_nontrivial_implementation": "-1 of the measured double cover",
                "selected": True,
            },
            "split_abstract_grading": {
                "requires": "a deck-side involution lifting with square +1",
                "refuted_by": "measured section obstruction over all five Klein four-subgroups",
                "excluded": True,
            },
        },
        "typing_controls": {
            "vec": "fails: trivial implementation contradicts the derived nontrivial grading",
            "svec_split": "fails: no source-realized independent central involution exists",
            "opposite_weyl": (
                "fails alone: reverses the measured orientation convention; composed "
                "with charge conjugation it is the recorded unordered-pair symmetry"
            ),
        },
        "forced_typing": "spin_odd_weyl_super",
        "scope": (
            "finite source-model scope: the forcing enumerates implementations by "
            "measured central structure; a continuum spin-statistics theorem and "
            "laboratory exchange measurement remain separate lanes"
        ),
    }


# ---------------------------------------------------------------------------
# The upstream current algebra, rebuilt from the pinned #566 source packet
# ---------------------------------------------------------------------------


class CurrentAlgebra:
    """The derived #566 current algebra: generators, rotations, structure."""

    def __init__(self, current_manifest: Mapping[str, Any], base_dir: Path) -> None:
        params = p566.validate_manifest(current_manifest, base_dir)
        require(
            params["model"] == "charged_double_triplet",
            "UPSTREAM_MODEL",
            "the matter lift requires the charged double-triplet response model",
        )
        carrier, _group_row, plus, carrier_manifest = p566.load_carrier(current_manifest, base_dir)
        verts = p566.standard_vertices()
        matched = p566.orientation_matched_assignments(carrier, verts)
        frame = p566.FrameRealization(carrier, matched[0], verts)
        model = p566.ChargedDoubleTripletModel(frame, params)

        self.carrier = carrier
        self.carrier_manifest = carrier_manifest
        self.frame = frame
        self.model = model
        self.plus = [tuple(g) for g in plus]
        self.generators = [model.generator(field) for field in p566.BASIS_FIELDS]
        p566.check_skew_adjoint(self.generators)
        self.rotations = {g: frame.rotation_of(g) for g in self.plus}
        # The kernel block transforms in the Galois-conjugate (3') action:
        # w_{g(p)} = conj(R_g) w_p, exactly as in the #566 covariance.
        self.kernel_rotations = {
            g: [[entry.conj() for entry in row] for row in rotation]
            for g, rotation in self.rotations.items()
        }
        self.orders = {g: e565.permutation_order(g) for g in self.plus}

        flats = [p566.flatten(blocks) for blocks in self.generators]
        require(rank([row[:] for row in flats]) == 12, "UPSTREAM_RANK", "current algebra image is not twelve dimensional")
        self.flats = flats

        self.structure_constants: dict[tuple[int, int], list[F5]] = {}
        for i in range(12):
            for j in range(i + 1, 12):
                bracket = tuple(
                    commutator(self.generators[i][b], self.generators[j][b])
                    for b in range(2)
                )
                coeffs = p566.solve_in_span(flats, p566.flatten(bracket))
                self.structure_constants[(i, j)] = coeffs

    def even_block(self, index: int) -> list[list[C5]]:
        return self.generators[index][0]

    def kernel_axis(self, index: int) -> tuple[F5, F5, F5]:
        """The kernel block is hat(w); recover w exactly from the skew entries."""

        block = self.generators[index][1]
        return (block[2][1].re, block[0][2].re, block[1][0].re)


# ---------------------------------------------------------------------------
# PORT-SPIN-LIFT: exact SU(2) lifts of the sixty proper implementers
# ---------------------------------------------------------------------------


PAULI = (
    ((IZERO, IONE), (IONE, IZERO)),
    ((IZERO, C5(ZERO, -ONE)), (C5(ZERO, ONE), IZERO)),
    ((IONE, IZERO), (IZERO, C5(-ONE, ZERO))),
)


def sigma_dot(w: Sequence[F5]) -> list[list[C5]]:
    out = czeros(2)
    for k in range(3):
        coeff = C5(w[k], ZERO)
        if w[k].is_zero():
            continue
        for i in range(2):
            for j in range(2):
                out[i][j] = out[i][j] + coeff * PAULI[k][i][j]
    return out


def rotation_axis_cosine(rotation: Sequence[Sequence[F5]]) -> tuple[list[F5] | None, F5]:
    """Exact (axis, cosine) of a proper rotation; axis is None for the identity."""

    identity = [[ONE if a == b else ZERO for b in range(3)] for a in range(3)]
    difference = [[rotation[i][j] - identity[i][j] for j in range(3)] for i in range(3)]
    trace = rotation[0][0] + rotation[1][1] + rotation[2][2]
    cosine = (trace - ONE) / F5(2)
    if all(entry.is_zero() for row in difference for entry in row):
        return None, F5(1)
    kernel = nullspace(difference)
    require(len(kernel) == 1, "SPIN_LIFT_AXIS", "rotation axis is not one-dimensional")
    return kernel[0], cosine


def spin_lift_of_rotation(rotation: Sequence[Sequence[F5]]) -> list[list[C5]]:
    """An exact U in SU(2) with U (sigma . v) U^dagger = sigma . (R v)."""

    axis, cosine = rotation_axis_cosine(rotation)
    if axis is None:
        return cidentity(2)
    half_sq = (ONE + cosine) / F5(2)
    u = sqrt_f5(half_sq)
    require(u is not None, "SPIN_LIFT_FIELD", "half-angle cosine has no square root in Q(sqrt5)")
    eta = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]
    require(not eta.is_zero(), "SPIN_LIFT_AXIS", "axis vector is null")
    c_sq = (ONE - u * u) / eta
    c = sqrt_f5(c_sq)
    require(c is not None, "SPIN_LIFT_FIELD", "scaled half-angle sine has no square root in Q(sqrt5)")
    for candidate_c in (c, -c):
        lift = [
            [
                C5(u, -(candidate_c * axis[2])),
                C5(-(candidate_c * axis[1]), -(candidate_c * axis[0])),
            ],
            [
                C5(candidate_c * axis[1], -(candidate_c * axis[0])),
                C5(u, candidate_c * axis[2]),
            ],
        ]
        if spin_lift_matches(lift, rotation):
            return lift
    raise CertificateError("SPIN_LIFT_ADJOINT", "no exact SU(2) lift reproduces the rotation")


def spin_lift_matches(lift: Sequence[Sequence[C5]], rotation: Sequence[Sequence[F5]]) -> bool:
    dagger = cdagger(lift)
    for k in range(3):
        conjugated = cmul(cmul([list(row) for row in lift], [list(row) for row in PAULI[k]]), dagger)
        expected = czeros(2)
        for j in range(3):
            if rotation[j][k].is_zero():
                continue
            coeff = C5(rotation[j][k], ZERO)
            for a in range(2):
                for b in range(2):
                    expected[a][b] = expected[a][b] + coeff * PAULI[j][a][b]
        if not c_is_zero(csub(conjugated, expected)):
            return False
    return True


def matrix_key(matrix: Sequence[Sequence[C5]]) -> tuple:
    return tuple(
        (entry.re.a, entry.re.b, entry.im.a, entry.im.b)
        for row in matrix
        for entry in row
    )


def spin_lift_certificate(algebra: CurrentAlgebra) -> dict[str, Any]:
    lifts: dict[tuple[int, ...], list[list[C5]]] = {}
    for g in algebra.plus:
        lift = spin_lift_of_rotation(algebra.kernel_rotations[g])
        product = cmul(cdagger(lift), lift)
        require(c_is_zero(csub(product, cidentity(2))), "SPIN_LIFT_UNITARY", "spin lift is not unitary")
        determinant = cdet(lift)
        require(
            (determinant - IONE).is_zero(),
            "SPIN_LIFT_SPECIAL",
            "spin lift determinant is not one",
        )
        lifts[g] = lift

    minus_identity = cscale(cidentity(2), C5(-ONE, ZERO))
    group: dict[tuple, list[list[C5]]] = {}
    for lift in lifts.values():
        group[matrix_key(lift)] = lift
        negated = cscale(lift, C5(-ONE, ZERO))
        group[matrix_key(negated)] = negated
    require(len(group) == 120, "SPIN_LIFT_COUNT", f"expected 120 lift elements, got {len(group)}")

    elements = list(group.values())
    keys = set(group.keys())
    for x in elements:
        for y in elements:
            require(
                matrix_key(cmul(x, y)) in keys,
                "SPIN_LIFT_CLOSURE",
                "the lift set is not closed under multiplication",
            )

    def element_order(matrix: list[list[C5]]) -> int:
        power = matrix
        identity = cidentity(2)
        for order in range(1, 121):
            if c_is_zero(csub(power, identity)):
                return order
            power = cmul(power, matrix)
        raise CertificateError("SPIN_LIFT_ORDER", "lift element order exceeded finite bound")

    order_profile: dict[int, int] = {}
    involutions = 0
    for matrix in elements:
        order = element_order(matrix)
        order_profile[order] = order_profile.get(order, 0) + 1
        if order == 2:
            involutions += 1
            require(
                c_is_zero(csub(matrix, minus_identity)),
                "SPIN_LIFT_SPLIT",
                "an involution other than -1 exists, so the extension would split",
            )
    require(
        involutions == 1,
        "SPIN_LIFT_SPLIT",
        f"the lift group must have a unique involution, got {involutions}",
    )
    expected_profile = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}
    require(
        order_profile == expected_profile,
        "SPIN_LIFT_PROFILE",
        f"expected the binary icosahedral order profile {expected_profile}, got {order_profile}",
    )

    for g in algebra.plus:
        if algebra.orders[g] == 2:
            square = cmul(lifts[g], lifts[g])
            require(
                c_is_zero(csub(square, minus_identity)),
                "SPIN_LIFT_SPLIT",
                "the lift of an involution must have order four",
            )

    irrational_spinor_traces = 0
    for g in algebra.plus:
        if algebra.orders[g] == 5:
            trace = lifts[g][0][0] + lifts[g][1][1]
            require(
                trace.im.is_zero() and trace.re.b != 0,
                "SPIN_LIFT_RELABELING",
                "an order-five spinor character is not irrational",
            )
            irrational_spinor_traces += 1

    centre = [
        matrix
        for matrix in elements
        if all(
            c_is_zero(csub(cmul(matrix, other), cmul(other, matrix)))
            for other in elements
        )
    ]
    require(
        len(centre) == 2
        and any(c_is_zero(csub(matrix, cidentity(2))) for matrix in centre)
        and any(c_is_zero(csub(matrix, minus_identity)) for matrix in centre),
        "SPIN_LIFT_CENTRE",
        "the lift group centre must be exactly {+1, -1}",
    )

    return {
        "lifts": lifts,
        "witness_count": len(lifts),
        "lift_group_order": len(group),
        "unique_involution": True,
        "order_profile": {str(k): v for k, v in sorted(order_profile.items())},
        "involution_lift_order": 4,
        "centre_order": len(centre),
        "irrational_order_five_spinor_traces": irrational_spinor_traces,
        "conclusion": (
            "the sixty proper implementers lift to SU(2) with a unique involution -1 in the "
            "120-element lift group, so the double cover is non-split (binary icosahedral) and "
            "PORT-SPIN-LIFT is realized on this branch; irrational spinor characters exclude any "
            "signed register-relabeling realization of the spinor sector"
        ),
    }


# ---------------------------------------------------------------------------
# Matter carrier V = C (+) W and the current transport
# ---------------------------------------------------------------------------


def matter_transport(
    algebra: CurrentAlgebra,
    y_color: Fraction,
    y_weak: Fraction,
) -> dict[str, Any]:
    """T(K_p) on V = C^3 (+) C^2: source even block, spin-lifted kernel block,
    declared trace-balanced redistribution of the central charge."""

    ratio = F5(y_weak / y_color) if y_color != 0 else None
    half = F5(Fraction(1, 2))
    images: list[list[list[C5]]] = []
    for p in range(12):
        even = algebra.even_block(p)
        trace = ctrace(even)
        require(trace.re.is_zero(), "MATTER_TRANSPORT", "even-block trace is not purely imaginary")
        tau_third = trace.im / F5(3)
        w = algebra.kernel_axis(p)
        sigma_w = sigma_dot(w)
        image = [[IZERO for _ in range(5)] for _ in range(5)]
        for i in range(3):
            for j in range(3):
                image[i][j] = even[i][j]
            if ratio is None:
                # A charge-dead contract annihilates the central lane.
                image[i][i] = image[i][i] - C5(ZERO, tau_third)
        for i in range(2):
            for j in range(2):
                image[3 + i][3 + j] = C5(ZERO, -half) * sigma_w[i][j]
            if ratio is not None:
                image[3 + i][3 + i] = image[3 + i][3 + i] + C5(ZERO, tau_third * ratio)
        images.append(image)

    for image in images:
        for i in range(5):
            for j in range(5):
                require(
                    (image[j][i].conj() + image[i][j]).is_zero(),
                    "MATTER_TRANSPORT",
                    "a transported generator is not skew-adjoint",
                )

    faithful_rank = span_dimension(images)
    require(
        faithful_rank == 12,
        "CURRENT_ACTION_NOT_FAITHFUL",
        f"the current algebra must act faithfully on the matter carrier, got rank {faithful_rank}",
    )

    homomorphism_checks = 0
    for (i, j), coeffs in algebra.structure_constants.items():
        bracket = commutator(images[i], images[j])
        expected = [[IZERO for _ in range(5)] for _ in range(5)]
        for r in range(12):
            if coeffs[r].is_zero():
                continue
            coeff = C5(coeffs[r], ZERO)
            for a in range(5):
                for b in range(5):
                    expected[a][b] = expected[a][b] + coeff * images[r][a][b]
        require(
            c_is_zero(csub(bracket, expected)),
            "MATTER_TRANSPORT",
            "the matter transport is not a Lie algebra homomorphism",
        )
        homomorphism_checks += 1

    return {
        "images": images,
        "faithful_rank_on_carrier": faithful_rank,
        "homomorphism_bracket_checks": homomorphism_checks,
    }


def carrier_conjugation_checks(
    algebra: CurrentAlgebra,
    spin_lifts: Mapping[tuple[int, ...], list[list[C5]]],
    images: Sequence[Sequence[Sequence[C5]]],
) -> int:
    """pi_V conjugation transports T(K_p) to T(K_{g(p)}) for every g and p."""

    checked = 0
    for g in algebra.plus:
        pi_v = carrier_implementer(algebra.rotations[g], spin_lifts[g])
        pi_dagger = cdagger(pi_v)
        for p in range(12):
            conjugated = cmul(cmul(pi_v, [list(row) for row in images[p]]), pi_dagger)
            require(
                c_is_zero(csub(conjugated, [list(row) for row in images[g[p]]])),
                "MATTER_COVARIANCE",
                "pi_V conjugation does not transport the matter generators covariantly",
            )
            checked += 1
    return checked


def carrier_implementer(rotation: Sequence[Sequence[F5]], lift: Sequence[Sequence[C5]]) -> list[list[C5]]:
    pi_v = [[IZERO for _ in range(5)] for _ in range(5)]
    for i in range(3):
        for j in range(3):
            pi_v[i][j] = C5(rotation[i][j], ZERO)
    for i in range(2):
        for j in range(2):
            pi_v[3 + i][3 + j] = lift[i][j]
    return pi_v


# ---------------------------------------------------------------------------
# Auxiliary CAR/Fock space over V
# ---------------------------------------------------------------------------


class FockSpace:
    """Lambda^* C^5 with exact CAR operators over Q(sqrt5)."""

    def __init__(self) -> None:
        self.subsets: list[tuple[int, ...]] = []
        for mask in range(32):
            self.subsets.append(tuple(i for i in range(5) if mask & (1 << i)))
        self.subsets.sort(key=lambda s: (len(s), s))
        self.index = {s: n for n, s in enumerate(self.subsets)}
        self.dim = 32

    def creation(self, mode: int) -> list[list[C5]]:
        out = czeros(self.dim)
        for s, subset in enumerate(self.subsets):
            if mode in subset:
                continue
            sign = (-1) ** sum(1 for j in subset if j < mode)
            new = tuple(sorted(subset + (mode,)))
            out[self.index[new]][s] = C5(F5(sign), ZERO)
        return out

    def dgamma(self, one_body: Sequence[Sequence[C5]]) -> list[list[C5]]:
        """Second quantization: dGamma(X) = sum_ij X_ij a_i^dagger a_j."""

        out = czeros(self.dim)
        for s, subset in enumerate(self.subsets):
            for j in subset:
                sign_j = (-1) ** sum(1 for m in subset if m < j)
                removed = tuple(m for m in subset if m != j)
                for i in range(5):
                    entry = one_body[i][j]
                    if entry.is_zero() or i in removed:
                        continue
                    sign_i = (-1) ** sum(1 for m in removed if m < i)
                    new = tuple(sorted(removed + (i,)))
                    coeff = C5(F5(sign_i * sign_j), ZERO) * entry
                    out[self.index[new]][s] = out[self.index[new]][s] + coeff
        return out

    def parity(self) -> list[list[C5]]:
        out = czeros(self.dim)
        for s, subset in enumerate(self.subsets):
            out[s][s] = C5(F5((-1) ** len(subset)), ZERO)
        return out

    def exterior_lift(self, one_body: Sequence[Sequence[C5]]) -> list[list[C5]]:
        """Gamma(g): the exterior-power (multiplicative) lift of a 5x5 matrix."""

        out = czeros(self.dim)
        for s, source in enumerate(self.subsets):
            k = len(source)
            for t, target in enumerate(self.subsets):
                if len(target) != k:
                    continue
                minor = [[one_body[r][c] for c in source] for r in target]
                out[t][s] = cdet(minor)
        return out


def car_certificate(fock: FockSpace) -> dict[str, Any]:
    creations = [fock.creation(i) for i in range(5)]
    annihilations = [cdagger(a) for a in creations]
    identity = cidentity(fock.dim)
    zero = czeros(fock.dim)
    checks = 0
    for i in range(5):
        for j in range(5):
            anti_cc = cadd(cmul(creations[i], creations[j]), cmul(creations[j], creations[i]))
            require(c_is_zero(anti_cc), "CAR_RELATIONS", "creation operators do not anticommute")
            anti_ca = cadd(cmul(annihilations[i], creations[j]), cmul(creations[j], annihilations[i]))
            expected = identity if i == j else zero
            require(
                c_is_zero(csub(anti_ca, expected)),
                "CAR_RELATIONS",
                "canonical anticommutation relations fail",
            )
            checks += 2
    # The Fock space is generated from the vacuum by the creation operators.
    vacuum_index = fock.index[()]
    generated: list[list[F5]] = []
    for subset in fock.subsets:
        vector = [IZERO] * fock.dim
        vector[vacuum_index] = IONE
        for mode in reversed(subset):
            matrix = creations[mode]
            new_vector = [IZERO] * fock.dim
            for col, amp in enumerate(vector):
                if amp.is_zero():
                    continue
                for row_index in range(fock.dim):
                    entry = matrix[row_index][col]
                    if not entry.is_zero():
                        new_vector[row_index] = new_vector[row_index] + entry * amp
            vector = new_vector
        flat: list[F5] = []
        for amp in vector:
            flat.append(amp.re)
            flat.append(amp.im)
        generated.append(flat)
    require(
        rank([row[:] for row in generated]) == 32,
        "CAR_CYCLIC",
        "the Fock space is not generated from the vacuum",
    )
    return {"creations": creations, "car_checks": checks, "cyclic_rank": 32}


# ---------------------------------------------------------------------------
# Certificate payload
# ---------------------------------------------------------------------------


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


# ---------------------------------------------------------------------------
# BLOCK-DETERMINANT-BALANCE and scalar/channel selection producers
# ---------------------------------------------------------------------------

# General-charge exterior package: field label -> (dimension, charge linear
# form (coeff_a, coeff_b) in the block charges a = y_color, b = y_weak,
# su(3) fundamental multiplicity, su(2) doublet multiplicity, triality,
# duality). The structure matches the projector-realized package verified
# operator-by-operator later in the payload.
GENERAL_PACKAGE: dict[str, dict[str, Any]] = {
    "Q": {"dim": 6, "charge": (1, 1), "su3_fund": 2, "su2_doublets": 3, "triality": 1, "duality": 1},
    "u_c": {"dim": 3, "charge": (2, 0), "su3_fund": 1, "su2_doublets": 0, "triality": 2, "duality": 0},
    "e_c": {"dim": 1, "charge": (0, 2), "su3_fund": 0, "su2_doublets": 0, "triality": 0, "duality": 0},
    "d_c": {"dim": 3, "charge": (2, 2), "su3_fund": 1, "su2_doublets": 0, "triality": 2, "duality": 0},
    "L": {"dim": 2, "charge": (3, 1), "su3_fund": 0, "su2_doublets": 1, "triality": 0, "duality": 1},
}


def _charge_at(field: str, a: Fraction, b: Fraction) -> Fraction:
    ca, cb = GENERAL_PACKAGE[field]["charge"]
    return ca * a + cb * b


def block_determinant_balance_certificate(
    y_color: Fraction, y_weak: Fraction
) -> dict[str, Any]:
    """Derive the trace-balanced charge pair from the realized package.

    On the realized exterior package with general block charges (a, b), the
    gravitational, SU(3)^2 U(1), and SU(2)^2 U(1) anomaly forms are exact
    linear polynomials proportional to 3a + 2b, and the U(1)^3 form vanishes
    identically on the balance line and only there among the sampled rays.
    Anomaly freedom of the realized package therefore forces the determinant
    balance. Primitive integral block characters solve
    ``3 q_C + 2 q_W = 0`` only as ``+/-(-2, 3)``. Writing ``Y=q/6`` therefore
    gives ``+/-(-1/3, 1/2)``. The common sign is charge conjugation and is not
    selected by the finite response artifact.
    """

    half = Fraction(1, 2)
    grav = (Fraction(0), Fraction(0))
    su3 = (Fraction(0), Fraction(0))
    su2 = (Fraction(0), Fraction(0))
    for field, row in GENERAL_PACKAGE.items():
        ca, cb = (Fraction(x) for x in row["charge"])
        grav = (grav[0] + row["dim"] * ca, grav[1] + row["dim"] * cb)
        su3 = (su3[0] + row["su3_fund"] * half * ca, su3[1] + row["su3_fund"] * half * cb)
        su2 = (
            su2[0] + row["su2_doublets"] * half * ca,
            su2[1] + row["su2_doublets"] * half * cb,
        )
    require(grav == (24, 16), "BLOCK_DETERMINANT_BALANCE", f"gravity^2 U1 form drifted: {grav}")
    require(su3 == (3, 2), "BLOCK_DETERMINANT_BALANCE", f"SU3^2 U1 form drifted: {su3}")
    require(su2 == (3, 2), "BLOCK_DETERMINANT_BALANCE", f"SU2^2 U1 form drifted: {su2}")

    def u1_cubed(a: Fraction, b: Fraction) -> Fraction:
        total = Fraction(0)
        for field, row in GENERAL_PACKAGE.items():
            total += row["dim"] * _charge_at(field, a, b) ** 3
        return total

    for probe in (Fraction(1), Fraction(-1, 3), Fraction(7, 5)):
        require(
            u1_cubed(probe, Fraction(-3, 2) * probe) == 0,
            "BLOCK_DETERMINANT_BALANCE",
            "U1^3 does not vanish on the balance line",
        )
    require(
        u1_cubed(Fraction(1), Fraction(0)) != 0
        and u1_cubed(Fraction(-1, 3), Fraction(1, 3)) != 0,
        "BLOCK_DETERMINANT_BALANCE",
        "U1^3 vanishes off the balance line; the balance is not selective",
    )

    # Primitive integral block characters on the balance line are exactly
    # +/-(-2, 3). The displayed Y normalization is q/6.
    unit = {
        field: 6 * _charge_at(field, Fraction(-1, 3), Fraction(1, 2))
        for field in GENERAL_PACKAGE
    }
    require(
        unit == {"Q": 1, "u_c": -4, "e_c": 6, "d_c": 2, "L": -3},
        "BLOCK_DETERMINANT_BALANCE",
        f"unit charge spectrum drifted: {unit}",
    )
    balance = 3 * y_color + 2 * y_weak
    require(
        balance == 0,
        "BLOCK_DETERMINANT_BALANCE",
        "the declared pair is off the derived balance line",
    )
    q_color = 6 * y_color
    q_weak = 6 * y_weak
    require(
        q_color.denominator == 1
        and q_weak.denominator == 1
        and 3 * q_color + 2 * q_weak == 0
        and math.gcd(abs(q_color.numerator), abs(q_weak.numerator)) == 1,
        "BLOCK_DETERMINANT_BALANCE",
        "the declared block characters are not a primitive integral solution "
        "of 3 q_C + 2 q_W = 0",
    )
    primitive_pairs = {
        (Fraction(-1, 3), Fraction(1, 2)),
        (Fraction(1, 3), Fraction(-1, 2)),
    }
    require(
        (y_color, y_weak) in primitive_pairs,
        "BLOCK_DETERMINANT_BALANCE",
        "the declared pair is not one of the two charge-conjugate primitive pairs",
    )
    declared_spectrum = {
        field: 6 * _charge_at(field, y_color, y_weak)
        for field in GENERAL_PACKAGE
    }
    return {
        "anomaly_forms_general_charges": {
            "gravity_squared_U1": "24 a + 16 b = 8 (3a + 2b)",
            "SU3_squared_U1": "3a + 2b",
            "SU2_squared_U1": "3a + 2b",
            "U1_cubed": "vanishes identically on the balance line 3a + 2b = 0; nonzero at off-line probes",
        },
        "conclusion_balance": (
            "anomaly freedom of the realized exterior package is equivalent to "
            "the determinant balance 3 y_C + 2 y_W = 0"
        ),
        "integral_spectrum_on_line": {"Q": 1, "u_c": -4, "e_c": 6, "d_c": 2, "L": -3},
        "declared_representative_integral_spectrum": {
            field: int(value) for field, value in declared_spectrum.items()
        },
        "primitivity": "the primitive integral block-character solutions are exactly (q_C,q_W)=+/-(-2,3)",
        "overall_sign_status": "charge-conjugation convention; not selected by the response artifact",
        "derived_pairs_up_to_charge_conjugation": [
            {"color_block": "-1/3", "weak_block": "1/2"},
            {"color_block": "1/3", "weak_block": "-1/2"},
        ],
        "calculation_convention": "one charge table is printed with positive weak-block charge; the source-derived physical object is the unordered charge-conjugate pair",
        "derived_pair": {"color_block": "-1/3", "weak_block": "1/2"},
        "declared_pair": {
            "color_block": frac_text(y_color),
            "weak_block": frac_text(y_weak),
        },
        "declared_matches_derived_pair_up_to_conjugation": True,
    }


def scalar_and_channel_selection_certificate(
    y_color: Fraction,
    y_weak: Fraction,
    channels: Sequence[tuple[str, str, str]],
) -> dict[str, Any]:
    """Derive scalar-charge/channel compatibility, not scalar economy.

    Conditional on a color-singlet weak-doublet scalar, enumerate every field
    pair that satisfies triality and duality. Charge conservation fixes the
    scalar charge for each pair, so no arbitrary bounded charge scan is used.
    The only compatible values are the conjugate pair q_S = +/-3, with three
    channels for either convention. Scalar existence, multiplicity, and the
    absence of additional light scalars remain assumptions outside this scan.
    """

    q = {field: 6 * _charge_at(field, y_color, y_weak) for field in GENERAL_PACKAGE}
    fields = sorted(GENERAL_PACKAGE)
    scan: dict[int, list[tuple[str, str, str]]] = {}
    for i, left in enumerate(fields):
        for right in fields[i:]:
            triality = (
                GENERAL_PACKAGE[left]["triality"]
                + GENERAL_PACKAGE[right]["triality"]
            ) % 3
            duality = (
                GENERAL_PACKAGE[left]["duality"]
                + 1
                + GENERAL_PACKAGE[right]["duality"]
            ) % 2
            if triality != 0 or duality != 0:
                continue
            pair_charge = q[left] + q[right]
            require(
                pair_charge.denominator == 1,
                "SCALAR_SELECTION",
                "a compatible matter pair has non-integral q=6Y charge",
            )
            for scalar, q_s in (("S", -pair_charge), ("Sbar", pair_charge)):
                if q_s == 0:
                    continue
                scan.setdefault(int(q_s), []).append((left, scalar, right))

    require(
        set(scan) == {3, -3},
        "SCALAR_SELECTION",
        f"admissible scalar charges are not the conjugate pair (3, -3): {sorted(scan)}",
    )
    q_weak = 6 * y_weak
    require(
        q_weak.denominator == 1 and int(q_weak) in scan,
        "SCALAR_SELECTION",
        "the weak block does not carry a compatible scalar charge",
    )
    representative_charge = int(q_weak)
    canonical = {
        (left, scalar, right)
        for left, scalar, right in scan[representative_charge]
    }
    declared = {
        (min(left, right), scalar, max(left, right)) for left, scalar, right in channels
    }
    require(
        canonical == declared,
        "SCALAR_SELECTION",
        f"the declared channel list does not equal the derived admissible set: "
        f"derived {sorted(canonical)}, declared {sorted(declared)}",
    )
    return {
        "enumeration_domain": "all unordered realized matter-field pairs; q_S is solved exactly from charge conservation",
        "selection_rules": [
            "charge sum zero in q = 6Y",
            "triality zero modulo three",
            "even total doublet count",
        ],
        "admissible_scalar_charges": [3, -3],
        "conjugation_relation": "q_S = -3 is the conjugate relabeling S <-> Sbar of q_S = +3",
        "representative_scalar_charge": representative_charge,
        "derived_channels_for_declared_representative": sorted(list(c) for c in canonical),
        "derived_channels_at_plus_three": sorted(
            list(c) for c in set(scan[3])
        ),
        "compatibility_conclusion": (
            "if a color-singlet weak-doublet scalar supports these trilinear "
            "channels, its primitive charge is +/-3 and the compatible channel "
            "set has three elements"
        ),
        "declared_channels_equal_compatible_set": True,
        "scalar_content_status": (
            "one weak-block scalar is a manifest assumption; existence, "
            "multiplicity, no-extra-light-scalar economy, and potential are not derived"
        ),
    }


def certificate_payload(
    manifest: Mapping[str, Any],
    base_dir: Path | None = None,
    *,
    allow_control_contracts: bool = False,
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    params = validate_manifest(manifest, allow_control_contracts=allow_control_contracts)
    upstream = load_upstream(manifest, base)

    # --- Semantic response artifact: physical refinement maps -----------------
    current_manifest = upstream["current_manifest"]
    artifact_ref = current_manifest.get("semantic_response_artifact")
    require(
        isinstance(artifact_ref, Mapping) and isinstance(artifact_ref.get("path"), str),
        "ARTIFACT_REFERENCE",
        "the upstream manifest does not bind a semantic response artifact",
    )
    artifact_path = Path(artifact_ref["path"])
    if not artifact_path.is_absolute():
        artifact_path = base / artifact_path
    semantic_artifact = load_json(artifact_path)
    require(
        semantic_artifact.get("artifact_sha256")
        == upstream["semantic_response_artifact_sha256"],
        "ARTIFACT_HASH",
        "the semantic artifact does not match the upstream receipt binding",
    )
    physical_maps_block = semantic_artifact.get("physical_refinement_maps", {})
    artifact_port_maps = (
        physical_maps_block.get("port_persistence_maps", [])
        if isinstance(physical_maps_block, Mapping)
        else []
    )

    balance_certificate = block_determinant_balance_certificate(
        params["y_color"], params["y_weak"]
    )
    scalar_certificate = scalar_and_channel_selection_certificate(
        params["y_color"], params["y_weak"], params["channels"]
    )

    # --- Source-produced statistics inputs -----------------------------------
    spin_artifact = load_spin_statistics_artifact(manifest, base, upstream)
    selection_scan = exterior_selection_scan(params["y_color"], params["y_weak"])

    algebra = CurrentAlgebra(upstream["current_manifest"], base)

    # --- PORT-SPIN-LIFT ------------------------------------------------------
    spin = spin_lift_certificate(algebra)
    require(
        spin["order_profile"] == spin_artifact["order_profile"]
        and spin["centre_order"] == spin_artifact["centre_order"],
        "SPIN_ARTIFACT",
        "the independently recomputed lift group does not match the measured artifact",
    )
    forcing = statistics_forcing_certificate(selection_scan, spin, spin_artifact)
    if params["category_typing"] == "vec":
        raise CertificateError(
            "VEC_TYPING",
            "a Vec-typed category cannot carry the derived nontrivial fermionic grading and non-split spin cocycle",
        )
    if not params["double_cover"]:
        raise CertificateError(
            "SPIN_LIFT_SPLIT",
            "the declared split (no-double-cover) lift contradicts the derived non-split binary icosahedral cover",
        )
    require(
        params["category_typing"] == "spin_odd_weyl_super",
        "CATEGORY_TYPING",
        "only the Spin/odd-Weyl super typing can carry the derived grading, braiding, and double cover together",
    )

    # --- Matter carrier and transport ----------------------------------------
    y_color = params["y_color"]
    y_weak = params["y_weak"]
    transport = matter_transport(algebra, y_color, y_weak)
    images = transport["images"]
    covariance_checks = carrier_conjugation_checks(algebra, spin["lifts"], images)

    hypercharge_v = [[IZERO for _ in range(5)] for _ in range(5)]
    for i in range(3):
        hypercharge_v[i][i] = C5(F5(y_color), ZERO)
    for i in range(2):
        hypercharge_v[3 + i][3 + i] = C5(F5(y_weak), ZERO)

    # --- Auxiliary CAR/Fock space --------------------------------------------
    fock = FockSpace()
    car = car_certificate(fock)
    creations = car["creations"]
    if params["matter_statistics"] != "fermionic_odd":
        raise CertificateError(
            "STATISTICS_TYPING",
            "the derived matter-building operators anticommute, so bosonic-even matter typing fails closed",
        )

    dgammas = [fock.dgamma(images[p]) for p in range(12)]
    parity = fock.parity()
    for p in range(12):
        require(
            c_is_zero(csub(cmul(parity, dgammas[p]), cmul(dgammas[p], parity))),
            "PARITY_TYPING",
            "second-quantized currents must be parity even",
        )
    for i in range(5):
        anti = cadd(cmul(parity, creations[i]), cmul(creations[i], parity))
        require(c_is_zero(anti), "PARITY_TYPING", "creation operators must be parity odd")

    # dGamma is a Lie homomorphism: it annihilates the vacuum and satisfies the
    # CAR-derivation (super-Leibniz) identity [dGamma(X), a^dagger(v)] =
    # a^dagger(Xv); with vacuum cyclicity this pins the bracket property.
    vacuum_index = fock.index[()]
    derivation_checks = 0
    for p in range(12):
        for row_index in range(fock.dim):
            if row_index == vacuum_index:
                continue
            require(
                dgammas[p][row_index][vacuum_index].is_zero(),
                "SECOND_QUANTIZATION",
                "dGamma must annihilate the Fock vacuum",
            )
        require(
            dgammas[p][vacuum_index][vacuum_index].is_zero(),
            "SECOND_QUANTIZATION",
            "dGamma must annihilate the Fock vacuum",
        )
        for mode in range(5):
            lhs = csub(cmul(dgammas[p], creations[mode]), cmul(creations[mode], dgammas[p]))
            rhs = czeros(fock.dim)
            for i in range(5):
                entry = images[p][i][mode]
                if entry.is_zero():
                    continue
                for a in range(fock.dim):
                    for b in range(fock.dim):
                        if not creations[i][a][b].is_zero():
                            rhs[a][b] = rhs[a][b] + entry * creations[i][a][b]
            require(
                c_is_zero(csub(lhs, rhs)),
                "SECOND_QUANTIZATION",
                "[dGamma(X), a^dagger(v)] = a^dagger(Xv) fails",
            )
            derivation_checks += 1

    dgamma_y = fock.dgamma(hypercharge_v)

    # --- Invariant sector and the selection projector -------------------------
    invariant_rows: list[list[F5]] = []
    for p in range(12):
        invariant_rows.extend(complex_constraint_rows(dgammas[p]))
    invariant_dimension = complex_nullspace_dimension(invariant_rows)
    require(
        invariant_dimension == 2,
        "INVARIANT_SECTOR",
        f"expected the derived invariant sector to be the vacuum and top lines, got complex dimension {invariant_dimension}",
    )
    top_index = fock.index[(0, 1, 2, 3, 4)]
    for p in range(12):
        for a in range(fock.dim):
            require(
                dgammas[p][a][top_index].is_zero(),
                "INVARIANT_SECTOR",
                "the top line is not exactly invariant, so trace balance fails at operator level",
            )

    # The source data determine a conjugate pair, not one physically preferred
    # projector. Each parity sector has rank sixteen and contains exactly one
    # of the two derived invariant lines. Removing that line gives the unique
    # rank-fifteen orthogonal complement within that parity sector. Subset
    # complementation (the exterior pairing) exchanges the two candidates.
    even_candidate = [
        n
        for n, subset in enumerate(fock.subsets)
        if len(subset) % 2 == 0 and n != vacuum_index
    ]
    odd_candidate = [
        n
        for n, subset in enumerate(fock.subsets)
        if len(subset) % 2 == 1 and n != top_index
    ]
    require(
        len(even_candidate) == len(odd_candidate) == 15,
        "SELECTION_PROJECTOR",
        "the two parity/invariant-line complements do not both have rank fifteen",
    )
    universe = set(range(5))
    complement_map = {
        index: fock.index[tuple(sorted(universe.difference(subset)))]
        for index, subset in enumerate(fock.subsets)
    }
    require(
        {complement_map[index] for index in even_candidate} == set(odd_candidate),
        "SELECTION_PROJECTOR",
        "exterior conjugation does not exchange the two rank-fifteen candidates",
    )
    derived_projector_rows: dict[str, dict[str, Any]] = {}
    for label, indices in (
        ("even_minus_vacuum", even_candidate),
        ("odd_minus_top", odd_candidate),
    ):
        candidate = czeros(fock.dim)
        for index in indices:
            candidate[index][index] = IONE
        require(
            c_is_zero(csub(cmul(candidate, candidate), candidate))
            and c_is_zero(csub(cdagger(candidate), candidate)),
            "SELECTION_PROJECTOR",
            f"the derived {label} candidate is not an orthogonal projector",
        )
        for p in range(12):
            require(
                c_is_zero(
                    csub(cmul(candidate, dgammas[p]), cmul(dgammas[p], candidate))
                ),
                "SELECTION_PROJECTOR",
                f"the derived {label} candidate is not current-equivariant",
            )
        derived_projector_rows[label] = {
            "rank": len(indices),
            "orthogonal_projector": True,
            "current_equivariant": True,
            "contains_invariant_line": False,
        }

    rule = params["selection_rule"]
    selected: list[int] = []
    for n, subset in enumerate(fock.subsets):
        size = len(subset)
        if rule == "parity_even_minus_derived_invariants":
            keep = size % 2 == 0 and n not in (vacuum_index, top_index)
        elif rule == "parity_odd_minus_derived_invariants":
            keep = size % 2 == 1 and n not in (vacuum_index, top_index)
        elif rule == "lambda2_only":
            keep = size == 2
        else:  # even_including_vacuum
            keep = size % 2 == 0
        if keep:
            selected.append(n)
    projector = czeros(fock.dim)
    for n in selected:
        projector[n][n] = IONE
    require(
        c_is_zero(csub(cmul(projector, projector), projector)),
        "SELECTION_PROJECTOR",
        "the selection projector is not idempotent",
    )
    require(
        c_is_zero(csub(cdagger(projector), projector)),
        "SELECTION_PROJECTOR",
        "the selection projector is not self-adjoint",
    )
    equivariance_checks = 0
    for p in range(12):
        require(
            c_is_zero(csub(cmul(projector, dgammas[p]), cmul(dgammas[p], projector))),
            "SELECTION_PROJECTOR",
            "the selection projector is not equivariant for the current action",
        )
        equivariance_checks += 1
    require(
        c_is_zero(csub(cmul(projector, parity), cmul(parity, projector))),
        "SELECTION_PROJECTOR",
        "the selection projector does not commute with fermionic parity",
    )

    matter_indices = selected
    matter_dimension = len(matter_indices)

    # No derived invariant line may survive inside the matter module.
    for n in matter_indices:
        require(
            n not in (vacuum_index, top_index),
            "TRIVIAL_LINE_IN_MATTER",
            "the selected matter module contains a derived invariant line",
        )

    # --- Realized module: faithfulness, blocks, multiplicity ------------------
    def restrict(matrix: Sequence[Sequence[C5]], indices: Sequence[int] | None = None) -> list[list[C5]]:
        chosen = matter_indices if indices is None else indices
        return [[matrix[a][b] for b in chosen] for a in chosen]

    matter_ops = [restrict(dgammas[p]) for p in range(12)]
    matter_rank = span_dimension(matter_ops)
    require(
        matter_rank == 12,
        "CURRENT_ACTION_NOT_FAITHFUL",
        f"the current algebra must act faithfully on the matter tensors, got rank {matter_rank}",
    )

    weight_of: dict[int, Fraction] = {}
    for n, subset in enumerate(fock.subsets):
        weight_of[n] = sum(
            (y_color if i < 3 else y_weak for i in subset),
            Fraction(0),
        )
    weight_multiset: dict[Fraction, int] = {}
    for n in matter_indices:
        weight_multiset[weight_of[n]] = weight_multiset.get(weight_of[n], 0) + 1

    def triality(subset: tuple[int, ...]) -> int:
        return sum(1 for i in subset if i < 3) % 3

    def duality(subset: tuple[int, ...]) -> int:
        return sum(1 for i in subset if i >= 3) % 2

    # su(3) and su(2) bases from the source transport.
    su3_candidates: list[list[list[C5]]] = []
    for p in range(12):
        even = algebra.even_block(p)
        trace = ctrace(even)
        third = C5(trace.re / F5(3), trace.im / F5(3))
        traceless = [[even[i][j] - (third if i == j else IZERO) for j in range(3)] for i in range(3)]
        su3_candidates.append(traceless)
    su3_basis = basis_by_pivots(su3_candidates)
    require(len(su3_basis) == 8, "MATTER_TRANSPORT", f"expected an eight-dimensional su(3) span, got {len(su3_basis)}")
    half = F5(Fraction(1, 2))
    su2_candidates: list[list[list[C5]]] = []
    for p in range(12):
        w = algebra.kernel_axis(p)
        su2_candidates.append(cscale(sigma_dot(w), C5(ZERO, -half)))
    su2_basis = basis_by_pivots(su2_candidates)
    require(len(su2_basis) == 3, "MATTER_TRANSPORT", f"expected a three-dimensional su(2) span, got {len(su2_basis)}")

    def embed_color(block: Sequence[Sequence[C5]]) -> list[list[C5]]:
        out = [[IZERO for _ in range(5)] for _ in range(5)]
        for i in range(3):
            for j in range(3):
                out[i][j] = block[i][j]
        return out

    def embed_weak(block: Sequence[Sequence[C5]]) -> list[list[C5]]:
        out = [[IZERO for _ in range(5)] for _ in range(5)]
        for i in range(2):
            for j in range(2):
                out[3 + i][3 + j] = block[i][j]
        return out

    su3_matter = [restrict(fock.dgamma(embed_color(b))) for b in su3_basis]
    su2_matter = [restrict(fock.dgamma(embed_weak(b))) for b in su2_basis]
    y_matter = restrict(dgamma_y)

    # Weak doublet count: a subset basis state lies in a doublet exactly when
    # every su(2) generator kills it or not; the diagonal of the summed squares
    # is -sum ||S v||^2, which vanishes exactly on su(2)-invariant states.
    su2_casimir = czeros(matter_dimension)
    for op in su2_matter:
        su2_casimir = cadd(su2_casimir, cmul(op, op))
    doublet_states = 0
    for n in range(matter_dimension):
        diagonal = su2_casimir[n][n]
        require(diagonal.im.is_zero(), "WITTEN_PARITY", "su(2) Casimir diagonal is not real")
        if not diagonal.re.is_zero():
            doublet_states += 1
    require(doublet_states % 2 == 0, "WITTEN_PARITY_STATES", "doublet state count must be even")
    weak_doublets = doublet_states // 2
    require(
        weak_doublets % 2 == 0,
        "WITTEN_PARITY",
        f"the Witten parity requires an even number of weak doublets, got {weak_doublets}",
    )

    # --- Realized perturbative anomaly traces ----------------------------------
    anomaly_traces: dict[str, str] = {}

    trace_y = ctrace(y_matter)
    require(trace_y.is_zero(), "PERTURBATIVE_ANOMALY", "the realized gravity^2 U1 trace does not vanish")
    anomaly_traces["gravity_squared_U1"] = "0"
    y_cubed = cmul(cmul(y_matter, y_matter), y_matter)
    require(ctrace(y_cubed).is_zero(), "PERTURBATIVE_ANOMALY", "the realized U1^3 trace does not vanish")
    anomaly_traces["U1_cubed"] = "0"
    for a, ta in enumerate(su3_matter):
        for b, tb in enumerate(su3_matter):
            value = ctrace(cmul(y_matter, cmul(ta, tb)))
            require(
                value.is_zero(),
                "PERTURBATIVE_ANOMALY",
                f"the realized SU3^2 U1 anomaly trace ({a},{b}) does not vanish",
            )
    anomaly_traces["SU3_squared_U1"] = "0"
    for a, sa in enumerate(su2_matter):
        for b, sb in enumerate(su2_matter):
            value = ctrace(cmul(y_matter, cmul(sa, sb)))
            require(
                value.is_zero(),
                "PERTURBATIVE_ANOMALY",
                f"the realized SU2^2 U1 anomaly trace ({a},{b}) does not vanish",
            )
    anomaly_traces["SU2_squared_U1"] = "0"
    dsymbol_checks = 0
    for a in range(8):
        for b in range(a, 8):
            for c in range(b, 8):
                anti = cadd(
                    cmul(su3_matter[b], su3_matter[c]),
                    cmul(su3_matter[c], su3_matter[b]),
                )
                value = ctrace(cmul(su3_matter[a], anti))
                require(
                    value.is_zero(),
                    "PERTURBATIVE_ANOMALY",
                    f"the realized SU3^3 d-symbol trace ({a},{b},{c}) does not vanish",
                )
                dsymbol_checks += 1
    anomaly_traces["SU3_cubed"] = "0"

    # --- Chirality and conjugation ---------------------------------------------
    dual_weights = {-value for value in weight_multiset}
    spectra_intersection = sorted(
        frac_text(v) for v in (set(weight_multiset) & dual_weights)
    )
    require(
        not spectra_intersection,
        "CHIRALITY",
        "the matter module and its dual share a charge eigenvalue, so chirality is not established",
    )

    complement_of = {
        n: fock.index[tuple(i for i in range(5) if i not in fock.subsets[n])]
        for n in range(fock.dim)
    }
    complement_indices = [complement_of[n] for n in matter_indices]

    def wedge_sign(subset: tuple[int, ...]) -> int:
        permutation = list(subset) + [i for i in range(5) if i not in subset]
        sign = 1
        for i in range(len(permutation)):
            for j in range(i + 1, len(permutation)):
                if permutation[i] > permutation[j]:
                    sign = -sign
        return sign

    pairing = czeros(matter_dimension)
    for row, n in enumerate(matter_indices):
        col = complement_indices.index(complement_of[n])
        pairing[row][col] = C5(F5(wedge_sign(fock.subsets[n])), ZERO)
    nondegenerate = all(
        sum(1 for col in range(matter_dimension) if not pairing[row][col].is_zero()) == 1
        for row in range(matter_dimension)
    )
    require(nondegenerate, "CONJUGATION", "the wedge pairing into the top line is degenerate")

    conjugation_invariance_checks = 0
    for p in range(12):
        left = restrict(dgammas[p])
        right = restrict(dgammas[p], complement_indices)
        transposed = [[left[b][a] for b in range(matter_dimension)] for a in range(matter_dimension)]
        residual = cadd(cmul(transposed, pairing), cmul(pairing, right))
        require(
            c_is_zero(residual),
            "CONJUGATION",
            "the wedge pairing into the top line is not invariant under the current action",
        )
        conjugation_invariance_checks += 1

    # --- Field blocks, multiplicity freeness, Yukawa channels -------------------
    blocks: dict[str, list[int]] = {}
    reference_signatures = {
        "Q": (y_color + y_weak, 1, 1),
        "u_c": (2 * y_color, 2, 0),
        "e_c": (2 * y_weak, 0, 0),
        "d_c": (2 * y_color + 2 * y_weak, 2, 0),
        "L": (3 * y_color + y_weak, 0, 1),
    }
    for n in matter_indices:
        subset = fock.subsets[n]
        signature = (weight_of[n], triality(subset), duality(subset))
        for label, expected in reference_signatures.items():
            if signature == expected:
                blocks.setdefault(label, []).append(n)
                break

    if set(blocks) != set(FIELD_LABELS):
        raise CertificateError(
            "YUKAWA_CHANNEL_EMPTY",
            "the selected module does not realize the declared exterior field blocks, so every declared channel is empty",
        )

    expected_dimensions = {"Q": 6, "u_c": 3, "e_c": 1, "d_c": 3, "L": 2}
    field_tables: dict[str, dict[str, Any]] = {}
    commutant_total = 0
    for label, indices in blocks.items():
        require(
            len(indices) == expected_dimensions[label],
            "PACKAGE_REALIZATION",
            f"field block {label} has dimension {len(indices)}, expected {expected_dimensions[label]}",
        )
        block_ops_here = [restrict(dgammas[p], indices) for p in range(12)]
        size = len(indices)
        rows: list[list[F5]] = []
        for op in block_ops_here:
            # S X = X S as real-linear constraints on S (size x size complex).
            for i in range(size):
                for j in range(size):
                    real_row = [ZERO] * (2 * size * size)
                    imag_row = [ZERO] * (2 * size * size)
                    for k in range(size):
                        # + S[i][k] X[k][j]
                        entry = op[k][j]
                        base = 2 * (i * size + k)
                        real_row[base] = real_row[base] + entry.re
                        real_row[base + 1] = real_row[base + 1] - entry.im
                        imag_row[base] = imag_row[base] + entry.im
                        imag_row[base + 1] = imag_row[base + 1] + entry.re
                        # - X[i][k] S[k][j]
                        entry = op[i][k]
                        base = 2 * (k * size + j)
                        real_row[base] = real_row[base] - entry.re
                        real_row[base + 1] = real_row[base + 1] + entry.im
                        imag_row[base] = imag_row[base] - entry.im
                        imag_row[base + 1] = imag_row[base + 1] - entry.re
                    rows.append(real_row)
                    rows.append(imag_row)
        commutant_dimension = complex_nullspace_dimension(rows)
        require(
            commutant_dimension == 1,
            "PACKAGE_REALIZATION",
            f"field block {label} is not irreducible: commutant dimension {commutant_dimension}",
        )
        commutant_total += commutant_dimension
        field_tables[label] = {
            "dimension": size,
            "charge": frac_text(reference_signatures[label][0]),
            "commutant_dimension": 1,
        }
    require(
        commutant_total == 5,
        "PACKAGE_REALIZATION",
        "the realized module is not multiplicity-free with five irreducible blocks",
    )

    def dual_ops(ops: Sequence[Sequence[Sequence[C5]]]) -> list[list[list[C5]]]:
        out = []
        for op in ops:
            n = len(op)
            out.append([[C5(-op[b][a].re, -op[b][a].im) for b in range(n)] for a in range(n)])
        return out

    def scalar_ops(kind: str) -> tuple[int, list[list[list[C5]]]]:
        base_ops: list[list[list[C5]]] = []
        for p in range(12):
            image = images[p]
            block = [[image[3 + i][3 + j] for j in range(2)] for i in range(2)]
            base_ops.append(block)
        if kind == "S":
            return 2, base_ops
        return 2, dual_ops(base_ops)

    def kron_sum_invariants(op_families: Sequence[tuple[int, Sequence[Sequence[Sequence[C5]]]]]) -> int:
        dims = [dim for dim, _ in op_families]
        total = 1
        for d in dims:
            total *= d
        rows: list[list[F5]] = []
        for p in range(12):
            matrix = [[IZERO for _ in range(total)] for _ in range(total)]
            for slot, (dim, ops) in enumerate(op_families):
                op = ops[p]
                stride_after = 1
                for d in dims[slot + 1 :]:
                    stride_after *= d
                stride_block = dim * stride_after
                for outer in range(total // stride_block):
                    for i in range(dim):
                        for j in range(dim):
                            entry = op[i][j]
                            if entry.is_zero():
                                continue
                            for inner in range(stride_after):
                                row_index = outer * stride_block + i * stride_after + inner
                                col_index = outer * stride_block + j * stride_after + inner
                                matrix[row_index][col_index] = matrix[row_index][col_index] + entry
            rows.extend(complex_constraint_rows(matrix))
        return complex_nullspace_dimension(rows)

    yukawa_rows: list[dict[str, Any]] = []
    for left, scalar, right in params["channels"]:
        families = [
            (len(blocks[left]), [restrict(dgammas[p], blocks[left]) for p in range(12)]),
            scalar_ops(scalar),
            (len(blocks[right]), [restrict(dgammas[p], blocks[right]) for p in range(12)]),
        ]
        dimension = kron_sum_invariants(families)
        require(
            dimension == 1,
            "YUKAWA_CHANNEL_EMPTY",
            f"the declared Yukawa channel {left} {scalar} {right} carries invariant dimension {dimension}, not one",
        )
        yukawa_rows.append({"channel": [left, scalar, right], "invariant_dimension": 1})
    forbidden_dimension = kron_sum_invariants(
        [
            (len(blocks["Q"]), [restrict(dgammas[p], blocks["Q"]) for p in range(12)]),
            scalar_ops("S"),
            (len(blocks["d_c"]), [restrict(dgammas[p], blocks["d_c"]) for p in range(12)]),
        ]
    )
    require(
        forbidden_dimension == 0,
        "YUKAWA_CONTROL",
        "the forbidden channel Q S d_c unexpectedly carries an invariant line",
    )

    # --- Common action kernel: emitted, never assumed ----------------------------
    weight_rows: list[tuple[str, Fraction, int, int]] = []
    for n in matter_indices:
        subset = fock.subsets[n]
        weight_rows.append((f"matter_state_{n}", weight_of[n], triality(subset), duality(subset)))
    for i in range(5):
        weight_rows.append(
            (f"carrier_mode_{i}", y_color if i < 3 else y_weak, 1 if i < 3 else 0, 0 if i < 3 else 1)
        )
    weight_rows.append(("scalar_S", y_weak, 0, 1))
    for scalar in params["extra_scalars"]:
        weight_rows.append(
            (f"scalar_{scalar['label']}", scalar["charge"], scalar["triality"], scalar["duality"])
        )

    normalization = 1
    for _, w, _, _ in weight_rows:
        normalization = normalization * w.denominator // math.gcd(normalization, w.denominator)
    integer_charges: list[tuple[str, int, int, int]] = []
    for name, w, t, d in weight_rows:
        q = w * normalization
        require(q.denominator == 1, "KERNEL_INTEGRALITY", "derived integral charges failed")
        integer_charges.append((name, int(q), t, d))
    nonzero_qs = [abs(q) for _, q, _, _ in integer_charges if q != 0]
    require(
        bool(nonzero_qs),
        "KERNEL_NOT_DISCRETE",
        "every realized charge vanished, so the whole central R sits inside the kernel; "
        "a non-discrete kernel cannot be emitted as data",
    )
    charge_gcd = 0
    for q in nonzero_qs:
        charge_gcd = math.gcd(charge_gcd, q)

    # The kernel lives on the simply connected cover R x SU(3) x SU(2): the
    # central factor is the non-compact R (in phase turns), NOT U(1); no
    # compactification quotient is chosen here.  An element (r, a, b) of
    # R x Z3 x Z2 acts trivially on a realized weight (q, t, d) iff
    # r q + a t/3 + b d/2 in Z.  Multiplying by 6 shows 6 r q in Z for every
    # integral charge q, and a Bezout combination gives
    # r in (1/(6 gcd)) Z.  Membership then depends on the numerator
    # k = 6 gcd r only through k mod (6 gcd) (adding 6 gcd to k shifts each
    # phase by the integer q), so enumerating one fundamental window of
    # residues determines the full kernel on the cover exactly.
    denominator_bound = 6 * charge_gcd
    kernel_residues: list[tuple[int, int, int]] = []
    for k in range(denominator_bound):
        r = Fraction(k, denominator_bound)
        for a in range(3):
            for b in range(2):
                trivial = True
                for _, q, t, d in integer_charges:
                    phase = r * q + Fraction(a * t, 3) + Fraction(b * d, 2)
                    if phase.denominator != 1:
                        trivial = False
                        break
                if trivial:
                    kernel_residues.append((k, a, b))
    require(
        (0, 0, 0) in kernel_residues,
        "KERNEL_CLOSURE",
        "the identity residue is missing from the enumerated kernel",
    )
    # The unit deck translation (one full central turn, trivial centers) is
    # always in the kernel: r = 1 times an integral charge is an integer.
    # Its residue is (0, 0, 0), so the kernel on the cover is INFINITE; the
    # finite object below is the residual modulo these pure translations.
    residual_order = len(kernel_residues)
    require(
        residual_order > 1,
        "KERNEL_TRIVIAL",
        "the kernel reduces to the pure full-turn deck translations; "
        "a trivial residual kernel cannot satisfy the packet",
    )

    # Residue-level closure: composition in R x Z3 x Z2 adds numerators in Z
    # (full-turn carries do not affect membership), so the kernel is a
    # subgroup iff the residue set is closed under residue addition.
    residue_set = set(kernel_residues)
    for x in kernel_residues:
        for y in kernel_residues:
            composed = ((x[0] + y[0]) % denominator_bound, (x[1] + y[1]) % 3, (x[2] + y[2]) % 2)
            require(
                composed in residue_set,
                "KERNEL_CLOSURE",
                "the emitted kernel is not closed under composition",
            )
    # Torsion-freeness on the cover: a torsion element must have r = 0
    # (R is torsion-free), so torsion shows up as a residue (0, a, b) with
    # (a, b) != (0, 0).  Its absence forces the su3/su2 components to be
    # determined by the r-numerator, so the kernel projects injectively to
    # (1/(6 gcd)) Z and is free of rank one: infinite cyclic.
    for k, a, b in kernel_residues:
        require(
            k != 0 or (a == 0 and b == 0),
            "KERNEL_CYCLIC",
            "the emitted kernel has torsion over the pure deck translations and is not infinite cyclic",
        )
    numerators = sorted({k for k, _, _ in kernel_residues})
    step = min(k for k in numerators if k > 0) if len(numerators) > 1 else denominator_bound
    require(
        denominator_bound % step == 0 and denominator_bound // step == residual_order,
        "KERNEL_CYCLIC",
        "the kernel residues are not generated by a single element",
    )
    generator_row = next((k, a, b) for k, a, b in sorted(kernel_residues) if k == step)
    generator = (Fraction(step, denominator_bound), generator_row[1], generator_row[2])
    # Verify cyclic generation and the deck relation exactly: the n-th power
    # of the generator has numerator n*step; at n = residual_order it equals
    # one full turn with trivial center components - the unit deck
    # translation, which is NOT the identity on the cover.
    for n in range(1, residual_order + 1):
        power = (
            (n * step) % denominator_bound,
            (n * generator[1]) % 3,
            (n * generator[2]) % 2,
        )
        require(
            power in residue_set,
            "KERNEL_CYCLIC",
            "a generator power leaves the enumerated kernel residues",
        )
    require(
        generator[0] * residual_order == 1
        and (generator[1] * residual_order) % 3 == 0
        and (generator[2] * residual_order) % 2 == 0,
        "KERNEL_CYCLIC",
        "the generator's residual-order power is not the unit deck translation",
    )

    kernel_payload = {
        "cover": "R x SU(3) x SU(2), the simply connected cover of the derived current group data",
        "central_factor": "the non-compact R (phase turns); no compactification of R to U(1) is chosen here",
        "integrality_normalization": normalization,
        "kernel_group_on_cover": "infinite cyclic (isomorphic to Z); the kernel is NOT finite on the cover",
        "kernel_generator": {
            "u1_phase_turns": frac_text(generator[0]),
            "su3_center_power": generator[1],
            "su2_center_power": generator[2],
        },
        "deck_relation": (
            f"generator^{residual_order} = (one full central turn, trivial centers), the unit deck "
            "translation, which is not the identity on the cover"
        ),
        "pure_deck_translation_subgroup": "generated by (1 turn, 0, 0) = generator^" + str(residual_order),
        "residual_order_modulo_deck_translations": residual_order,
        "kernel_residues_modulo_deck_translations": [
            {
                "u1_phase_turns": frac_text(Fraction(k, denominator_bound)),
                "su3_center_power": a,
                "su2_center_power": b,
            }
            for k, a, b in sorted(kernel_residues)
        ],
        "verified_trivial_on": "every realized matter state, every carrier mode, and the declared scalar",
        "tensor_additivity": (
            "charge, triality, and duality are additive over wedge and tensor factors by construction, "
            "so triviality on the verified weight list extends to every realized matter tensor"
        ),
        "global_quotient_assumed": False,
        "downstream_consumer": (
            "AXIS-CENTER-DESCENT (global-form descent); neither the compactification of the central R "
            "nor any global quotient is formed here - the emitted generator, deck relation, and residual "
            "determine the kernel image in every candidate quotient"
        ),
    }

    # --- Refinement descent --------------------------------------------------------
    def check_refinement_permutation(permutation: tuple[int, ...], error_code: str) -> None:
        rotation = algebra.frame.rotation_of(permutation)
        kernel_rotation = [[entry.conj() for entry in row] for row in rotation]
        lift = spin["lifts"].get(tuple(permutation)) or spin_lift_of_rotation(kernel_rotation)
        pi_v = carrier_implementer(rotation, lift)
        gamma = fock.exterior_lift(pi_v)
        require(
            c_is_zero(csub(cmul(gamma, projector), cmul(projector, gamma))),
            error_code,
            "a refinement map does not commute with the matter selection projector",
        )
        pi_dagger = cdagger(pi_v)
        gamma_dagger = cdagger(gamma)
        for p in range(12):
            target_index = permutation[p]
            conjugated_v = cmul(cmul(pi_v, [list(row) for row in images[p]]), pi_dagger)
            require(
                c_is_zero(csub(conjugated_v, [list(row) for row in images[target_index]])),
                error_code,
                "a refinement map is not intertwined on the matter carrier",
            )
            conjugated_f = cmul(cmul(gamma, dgammas[p]), gamma_dagger)
            require(
                c_is_zero(csub(conjugated_f, dgammas[target_index])),
                error_code,
                "a refinement map is not intertwined on the Fock realization",
            )

    tower = algebra.carrier_manifest["refinement_tower"]
    refinement_rows = []
    for item in tower["maps"]:
        permutation = e565.parse_port_permutation(item["port_map"], algebra.carrier)
        check_refinement_permutation(tuple(permutation), "REFINEMENT_DESCENT")
        refinement_rows.append({"source": item["source"], "target": item["target"], "intertwined": True})

    # Physical refinement maps from the bound semantic artifact: the
    # defect-port persistence maps of the geodesic tower, each intertwined on
    # the matter carrier and the Fock realization.
    physical_refinement_rows = []
    require(
        len(artifact_port_maps) > 0,
        "PHYSICAL_REFINEMENT",
        "the bound artifact carries no physical refinement maps",
    )
    for row in artifact_port_maps:
        port_map = row.get("port_map")
        require(
            isinstance(port_map, list) and sorted(port_map) == list(range(12)),
            "PHYSICAL_REFINEMENT",
            "a physical refinement map is not a port permutation",
        )
        check_refinement_permutation(tuple(int(v) for v in port_map), "PHYSICAL_REFINEMENT")
        physical_refinement_rows.append(
            {
                "source_level": row.get("source_level"),
                "target_level": row.get("target_level"),
                "origin": row.get("origin"),
                "intertwined": True,
            }
        )

    # --- Gate ------------------------------------------------------------------------
    require(matter_dimension == 15, "PACKAGE_REALIZATION", f"expected fifteen matter states, got {matter_dimension}")
    gate = {
        "port_spin_lift_nonsplit_double_cover": True,
        "faithful_current_action_on_matter_tensors": True,
        "fifteen_state_module_from_equivariant_projector": True,
        "exterior_package_realized_on_cover": True,
        "listed_anomalies_and_witten_parity_checked": True,
        "chiral_no_common_summand_with_dual": True,
        "conjugation_and_super_tensor_derived": True,
        "nonzero_invariant_sector": True,
        "common_action_kernel_emitted_not_assumed": True,
        "candidate_matter_class_nonempty_witnessed": True,
        "declared_tower_descent": True,
        "family_and_potential_and_mass_firewalled": True,
    }
    require(all(gate.values()), "GATE", "conditional algebraic matter-lift gate did not pass")

    weight_spectrum = {
        frac_text(value): count for value, count in sorted(weight_multiset.items())
    }

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 314,
        "manifest_sha256": sha256_json(manifest),
        "upstream": {
            "current_manifest_sha256": upstream["current_manifest_sha256"],
            "current_receipt_sha256": upstream["current_receipt_sha256"],
            "carrier_manifest_sha256": upstream["carrier_manifest_sha256"],
            "semantic_response_artifact_sha256": upstream["semantic_response_artifact_sha256"],
            "dependencies": [
                "#565 (carrier packet)",
                "#566 (finite port-current algebra)",
                "#599 (source-bound target-blind carrier impulse/readback response)",
            ],
            "inherited_scope": (
                "the pinned #566 response representation passes its source gate by "
                "an exact target-blind impulse/readback protocol on the carrier "
                "recurrence; this is a finite carrier-dynamics result, not a "
                "laboratory gauge-current attachment"
            ),
        },
        "source_firewall": {
            "forbidden_dependency_hits": [],
            "uses_only": [
                "hash-pinned #566 conditional finite current packet (manifest and receipt)",
                "the bound semantic response artifact (finite response definition and physical refinement maps)",
                "the trace-balanced charge pair up to overall charge conjugation, derived by BLOCK-DETERMINANT-BALANCE",
                "scalar-charge and Yukawa-channel compatibility conditional on the declared weak-doublet scalar content",
            "typed statistics and Spin/odd-Weyl category contracts; the controls test consistency but do not source those physical inputs",
                "kernel emission contract and candidate matter class declaration",
            ],
        },
        "block_determinant_balance": balance_certificate,
        "scalar_and_channel_selection": scalar_certificate,
        "category_typing_source_binding": {
            "double_cover_forced": "the derived lift group has a unique nontrivial involution and does not factor through the sixty rotations, so a split lift is impossible",
            "vec_control": "the Vec typing fails closed on the realized module (negative control vec_typing)",
            "svec_control": "the split-spin sVec control fails closed on the realized module",
            "weyl_pair": "the even-minus-vacuum and odd-minus-top projectors form a source-derived unordered charge-conjugate pair; the source does not select either representative",
            "conclusion": "the exact lift and projector-pair calculations are derived, but the manifest still types the physical matter object as fermionic Spin/odd-Weyl; failing alternative controls do not derive that premise",
        },
        "port_spin_lift": {
            "witness_count": spin["witness_count"],
            "lift_group_order": spin["lift_group_order"],
            "unique_involution": spin["unique_involution"],
            "order_profile": spin["order_profile"],
            "involution_lift_order": spin["involution_lift_order"],
            "irrational_order_five_spinor_traces": spin["irrational_order_five_spinor_traces"],
            "conclusion": spin["conclusion"],
        },
        "matter_carrier": {
            "carrier": "V = C (+) W = (even response sector C^3) (+) (spin-lifted weak doublet C^2)",
            "block_trace_charges": {"color_block": frac_text(y_color), "weak_block": frac_text(y_weak)},
            "trace_balance": "3 y_C + 2 y_W = 0, checked exactly; the top line is exactly invariant because of it",
            "raw_source_central_charge": "i on the even response sector, 0 on the kernel sector (from the #566 packet)",
            "central_charge_provenance": (
                "the trace-balanced charge pair is derived up to charge conjugation by "
                "BLOCK-DETERMINANT-BALANCE: anomaly freedom of the realized package forces the balance "
                "line and primitive integral block characters give +/-(-2,3); the printed "
                "(y_C,y_W)=(-1/3,1/2) representative is a convention, not an artifact measurement"
            ),
            "transport_homomorphism_bracket_checks": transport["homomorphism_bracket_checks"],
            "transport_covariance_checks": covariance_checks,
            "faithful_rank_on_carrier": transport["faithful_rank_on_carrier"],
        },
        "auxiliary_car_fock": {
            "dimension": 32,
            "car_relation_checks": car["car_checks"],
            "vacuum_cyclic_rank": car["cyclic_rank"],
            "fermionic_parity": "(-1)^N, derived from the CAR grading; currents are parity even, creations parity odd",
            "second_quantization_derivation_checks": derivation_checks,
            "koszul_braiding": "creation operators anticommute exactly (CAR), giving the super braiding and the wedge Leibniz rule",
        },
        "selection": {
            "rule": rule,
            "derived_invariant_sector": "complex dimension 2: the Fock vacuum and the top line",
            "derived_conjugate_projector_pair": derived_projector_rows,
            "uniqueness_scope": (
                "within projectors obtained by choosing one fermionic-parity "
                "eigenspace and removing its complete invariant line, these are "
                "the only two rank-fifteen orthogonal equivariant complements"
            ),
            "conjugation_exchange_verified": True,
            "physical_object": "unordered conjugacy class {even parity minus vacuum, odd parity minus top}",
            "calculation_representative": "even parity minus vacuum, used only to print one charge table; every physical conclusion is stated modulo charge conjugation",
            "projector_rank": matter_dimension,
            "equivariance_checks": equivariance_checks,
            "commutes_with_parity": True,
            "realization": "operator projector on the auxiliary CAR/Fock space, not representation arithmetic",
        },
        "realized_package": {
            "module": "M1 = Lambda^2 V (+) Lambda^4 V",
            "dimension": matter_dimension,
            "charge_spectrum": weight_spectrum,
            "integrality_normalization": normalization,
            "fields": field_tables,
            "multiplicity_free": True,
            "irreducible_block_commutants": commutant_total,
            "contains_no_invariant_line": True,
            "faithful_rank_on_matter": matter_rank,
        },
        "chirality": {
            "matter_spectrum_disjoint_from_dual": True,
            "argument": (
                "any intertwiner S with the dual module satisfies S dGamma(Y) = -dGamma(Y)^T S, so S maps each "
                "charge eigenspace into the negated eigenvalue; the exact spectra are disjoint, hence S = 0"
            ),
            "hom_dimension_with_dual": 0,
        },
        "conjugation": {
            "pairing": "wedge into the top line Lambda^5 V, which is exactly invariant by trace balance",
            "invariance_checks": conjugation_invariance_checks,
            "nondegenerate": True,
            "statement": "the parity-complement selection realizes the conjugate module inside the same Fock space",
        },
        "anomalies": {
            "scope": "listed four-dimensional perturbative anomalies, realized as exact operator traces on the matter module",
            "traces": anomaly_traces,
            "su3_d_symbol_checks": dsymbol_checks,
            "witten_parity": {
                "weak_doublets": weak_doublets,
                "even": True,
                "scope": "finite Witten-parity surrogate: the mod-2 count of realized weak doublets",
            },
        },
        "yukawa_sector": {
            "channels": yukawa_rows,
            "forbidden_channel_control": {"channel": ["Q", "S", "d_c"], "invariant_dimension": 0},
            "invariant_sector_dimension": len(yukawa_rows),
            "statement": "the invariant Gauss sector is nonzero: one exact invariant line per declared channel",
        },
        "kernel_emission": kernel_payload,
        "refinement": {
            "natural": True,
            "maps": refinement_rows,
            "physical_maps": physical_refinement_rows,
            "scope": (
                "naturality along the declared algebraic tower maps and the "
                "artifact-bound physical defect-port persistence maps, each "
                "intertwined on the matter carrier and the Fock realization"
            ),
        },
        "candidate_matter_class": {
            "declared": "one_generation_one_scalar_chiral_anomaly_free",
            "nonempty": True,
            "witness": (
                "the realized packet: fifteen multiplicity-free states (one generation), one declared scalar, "
                "exact chirality, vanishing realized anomaly traces"
            ),
            "uniqueness_promoted": False,
            "note": "this packet discharges the nonemptiness precondition; class uniqueness stays in its own lane",
        },
        "category": {
            "typing": "spin_odd_weyl_super",
            "objects": "subquotients of tensor powers of V and its dual, realized on the auxiliary CAR/Fock space",
            "parity": "derived from the CAR grading",
            "spin": "the derived algebraic PORT-SPIN-LIFT target with non-split double cover",
            "tensor": "wedge product with Koszul braiding; the CAR-derivation identity is the exact Leibniz rule",
            "conjugation": "wedge pairing into the invariant top line",
            "nonempty": "the odd matter object is nonzero (fifteen states) with a faithful current action",
        },
        "matter_selection_scan": selection_scan,
        "statistics_forcing": forcing,
        "source_statistics_binding": {
            "spin_statistics_artifact_sha256": spin_artifact["artifact_sha256"],
            "measured_order_profile": spin_artifact["order_profile"],
            "measured_centre_order": spin_artifact["centre_order"],
            "measured_no_section_over_klein_four": spin_artifact["no_section"],
            "measured_spin_structure_count": spin_artifact["spin_structure_count"],
            "recomputed_lift_matches_measured": True,
            "conclusion": (
                "the fermionic-parity grading is derived by the exhaustive "
                "selection scan and its unique source implementation is the "
                "measured central -1 of the non-split transport double cover; "
                "the declared statistics and category contracts are validated "
                "against, not substituted for, this derivation"
            ),
        },
        "conditional_algebraic_gate": {**gate, "passed": True},
        "physical_source_gate": {
            "charge_pair_derived_up_to_charge_conjugation": bool(
                balance_certificate["declared_matches_derived_pair_up_to_conjugation"]
            ),
            "conjugate_projector_pair_source_derived": True,
            "unordered_pair_selected_by_exhaustive_anomaly_scan": bool(
                selection_scan["survivor_count"] == 2
                and selection_scan["survivors_equal_parity_sectors"]
            ),
            "single_projector_representative_source_selected": False,
            "fermionic_statistics_source_derived": bool(
                selection_scan["survivors_equal_parity_sectors"]
                and forcing["forced_typing"] == "spin_odd_weyl_super"
            ),
            "spin_odd_weyl_category_source_derived": bool(
                forcing["forced_typing"] == "spin_odd_weyl_super"
            ),
            "spin_statistics_artifact_source_bound": True,
            "declared_scalar_content_source_bound": False,
            "scalar_economy_source_bound": False,
            "upstream_response_representation_source_bound": bool(
                upstream["upstream_physical_source_gate_passed"]
            ),
            "physical_refinement_intertwining_source_bound": bool(
                physical_refinement_rows
            ),
            "passed": bool(
                balance_certificate["declared_matches_derived_pair_up_to_conjugation"]
                and selection_scan["survivor_count"] == 2
                and selection_scan["survivors_equal_parity_sectors"]
                and forcing["forced_typing"] == "spin_odd_weyl_super"
                and upstream["upstream_physical_source_gate_passed"]
                and physical_refinement_rows
            ),
            "composition": {
                "passed_over": [
                    "charge_pair_derived_up_to_charge_conjugation",
                    "conjugate_projector_pair_source_derived",
                    "unordered_pair_selected_by_exhaustive_anomaly_scan",
                    "fermionic_statistics_source_derived",
                    "spin_odd_weyl_category_source_derived",
                    "spin_statistics_artifact_source_bound",
                    "upstream_response_representation_source_bound",
                    "physical_refinement_intertwining_source_bound",
                ],
                "false_by_design": {
                    "single_projector_representative_source_selected": (
                        "charge conjugation is a physical symmetry of the source; "
                        "the unordered pair is the physical object and the SELECTION_RULE "
                        "guard forbids preferring a member"
                    )
                },
                "deferred": {
                    "declared_scalar_content_source_bound": "issue 609",
                    "scalar_economy_source_bound": "issue 609",
                },
            },
        },
        "derivation_chain": [
            {
                "step": 1,
                "premise": "matter-lift manifest with typed contracts",
                "uses": ["schema check", "matter firewall", "typed contracts"],
                "source_artifact": "validate_manifest",
                "conclusion": "the source packet is admissible: trace-balanced exterior contract, fermionic statistics, spin typing, kernel emission, candidate class declaration; the contract values are matched against the derivations of steps 2a and 2b",
            },
            {
                "step": 2,
                "premise": "hash-pinned #566 manifest and receipt with the bound semantic response artifact",
                "uses": ["sha256 pins", "gate check on the stored receipt", "recomputed semantic binding requirement"],
                "source_artifact": "load_upstream",
                "conclusion": "the finite current algebra u(3) (+) so(3) is strictly upstream and its #599 exact carrier-dynamics source gate passes; laboratory gauge-current attachment is outside this packet",
            },
            {
                "step": "2a",
                "premise": "the realized exterior package with general block charges (a, b)",
                "uses": ["exact anomaly polynomials", "primitive integral block-character equation"],
                "source_artifact": "block_determinant_balance_certificate",
                "conclusion": "anomaly freedom forces 3a + 2b = 0 and primitive integral block characters give the charge-conjugate pair +/-(-1/3,1/2); the receipt's sign is conventional",
            },
            {
                "step": "2b",
                "premise": "the integer charge, triality, and duality arithmetic of the realized fields",
                "uses": ["exact solution of q_S from every admissible field pair", "triality and duality rules"],
                "source_artifact": "scalar_and_channel_selection_certificate",
                "conclusion": "conditional on a color-singlet weak-doublet scalar, compatibility forces q_S=+-3 and the three displayed channels; scalar existence, multiplicity, and economy are not derived",
            },
            {
                "step": 3,
                "premise": "the #566 source packet",
                "uses": ["frame realization", "generator rebuild", "66 structure-constant solves"],
                "source_artifact": "CurrentAlgebra",
                "conclusion": "the twelve current generators, sixty rotations, and exact structure constants are rebuilt from source",
            },
            {
                "step": 4,
                "premise": "sixty exact kernel-block rotations",
                "uses": ["exact Q(sqrt5) square roots", "adjoint transport checks", "group closure", "unique involution"],
                "source_artifact": "spin_lift_certificate",
                "conclusion": "PORT-SPIN-LIFT: a non-split SU(2) double cover (binary icosahedral) with irrational spinor characters",
            },
            {
                "step": 5,
                "premise": "source even block, spin lift, declared trace-balanced charges",
                "uses": ["66 bracket checks", "720 conjugation transports", "rank 12"],
                "source_artifact": "matter_transport",
                "conclusion": "a faithful skew-adjoint Lie algebra homomorphism onto the matter carrier V = C (+) W",
            },
            {
                "step": 6,
                "premise": "the auxiliary CAR algebra over V",
                "uses": ["50 CAR relation checks", "vacuum cyclicity rank 32", "parity grading"],
                "source_artifact": "FockSpace",
                "conclusion": "the 32-state auxiliary Fock space with derived fermionic parity and Koszul super structure",
            },
            {
                "step": 7,
                "premise": "second quantization of the matter transport",
                "uses": ["vacuum annihilation", "60 CAR-derivation checks", "parity evenness"],
                "source_artifact": "FockSpace.dgamma",
                "conclusion": "the current algebra acts on the auxiliary Fock space by parity-even super derivations",
            },
            {
                "step": 8,
                "premise": "the joint kernel of the twelve second-quantized currents",
                "uses": ["exact nullspace, complex dimension 2"],
                "source_artifact": "certificate_payload",
                "conclusion": "the derived invariant sector is exactly the vacuum and top lines; the top line is invariant by trace balance",
            },
            {
                "step": 9,
                "premise": "the two parity eigenspaces and the complete two-line invariant sector",
                "uses": ["idempotence", "self-adjointness", "equivariance", "rank 15", "exterior-complement exchange"],
                "source_artifact": "certificate_payload",
                "conclusion": "exactly two parity/invariant-line-complement projectors form a conjugate rank-fifteen pair; the even representative is printed by convention",
            },
            {
                "step": 10,
                "premise": "the realized fifteen-state module",
                "uses": ["exact charge spectrum", "faithful rank 12", "five irreducible blocks with scalar commutants", "no invariant line"],
                "source_artifact": "certificate_payload",
                "conclusion": "the exterior package is realized with fields Q, u_c, e_c, d_c, L and derived integrality normalization six",
            },
            {
                "step": 11,
                "premise": "realized operator traces on the matter module",
                "uses": ["U1^3", "grav^2 U1", "SU3^2 U1", "SU2^2 U1", "SU3^3 d-symbol", "doublet count"],
                "source_artifact": "certificate_payload",
                "conclusion": "the listed perturbative anomalies vanish exactly and the Witten parity is even (four weak doublets)",
            },
            {
                "step": 12,
                "premise": "disjoint charge spectra and the wedge pairing",
                "uses": ["exact spectral disjointness", "12 pairing invariance checks", "nondegeneracy"],
                "source_artifact": "certificate_payload",
                "conclusion": "the matter module is exactly chiral and its conjugate is realized in the opposite parity sector",
            },
            {
                "step": 13,
                "premise": "declared Yukawa channels on realized blocks",
                "uses": ["exact joint invariants per channel", "forbidden-channel control"],
                "source_artifact": "kron_sum_invariants",
                "conclusion": "the invariant Gauss sector is nonzero: exactly one invariant line per declared channel",
            },
            {
                "step": 14,
                "premise": "realized weights with triality and duality",
                "uses": ["derived integrality normalization", "exact congruence enumeration", "closure and cyclicity"],
                "source_artifact": "certificate_payload",
                "conclusion": "the common action kernel on the cover (infinite cyclic, generator (zeta_6, omega, -1), sixth power the unit deck translation, residual order six modulo pure deck translations) is emitted as data; neither the central compactification nor the global quotient is formed",
            },
            {
                "step": 15,
                "premise": "declared refinement tower",
                "uses": ["exterior lifts", "projector commutation", "generator intertwining"],
                "source_artifact": "certificate_payload",
                "conclusion": "the matter lift descends naturally along the declared algebraic tower maps",
            },
            {
                "step": 16,
                "premise": "the hash-pinned measured spin statistics artifact",
                "uses": [
                    "self-hash recompute",
                    "carrier binding equality with the response packet",
                    "order profile, centre, Klein-four section obstruction, spin structure count",
                ],
                "source_artifact": "load_spin_statistics_artifact",
                "conclusion": (
                    "the measured transport double cover is binary icosahedral, non-split over "
                    "every Klein four-subgroup, with centre {+1,-1} and a unique spin structure "
                    "on the oriented support; the internal lift recomputation matches it exactly"
                ),
            },
            {
                "step": 17,
                "premise": "the derived primitive block charges on the exterior module",
                "uses": ["exhaustive 1024-subset enumeration", "chirality and four anomaly constraints"],
                "source_artifact": "exterior_selection_scan",
                "conclusion": (
                    "exactly two subsets survive, they are the two fermionic-parity sectors and "
                    "a charge-conjugate pair; the matter grading is an output of the scan"
                ),
            },
            {
                "step": 18,
                "premise": "the derived grading and the measured central structure",
                "uses": [
                    "36-candidate gauge-centre exclusion (Lean no_universal_fermion_minus_one)",
                    "measured centre {+1,-1}",
                    "measured section obstruction",
                ],
                "source_artifact": "statistics_forcing_certificate",
                "conclusion": (
                    "the unique source implementation of the derived grading is the measured "
                    "central -1 of the non-split double cover, forcing the Spin/odd-Weyl super "
                    "typing; Vec, split sVec, and lone opposite-Weyl relabelings fail for derived reasons"
                ),
            },
            {
                "step": 19,
                "premise": "gate aggregation and finite countermodels",
                "uses": ["typed negative controls"],
                "source_artifact": "negative_controls/issue_314_negative_controls.json",
                "conclusion": (
                    "the conditional algebraic gate and the physical source gate both pass on "
                    "the reference packet and fail closed on every countermodel; scalar content "
                    "stays deferred to #609"
                ),
            },
        ],
        "factor_origins": {
            "dimensions_32_16_15": "2^5 auxiliary Fock states; the even parity sector; even minus the derived vacuum line",
            "order_120_unique_involution": "the exact SU(2) lift set of the sixty proper rotations with its single order-two element -1",
            "order_profile_1_1_20_30_24_20_24": "exact element orders of the binary icosahedral lift group (orders 1,2,3,4,5,6,10)",
            "charges_1/6_-2/3_1_1/3_-1/2": "additive weights of the declared (y_C, y_W) = (-1/3, 1/2) on the realized exterior blocks",
            "normalization_6": "least common multiple of the realized charge denominators",
            "kernel_residual_order_6": "exact count of kernel residues modulo the pure full-turn deck translations; on the cover itself the kernel is infinite cyclic, since the unit deck translation acts trivially on every integral weight without being the identity",
            "weak_doublets_4": "half the number of matter states with nonzero su(2) Casimir diagonal",
            "yukawa_lines_3": "exact joint-invariant dimensions of the three declared channels",
        },
        "branch_scope": {
            "branch": "source-bound echosahedral response branch",
            "upstream_packets": (
                "the certified #565 carrier, the source-bound #566/#599 current packet, and the "
                "measured #314 spin statistics artifact, all hash-pinned"
            ),
            "declared_branch_premises": (
                "one weak-block scalar (owned by #609), the kernel emission contract, and the candidate matter "
                "class declaration; the statistics and category contracts are validated against the "
                "scan and forcing derivations rather than accepted as premises"
            ),
            "not_claimed": (
                "no scalar existence/economy source binding, no preferred overall charge sign, "
                "no physically attached global-form choice (the #567 packet carries that), no family "
                "attachment, no class-uniqueness promotion, no scalar potential, no pole mass, no continuum "
                "spin-statistics theorem, no laboratory exchange measurement, no identification "
                "with physical particle content"
            ),
        },
        "acceptance_criteria_status": {
            "fermionic_parity_spin_lift_chirality_conjugation_tensor_product_source_derived": True,
            "current_algebra_acts_faithfully_on_matter_tensors": True,
            "exterior_package_realized_on_cover_with_anomalies_and_witten_checked": True,
            "common_action_kernel_emitted_not_assumed_as_z6_quotient": True,
            "candidate_class_proved_nonempty_before_uniqueness_promoted": True,
            "family_attachment_scalar_potential_pole_mass_outside_packet": True,
            "spin_odd_weyl_nonempty_and_vec_svec_opposite_weyl_controls_fail": True,
            "faithful_action_and_nonzero_invariant_sector_with_gauss_and_kernel_gates": True,
            "fifteen_state_module_selected_by_derived_equivariant_projector": True,
            "fifteen_state_module_selected_by_exhaustive_anomaly_scan": True,
            "projector_pair_not_single_physical_weyl_selection": True,
            "scalar_economy_explicitly_not_derived": True,
        },
        "issue_closure_condition": {
            "produced_locally": (
                "the exact matter lift on the finite current packet at source scope: the non-split "
                "PORT-SPIN-LIFT cross-checked against the measured transport artifact, the faithful "
                "current action, the unordered conjugate rank-fifteen pair selected by the exhaustive "
                "1024-subset anomaly scan with the parity grading as an output, the statistics/category "
                "typing forced by the measured centre and section obstruction, realized anomaly and "
                "Witten checks, chirality, conjugation, Yukawa invariant lines, the emitted action "
                "kernel (infinite cyclic on the cover, residual order six), declared-tower descent, the "
                "artifact-bound physical refinement descent, the derived BLOCK-DETERMINANT-BALANCE "
                "charge pair up to conjugation, and scalar/channel compatibility"
            ),
            "branch_premises": (
                "the hash-pinned #565/#566 packets with the source-bound #599 semantic response "
                "artifact and the measured #314 spin statistics artifact; the matter charges, the "
                "unordered projector pair, the parity grading, and the Spin/odd-Weyl typing are "
                "derived; declared scalar content remains a typed branch premise owned by #609"
            ),
            "conditional_algebraic_gate_passed": True,
            "physical_source_realization_gate_passed": True,
            "met_locally": True,
            "remaining_open_lanes": {
                "scalar_content_and_economy": "issue 609 (no-extra-light-sector theorem)",
                "family_and_laboratory_attachment": "issue 569",
                "physical_global_form": "issue 567 (own packet)",
                "continuum_spin_statistics_theorem": "continuum lane, not claimed",
            },
        },
        "dependency_acyclicity_note": {
            "upstream": [
                "manifests/echosahedral_federation_reference.json (#565 carrier packet)",
                "manifests/port_current_response_reference.json and receipts/port_current_inner_reference.receipt.json (#566)",
            ],
            "downstream": [
                "AXIS-CENTER-DESCENT consumes the emitted action kernel; a5_screen_sm_closure.py and exterior_sm_completion.py reference this closure in their gate ledgers",
            ],
            "summary": "carrier packet -> current packet -> matter lift receipt -> ledger references; the graph is acyclic",
        },
        "verifier_command": (
            "python3 code/a5_closure/super_tannakian_matter_lift_certificate.py verify "
            "--manifest code/a5_closure/manifests/super_tannakian_matter_reference.json "
            "--receipt code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json"
        ),
        "claim_boundary": {
            "proves": (
                "the source-bound super-Tannakian matter lift at finite source-model scope: "
                "PORT-SPIN-LIFT cross-checked against the measured transport artifact, the "
                "unordered conjugate rank-fifteen pair selected by the exhaustive anomaly scan "
                "with the parity grading as an output, the statistics/category typing forced by "
                "the measured centre and section obstruction, anomaly balance up to charge "
                "conjugation, compatible scalar charges/channels, and refinement intertwining"
            ),
            "status": "source_bound_at_finite_scope_conditional_on_declared_scalar_content",
            "contract_provenance": (
                "BLOCK-DETERMINANT-BALANCE derives the primitive charge-conjugate pair; the "
                "exhaustive 1024-subset scan selects the unordered pair and outputs the parity "
                "grading; the measured lift centre and section obstruction force the "
                "Spin/odd-Weyl implementation, with the gauge-centre branch excluded by the "
                "Lean fermion-parity no-go; scalar-charge and channel compatibility are derived "
                "only conditional on a declared weak-doublet scalar, whose existence and economy "
                "are owned by #609"
            ),
            "does_not_close": [
                "physical AXIS-CENTER-DESCENT/global-form attachment (the #567 packet carries that)",
                "class uniqueness (only nonemptiness is discharged here)",
                "A5-FAMILY-ATTACHMENT, family structure, and any three-family claim (#569)",
                "exclusion of other anomaly-free light sectors beyond the exterior module (#609)",
                "scalar existence, one-scalar economy, or absence of additional light scalars (#609)",
                "laboratory measurement of any matter observable",
                "scalar potential, pole masses, measured couplings, continuum spin-statistics, or quantum field theory",
            ],
        },
    }


# ---------------------------------------------------------------------------
# Negative controls
# ---------------------------------------------------------------------------


def negative_control_cases(manifest: Mapping[str, Any]) -> list[tuple[str, dict[str, Any], str]]:
    cases: list[tuple[str, dict[str, Any], str]] = []

    vec = copy.deepcopy(manifest)
    vec["category_contract"]["typing"] = "vec"
    cases.append(("vec_typing", vec, "VEC_TYPING"))

    svec = copy.deepcopy(manifest)
    svec["category_contract"]["typing"] = "svec"
    svec["category_contract"]["spin_lift"]["double_cover"] = False
    cases.append(("svec_split_spin", svec, "SPIN_LIFT_SPLIT"))

    opposite = copy.deepcopy(manifest)
    opposite["category_contract"]["selection_rule"] = "parity_odd_minus_derived_invariants"
    cases.append(("mixed_weyl_representative", opposite, "YUKAWA_CHANNEL_EMPTY"))

    bosonic = copy.deepcopy(manifest)
    bosonic["statistics_contract"]["matter_statistics"] = "bosonic_even"
    cases.append(("bosonic_matter_statistics", bosonic, "STATISTICS_TYPING"))

    truncated = copy.deepcopy(manifest)
    truncated["category_contract"]["selection_rule"] = "lambda2_only"
    cases.append(("truncated_lambda2_selection", truncated, "WITTEN_PARITY"))

    full_even = copy.deepcopy(manifest)
    full_even["category_contract"]["selection_rule"] = "even_including_vacuum"
    cases.append(("full_even_clifford_module", full_even, "TRIVIAL_LINE_IN_MATTER"))

    empty_gauss = copy.deepcopy(manifest)
    empty_gauss["exterior_matter_contract"]["yukawa_channels"] = []
    cases.append(("empty_gauss_data", empty_gauss, "GAUSS_DATA_EMPTY"))

    assumed = copy.deepcopy(manifest)
    assumed["kernel_emission_contract"]["assume_global_quotient"] = True
    cases.append(("assumed_global_quotient", assumed, "KERNEL_EMISSION_CONTRACT"))

    killing = copy.deepcopy(manifest)
    killing["exterior_matter_contract"]["extra_scalars"] = [
        {"label": "kernel_killing_singlet", "charge": "1/6", "triality": 0, "duality": 0}
    ]
    cases.append(("kernel_killing_extra_scalar", killing, "KERNEL_TRIVIAL"))

    arithmetic = copy.deepcopy(manifest)
    arithmetic["category_contract"]["realization"] = "representation_arithmetic"
    cases.append(("representation_arithmetic_only", arithmetic, "REALIZATION_NOT_OPERATOR"))

    charge_dead = copy.deepcopy(manifest)
    charge_dead["exterior_matter_contract"]["block_trace_charges"] = {
        "color_block": "0",
        "weak_block": "0",
    }
    # The zero pair sits on the balance line at t = 0, so the derivation
    # rejects it as non-primitive before the faithfulness check runs.
    cases.append(("charge_dead_package", charge_dead, "BLOCK_DETERMINANT_BALANCE"))

    unbalanced = copy.deepcopy(manifest)
    unbalanced["exterior_matter_contract"]["block_trace_charges"]["weak_block"] = "1/3"
    cases.append(("unbalanced_trace_charges", unbalanced, "TRACE_BALANCE"))

    non_primitive = copy.deepcopy(manifest)
    non_primitive["exterior_matter_contract"]["block_trace_charges"] = {
        "color_block": "-2/3",
        "weak_block": "1",
    }
    cases.append(
        ("non_primitive_balanced_pair", non_primitive, "BLOCK_DETERMINANT_BALANCE")
    )

    preferred_projector = copy.deepcopy(manifest)
    preferred_projector["category_contract"][
        "projector_representative_convention"
    ] = "even_minus_vacuum"
    cases.append(
        ("preferred_projector_representative", preferred_projector, "SELECTION_RULE")
    )

    extra_channel = copy.deepcopy(manifest)
    extra_channel["exterior_matter_contract"]["yukawa_channels"] = [
        ["Q", "S", "u_c"],
        ["Q", "Sbar", "d_c"],
        ["L", "Sbar", "e_c"],
        ["Q", "S", "d_c"],
    ]
    cases.append(("undeclared_forbidden_channel", extra_channel, "SCALAR_SELECTION"))

    promoted = copy.deepcopy(manifest)
    promoted["candidate_matter_class"]["promote_uniqueness"] = True
    cases.append(("class_uniqueness_promoted", promoted, "CLASS_UNIQUENESS_PROMOTION"))

    family = copy.deepcopy(manifest)
    family["downstream_hint"] = {"attachment_target": "three family attachment"}
    cases.append(("family_attachment_injection", family, "FORBIDDEN_DEPENDENCY"))

    potential = copy.deepcopy(manifest)
    potential["downstream_hint"] = {"scalar_sector": "scalar potential quartic and pole mass"}
    cases.append(("scalar_potential_injection", potential, "FORBIDDEN_DEPENDENCY"))

    artifact_drift = copy.deepcopy(manifest)
    artifact_drift["spin_statistics_artifact_sha256"] = "sha256:" + "0" * 64
    cases.append(("spin_artifact_pin_drift", artifact_drift, "UPSTREAM_HASH"))

    wrong_artifact = copy.deepcopy(manifest)
    wrong_artifact["spin_statistics_artifact_path"] = str(
        manifest.get("current_manifest_path")
    )
    wrong_artifact["spin_statistics_artifact_sha256"] = str(
        manifest.get("current_manifest_sha256")
    )
    cases.append(("wrong_spin_artifact_schema", wrong_artifact, "UPSTREAM_HASH"))

    return cases


def negative_control_payload(manifest: Mapping[str, Any], base_dir: Path | None = None) -> dict[str, Any]:
    results: list[dict[str, Any]] = []
    for name, mutant, expected_code in negative_control_cases(manifest):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir, allow_control_contracts=True)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append({"name": name, "expected_error": expected_code, "actual_error": actual_code, "passed": True})
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 314,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "vec_and_svec": {
                "vec": "an ungraded category cannot carry the derived nontrivial parity and non-split spin cocycle",
                "svec_split": "a split lift would give the 120-element group more than one involution; the derived group has exactly one",
                "conclusion": "the Vec and sVec same-reduct typings fail against derived facts, not by fiat",
            },
            "mixed_weyl_convention": {
                "valid_pair": "the odd-minus-top projector is the valid conjugate of even-minus-vacuum",
                "failure": "changing only the projector while retaining the even-representative charges and channel labels is an inconsistent mixed convention",
            },
            "selection": {
                "lambda2_only": "ten states with three weak doublets: odd Witten parity",
                "even_including_vacuum": "the full even Clifford module keeps the trivial vacuum line inside matter",
            },
            "kernel": {
                "assumed_quotient": "assuming the global quotient violates the emission contract",
                "kernel_killing_scalar": "an extra integral-charge singlet collapses the kernel to the pure full-turn deck translations; a trivial residual cannot satisfy the packet",
            },
            "typing": {
                "bosonic_matter": "the derived matter-building operators anticommute; bosonic typing fails closed",
                "representation_arithmetic": "representation arithmetic without the operator projector is not physical realization",
                "uniqueness_promotion": "promoting class uniqueness inside this packet is rejected; only nonemptiness is discharged",
                "charge_dead": "a charge-dead contract annihilates the central lane, so the current action is unfaithful",
            },
        },
    }


def verify_receipt(manifest: Mapping[str, Any], receipt: Mapping[str, Any], base_dir: Path | None = None) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json",
        MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_314_negative_controls.json",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify", help="create the deterministic exact receipt")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify", help="recompute and compare a receipt")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls", help="run and write the finite countermodel bundle")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all", help="regenerate receipt and negative controls at repository-default paths")
    all_cmd.add_argument("--manifest", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output), "sha256": sha256_json(receipt)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS", "negative_controls": str(args.output)}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {"status": "PASS", "receipt": str(default_receipt), "negative_controls": str(default_negative)},
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
