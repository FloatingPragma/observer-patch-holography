#!/usr/bin/env python3
"""Exact certificate for GitHub issue #314: the conditional super-Tannakian matter lift.

The input is a matter-lift manifest.  It declares only:

* the conditional port-current response manifest of issue #566 (by path
  and hash) together with its stored receipt hash - the conditional
  current algebra u(3) (+) so(3) on the charged response space is
  strictly upstream, and its physical source binding is open in #599;
* the trace-balanced exterior matter contract: the exact charge pair on
  the color and weak blocks, the one-scalar choice (the weak block
  itself), and the declared invariant Yukawa channel list - a typed
  branch premise, not a physical measurement;
* the fermionic statistics contract, the Spin/odd-Weyl category typing
  with a genuine double-cover spin lift, the operator-projector
  realization contract, the kernel emission contract, and the declared
  MAR class with uniqueness promotion switched off.

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
* the derived equivariant projector (parity projector minus the derived
  invariant line) of exact rank fifteen whose image realizes the exterior
  package Lambda^2 V (+) Lambda^4 V;
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

This is a conditional exact algebraic theorem.  It does not source-bind
the matter contract or the upstream response representation to physical
carrier response; those remain separate source-binding gates.  The
receipt therefore separates the verified conditional algebraic gate
from the open physical source-realization gate, and issue #314 remains
open at that physical boundary.

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

SCHEMA = "oph.super_tannakian_matter_manifest.v2"
RECEIPT_SCHEMA = "oph.super_tannakian_matter_receipt.v2"
NEGATIVE_SCHEMA = "oph.super_tannakian_matter_negative_controls.v2"

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
    selection_rule = category.get("selection_rule")
    allowed_rules = (
        SELECTION_RULES if allow_control_contracts else ("parity_even_minus_derived_invariants",)
    )
    require(
        selection_rule in allowed_rules,
        "SELECTION_RULE",
        "the production packet must declare the parity_even_minus_derived_invariants selection rule",
    )
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

    mar = manifest.get("mar_class")
    require(isinstance(mar, Mapping), "MAR_CONTRACT", "mar_class is missing")
    require(
        mar.get("declared") == "one_generation_one_scalar_chiral_anomaly_free",
        "MAR_CONTRACT",
        "the declared MAR class must be the one-generation one-scalar chiral anomaly-free class",
    )
    require(
        mar.get("promote_uniqueness") is False,
        "MAR_UNIQUENESS_PROMOTION",
        "MAR uniqueness may not be promoted before the class is proved nonempty and the uniqueness lane runs",
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
        isinstance(physical_gate, Mapping) and physical_gate.get("passed") is False,
        "UPSTREAM_RECEIPT",
        "the upstream premise-binding status must be recorded honestly: the #566 response premises are "
        "declared, and their source binding is tracked in #599, so the stored gate must read open",
    )
    return {
        "current_manifest": current_manifest,
        "current_manifest_sha256": digest,
        "current_receipt_sha256": sha256_json(current_receipt),
        "carrier_manifest_sha256": current_receipt.get("carrier_manifest_sha256"),
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

    return {
        "lifts": lifts,
        "witness_count": len(lifts),
        "lift_group_order": len(group),
        "unique_involution": True,
        "order_profile": {str(k): v for k, v in sorted(order_profile.items())},
        "involution_lift_order": 4,
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


def certificate_payload(
    manifest: Mapping[str, Any],
    base_dir: Path | None = None,
    *,
    allow_control_contracts: bool = False,
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    params = validate_manifest(manifest, allow_control_contracts=allow_control_contracts)
    upstream = load_upstream(manifest, base)

    algebra = CurrentAlgebra(upstream["current_manifest"], base)

    # --- PORT-SPIN-LIFT ------------------------------------------------------
    spin = spin_lift_certificate(algebra)
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
    tower = algebra.carrier_manifest["refinement_tower"]
    refinement_rows = []
    for item in tower["maps"]:
        permutation = e565.parse_port_permutation(item["port_map"], algebra.carrier)
        rotation = algebra.frame.rotation_of(permutation)
        kernel_rotation = [[entry.conj() for entry in row] for row in rotation]
        lift = spin["lifts"].get(tuple(permutation)) or spin_lift_of_rotation(kernel_rotation)
        pi_v = carrier_implementer(rotation, lift)
        gamma = fock.exterior_lift(pi_v)
        require(
            c_is_zero(csub(cmul(gamma, projector), cmul(projector, gamma))),
            "REFINEMENT_DESCENT",
            "a refinement map does not commute with the matter selection projector",
        )
        pi_dagger = cdagger(pi_v)
        gamma_dagger = cdagger(gamma)
        for p in range(12):
            target_index = permutation[p]
            conjugated_v = cmul(cmul(pi_v, [list(row) for row in images[p]]), pi_dagger)
            require(
                c_is_zero(csub(conjugated_v, [list(row) for row in images[target_index]])),
                "REFINEMENT_DESCENT",
                "a refinement map is not intertwined on the matter carrier",
            )
            conjugated_f = cmul(cmul(gamma, dgammas[p]), gamma_dagger)
            require(
                c_is_zero(csub(conjugated_f, dgammas[target_index])),
                "REFINEMENT_DESCENT",
                "a refinement map is not intertwined on the Fock realization",
            )
        refinement_rows.append({"source": item["source"], "target": item["target"], "intertwined": True})

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
        "mar_class_nonempty_witnessed": True,
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
            "dependencies": ["#565 (carrier packet)", "#566 (conditional port-current algebra)"],
            "inherited_scope": (
                "the upstream response premises are declared in the closed #566 packet; their physical "
                "source binding is tracked in #599, which is not a dependency of this issue"
            ),
        },
        "source_firewall": {
            "forbidden_dependency_hits": [],
            "uses_only": [
                "hash-pinned #566 conditional current packet (manifest and receipt)",
                "declared trace-balanced exterior matter contract (branch premise)",
                "declared fermionic statistics and Spin/odd-Weyl category contracts (branch premises)",
                "declared kernel emission contract and MAR class declaration (branch premises)",
            ],
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
                "the trace-balanced redistribution of the central charge onto (y_C, y_W) is the declared "
                "conditional exterior matter contract; BLOCK-DETERMINANT-BALANCE (the physical selection of "
                "trace balance) is not closed here"
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
            "projector": "parity projector minus the derived invariant line, diagonal in the subset basis, exact rank 15",
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
            "scope": "naturality along the declared algebraic tower maps; physical refinement intertwining is not source-bound here",
        },
        "mar_class": {
            "declared": "one_generation_one_scalar_chiral_anomaly_free",
            "nonempty": True,
            "witness": (
                "the realized packet: fifteen multiplicity-free states (one generation), one declared scalar, "
                "exact chirality, vanishing realized anomaly traces"
            ),
            "uniqueness_promoted": False,
            "note": "this packet discharges the nonemptiness precondition; MAR uniqueness stays in its own lane",
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
        "conditional_algebraic_gate": {**gate, "passed": True},
        "physical_source_gate": {
            "matter_contract_source_bound": False,
            "upstream_response_representation_source_bound": False,
            "physical_refinement_intertwining_source_bound": False,
            "passed": False,
        },
        "derivation_chain": [
            {
                "step": 1,
                "premise": "declared matter-lift manifest",
                "uses": ["schema check", "matter firewall", "typed contracts"],
                "source_artifact": "validate_manifest",
                "conclusion": "the source packet is admissible: trace-balanced exterior contract, fermionic statistics, spin typing, kernel emission, MAR declaration",
            },
            {
                "step": 2,
                "premise": "hash-pinned #566 manifest and receipt",
                "uses": ["sha256 pins", "gate check on the stored receipt"],
                "source_artifact": "load_upstream",
                "conclusion": "the conditional current algebra u(3) (+) so(3) with its charged response space is strictly upstream; its physical source gate is recorded open in #599",
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
                "premise": "parity projector minus the derived invariant line",
                "uses": ["idempotence", "self-adjointness", "12 equivariance checks", "rank 15"],
                "source_artifact": "certificate_payload",
                "conclusion": "the derived equivariant projector selects the fifteen-state matter module from the full auxiliary Fock space",
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
                "premise": "gate aggregation and finite countermodels",
                "uses": ["typed negative controls"],
                "source_artifact": "negative_controls/issue_314_negative_controls.json",
                "conclusion": "the conditional algebraic gate passes on the reference packet and fails closed on every countermodel; the physical source gate is recorded open",
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
            "branch": "declared echosahedral response branch",
            "upstream_packets": "the certified #565 carrier and the conditional #566 current algebra, both hash-pinned",
            "declared_branch_premises": (
                "the trace-balanced exterior matter contract (block trace charges, the one-scalar choice, the "
                "Yukawa channel list), the fermionic statistics contract, the Spin/odd-Weyl category typing, "
                "the kernel emission contract, and the MAR class declaration - typed premises, not measurements"
            ),
            "not_claimed": (
                "no physical source binding of the matter contract or the upstream response representation, "
                "no derivation of trace balance from source dynamics (BLOCK-DETERMINANT-BALANCE open), no global-form "
                "choice (AXIS-CENTER-DESCENT open), no family attachment, no MAR uniqueness, no scalar potential, no "
                "pole mass, no continuum spin-statistics theorem, no identification with physical particle content"
            ),
        },
        "acceptance_criteria_status": {
            "fermionic_parity_spin_lift_chirality_conjugation_tensor_product_source_derived": False,
            "current_algebra_acts_faithfully_on_matter_tensors": True,
            "exterior_package_realized_on_cover_with_anomalies_and_witten_checked": True,
            "common_action_kernel_emitted_not_assumed_as_z6_quotient": True,
            "mar_class_proved_nonempty_before_uniqueness_promoted": True,
            "family_attachment_scalar_potential_pole_mass_outside_packet": True,
            "spin_odd_weyl_nonempty_and_vec_svec_opposite_weyl_controls_fail": True,
            "faithful_action_and_nonzero_invariant_sector_with_gauss_and_kernel_gates": True,
            "fifteen_state_module_selected_by_derived_equivariant_projector": True,
        },
        "issue_closure_condition": {
            "produced_locally": (
                "the conditional exact matter lift: the non-split algebraic PORT-SPIN-LIFT target, faithful "
                "current action, derived equivariant selection of the fifteen-state module, realized anomaly "
                "and Witten checks, chirality, conjugation, Yukawa invariant lines, the emitted action "
                "kernel (infinite cyclic on the cover, residual order six), and declared-tower descent"
            ),
            "branch_premises": (
                "the hash-pinned #565/#566 packets plus the declared matter-lift contracts (trace-balanced "
                "block charges, one scalar, channel list, statistics, category typing, kernel emission, MAR class)"
            ),
            "conditional_algebraic_gate_passed": True,
            "physical_source_realization_gate_passed": False,
            "met_locally": False,
            "remaining_producer": (
                "physical source binding: #599 must source-bind the upstream response representation, "
                "the trace-balanced matter contract must be physically selected, and the physical refinement "
                "maps must be derived from carrier response before the issue's source-derived matter category closes"
            ),
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
            "proves": "the conditional exact super-Tannakian matter lift for the declared matter contracts over the pinned conditional current packet, including the algebraic PORT-SPIN-LIFT target",
            "status": "proved_conditional_on_declared_matter_contracts",
            "declared_branch_premises": (
                "the trace-balanced exterior matter contract, statistics and category typing, kernel emission "
                "contract, and MAR class declaration enter as typed branch premises, not as physical measurements"
            ),
            "does_not_close": [
                "source binding of the upstream response representation and coefficients (tracked in #599)",
                "PORT-SPIN-LIFT beyond its algebraic target (inherits the #599 premise binding)",
                "BLOCK-DETERMINANT-BALANCE (physical selection of the trace-balanced charge pair)",
                "physical refinement intertwining beyond the declared algebraic tower maps",
                "AXIS-CENTER-DESCENT and the global form (the kernel is emitted, the quotient is not chosen)",
                "MAR uniqueness (only nonemptiness is discharged here)",
                "A5-FAMILY-ATTACHMENT, family structure, and any three-family claim",
                "exclusion of other anomaly-free light sectors (MGFC-grade no-extra-sector)",
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
    cases.append(("opposite_weyl_selection", opposite, "YUKAWA_CHANNEL_EMPTY"))

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
    cases.append(("charge_dead_package", charge_dead, "CURRENT_ACTION_NOT_FAITHFUL"))

    unbalanced = copy.deepcopy(manifest)
    unbalanced["exterior_matter_contract"]["block_trace_charges"]["weak_block"] = "1/3"
    cases.append(("unbalanced_trace_charges", unbalanced, "TRACE_BALANCE"))

    promoted = copy.deepcopy(manifest)
    promoted["mar_class"]["promote_uniqueness"] = True
    cases.append(("mar_uniqueness_promoted", promoted, "MAR_UNIQUENESS_PROMOTION"))

    family = copy.deepcopy(manifest)
    family["downstream_hint"] = {"attachment_target": "three family attachment"}
    cases.append(("family_attachment_injection", family, "FORBIDDEN_DEPENDENCY"))

    potential = copy.deepcopy(manifest)
    potential["downstream_hint"] = {"scalar_sector": "scalar potential quartic and pole mass"}
    cases.append(("scalar_potential_injection", potential, "FORBIDDEN_DEPENDENCY"))

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
            "opposite_weyl": {
                "same_reduct": "the odd-parity selection realizes the conjugate module",
                "failure": "every declared one-scalar Yukawa channel then has charge sum away from zero, so the Gauss sector is empty",
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
                "uniqueness_promotion": "promoting MAR uniqueness inside this packet is rejected; only nonemptiness is discharged",
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
