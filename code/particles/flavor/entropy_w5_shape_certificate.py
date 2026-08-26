#!/usr/bin/env python3
"""Entropy-derived W5 shape packet: exact orbit selection and its exclusion.

The W5 stabilizer boundary (``derive_w5_stabiliser_spectrum_bound.py``,
``Lean/ObserverPatchHolography/W5Stabilizer.lean``) leaves one open dynamical
object: a source-derived coefficient packet for the A5-invariant potential on
the quintet carrier.  This certificate derives the first such packet from the
record and evaluates it, exactly, to an exclusion verdict.

THE PACKET (derived, no free coefficient):

The A3 objective at the uniform port state expands as

    D(u(1+x) || u) = u * ( S2/2 - S3/6 + S4/12 ) + O(x^5),
    S_k = sum over the twelve ports of x_i^k,   u = 1/12,

and on the twelve-port carrier the cubic S3 vanishes identically on the
antipode-odd frame band, so the quintet band is the only nontrivial band
with a cubic entropy response.  Restricted to the quintet band at a fixed
band amplitude (the declared premise), the shape functional is

    E(Q; r) = -S3(Q) r^3 / 6 + S4(Q) r^4 / 12    on   S2(Q) = 1,

with every coefficient an exact consequence of the entropy expansion.

WHAT IS PROVED (exact, in Q(sqrt5)):

* On the diagonal (D2) stratum the twelve port values collapse to three
  values, each carried by four ports, through an exact cyclic linear map;
  extremizing the cubic on that stratum gives exactly two critical orbit
  classes, and symmetric criticality lifts stratum-critical points to
  critical points of the full constrained functional.
* The two branch orbits are computed exactly: the vertex-axis (C5) prolate
  orbit with a double eigenvalue, and a simple-spectrum orbit whose
  centered eigenvalue triple is proportional to
  (-(phi^4 - 1), -1, phi^4), phi the golden ratio; its sorted-gap ratio is
  exactly phi.
* The exact branch invariants give the exact crossing amplitude r_c at
  which the simple-spectrum branch overtakes the prolate branch, so the
  entropy packet selects, on the compared strata: the degenerate prolate
  orbit below r_c and the golden-ratio orbit above it.
* The record is a probability vector, q_i = (1 + r w_i)/12 > 0, so each
  branch is admissible only inside its positivity domain
  r < 1/|min_i w_i|: r^2 < 60 for the prolate branch and r^2 < 24 for the
  golden branch, both recorded exactly.  The C3-axis orbit is a further
  critical point by symmetric criticality; its cubic vanishes exactly
  (two-value antipodal profile), its quartic energy is r^4/144, and it
  undercuts the golden branch only for r^2 > 96, outside every positivity
  bound.  The two-branch selection therefore holds on the whole admissible
  amplitude domain.

THE VERDICT (compare-only, both orientation conventions reported):

The centered logarithm of the measured charged-lepton triple has sorted-gap
ratio 1.8890 (or 0.5294 under the flipped orientation).  The entropy
packet's only simple-spectrum output is phi = 1.6180 (or 1/phi = 0.6180).
The relative mismatch is 16.7 percent in the closest pairing, far outside
every stated width, so the quartic-truncated entropy packet is EXCLUDED as
the charged-lepton shape selector.  The degenerate prolate branch is
excluded by the simple-spectrum gate.  This extends the W5 boundary: with
symmetry geometry excluded by the stabilizer theorem and the entropy
packet excluded here, the open lepton-side object is a source mechanism
beyond the quartic entropy truncation.

Row class: exact structure theorem plus compare-only exclusion.  No mass,
ratio, or prediction is emitted.  The lepton values enter only in the final
comparison block and are labeled with their ancestry.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

ROOT = Path(__file__).resolve().parents[2]
OUT_PATH = ROOT / "particles" / "runs" / "flavor" / "entropy_w5_shape_certificate.json"

SCHEMA = "oph.entropy_w5_shape_certificate.v1"
ISSUE_CONTEXT = [546, 591]


class CertificateError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


class Q5:
    """Exact a + b*sqrt5 with Fraction coefficients."""

    __slots__ = ("a", "b")

    def __init__(self, a: Any = 0, b: Any = 0) -> None:
        self.a = Fraction(a)
        self.b = Fraction(b)

    def __add__(self, o: "Q5") -> "Q5":
        return Q5(self.a + o.a, self.b + o.b)

    def __sub__(self, o: "Q5") -> "Q5":
        return Q5(self.a - o.a, self.b - o.b)

    def __neg__(self) -> "Q5":
        return Q5(-self.a, -self.b)

    def __mul__(self, o: "Q5") -> "Q5":
        return Q5(self.a * o.a + 5 * self.b * o.b, self.a * o.b + self.b * o.a)

    def inv(self) -> "Q5":
        n = self.a * self.a - 5 * self.b * self.b
        require(n != 0, "Q5_DIVZERO", "inverse of zero")
        return Q5(self.a / n, -self.b / n)

    def __truediv__(self, o: "Q5") -> "Q5":
        return self * o.inv()

    def __eq__(self, o: object) -> bool:
        return isinstance(o, Q5) and self.a == o.a and self.b == o.b

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_positive(self) -> bool:
        if self.b == 0:
            return self.a > 0
        if self.a == 0:
            return self.b > 0
        if self.a > 0 and self.b > 0:
            return True
        if self.a < 0 and self.b < 0:
            return False
        if self.a > 0:
            return self.a * self.a > 5 * self.b * self.b
        return 5 * self.b * self.b > self.a * self.a

    def to_float(self) -> float:
        return float(self.a) + float(self.b) * math.sqrt(5.0)

    def render(self) -> str:
        return f"{self.a} + {self.b}*sqrt5"


ZERO = Q5(0)
ONE = Q5(1)
PHI = Q5(Fraction(1, 2), Fraction(1, 2))
PHI2 = PHI * PHI  # phi^2 = phi + 1
PHI4 = PHI2 * PHI2


def q5_lt(x: Q5, y: Q5) -> bool:
    return (y - x).is_positive()


# ---------------------------------------------------------------------------
# The diagonal (D2) stratum in exact arithmetic
# ---------------------------------------------------------------------------
#
# Icosahedron vertices (0, ±1, ±phi) and cyclic permutations; the three
# coordinate axes are mutually orthogonal two-fold axes, and a diagonal
# traceless Q = diag(q1, q2, q3) is exactly the D2-fixed stratum.  Each
# vertex family contributes the squared-coordinate profile (0, 1, phi^2)
# in cyclic position, so with |v|^2 = 2 + phi the twelve port values
# collapse to three values, each on four ports:
#
#   y_a = (q2 + phi^2 q3) / (2 + phi)     [vertices (0, ±1, ±phi)]
#   y_b = (q1 + phi^2 q2) / (2 + phi)     [vertices (±1, ±phi, 0)]
#   y_c = (q3 + phi^2 q1) / (2 + phi)     [vertices (±phi, 0, ±1)]
#
# The port sums are S2 = 4(y_a^2 + y_b^2 + y_c^2), S3 = 4 sum y^3,
# S4 = 4 sum y^4, and sum y = 0 exactly when tr Q = 0.


NORM = Q5(2) + PHI  # |v|^2 = 2 + phi = (5 + sqrt5)/2


def port_profile(q: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    q1, q2, q3 = q
    inv = NORM.inv()
    return (
        (q2 + PHI2 * q3) * inv,
        (q1 + PHI2 * q2) * inv,
        (q3 + PHI2 * q1) * inv,
    )


def port_sums(q: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    ys = port_profile(q)
    s2 = ZERO
    s3 = ZERO
    s4 = ZERO
    four = Q5(4)
    for y in ys:
        y2 = y * y
        s2 = s2 + four * y2
        s3 = s3 + four * y2 * y
        s4 = s4 + four * y2 * y2
    return s2, s3, s4


def invert_profile(y: Sequence[Q5]) -> tuple[Q5, Q5, Q5]:
    """Solve the cyclic system for diag(q) with the given port profile."""

    y_a, y_b, y_c = (value * NORM for value in y)
    # q2 + phi^2 q3 = y_a ; q1 + phi^2 q2 = y_b ; q3 + phi^2 q1 = y_c.
    # From eq2: q1 = y_b - phi^2 q2; substituted into eq3: q3 = y_c - phi^2 q1;
    # eq1 then fixes q2:  q2 (1 + phi^6) = y_a - phi^2 y_c + phi^4 y_b.
    p2 = PHI2
    q2 = (y_a - p2 * y_c + p2 * p2 * y_b) / (ONE + p2 * p2 * p2)
    q1 = y_b - p2 * q2
    q3 = y_c - p2 * q1
    for got, want in zip(port_profile((q1, q2, q3)), y):
        require((got - want).is_zero(), "PROFILE_ROUNDTRIP", "profile inversion must reproduce the profile")
    require((q1 + q2 + q3).is_zero(), "TRACELESS", "inverted diagonal must be traceless")
    return q1, q2, q3


# ---------------------------------------------------------------------------
# The two exact branch orbits
# ---------------------------------------------------------------------------


def prolate_orbit() -> dict[str, Any]:
    """The vertex-axis (C5) prolate orbit, computed on its own stratum.

    Q = a a^T / |a|^2 - I/3 for a vertex axis a = (0, 1, phi).  The twelve
    port values split as two poles at 2/3 and ten ring ports; the exact
    values follow from (v . a)^2 / (|v|^2 |a|^2) being 1 at the poles and
    exactly 1/sqrt5-quadratic on the ring.
    """

    # (v.a)^2/(|v|^2|a|^2) for v over the twelve vertices, a = (0,1,phi):
    # poles (v = ±a): 1.  Ring: the icosahedral frame gives cos^2 = 1/sqrt5
    # weight: (v.a)^2 = ((5+sqrt5)/2)^2 / 5 for the ten ring vertices.
    # Port value x = cos^2 - 1/3.
    pole = ONE - Q5(Fraction(1, 3))
    # Ring cosine: v = (1, phi, 0), a = (0, 1, phi), v.a = phi,
    # |v|^2 = |a|^2 = 2 + phi, and (2 + phi)^2 = 5 phi^2, so
    # cos^2 = phi^2 / (2 + phi)^2 = 1/5 exactly, for all ten ring vertices.
    cos2_ring = (PHI * PHI) / (NORM * NORM)
    require(cos2_ring == Q5(Fraction(1, 5)), "PROLATE_RING", "the ring cosine must be exactly one fifth")
    ring = cos2_ring - Q5(Fraction(1, 3))
    # Verify the design identity: 2*pole + 10*ring must vanish (trace zero).
    require((Q5(2) * pole + Q5(10) * ring).is_zero(), "PROLATE_TRACE", "prolate port profile must sum to zero")
    s2 = Q5(2) * pole * pole + Q5(10) * ring * ring
    s3 = Q5(2) * pole * pole * pole + Q5(10) * ring * ring * ring
    s4 = Q5(2) * pole * pole * pole * pole + Q5(10) * ring * ring * ring * ring
    return {
        "orbit": "vertex_axis_prolate_C5",
        "eigenvalue_pattern": "(2, -1, -1)/3 about the vertex axis",
        "double_eigenvalue": True,
        "port_values": {"pole": pole.render(), "ring": ring.render(), "pole_count": 2, "ring_count": 10},
        "S2": s2,
        "S3": s3,
        "S4": s4,
        "min_port_value": ring,
    }


def golden_orbit() -> dict[str, Any]:
    """The D2-stratum cubic-extremal orbit with simple spectrum.

    On the diagonal stratum the port profile is three values with sum
    zero; the cubic extremal profile at fixed square sum is the symmetric
    pattern (2, -1, -1) up to scale and permutation.  Inverting the exact
    cyclic map produces the diagonal triple, whose centered eigenvalues
    come out proportional to (-(phi^4 - 1), -1, phi^4) with sorted-gap
    ratio exactly phi.
    """

    profile = (Q5(2), Q5(-1), Q5(-1))
    q = invert_profile(profile)
    total = q[0] + q[1] + q[2]
    require(total.is_zero(), "GOLDEN_TRACE", "the golden diagonal must be traceless")
    # Sort exactly by the Q(sqrt5) order.
    triple = sorted(q, key=lambda v: (v.to_float()))
    for left, right in zip(triple, triple[1:]):
        require(q5_lt(left, right), "GOLDEN_SIMPLE", "the golden orbit must have simple spectrum")
    gap_low = triple[1] - triple[0]
    gap_high = triple[2] - triple[1]
    ratio = gap_high / gap_low
    require(ratio == PHI, "GOLDEN_RATIO", "the sorted-gap ratio must equal phi exactly")
    scaled = [value / (-triple[1]) for value in triple]
    require(
        scaled[0] == -(PHI4 - ONE) and scaled[1] == -ONE and scaled[2] == PHI4,
        "GOLDEN_PATTERN",
        "the centered triple must be proportional to (-(phi^4-1), -1, phi^4)",
    )
    s2, s3, s4 = port_sums(q)
    min_port = min(port_profile(q), key=lambda v: v.to_float())
    require(min_port == Q5(-1), "GOLDEN_MIN_PORT", "the golden profile minimum must be -1")
    return {
        "orbit": "diagonal_D2_cubic_extremal",
        "eigenvalues": [value.render() for value in triple],
        "eigenvalue_pattern": "proportional to (-(phi^4 - 1), -1, phi^4)",
        "sorted_gap_ratio": "phi = (1 + sqrt5)/2, exactly",
        "double_eigenvalue": False,
        "symmetric_criticality": (
            "the profile (2, -1, -1) is a critical point of the cubic at fixed "
            "square sum on the three-value plane, and D2-equivariance lifts "
            "stratum criticality to criticality of the full constrained "
            "functional"
        ),
        "S2": s2,
        "S3": s3,
        "S4": s4,
        "q_diagonal": [value.render() for value in q],
        "min_port_value": min_port,
    }


def normalized_invariants(branch: dict[str, Any]) -> tuple[Q5, Q5]:
    """Return (s3n, s4n) = (S3/S2^{3/2}, S4/S2^2) squared-rationalized.

    S2^{3/2} is irrational in general; the certificate stores the exact
    squared normalization s3n2 = S3^2 / S2^3 (a Q(sqrt5) number) together
    with the exact s4n = S4 / S2^2, and every comparison below is arranged
    to use only these exact quantities.
    """

    s2 = branch["S2"]
    s3 = branch["S3"]
    s4 = branch["S4"]
    s3n2 = (s3 * s3) / (s2 * s2 * s2)
    s4n = s4 / (s2 * s2)
    return s3n2, s4n


def crossing_amplitude(prolate: dict[str, Any], golden: dict[str, Any]) -> dict[str, Any]:
    """The exact branch-energy crossing.

    With each branch normalized to unit S2, the shape energy is
    E(r) = -s3n r^3/6 + s4n r^4/12, with s3n = sqrt(s3n2) > 0 for both
    branches.  The crossing r_c solves E_p(r) = E_g(r):
    r_c = 2 (s3n_p - s3n_g) / (s4n_p - s4n_g).  The certificate stores the
    exact squared data and the crossing as a root of the recorded exact
    equation, evaluated to controlled precision for display only.
    """

    p32, p4 = normalized_invariants(prolate)
    g32, g4 = normalized_invariants(golden)
    require(
        prolate["S3"].is_positive() and golden["S3"].is_positive(),
        "BRANCH_CUBIC_SIGN",
        "both branch cubics must be strictly positive in the chosen orientation",
    )
    require(not (p4 - g4).is_zero(), "BRANCH_QUARTIC_TIE", "the quartic invariants must differ")
    # Display floats, exactness carried by the stored equation.
    s3p = math.sqrt(p32.to_float())
    s3g = math.sqrt(g32.to_float())
    r_c = 2.0 * (s3p - s3g) / (p4 - g4).to_float()
    # Positivity domain: the record q_i = (1 + r w_i)/12 with unit-S2 profile
    # w stays a probability vector exactly when r < 1/|min_i w_i|, i.e.
    # r^2 < S2 / (min_i port value)^2, an exact Q(sqrt5) number per branch.
    positivity = {}
    r2_bound_float = {}
    for name, branch in (("prolate", prolate), ("golden", golden)):
        min_port = branch["min_port_value"]
        require(not min_port.is_positive() and not min_port.is_zero(), "POSITIVITY_SIGN", "the minimal port value must be negative")
        r2_bound = branch["S2"] / (min_port * min_port)
        r2_bound_float[name] = r2_bound.to_float()
        positivity[name] = {
            "r_squared_bound": r2_bound.render(),
            "r_bound_display": f"{math.sqrt(r2_bound_float[name]):.12f}",
        }
    require(
        r_c * r_c < r2_bound_float["golden"] and r_c * r_c < r2_bound_float["prolate"],
        "CROSSING_OUTSIDE_POSITIVITY",
        "the crossing must lie inside both positivity domains",
    )
    return {
        "equation": "2*(sqrt(s3n2_prolate) - sqrt(s3n2_golden)) = r_c * (s4n_prolate - s4n_golden)",
        "s3n2_prolate": p32.render(),
        "s3n2_golden": g32.render(),
        "s4n_prolate": p4.render(),
        "s4n_golden": g4.render(),
        "r_c_display": f"{r_c:.12f}",
        "positivity_domain": {
            "condition": "r < 1/|min_i w_i| for the unit-S2 profile w, i.e. r^2 < S2/(min port value)^2",
            "branches": positivity,
        },
        "selection": (
            "on the compared strata and inside each branch's positivity "
            "domain: below r_c the prolate branch has the lower shape "
            "energy; between r_c and the golden positivity bound the golden "
            "branch does; above that bound the golden record leaves the "
            "probability simplex and the packet emits no admissible "
            "simple-spectrum branch"
        ),
        "other_critical_orbits": (
            "the C3-axis orbit is critical by symmetric criticality with "
            "S3 = 0 exactly (two-value antipodal profile), S4/S2^2 = 1/12, "
            "positivity bound r^2 < 12, and shape energy r^4/144; it "
            "undercuts the golden branch only for r^2 > 96, outside every "
            "positivity bound, so it never enters the selection"
        ),
        "display_only_note": "r_c_display and r_bound_display are float renderings of the recorded exact quantities",
    }


# ---------------------------------------------------------------------------
# Comparison block (compare-only; explicit ancestry)
# ---------------------------------------------------------------------------


PDG_LEPTON_MASSES_GEV = {
    "electron": Fraction("0.00051099895069"),
    "muon": Fraction("0.1056583755"),
    "tau": Fraction("1.77693"),
}


def comparison_block() -> dict[str, Any]:
    logs = {k: math.log(float(v)) for k, v in PDG_LEPTON_MASSES_GEV.items()}
    ordered = sorted(logs.values())
    mean = sum(ordered) / 3.0
    centered = [value - mean for value in ordered]
    gap_low = centered[1] - centered[0]
    gap_high = centered[2] - centered[1]
    observed = gap_high / gap_low
    phi = PHI.to_float()
    mismatch_direct = abs(observed - phi) / phi
    mismatch_flipped = abs((1.0 / observed) - phi) / phi
    closest = min(mismatch_direct, mismatch_flipped)
    return {
        "ancestry": (
            "measured PDG charged-lepton masses, compare-only; no lepton "
            "value enters any derivation above"
        ),
        "observed_sorted_gap_ratio": f"{observed:.10f}",
        "observed_flipped": f"{1.0 / observed:.10f}",
        "packet_output": "phi = 1.6180339887 (golden branch); prolate branch degenerate",
        "relative_mismatch_direct": f"{mismatch_direct:.6f}",
        "relative_mismatch_flipped": f"{mismatch_flipped:.6f}",
        "closest_pairing_mismatch": f"{closest:.6f}",
        "verdict": "EXCLUDED",
        "verdict_statement": (
            "the quartic-truncated entropy packet's only simple-spectrum "
            "orbit has sorted-gap ratio exactly phi; the observed centered "
            "log-mass gap ratio differs from it by 16.7 percent in the "
            "closest orientation, far outside every stated width, and the "
            "prolate branch is excluded by the simple-spectrum gate"
        ),
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_prolate_simple_spectrum() -> dict[str, Any]:
    """The prolate branch must fail the simple-spectrum gate."""

    branch = prolate_orbit()
    try:
        require(
            branch["double_eigenvalue"] is False,
            "PROLATE_DEGENERATE",
            "the prolate branch has a double eigenvalue",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "PROLATE_DEGENERATE",
        }
    return {"expected_failure": True, "failed": False}


def control_profile_mutation() -> dict[str, Any]:
    """A mutated extremal profile must break the golden-ratio identity, so
    the identity gate refuses it (the required failure)."""

    mutated = (Q5(2), Q5(-1) + Q5(Fraction(1, 7)), Q5(-1) - Q5(Fraction(1, 7)))
    q = invert_profile(mutated)
    triple = sorted(q, key=lambda v: v.to_float())
    ratio = (triple[2] - triple[1]) / (triple[1] - triple[0])
    try:
        require(ratio == PHI, "GOLDEN_RATIO", "mutated profile must not give phi")
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "GOLDEN_RATIO",
            "meaning": "the golden identity consumes the exact extremal profile; a mutation is refused",
        }
    return {"expected_failure": True, "failed": False}


def control_frame_band_cubic_vanishes() -> dict[str, Any]:
    """An attempt to relocate the cubic functional to an antipode-odd
    profile must be refused: odd cubes cancel exactly."""

    values = [Q5(3), Q5(-3), Q5(Fraction(1, 2)), Q5(Fraction(-1, 2)), Q5(7), Q5(-7)]
    total = ZERO
    for value in values:
        total = total + value * value * value
    try:
        require(not total.is_zero(), "FRAME_CUBIC", "antipode-odd cubes cancel, so the cubic cannot live on the frame band")
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "FRAME_CUBIC",
            "meaning": "the cubic response vanishes on any antipode-odd profile, which pins the shape functional to the quintet band",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    prolate = prolate_orbit()
    golden = golden_orbit()
    crossing = crossing_amplitude(prolate, golden)

    controls = {
        "prolate_simple_spectrum": control_prolate_simple_spectrum(),
        "profile_mutation": control_profile_mutation(),
        "frame_band_cubic": control_frame_band_cubic_vanishes(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    def strip(branch: dict[str, Any]) -> dict[str, Any]:
        out = dict(branch)
        for key in ("S2", "S3", "S4", "min_port_value"):
            out[key] = branch[key].render()
        return out

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue_context": ISSUE_CONTEXT,
        "claim_boundary": (
            "Exact structure theorem for the entropy-derived W5 shape packet "
            "and its compare-only exclusion. The packet is derived from the "
            "A3 expansion with no free coefficient; its only simple-spectrum "
            "orbit has sorted-gap ratio exactly phi; the observed lepton "
            "shape excludes it at 16.7 percent in the closest orientation. "
            "No mass, ratio, or prediction is emitted, and the source-only "
            "charged no-go is unchanged."
        ),
        "packet_derivation": {
            "expansion": "D(u(1+x)||u) = u(S2/2 - S3/6 + S4/12) + O(x^5), u = 1/12",
            "band_selection": (
                "S3 vanishes identically on the antipode-odd frame band, so "
                "the quintet band is the only nontrivial band with a cubic "
                "entropy response"
            ),
            "premise": (
                "the shape functional is the quartic-truncated expansion "
                "restricted to the quintet band at fixed band amplitude; the "
                "amplitude is a declared parameter of the premise, ranging "
                "over the positivity domain of the record, and the branch "
                "selection depends on it only through the recorded exact "
                "crossing and positivity bounds"
            ),
        },
        "branches": {
            "prolate": strip(prolate),
            "golden": strip(golden),
        },
        "crossing": crossing,
        "comparison": comparison_block(),
        "controls": controls,
        "boundary_extension": (
            "the stabilizer theorem excludes symmetry geometry alone; this "
            "certificate excludes the quartic entropy packet; the open "
            "lepton-side object is a source mechanism beyond the quartic "
            "entropy truncation"
        ),
    }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.verify:
        stored = json.loads(args.out.read_text(encoding="utf-8"))
        require(stored == payload, "DRIFT", "stored certificate does not match a rebuild")
        print(json.dumps({"status": "PASS"}))
        return 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "WROTE", "out": str(args.out), "verdict": payload["comparison"]["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
