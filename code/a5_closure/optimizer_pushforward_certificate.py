#!/usr/bin/env python3
"""Exact certificate for GitHub issue #619: optimizer pushforward across refinement.

The question is whether independently optimized A3 states push forward across
refinement.  The setting is classical and exact: the fine level is the
probability simplex on the twelve carrier ports, the coarse level is the
probability simplex on the six incidence axes, and the refinement channel T
sums each antipodal port pair into its axis.  The pairing is A1-sourced: it is
read from the hash-pinned #565 selector receipt, whose inverse pairing is the
unique fixed-point-free graph-distance-three involution of the carrier, and
the six axes are the derived incidence axes.  References are the uniform fine
state (1/12 per port) and the uniform coarse state (1/6 per axis); the coarse
reference equals the pushforward of the fine reference, verified exactly.  The
objective at each level is relative entropy to the level's reference, and the
optimizers are information projections onto feasible sets.

Every projection in this file is certified by an exact stationarity check that
needs no logarithms.  For a feasible set that fixes the total mass of each
block of a partition, the chain rule splits the relative entropy into the
block-marginal term plus a nonnegative mixture of block-conditional
divergences.  The marginal term is constant on the feasible set, and the
conditional term vanishes exactly when every block-conditional distribution
equals the block-conditional reference (Gibbs' inequality, with uniqueness
from strict convexity).  Matching block conditionals is a rational identity,
so the Karush-Kuhn-Tucker stationarity of each claimed projection is verified
by exact fraction arithmetic.  For a feasible set that contains the reference
itself, the projection is the reference, because the divergence is zero there
and positive elsewhere.

The certificate derives, rather than assumes:

* the positive sufficiency family: for every rational coarse target m in the
  open six-simplex, the fine feasible set with axis marginals m has the
  within-pair equal split as its information projection, the coarse feasible
  set is the singleton {m}, and the two optimizers commute with T exactly;
  a battery of rational targets, including asymmetric ones, witnesses the
  one-parameter-vector family;
* the quantitative nonclosure countermodel: pinning a single port to 1/6 is a
  within-axis constraint that T cannot express; fine projection then
  pushforward gives axis weight 8/33, pushforward then projection gives 1/6,
  and the nonclosure defect is exactly 5/66;
* the exact image of the pinned-port feasible set under T, an axis-mass
  inequality set, with a preimage construction battery and an exact
  no-preimage witness below the threshold;
* the control battery: wrong-reference, changed-weight, incomplete-constraint,
  and non-sufficient-channel controls each break commutation with an exact
  rational witness.  On the full-marginal family the coarse feasible set is a
  single point and its projection ignores the reference, so the two
  reference-mismatch controls exhibit their breakage on the smallest
  channel-expressible relaxation, a single-axis mass constraint, and the
  masking is recorded;
* the verdict: cross-refinement optimizer compatibility is a conditional open
  interface.  It is exact under the named sufficiency conditions (compatible
  references and channel-expressible constraints) and fails by the exact
  defect 5/66 otherwise, so it is a premise to be discharged per constraint
  family, and it is an optimizer-interface statement rather than an instance
  of A2 naturality.

All arithmetic is exact rational arithmetic on classical probability
simplices; no floating point number appears anywhere in this certificate.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

SCHEMA = "oph.optimizer_pushforward_certificate.v1"
RECEIPT_SCHEMA = "oph.optimizer_pushforward_receipt.v1"
NEGATIVE_SCHEMA = "oph.optimizer_pushforward_negative_controls.v1"
AXIS_RECEIPT_SCHEMA = "oph.echosahedral_selector_receipt.v1"

FINE_PORTS = 12
COARSE_AXES = 6
ARITHMETIC_CONVENTION = "exact_rational_simplex"

# Manifest keys that would smuggle in approximation, declared answers, or a
# quantum extension.  The certificate is classical and exact by construction.
FORBIDDEN_MANIFEST_KEYS = (
    "float_tolerance",
    "numeric_epsilon",
    "approximate_projection",
    "quantum_state_space",
    "density_matrix",
    "nonclosure_defect",
    "measured_coupling",
)

Vector = tuple[Fraction, ...]


def frac(text: str | int) -> Fraction:
    return Fraction(text)


def vec(entries: Sequence[str | int]) -> Vector:
    return tuple(Fraction(entry) for entry in entries)


def vec_str(vector: Sequence[Fraction]) -> list[str]:
    return [str(entry) for entry in vector]


def require_simplex(
    vector: Sequence[Fraction], size: int, open_interior: bool, code: str, label: str
) -> None:
    require(len(vector) == size, code, f"{label}: expected {size} coordinates")
    require(sum(vector) == 1, code, f"{label}: coordinates do not sum to one")
    for entry in vector:
        if open_interior:
            require(entry > 0, code, f"{label}: coordinate {entry} is not interior")
        else:
            require(entry >= 0, code, f"{label}: coordinate {entry} is negative")


def uniform(size: int) -> Vector:
    return tuple(Fraction(1, size) for _ in range(size))


def pushforward(rho: Sequence[Fraction], axes: Sequence[tuple[int, int]]) -> Vector:
    return tuple(rho[i] + rho[j] for i, j in axes)


def block_mass_projection(
    reference: Sequence[Fraction],
    blocks: Sequence[Sequence[int]],
    masses: Sequence[Fraction],
) -> Vector:
    """Information projection onto {rho : sum over block b equals masses[b]}.

    The blocks partition the coordinate set.  By the relative-entropy chain
    rule over the block partition, D(rho || ref) equals the divergence of the
    block-mass vectors plus the mass-weighted sum of block-conditional
    divergences.  The first term is constant on the feasible set and the
    second is nonnegative, vanishing exactly when every block-conditional
    distribution equals the block-conditional reference.  The unique minimizer
    therefore distributes each prescribed block mass over its block in
    proportion to the reference, which is the exact formula used here; its
    stationarity is re-verified by ``conditional_matches_reference``.
    """

    covered = sorted(index for block in blocks for index in block)
    require(
        covered == list(range(len(reference))),
        "PROJECTION",
        "the blocks do not partition the coordinate set",
    )
    require(
        sum(masses) == 1 and all(mass >= 0 for mass in masses),
        "PROJECTION",
        "the prescribed block masses are not a probability vector",
    )
    result = [Fraction(0)] * len(reference)
    for block, mass in zip(blocks, masses):
        block_reference_mass = sum(reference[index] for index in block)
        require(
            block_reference_mass > 0,
            "PROJECTION",
            "a block carries zero reference mass",
        )
        for index in block:
            result[index] = mass * reference[index] / block_reference_mass
    return tuple(result)


def conditional_matches_reference(
    candidate: Sequence[Fraction],
    reference: Sequence[Fraction],
    block: Sequence[int],
) -> bool:
    """Exact stationarity check: block conditionals of candidate and reference agree.

    The identity candidate[i] * ref_mass == reference[i] * candidate_mass for
    every i in the block states that the two block-conditional distributions
    coincide, which is the exact first-order condition for the chain-rule
    argument in ``block_mass_projection``.  It is a rational identity and
    involves no logarithm.
    """

    candidate_mass = sum(candidate[index] for index in block)
    reference_mass = sum(reference[index] for index in block)
    return all(
        candidate[index] * reference_mass == reference[index] * candidate_mass
        for index in block
    )


def load_axis_pairing(
    manifest: Mapping[str, Any], base_dir: Path
) -> tuple[tuple[tuple[int, int], ...], str]:
    """Load the hash-pinned #565 receipt and derive the six antipodal axes.

    The channel is A1-sourced: the pairing is the receipt's inverse pairing,
    the unique fixed-point-free graph-distance-three involution of the
    certified carrier, and the axes are its orbits.
    """

    path_raw = manifest.get("axis_receipt_path")
    require(isinstance(path_raw, str), "UPSTREAM_REFERENCE", "axis_receipt_path is missing")
    path = Path(path_raw)
    if not path.is_absolute():
        path = base_dir / path
    receipt = load_json(path)
    require(
        manifest.get("axis_receipt_sha256") == sha256_json(receipt),
        "UPSTREAM_HASH",
        "the #565 axis receipt hash does not match the declared pin",
    )
    require(
        receipt.get("issue") == 565 and receipt.get("schema") == AXIS_RECEIPT_SCHEMA,
        "UPSTREAM_RECEIPT",
        "the pinned receipt is not a #565 selector receipt",
    )
    pairing = receipt.get("inverse_pairing", {})
    require(
        pairing.get("fixed_point_free") is True
        and pairing.get("involutive") is True
        and pairing.get("unique") is True,
        "UPSTREAM_RECEIPT",
        "the pinned receipt does not certify a unique fixed-point-free involution",
    )
    raw_pairs = pairing.get("pairs")
    require(
        isinstance(raw_pairs, list) and len(raw_pairs) == COARSE_AXES,
        "AXIS_PAIRING",
        "the receipt does not carry six antipodal pairs",
    )
    pairs: list[tuple[int, int]] = []
    for row in raw_pairs:
        require(
            isinstance(row, list) and len(row) == 2,
            "AXIS_PAIRING",
            "a pairing row is malformed",
        )
        indices = []
        for name in row:
            require(
                isinstance(name, str) and name.startswith("p") and name[1:].isdigit(),
                "AXIS_PAIRING",
                f"port name {name!r} is not in the pNN convention",
            )
            indices.append(int(name[1:]))
        low, high = sorted(indices)
        require(low != high, "AXIS_PAIRING", "a pair repeats a port")
        pairs.append((low, high))
    pairs.sort()
    covered = sorted(index for pair in pairs for index in pair)
    require(
        covered == list(range(FINE_PORTS)),
        "AXIS_PAIRING",
        "the six pairs do not partition the twelve ports",
    )
    return tuple(pairs), sha256_json(receipt)


def validate_manifest(manifest: Mapping[str, Any]) -> None:
    require(manifest.get("schema") == SCHEMA, "SCHEMA", f"expected {SCHEMA}")
    require(
        manifest.get("arithmetic_convention") == ARITHMETIC_CONVENTION,
        "CONVENTION",
        "the manifest must declare exact rational arithmetic on classical simplices",
    )
    for key in FORBIDDEN_MANIFEST_KEYS:
        require(key not in manifest, "FORBIDDEN_DEPENDENCY", f"forbidden manifest key {key}")


# Rational coarse targets for the sufficiency family, including asymmetric
# ones; each lies in the open six-simplex.
POSITIVE_BATTERY: tuple[tuple[str, ...], ...] = (
    ("1/6", "1/6", "1/6", "1/6", "1/6", "1/6"),
    ("1/21", "2/21", "3/21", "4/21", "5/21", "6/21"),
    ("1/2", "1/4", "1/8", "1/16", "1/32", "1/32"),
    ("1/3", "1/4", "1/6", "1/8", "1/12", "1/24"),
    ("2/7", "1/7", "1/7", "1/7", "1/7", "1/7"),
    ("5/12", "1/12", "1/6", "1/12", "1/6", "1/12"),
)


def positive_sufficiency_family(
    axes: Sequence[tuple[int, int]],
    fine_reference: Vector,
    coarse_reference: Vector,
) -> dict[str, Any]:
    """The exact one-parameter-vector family where optimizers commute.

    For a rational target m in the open six-simplex the fine feasible set
    K_m = {rho : axis marginals equal m} fixes the mass of every antipodal
    pair, so the chain-rule argument applies with the pairs as blocks: the
    unique information projection splits each axis mass within its pair in
    proportion to the fine reference, which is the equal split for the
    uniform reference.  The coarse feasible set T(K_m) is the singleton {m}
    and its projection is m.  Commutation T(rho*) = m is then an exact
    rational identity, checked per battery member.
    """

    battery_rows = []
    for target_raw in POSITIVE_BATTERY:
        target = vec(target_raw)
        require_simplex(target, COARSE_AXES, True, "POSITIVE_FAMILY", "coarse target")
        fine_projection = block_mass_projection(fine_reference, axes, target)
        require_simplex(
            fine_projection, FINE_PORTS, True, "POSITIVE_FAMILY", "fine projection"
        )
        for pair in axes:
            require(
                conditional_matches_reference(fine_projection, fine_reference, pair),
                "POSITIVE_FAMILY",
                "the within-pair conditional does not match the reference",
            )
            i, j = pair
            require(
                fine_projection[i] == fine_projection[j],
                "POSITIVE_FAMILY",
                "the uniform-reference projection is not the equal split",
            )
        pushed = pushforward(fine_projection, axes)
        require(
            pushed == target,
            "POSITIVE_FAMILY",
            "the pushforward of the fine optimizer misses the coarse target",
        )
        # The coarse feasible set is the singleton {m}: the axis-marginal
        # constraints determine every coarse coordinate, so the coarse
        # projection is m for any coarse reference.
        coarse_projection = target
        require(
            pushed == coarse_projection,
            "POSITIVE_FAMILY",
            "fine-then-push and push-then-optimize disagree on the family",
        )
        battery_rows.append(
            {
                "coarse_target": vec_str(target),
                "fine_projection": vec_str(fine_projection),
                "pushforward_of_fine_projection": vec_str(pushed),
                "coarse_projection": vec_str(coarse_projection),
                "within_pair_split": "equal",
                "commutes": True,
            }
        )
    asymmetric = sum(
        1 for row in battery_rows if len(set(row["coarse_target"])) == COARSE_AXES
    )
    require(
        len(battery_rows) >= 5 and asymmetric >= 2,
        "POSITIVE_FAMILY",
        "the battery must carry at least five targets including asymmetric ones",
    )
    require(
        pushforward(fine_reference, axes) == coarse_reference,
        "POSITIVE_FAMILY",
        "the coarse reference is not the pushforward of the fine reference",
    )
    return {
        "family": (
            "for every rational m in the open six-simplex, the fine feasible "
            "set fixes the axis marginals at m and the coarse feasible set is "
            "the singleton {m}"
        ),
        "free_rational_parameters": COARSE_AXES - 1,
        "kkt_argument": (
            "chain rule over the antipodal pairs: the axis-marginal divergence "
            "is constant on the feasible set and the pair-conditional term is "
            "nonnegative, vanishing exactly at the reference-proportional "
            "within-pair split; uniqueness follows from strict convexity, and "
            "stationarity is verified as the rational identity that pair "
            "conditionals of projection and reference coincide"
        ),
        "battery": battery_rows,
        "battery_size": len(battery_rows),
        "fully_asymmetric_targets": asymmetric,
        "conclusion": (
            "when the constraints are channel-expressible and the coarse "
            "reference is the pushforward of the fine reference, the "
            "information projections commute with refinement, exactly"
        ),
    }


def nonclosure_countermodel(
    axes: Sequence[tuple[int, int]],
    fine_reference: Vector,
    coarse_reference: Vector,
) -> dict[str, Any]:
    """The exact quantitative countermodel: a within-axis constraint.

    The fine feasible set pins port 0 to 1/6, a constraint T cannot express
    because it separates states with equal pushforward.  The fine projection
    pins port 0 and, by the chain rule over the partition {port 0, free face},
    distributes the remaining 5/6 uniformly over the eleven free ports.  Its
    pushforward puts 8/33 on the axis of port 0.  The pushed feasible set is
    the axis-mass inequality set {q : q_0 >= 1/6}, which contains the uniform
    coarse reference, so the coarse projection is the reference itself.  The
    nonclosure defect is 8/33 - 1/6 = 5/66, exactly.
    """

    pinned_port = 0
    pinned_value = Fraction(1, 6)
    axis_of_pinned = next(
        index for index, pair in enumerate(axes) if pinned_port in pair
    )
    partner_port = next(
        port for port in axes[axis_of_pinned] if port != pinned_port
    )
    free_ports = [index for index in range(FINE_PORTS) if index != pinned_port]

    # Fine projection: chain rule over the partition {pinned port, free face}.
    fine_projection = block_mass_projection(
        fine_reference,
        ([pinned_port], free_ports),
        (pinned_value, 1 - pinned_value),
    )
    require_simplex(
        fine_projection, FINE_PORTS, True, "COUNTERMODEL", "fine projection"
    )
    require(
        fine_projection[pinned_port] == pinned_value,
        "COUNTERMODEL",
        "the pinned coordinate drifted",
    )
    require(
        all(fine_projection[index] == Fraction(5, 66) for index in free_ports),
        "COUNTERMODEL",
        "the free-face projection is not uniform at 5/66",
    )
    require(
        conditional_matches_reference(fine_projection, fine_reference, free_ports),
        "COUNTERMODEL",
        "the free-face conditional does not match the reference",
    )

    pushed = pushforward(fine_projection, axes)
    require(
        pushed[axis_of_pinned] == Fraction(8, 33)
        and pushed[axis_of_pinned] == pinned_value + Fraction(5, 66),
        "COUNTERMODEL",
        "the pushed axis weight is not 1/6 + 5/66 = 8/33",
    )
    require(
        all(
            pushed[index] == Fraction(5, 33)
            for index in range(COARSE_AXES)
            if index != axis_of_pinned
        ),
        "COUNTERMODEL",
        "the pushed free axes are not uniform at 5/33",
    )

    # Exact image of the feasible set under T.  Any feasible rho has axis
    # weight 1/6 + rho[partner] with rho[partner] in [0, 5/6], so the axis
    # weight ranges over [1/6, 1]; the other axes are unconstrained.  The
    # construction below realizes any such coarse state exactly.
    image_battery = []
    for q_raw in (
        ("1/6", "1/6", "1/6", "1/6", "1/6", "1/6"),
        ("8/33", "5/33", "5/33", "5/33", "5/33", "5/33"),
        ("1/3", "2/15", "2/15", "2/15", "2/15", "2/15"),
        ("1/2", "1/10", "1/10", "1/10", "1/10", "1/10"),
        ("5/6", "1/30", "1/30", "1/30", "1/30", "1/30"),
        ("1", "0", "0", "0", "0", "0"),
    ):
        q = vec(q_raw)
        require_simplex(q, COARSE_AXES, False, "COUNTERMODEL", "image point")
        require(
            q[axis_of_pinned] >= pinned_value,
            "COUNTERMODEL",
            "an image battery point violates the axis threshold",
        )
        preimage = [Fraction(0)] * FINE_PORTS
        preimage[pinned_port] = pinned_value
        preimage[partner_port] = q[axis_of_pinned] - pinned_value
        for index, pair in enumerate(axes):
            if index == axis_of_pinned:
                continue
            i, j = pair
            preimage[i] = q[index] / 2
            preimage[j] = q[index] / 2
        require_simplex(preimage, FINE_PORTS, False, "COUNTERMODEL", "preimage")
        require(
            preimage[pinned_port] == pinned_value
            and pushforward(preimage, axes) == q,
            "COUNTERMODEL",
            "the preimage construction fails",
        )
        image_battery.append(
            {"coarse_state": vec_str(q), "preimage": vec_str(preimage)}
        )
    below_threshold = vec(("1/12", "11/60", "11/60", "11/60", "11/60", "11/60"))
    require_simplex(below_threshold, COARSE_AXES, False, "COUNTERMODEL", "witness")
    require(
        below_threshold[axis_of_pinned] - pinned_value < 0,
        "COUNTERMODEL",
        "the no-preimage witness is not below the threshold",
    )

    # Coarse projection: the inequality set contains the coarse reference
    # (1/6 >= 1/6 exactly), the divergence vanishes there and is positive
    # elsewhere, so the projection is the reference itself.
    require(
        coarse_reference[axis_of_pinned] >= pinned_value,
        "COUNTERMODEL",
        "the coarse reference is not feasible for the pushed constraint",
    )
    coarse_projection = coarse_reference

    defect = pushed[axis_of_pinned] - coarse_projection[axis_of_pinned]
    require(
        defect == Fraction(5, 66)
        and defect == Fraction(16, 66) - Fraction(11, 66),
        "COUNTERMODEL",
        f"the nonclosure defect drifted: {defect}",
    )

    # Exact non-expressibility witness: two fine states with the same
    # pushforward, exactly one of which is feasible.
    swapped = list(fine_projection)
    swapped[pinned_port], swapped[partner_port] = (
        swapped[partner_port],
        swapped[pinned_port],
    )
    require(
        pushforward(swapped, axes) == pushed
        and swapped[pinned_port] != pinned_value,
        "COUNTERMODEL",
        "the non-expressibility witness fails",
    )

    return {
        "fine_constraint": "rho[port 0] = 1/6, a within-axis constraint",
        "pinned_port": pinned_port,
        "partner_port": partner_port,
        "axis_of_pinned_port": axis_of_pinned,
        "fine_projection": vec_str(fine_projection),
        "free_face_value": "5/66",
        "pushforward_of_fine_projection": vec_str(pushed),
        "pushed_axis_weight": "8/33",
        "pushed_feasible_set": (
            "{q in the six-simplex : q[axis of port 0] >= 1/6}; the pinned "
            "port contributes 1/6 and the free partner port sweeps [0, 5/6]"
        ),
        "image_surjectivity_battery": image_battery,
        "no_preimage_witness": {
            "coarse_state": vec_str(below_threshold),
            "reason": (
                "a preimage would need the partner port at 1/12 - 1/6 = -1/12, "
                "which is negative"
            ),
        },
        "coarse_projection": vec_str(coarse_projection),
        "coarse_projection_argument": (
            "the reference lies in the pushed feasible set (1/6 >= 1/6), the "
            "divergence to the reference is zero there and positive elsewhere, "
            "so the projection is the reference"
        ),
        "nonclosure_defect": str(defect),
        "defect_arithmetic": "8/33 - 1/6 = 16/66 - 11/66 = 5/66",
        "non_expressibility_witness": {
            "feasible_state": vec_str(fine_projection),
            "infeasible_state_with_equal_pushforward": vec_str(tuple(swapped)),
            "conclusion": (
                "the constraint separates states with equal pushforward, so no "
                "constraint on T(rho) expresses it"
            ),
        },
        "conclusion": (
            "fine-then-push and push-then-optimize disagree by exactly 5/66 on "
            "the constrained axis; optimizer compatibility fails quantitatively "
            "for constraints the channel cannot express"
        ),
    }


def control_battery(
    axes: Sequence[tuple[int, int]],
    fine_reference: Vector,
    coarse_reference: Vector,
) -> dict[str, Any]:
    """Four controls, each breaking commutation with an exact witness.

    On the full-marginal family the coarse feasible set is a single point and
    its projection ignores the coarse reference, so a reference tilt is masked
    there.  The two reference-mismatch controls therefore run on the smallest
    channel-expressible relaxation, the single constraint that the axis of
    ports 0 and 3 carries mass 1/3, where the coarse feasible set has free
    directions and the reference matters.
    """

    arena_axis = 0
    arena_mass = Fraction(1, 3)
    arena_pair = axes[arena_axis]
    free_axes = [index for index in range(COARSE_AXES) if index != arena_axis]
    free_ports = [
        index for index in range(FINE_PORTS) if index not in arena_pair
    ]
    arena_blocks = (list(arena_pair), free_ports)
    arena_masses = (arena_mass, 1 - arena_mass)
    coarse_blocks = ([arena_axis], free_axes)

    # --- wrong-reference: tilted coarse reference, uniform fine reference ----
    tilted_coarse = vec(("1/6", "1/3", "1/12", "1/12", "1/6", "1/6"))
    require_simplex(tilted_coarse, COARSE_AXES, True, "CONTROL", "tilted coarse")
    require(
        tilted_coarse != pushforward(fine_reference, axes),
        "CONTROL",
        "the tilted coarse reference must differ from the pushforward",
    )
    fine_opt = block_mass_projection(fine_reference, arena_blocks, arena_masses)
    fine_then_push = pushforward(fine_opt, axes)
    require(
        fine_then_push
        == vec(("1/3", "2/15", "2/15", "2/15", "2/15", "2/15")),
        "CONTROL",
        "the wrong-reference fine-then-push lane drifted",
    )
    coarse_opt_tilted = block_mass_projection(
        tilted_coarse, coarse_blocks, (arena_mass, 1 - arena_mass)
    )
    require(
        coarse_opt_tilted
        == vec(("1/3", "4/15", "1/15", "1/15", "2/15", "2/15")),
        "CONTROL",
        "the wrong-reference coarse projection drifted",
    )
    require(
        conditional_matches_reference(coarse_opt_tilted, tilted_coarse, free_axes),
        "CONTROL",
        "the tilted coarse projection fails its stationarity check",
    )
    wrong_reference_defect = coarse_opt_tilted[1] - fine_then_push[1]
    require(
        fine_then_push != coarse_opt_tilted
        and wrong_reference_defect == Fraction(2, 15),
        "CONTROL",
        "the wrong-reference control did not fail closed",
    )
    wrong_reference = {
        "coarse_reference": vec_str(tilted_coarse),
        "fine_reference": vec_str(fine_reference),
        "arena": "axis 0 mass pinned to 1/3, a channel-expressible constraint",
        "fine_then_push": vec_str(fine_then_push),
        "push_then_optimize": vec_str(coarse_opt_tilted),
        "defect_on_axis_1": str(wrong_reference_defect),
        "masking_note": (
            "on the full-marginal family the coarse feasible set is the "
            "singleton {m}, its projection is m for any coarse reference, and "
            "the tilt is invisible; the breakage appears exactly on the "
            "smallest channel-expressible relaxation"
        ),
        "fails_closed": True,
    }

    # --- changed-weight ------------------------------------------------------
    # Part one: a within-pair tilt of the fine reference (2/18 against 1/18 on
    # ports 0 and 3) keeps every axis mass at 1/6, so the coarse reference
    # stays the exact pushforward and commutation survives; what breaks is the
    # equal-split formula, since the split follows the reference ratio 2:1.
    within_pair_tilt = [Fraction(1, 12)] * FINE_PORTS
    within_pair_tilt[arena_pair[0]] = Fraction(2, 18)
    within_pair_tilt[arena_pair[1]] = Fraction(1, 18)
    within_pair_tilt = tuple(within_pair_tilt)
    require_simplex(
        within_pair_tilt, FINE_PORTS, True, "CONTROL", "within-pair tilt"
    )
    require(
        pushforward(within_pair_tilt, axes) == coarse_reference,
        "CONTROL",
        "the within-pair tilt must leave the pushforward uniform",
    )
    target = vec(("1/3", "1/4", "1/6", "1/8", "1/12", "1/24"))
    tilted_projection = block_mass_projection(within_pair_tilt, axes, target)
    require(
        tilted_projection[arena_pair[0]] == Fraction(2, 9)
        and tilted_projection[arena_pair[1]] == Fraction(1, 9)
        and tilted_projection[arena_pair[0]] != target[arena_axis] / 2,
        "CONTROL",
        "the changed-weight split does not follow the reference ratio 2:1",
    )
    require(
        pushforward(tilted_projection, axes) == target,
        "CONTROL",
        "the T-invisible tilt must leave commutation intact",
    )
    # Part two: the mismatched pair.  A cross-axis tilt of the fine reference
    # with the coarse reference left uniform makes the coarse reference differ
    # from the pushforward, and commutation fails on the arena constraint.
    cross_axis_tilt_axis_masses = ("1/6", "1/3", "1/12", "1/12", "1/6", "1/6")
    cross_axis_tilt = [Fraction(0)] * FINE_PORTS
    for index, pair in enumerate(axes):
        for port in pair:
            cross_axis_tilt[port] = Fraction(cross_axis_tilt_axis_masses[index]) / 2
    cross_axis_tilt = tuple(cross_axis_tilt)
    require_simplex(cross_axis_tilt, FINE_PORTS, True, "CONTROL", "cross-axis tilt")
    require(
        pushforward(cross_axis_tilt, axes) != coarse_reference,
        "CONTROL",
        "the mismatched pair must break reference compatibility",
    )
    mismatched_fine_opt = block_mass_projection(
        cross_axis_tilt, arena_blocks, arena_masses
    )
    mismatched_push = pushforward(mismatched_fine_opt, axes)
    require(
        mismatched_push
        == vec(("1/3", "4/15", "1/15", "1/15", "2/15", "2/15")),
        "CONTROL",
        "the mismatched-pair fine-then-push lane drifted",
    )
    coarse_opt_uniform = block_mass_projection(
        coarse_reference, coarse_blocks, (arena_mass, 1 - arena_mass)
    )
    require(
        coarse_opt_uniform
        == vec(("1/3", "2/15", "2/15", "2/15", "2/15", "2/15")),
        "CONTROL",
        "the uniform coarse projection drifted",
    )
    changed_weight_defect = mismatched_push[1] - coarse_opt_uniform[1]
    require(
        mismatched_push != coarse_opt_uniform
        and changed_weight_defect == Fraction(2, 15),
        "CONTROL",
        "the changed-weight control did not fail closed",
    )
    changed_weight = {
        "t_invisible_tilt": {
            "fine_reference": vec_str(within_pair_tilt),
            "coarse_target": vec_str(target),
            "fine_projection": vec_str(tilted_projection),
            "within_pair_split_ratio": "2:1, the reference ratio",
            "equal_split_formula_breaks": True,
            "pushforward_of_fine_reference_stays_uniform": True,
            "commutation_survives": True,
            "finding": (
                "an unequal within-pair reference changes the optimizer while "
                "leaving its pushforward fixed, because the tilt is invisible "
                "to T; commutation is unharmed when both references stay "
                "compatible"
            ),
        },
        "mismatched_pair": {
            "fine_reference": vec_str(cross_axis_tilt),
            "coarse_reference": vec_str(coarse_reference),
            "coarse_reference_is_pushforward": False,
            "arena": "axis 0 mass pinned to 1/3",
            "fine_then_push": vec_str(mismatched_push),
            "push_then_optimize": vec_str(coarse_opt_uniform),
            "defect_on_axis_1": str(changed_weight_defect),
            "finding": (
                "a cross-axis fine tilt with the coarse reference left uniform "
                "destroys reference compatibility and commutation fails by "
                "exactly 2/15 on axis 1"
            ),
        },
        "fails_closed": True,
    }

    # --- incomplete-constraint ----------------------------------------------
    dropped_target = vec(("1/21", "2/21", "3/21", "4/21", "5/21", "6/21"))
    dropped_fine = block_mass_projection(fine_reference, axes, dropped_target)
    dropped_push = pushforward(dropped_fine, axes)
    require(
        dropped_push == dropped_target,
        "CONTROL",
        "the incomplete-constraint fine lane drifted",
    )
    # The coarse problem drops the marginal constraint: the feasible set is
    # the whole simplex, which contains the reference, so the projection is
    # the reference.
    dropped_coarse = coarse_reference
    incomplete_defect = dropped_push[COARSE_AXES - 1] - dropped_coarse[COARSE_AXES - 1]
    require(
        dropped_push != dropped_coarse and incomplete_defect == Fraction(5, 42),
        "CONTROL",
        "the incomplete-constraint control did not fail closed",
    )
    incomplete_constraint = {
        "fine_constraint": "axis marginals fixed at (1/21, ..., 6/21)",
        "coarse_constraint": "dropped; the feasible set is the whole simplex",
        "fine_then_push": vec_str(dropped_push),
        "push_then_optimize": vec_str(dropped_coarse),
        "defect_on_axis_5": str(incomplete_defect),
        "defect_arithmetic": "6/21 - 1/6 = 12/42 - 7/42 = 5/42",
        "fails_closed": True,
    }

    non_sufficient_channel = {
        "cross_reference": "nonclosure_countermodel",
        "statement": (
            "the pinned-port constraint is the non-sufficient-channel control "
            "itself: it separates states with equal pushforward, and the "
            "recorded defect is 5/66"
        ),
        "defect": "5/66",
        "fails_closed": True,
    }

    return {
        "wrong_reference": wrong_reference,
        "changed_weight": changed_weight,
        "incomplete_constraint": incomplete_constraint,
        "non_sufficient_channel": non_sufficient_channel,
    }


def certificate_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    base = base_dir or MODULE_DIR
    validate_manifest(manifest)
    axes, axis_receipt_sha = load_axis_pairing(manifest, base)

    fine_reference = uniform(FINE_PORTS)
    coarse_reference = uniform(COARSE_AXES)
    require_simplex(fine_reference, FINE_PORTS, True, "REFERENCE", "fine reference")
    require_simplex(coarse_reference, COARSE_AXES, True, "REFERENCE", "coarse reference")
    require(
        pushforward(fine_reference, axes) == coarse_reference,
        "REFERENCE",
        "the coarse reference must equal the pushforward of the fine reference",
    )

    positive = positive_sufficiency_family(axes, fine_reference, coarse_reference)
    countermodel = nonclosure_countermodel(axes, fine_reference, coarse_reference)
    controls = control_battery(axes, fine_reference, coarse_reference)

    return {
        "schema": RECEIPT_SCHEMA,
        "issue": 619,
        "manifest_sha256": sha256_json(manifest),
        "axis_receipt_sha256": axis_receipt_sha,
        "arithmetic": {
            "convention": ARITHMETIC_CONVENTION,
            "statement": (
                "every quantity is an exact rational number on a classical "
                "probability simplex; stationarity is checked through the "
                "chain-rule conditional identity, which needs no logarithm"
            ),
        },
        "refinement_channel": {
            "fine_level": "the twelve-port probability simplex",
            "coarse_level": "the six-axis probability simplex",
            "channel": "T(rho)[axis] = rho[i] + rho[j] over the antipodal pair (i, j)",
            "axis_pairs": [list(pair) for pair in axes],
            "pairing_source": (
                "the hash-pinned #565 selector receipt: the pairs are the "
                "orbits of the unique fixed-point-free graph-distance-three "
                "involution of the certified carrier, so the six axes are the "
                "derived incidence axes and the channel is A1-sourced"
            ),
        },
        "references": {
            "fine": vec_str(fine_reference),
            "coarse": vec_str(coarse_reference),
            "coarse_equals_pushforward_of_fine": True,
            "objective": "relative entropy to the reference at each level",
        },
        "positive_sufficiency_family": positive,
        "nonclosure_countermodel": countermodel,
        "control_battery": controls,
        "verdict": {
            "positive_lane": "exact_sufficiency_family",
            "negative_lane": "exact_nonclosure_defect_5_66",
            "axiom_status": "conditional_open_interface",
            "bounded_exit": "conditional_open_interface",
            "sufficiency_conditions": [
                "the coarse reference equals the pushforward of the fine reference",
                (
                    "the fine feasible set is channel-expressible: it is the "
                    "full T-preimage of its own image"
                ),
            ],
            "statement": (
                "cross-refinement optimizer compatibility is a premise, "
                "discharged per constraint family: A1, A2, and A3 do not force "
                "it, it holds exactly under the named sufficiency conditions, "
                "and it fails by the exact defect 5/66 for the pinned-port "
                "countermodel"
            ),
            "controls": {
                "wrong_reference": "fails closed with defect 2/15 on axis 1",
                "changed_weight": (
                    "the T-invisible tilt leaves commutation intact; the "
                    "mismatched pair fails closed with defect 2/15 on axis 1"
                ),
                "incomplete_constraint": "fails closed with defect 5/42 on axis 5",
                "non_sufficient_channel": "fails closed with defect 5/66 on axis 0",
            },
            "gravity_ladder_composition": (
                "the result is an optimizer-interface statement about A3 "
                "projections across the A1-sourced channel; it is consumed by "
                "the #618 ladder as a named premise and is not an instance of "
                "A2 naturality of meanings"
            ),
        },
        "verifier_command": (
            "python3 code/a5_closure/optimizer_pushforward_certificate.py verify "
            "--manifest code/a5_closure/manifests/optimizer_pushforward_reference.json "
            "--receipt code/a5_closure/receipts/optimizer_pushforward_reference.receipt.json"
        ),
    }


def negative_control_cases(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> list[tuple[str, dict[str, Any], str]]:
    base = base_dir or MODULE_DIR
    cases: list[tuple[str, dict[str, Any], str]] = []

    wrong_schema = copy.deepcopy(dict(manifest))
    wrong_schema["schema"] = "oph.optimizer_pushforward_certificate.v0"
    cases.append(("wrong_schema", wrong_schema, "SCHEMA"))

    wrong_arithmetic = copy.deepcopy(dict(manifest))
    wrong_arithmetic["arithmetic_convention"] = "binary64"
    cases.append(("float_arithmetic_declaration", wrong_arithmetic, "CONVENTION"))

    wrong_pin = copy.deepcopy(dict(manifest))
    wrong_pin["axis_receipt_sha256"] = "0" * 64
    cases.append(("wrong_axis_receipt_pin", wrong_pin, "UPSTREAM_HASH"))

    wrong_target = copy.deepcopy(dict(manifest))
    other_receipt_path = "receipts/axis_center_descent_reference.receipt.json"
    wrong_target["axis_receipt_path"] = other_receipt_path
    wrong_target["axis_receipt_sha256"] = sha256_json(
        load_json(base / other_receipt_path)
    )
    cases.append(("wrong_axis_receipt_target", wrong_target, "UPSTREAM_RECEIPT"))

    for name, key in (
        ("float_tolerance_injection", "float_tolerance"),
        ("approximate_projection_injection", "approximate_projection"),
        ("quantum_state_space_injection", "quantum_state_space"),
        ("declared_defect_injection", "nonclosure_defect"),
    ):
        mutant = copy.deepcopy(dict(manifest))
        mutant[key] = {"declared_without_derivation": True}
        cases.append((name, mutant, "FORBIDDEN_DEPENDENCY"))

    return cases


def negative_control_payload(
    manifest: Mapping[str, Any], base_dir: Path | None = None
) -> dict[str, Any]:
    results = []
    for name, mutant, expected_code in negative_control_cases(manifest, base_dir):
        actual_code = "ACCEPTED"
        try:
            certificate_payload(mutant, base_dir)
        except CertificateError as exc:
            actual_code = exc.code
        require(
            actual_code == expected_code,
            "NEGATIVE_CONTROL_FAILED",
            f"{name}: expected {expected_code}, got {actual_code}",
        )
        results.append(
            {
                "name": name,
                "expected_error": expected_code,
                "actual_error": actual_code,
                "passed": True,
            }
        )
    return {
        "schema": NEGATIVE_SCHEMA,
        "issue": 619,
        "manifest_sha256": sha256_json(manifest),
        "finite_controls": results,
        "countermodel_witnesses": {
            "wrong_reference": (
                "a tilted coarse reference on the single-axis arena moves the "
                "coarse projection to (1/3, 4/15, 1/15, 1/15, 2/15, 2/15) "
                "against the pushed fine optimizer (1/3, 2/15, ..., 2/15)"
            ),
            "changed_weight": (
                "a within-pair tilt changes the split to the reference ratio "
                "2:1 without harming commutation; a cross-axis tilt with the "
                "coarse reference left uniform breaks commutation by 2/15"
            ),
            "incomplete_constraint": (
                "dropping the coarse marginal constraint returns the uniform "
                "reference against the pushed target (1/21, ..., 6/21), a "
                "defect of 5/42 on axis 5"
            ),
            "non_sufficient_channel": (
                "the pinned-port constraint separates states with equal "
                "pushforward and leaves the exact nonclosure defect 5/66"
            ),
        },
    }


def verify_receipt(
    manifest: Mapping[str, Any],
    receipt: Mapping[str, Any],
    base_dir: Path | None = None,
) -> None:
    expected = certificate_payload(manifest, base_dir)
    require(receipt == expected, "RECEIPT_MISMATCH", "receipt is stale, malformed, or tampered")


def default_paths() -> tuple[Path, Path, Path]:
    return (
        MODULE_DIR / "manifests" / "optimizer_pushforward_reference.json",
        MODULE_DIR / "receipts" / "optimizer_pushforward_reference.receipt.json",
        MODULE_DIR / "negative_controls" / "issue_619_negative_controls.json",
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    certify = sub.add_parser("certify")
    certify.add_argument("--manifest", type=Path, required=True)
    certify.add_argument("--output", type=Path, required=True)
    verify = sub.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--receipt", type=Path, required=True)
    negative = sub.add_parser("negative-controls")
    negative.add_argument("--manifest", type=Path, required=True)
    negative.add_argument("--output", type=Path, required=True)
    all_cmd = sub.add_parser("all")
    all_cmd.add_argument("--manifest", type=Path)
    args = parser.parse_args(argv)
    if args.command == "certify":
        manifest = load_json(args.manifest)
        receipt = certificate_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, receipt)
        print(json.dumps({"status": "PASS", "receipt": str(args.output)}, indent=2))
    elif args.command == "verify":
        manifest = load_json(args.manifest)
        receipt = load_json(args.receipt)
        verify_receipt(manifest, receipt, args.manifest.resolve().parent.parent)
        print(json.dumps({"status": "PASS", "receipt": str(args.receipt)}, indent=2))
    elif args.command == "negative-controls":
        manifest = load_json(args.manifest)
        payload = negative_control_payload(manifest, args.manifest.resolve().parent.parent)
        write_json(args.output, payload)
        print(json.dumps({"status": "PASS"}, indent=2))
    else:
        default_manifest, default_receipt, default_negative = default_paths()
        manifest_path = args.manifest or default_manifest
        manifest = load_json(manifest_path)
        write_json(default_receipt, certificate_payload(manifest))
        write_json(default_negative, negative_control_payload(manifest))
        print(
            json.dumps(
                {
                    "status": "PASS",
                    "receipt": str(default_receipt),
                    "negative_controls": str(default_negative),
                },
                indent=2,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
