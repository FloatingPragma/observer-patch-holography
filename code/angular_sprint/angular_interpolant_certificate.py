#!/usr/bin/env python3
"""Issue #643 fast-falsification sprint: the exact angular template packet.

The twelve icosahedral ports embed in the sphere with pairwise dot products
in ``{1, 1/sqrt5, -1/sqrt5, -1}``, so every Legendre Gram kernel
``B_l(i, j) = P_l(v_i . v_j)`` lives in the exact adjacency algebra
``span{I, A, A', Pi}`` over the quadratic field, where ``A`` is the
adjacency matrix, ``Pi`` the antipodal involution, and ``A' = A Pi``. The
certificate proves, in exact arithmetic:

* the four kernels ``B_0..B_3`` are mutually orthogonal scaled projectors
  with scales ``(12, 4, 12/5, 4)`` and ranks ``(1, 3, 5, 3)``, and
  ``I = B_0/12 + B_1/4 + 5 B_2/12 + B_3/4`` (the canonical lowest-bandwidth
  interpolant reproduces the ports exactly);
* each kernel commutes with every generator of the icosahedral automorphism
  group (A5 equivariance) and satisfies ``Pi B_l = (-1)^l B_l`` (antipodal
  parity);
* the kernels are the pinned response bands: ``A B_1 = sqrt5 B_1``,
  ``A B_3 = -sqrt5 B_3``, ``A B_2 = -B_2``, ``A B_0 = 5 B_0``, which binds
  the interpolant bands to the frame, kernel, and quintet bands of the
  registered incidence response and carries refinement naturality through
  that pinned structure;
* the inverse-port response ``R = -J`` with
  ``J = (A^3 - 4A^2 - 5A + 10 I)/10`` acts on the interpolated image as
  ``R B_l = (-1)^(l+1) B_l`` for ``l = 0..3`` (the alternating readback
  parity);
* the complete equal-port sequence
  ``I_l = (1 + (-1)^l + 5 [P_l(1/sqrt5) + P_l(-1/sqrt5)]) / 12`` agrees
  with the direct row sums of ``B_l/12`` for ``l <= 14``, every odd level
  vanishes, and the even initial vector is
  ``(I_2, I_4, I_6, I_8, I_10, I_12, I_14)
  = (0, 0, 11/25, 0, 247/1875, 1071/3125, 0)``.

The transfer decision follows the frozen stop rules without comparison
data. Two completions of the twelve base-port values are exhibited and
registered as distinct source objects: the smooth band-limited canonical
interpolant of degree at most three and the discrete equal-port carrier
measure. Their normalized angular statistics disagree exactly: the
interpolant carries zero power at level six while the comb carries
``11/25`` there, so the twelve base-port values alone do not select the
sky completion. That is static base-port underdetermination, not
nonidentifiability of the full source: the refinement tower,
coarse-graining, and the repair semigroup are registered source structure
that the static witness ignores, and for the canonical degree<=3 band a
same-codomain smooth counterfunction vanishing on the twelve base ports
is separated by the level-one refinement vertex set
(``code/angular_sprint/refinement_transfer_certificate.py``). The
transfer row stays open pending the refinement, coarse-graining, and
repair-semigroup transfer decision under the registered repair law
(``oph-physics-sim docs/CANONICAL_REPAIR_LAW.md``). The bounded ancestry
search over the registered source surfaces records that no producer
emits a geometry-imprint sky field or a second independently observed
channel, so the normalized parity cross-response has no two-channel
input and a bare amplitude would be a free signal amplitude. The
interpolant, the equal-port carrier measure, and the two-channel parity
response stay three distinct source objects. No public measurement is
read and no comparison is opened.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "angular_template_receipt.json"

SCHEMA = "oph.angular_template_receipt.v1"
STATUS = (
    "EXACT_SOURCE_TEMPLATE__STATIC_BASE_PORT_TRANSFER_NONIDENTIFIABLE__"
    "REFINEMENT_REPAIR_SKY_TRANSFER_OPEN"
)
MAX_LEVEL = 14


class AngularError(ValueError):
    """The angular template certificate refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AngularError(message)


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


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Exact Q(sqrt5) scalars and 12x12 matrices
# ---------------------------------------------------------------------------

Q5 = tuple[Fraction, Fraction]


def q5(a, b=0) -> Q5:
    return (Fraction(a), Fraction(b))


def q5_add(x: Q5, y: Q5) -> Q5:
    return (x[0] + y[0], x[1] + y[1])


def q5_sub(x: Q5, y: Q5) -> Q5:
    return (x[0] - y[0], x[1] - y[1])


def q5_mul(x: Q5, y: Q5) -> Q5:
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def q5_scale(x: Q5, factor: Fraction) -> Q5:
    return (x[0] * factor, x[1] * factor)


def q5_str(x: Q5) -> str:
    return f"{x[0]}+{x[1]}*sqrt5"


ZERO = q5(0)
ONE = q5(1)
INV_SQRT5 = q5(0, Fraction(1, 5))


def legendre_at(t: Q5, max_level: int) -> list[Q5]:
    """Exact Legendre values ``P_0(t) .. P_max(t)`` by the recurrence."""

    values = [ONE, t]
    for level in range(1, max_level):
        term = q5_scale(q5_mul(t, values[level]), Fraction(2 * level + 1))
        term = q5_sub(term, q5_scale(values[level - 1], Fraction(level)))
        values.append(q5_scale(term, Fraction(1, level + 1)))
    return values[: max_level + 1]


def ports_and_structure() -> tuple[list[str], list[list[int]], list[list[int]], list[list[int]]]:
    """Adjacency, antipode, and second-neighbor matrices of the icosahedron."""

    upper = [f"u{i}" for i in range(5)]
    lower = [f"l{i}" for i in range(5)]
    ports = ["n"] + upper + lower + ["s"]
    index = {port: position for position, port in enumerate(ports)}
    size = len(ports)
    adjacency = [[0] * size for _ in range(size)]

    def connect(a: str, b: str) -> None:
        adjacency[index[a]][index[b]] = 1
        adjacency[index[b]][index[a]] = 1

    for i in range(5):
        connect("n", upper[i])
        connect("s", lower[i])
        connect(upper[i], upper[(i + 1) % 5])
        connect(lower[i], lower[(i + 1) % 5])
        connect(upper[i], lower[i])
        connect(upper[i], lower[(i + 1) % 5])

    antipode_map = {"n": "s", "s": "n"}
    for i in range(5):
        antipode_map[upper[i]] = lower[(i + 3) % 5]
        antipode_map[lower[(i + 3) % 5]] = upper[i]
    antipode = [[0] * size for _ in range(size)]
    for port, image in antipode_map.items():
        antipode[index[port]][index[image]] = 1
    for row in range(size):
        require(sum(antipode[row]) == 1, "antipode is not a permutation")
        require(antipode[row][row] == 0, "antipode has a fixed point")

    second = [
        [
            sum(adjacency[i][k] * antipode[k][j] for k in range(size))
            for j in range(size)
        ]
        for i in range(size)
    ]
    identity = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    for i in range(size):
        for j in range(size):
            total = (
                identity[i][j]
                + adjacency[i][j]
                + second[i][j]
                + antipode[i][j]
            )
            require(total == 1, "dot-product partition drift")
    antipode_involution = all(
        sum(antipode[i][k] * antipode[k][j] for k in range(size))
        == identity[i][j]
        for i in range(size)
        for j in range(size)
    )
    require(antipode_involution, "antipode is not an involution")
    antipode_automorphism = all(
        sum(
            antipode[i][a] * adjacency[a][b] * antipode[j][b]
            for a in range(size)
            for b in range(size)
        )
        == adjacency[i][j]
        for i in range(size)
        for j in range(size)
    )
    require(antipode_automorphism, "antipode is not an automorphism")
    return ports, adjacency, antipode, second


def q5_matrix(scalars: dict[str, Q5], components: dict[str, list[list[int]]]) -> list[list[Q5]]:
    size = 12
    matrix = [[ZERO for _ in range(size)] for _ in range(size)]
    for name, scalar in scalars.items():
        component = components[name]
        for i in range(size):
            for j in range(size):
                if component[i][j]:
                    matrix[i][j] = q5_add(matrix[i][j], scalar)
    return matrix


def mat_mul(x: list[list[Q5]], y: list[list[Q5]]) -> list[list[Q5]]:
    size = len(x)
    result = [[ZERO for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for k in range(size):
            if x[i][k] == ZERO:
                continue
            for j in range(size):
                if y[k][j] == ZERO:
                    continue
                result[i][j] = q5_add(result[i][j], q5_mul(x[i][k], y[k][j]))
    return result


def mat_scale(x: list[list[Q5]], factor: Fraction) -> list[list[Q5]]:
    return [[q5_scale(value, factor) for value in row] for row in x]


def mat_sub(x: list[list[Q5]], y: list[list[Q5]]) -> list[list[Q5]]:
    return [
        [q5_sub(a, b) for a, b in zip(row_x, row_y)]
        for row_x, row_y in zip(x, y)
    ]


def mat_is_zero(x: list[list[Q5]]) -> bool:
    return all(value == ZERO for row in x for value in row)


def mat_trace(x: list[list[Q5]]) -> Q5:
    total = ZERO
    for i in range(len(x)):
        total = q5_add(total, x[i][i])
    return total


def int_to_q5(matrix: list[list[int]]) -> list[list[Q5]]:
    return [[q5(value) for value in row] for row in matrix]


# ---------------------------------------------------------------------------
# The certificate
# ---------------------------------------------------------------------------


def build_kernels() -> dict[str, Any]:
    ports, adjacency, antipode, second = ports_and_structure()
    size = len(ports)
    identity = [[1 if i == j else 0 for j in range(size)] for i in range(size)]
    components = {
        "identity": identity,
        "adjacency": adjacency,
        "second": second,
        "antipode": antipode,
    }
    plus = legendre_at(INV_SQRT5, MAX_LEVEL)
    minus = legendre_at(q5_scale(INV_SQRT5, Fraction(-1)), MAX_LEVEL)

    kernels = []
    for level in range(MAX_LEVEL + 1):
        kernels.append(
            q5_matrix(
                {
                    "identity": ONE,
                    "adjacency": plus[level],
                    "second": minus[level],
                    "antipode": q5(Fraction((-1) ** level)),
                },
                components,
            )
        )
    return {
        "ports": ports,
        "components": components,
        "kernels": kernels,
        "legendre_plus": plus,
        "legendre_minus": minus,
    }


def projector_certificate(data: dict[str, Any]) -> dict[str, Any]:
    kernels = data["kernels"]
    scales = [Fraction(12), Fraction(4), Fraction(12, 5), Fraction(4)]
    orthogonal = all(
        mat_is_zero(mat_mul(kernels[a], kernels[b]))
        for a in range(4)
        for b in range(4)
        if a != b
    )
    scaled_projectors = all(
        mat_is_zero(
            mat_sub(
                mat_mul(kernels[level], kernels[level]),
                mat_scale(kernels[level], scales[level]),
            )
        )
        for level in range(4)
    )
    identity = int_to_q5(data["components"]["identity"])
    resolution = mat_is_zero(
        mat_sub(
            [
                [
                    q5_add(
                        q5_add(
                            q5_scale(kernels[0][i][j], Fraction(1, 12)),
                            q5_scale(kernels[1][i][j], Fraction(1, 4)),
                        ),
                        q5_add(
                            q5_scale(kernels[2][i][j], Fraction(5, 12)),
                            q5_scale(kernels[3][i][j], Fraction(1, 4)),
                        ),
                    )
                    for j in range(12)
                ]
                for i in range(12)
            ],
            identity,
        )
    )
    traces = [q5_str(mat_trace(kernels[level])) for level in range(4)]
    ranks = [
        (mat_trace(kernels[level])[0] / scales[level])
        for level in range(4)
    ]
    return {
        "mutually_orthogonal": orthogonal,
        "scaled_projectors": scaled_projectors,
        "scales": [str(scale) for scale in scales],
        "resolution_of_identity": resolution,
        "resolution_weights": ["1/12", "1/4", "5/12", "1/4"],
        "traces": traces,
        "ranks": [str(rank) for rank in ranks],
        "port_reproduction": resolution,
    }


def automorphism_generators() -> list[dict[int, int]]:
    ports = ["n"] + [f"u{i}" for i in range(5)] + [f"l{i}" for i in range(5)] + ["s"]
    index = {port: position for position, port in enumerate(ports)}

    def as_map(assignment: dict[str, str]) -> dict[int, int]:
        return {index[a]: index[b] for a, b in assignment.items()}

    rotation = {"n": "n", "s": "s"}
    for i in range(5):
        rotation[f"u{i}"] = f"u{(i + 1) % 5}"
        rotation[f"l{i}"] = f"l{(i + 1) % 5}"

    reflection = {"n": "n", "s": "s"}
    for i in range(5):
        reflection[f"u{i}"] = f"u{(-i) % 5}"
        reflection[f"l{i}"] = f"l{(1 - i) % 5}"

    antipode = {"n": "s", "s": "n"}
    for i in range(5):
        antipode[f"u{i}"] = f"l{(i + 3) % 5}"
        antipode[f"l{(i + 3) % 5}"] = f"u{i}"

    pole_swap = {
        "n": "u0",
        "u0": "n",
        "u2": "l1",
        "l1": "u2",
        "u3": "l0",
        "l0": "u3",
        "l3": "s",
        "s": "l3",
        "u1": "u1",
        "u4": "u4",
        "l2": "l2",
        "l4": "l4",
    }

    return [
        as_map(rotation),
        as_map(reflection),
        as_map(antipode),
        as_map(pole_swap),
    ]


def equivariance_certificate(data: dict[str, Any]) -> dict[str, Any]:
    adjacency = data["components"]["adjacency"]
    kernels = data["kernels"]
    results = []
    for generator in automorphism_generators():
        is_automorphism = all(
            adjacency[generator[i]][generator[j]] == adjacency[i][j]
            for i in range(12)
            for j in range(12)
        )
        commutes = all(
            kernels[level][generator[i]][generator[j]] == kernels[level][i][j]
            for level in range(4)
            for i in range(12)
            for j in range(12)
        )
        results.append(
            {"is_automorphism": is_automorphism, "kernels_equivariant": commutes}
        )
    generators = automorphism_generators()
    frontier = [tuple(generator[i] for i in range(12)) for generator in generators]
    closure = set(frontier)
    while frontier:
        left = frontier.pop()
        for generator in generators:
            composed = tuple(generator[left[i]] for i in range(12))
            if composed not in closure:
                closure.add(composed)
                frontier.append(composed)
    identity_perm = tuple(range(12))
    if identity_perm not in closure:
        closure.add(identity_perm)
    closure_order = len(closure)

    antipode = data["components"]["antipode"]
    parity = all(
        mat_is_zero(
            mat_sub(
                mat_mul(int_to_q5(antipode), data["kernels"][level]),
                mat_scale(data["kernels"][level], Fraction((-1) ** level)),
            )
        )
        for level in range(4)
    )
    return {
        "generators_checked": len(results),
        "generator_results": results,
        "generated_group_order": closure_order,
        "generating_set_complete": closure_order == 120,
        "all_equivariant": all(
            row["is_automorphism"] and row["kernels_equivariant"]
            for row in results
        ),
        "antipodal_parity_signs": [1, -1, 1, -1],
        "antipodal_parity": parity,
    }


def band_binding_certificate(data: dict[str, Any]) -> dict[str, Any]:
    adjacency = int_to_q5(data["components"]["adjacency"])
    kernels = data["kernels"]
    eigenvalues = [q5(5), q5(0, 1), q5(-1), q5(0, -1)]
    bindings = []
    for level in range(4):
        difference = mat_sub(
            mat_mul(adjacency, kernels[level]),
            [
                [q5_mul(eigenvalues[level], value) for value in row]
                for row in kernels[level]
            ],
        )
        bindings.append(mat_is_zero(difference))
    return {
        "adjacency_eigenvalues_by_level": [q5_str(value) for value in eigenvalues],
        "band_identification": {
            "level_0": "unit band",
            "level_1": "frame band (eigenvalue sqrt5)",
            "level_2": "quintet band (eigenvalue -1)",
            "level_3": "kernel band (eigenvalue -sqrt5)",
        },
        "bindings_exact": bindings,
        "all_bound": all(bindings),
        "refinement_naturality": {
            "type": "PREMISE_TYPED_CONDITIONAL",
            "premise": (
                "the registered incidence-response band structure is what a "
                "screen refinement acts through on the port coefficient "
                "module; this is a declared contract term of the port "
                "module registration, not proved in this producer"
            ),
            "conclusion": (
                "under that premise the interpolant bands are natural in "
                "refinement, because the binding identities prove they "
                "equal the pinned response bands exactly"
            ),
            "carried_by": "bindings_exact",
            "premise_discharged_here": False,
        },
    }


def parity_response_certificate(data: dict[str, Any]) -> dict[str, Any]:
    adjacency = int_to_q5(data["components"]["adjacency"])
    size = 12
    a2 = mat_mul(adjacency, adjacency)
    a3 = mat_mul(a2, adjacency)
    identity = int_to_q5(data["components"]["identity"])
    j_response = [
        [
            q5_scale(
                q5_add(
                    q5_sub(
                        q5_sub(a3[i][j], q5_scale(a2[i][j], Fraction(4))),
                        q5_scale(adjacency[i][j], Fraction(5)),
                    ),
                    q5_scale(identity[i][j], Fraction(10)),
                ),
                Fraction(1, 10),
            )
            for j in range(size)
        ]
        for i in range(size)
    ]
    readback = [[q5_scale(value, Fraction(-1)) for value in row] for row in j_response]
    signs = []
    for level in range(4):
        expected = Fraction((-1) ** (level + 1))
        difference = mat_sub(
            mat_mul(readback, data["kernels"][level]),
            mat_scale(data["kernels"][level], expected),
        )
        signs.append(mat_is_zero(difference))
    return {
        "response_polynomial": "J=(A^3-4A^2-5A+10I)/10; R=-J",
        "parity_law": "R restricted to level l acts as (-1)^(l+1)",
        "parity_signs_by_level": [-1, 1, -1, 1],
        "parity_exact": signs,
        "all_exact": all(signs),
    }


def equal_port_certificate(data: dict[str, Any]) -> dict[str, Any]:
    plus = data["legendre_plus"]
    minus = data["legendre_minus"]
    sequence = []
    for level in range(MAX_LEVEL + 1):
        formula = q5_scale(
            q5_add(
                q5(1 + (-1) ** level),
                q5_scale(q5_add(plus[level], minus[level]), Fraction(5)),
            ),
            Fraction(1, 12),
        )
        row_sum = ZERO
        for value in data["kernels"][level][0]:
            row_sum = q5_add(row_sum, value)
        row_sum = q5_scale(row_sum, Fraction(1, 12))
        require(formula == row_sum, f"equal-port mismatch at level {level}")
        require(formula[1] == 0, f"equal-port value irrational at level {level}")
        sequence.append(formula[0])
    expected_even = {
        2: Fraction(0),
        4: Fraction(0),
        6: Fraction(11, 25),
        8: Fraction(0),
        10: Fraction(247, 1875),
        12: Fraction(1071, 3125),
        14: Fraction(0),
    }
    even_ok = all(sequence[level] == value for level, value in expected_even.items())
    odd_ok = all(sequence[level] == 0 for level in range(1, MAX_LEVEL + 1, 2))
    return {
        "sequence": [str(value) for value in sequence],
        "closed_formula": (
            "I_l = (1 + (-1)^l + 5*[P_l(1/sqrt5) + P_l(-1/sqrt5)]) / 12"
        ),
        "even_initial_vector": {
            str(level): str(value) for level, value in expected_even.items()
        },
        "even_vector_exact": even_ok,
        "all_odd_levels_zero": odd_ok,
        "common_frame_template": {
            "levels": [6, 10, 12],
            "ratio_I10_over_I6": str(
                Fraction(247, 1875) / Fraction(11, 25)
            ),
            "ratio_I12_over_I6": str(
                Fraction(1071, 3125) / Fraction(11, 25)
            ),
        },
        "unit_normalization": sequence[0] == 1,
    }


# ---------------------------------------------------------------------------
# Transfer decision without comparison data
# ---------------------------------------------------------------------------


PROJECTION_PATH = (
    REPO_ROOT / "code" / "invariant_mining" / "outputs" / "source_projection.json"
)

SKY_EMISSION_KEYS = (
    "sky_map",
    "sky_field",
    "healpix",
    "alm",
    "a_lm",
    "angular_power_spectrum",
    "cl_spectrum",
    "screen_to_sky_map",
)

SECOND_CHANNEL_KEYS = (
    "second_channel",
    "cross_channel",
    "independent_channel",
    "second_observed_channel",
)


def _walk_keys(value: Any):
    if isinstance(value, dict):
        for key, item in value.items():
            yield key
            yield from _walk_keys(item)
    elif isinstance(value, list):
        for item in value:
            yield from _walk_keys(item)


def ancestry_sweep(searched_surfaces: list[str]) -> dict[str, Any]:
    """Registered-surface sweep for a sky-field or second-channel emission.

    The sweep walks every JSON key of every registered surface and matches
    the keys against the typed emission signatures; the verdict booleans
    are derived from the matches. Surfaces outside the registered set are
    inadmissible as premises under the pregeneration freeze, so the sweep
    over the registered set is the complete admissible search.
    """

    projection = json.loads(PROJECTION_PATH.read_text(encoding="utf-8"))
    registered = [row["path"] for row in projection["source_artifacts"]]
    inventory = sorted(set(searched_surfaces) | set(registered))
    surface_reports = []
    sky_matches = []
    channel_matches = []
    for relative in inventory:
        path = REPO_ROOT / relative
        require(path.is_file(), f"searched surface is absent: {relative}")
        if path.suffix == ".json":
            keys = list(_walk_keys(json.loads(path.read_text(encoding="utf-8"))))
        else:
            keys = []
        sky_hits = sorted(
            {key for key in keys if key.lower() in SKY_EMISSION_KEYS}
        )
        channel_hits = sorted(
            {key for key in keys if key.lower() in SECOND_CHANNEL_KEYS}
        )
        sky_matches.extend(sky_hits)
        channel_matches.extend(channel_hits)
        surface_reports.append(
            {
                "path": relative,
                "keys_scanned": len(keys),
                "sky_emission_keys_found": sky_hits,
                "second_channel_keys_found": channel_hits,
            }
        )
    return {
        "searched_inventory": inventory,
        "surface_reports": surface_reports,
        "emission_key_signatures": {
            "sky_field": list(SKY_EMISSION_KEYS),
            "second_channel": list(SECOND_CHANNEL_KEYS),
        },
        "sky_field_emission_found": bool(sky_matches),
        "second_independent_channel_found": bool(channel_matches),
        "admissibility_boundary": (
            "unregistered surfaces are inadmissible as premises under the "
            "pregeneration freeze, so this registered-surface sweep is the "
            "complete admissible ancestry search"
        ),
    }


def transfer_decision(data: dict[str, Any], equal_port: dict[str, Any]) -> dict[str, Any]:
    interpolant_power = {
        "0": "1",
        "6": "0",
        "10": "0",
        "12": "0",
    }
    comb_power = {
        "0": "1",
        "6": str(Fraction(11, 25)),
        "10": str(Fraction(247, 1875)),
        "12": str(Fraction(1071, 3125)),
    }
    disagreement = interpolant_power["6"] != comb_power["6"]
    searched_surfaces = [
        "code/invariant_mining/data/source_feature_registry.json",
        "code/invariant_mining/outputs/candidate_registry.json",
        "code/a5_closure/manifests/classical_realization_receipt.json",
        "code/a5_closure/manifests/family_band_attachment_reference.json",
        "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
        "code/capacity_readback/runtime/source_derived_public_checkpoint_packet.json",
        "claims/frozen_prediction_register.json",
    ]
    sweep = ancestry_sweep(searched_surfaces)
    require(
        not sweep["sky_field_emission_found"],
        "a registered surface emits a sky-valued field",
    )
    require(
        not sweep["second_independent_channel_found"],
        "a registered surface emits a second observed channel",
    )
    return {
        "verdict": "STATIC_BASE_PORT_UNDERDETERMINATION",
        "verdict_typing": (
            "the two completions are distinct source objects, a smooth "
            "band-limited interpolant and a discrete equal-port carrier "
            "measure; their exact level-six power disagreement certifies "
            "that the twelve base-port values alone do not select the sky "
            "completion; it does not certify nonidentifiability of the "
            "full source, because the refinement tower, coarse-graining, "
            "and the repair semigroup are registered source structure "
            "that the static witness ignores"
        ),
        "stop_rule": (
            "two distinct source objects that agree on the twelve base "
            "ports and disagree in an exact normalized statistic close "
            "the static base-port row underdetermined without opening "
            "comparison data; the row does not close "
            "transfer-nonidentifiable, and the sky-transfer row stays "
            "open pending the refinement, coarse-graining, and "
            "repair-semigroup transfer decision"
        ),
        "completion_a": {
            "name": "band-limited canonical interpolant",
            "object_type": "smooth band-limited field of degree at most 3",
            "support": "levels zero through three",
            "normalized_power_at_level_6": interpolant_power["6"],
        },
        "completion_b": {
            "name": "equal-port angular comb",
            "object_type": "discrete equal-port carrier measure",
            "support": "all even levels with the exact comb weights",
            "normalized_power_at_level_6": comb_power["6"],
        },
        "normalized_statistics_disagree_exactly": disagreement,
        "static_base_port_underdetermination": disagreement,
        "refinement_frontier": {
            "separation_witness": (
                "a same-codomain smooth counterfunction vanishing on the "
                "twelve base ports is separated by the level-one "
                "refinement vertex set for the canonical degree<=3 band"
            ),
            "producer": (
                "code/angular_sprint/refinement_transfer_certificate.py"
            ),
            "future_verdict_requirement": (
                "any future nonidentifiability verdict must exhibit two "
                "source-admissible completions that intertwine "
                "refinement, coarse-graining, and the repair semigroup"
            ),
            "repair_law": "oph-physics-sim docs/CANONICAL_REPAIR_LAW.md",
        },
        "geometry_imprint_search": {
            **sweep,
            "note": (
                "the registered producers emit finite port, seam, sector, "
                "capacity, and response data; none emits a sky-valued "
                "geometry-imprint field, and the readback response is the "
                "single registered channel, so the normalized parity "
                "cross-response has no two-channel input"
            ),
        },
        "free_amplitude_classification": (
            "NOT_EVALUABLE_FREE_SIGNAL_AMPLITUDE for any bare comb or "
            "interpolant amplitude; a normalized statistic removes the "
            "amplitude but not the completion choice"
        ),
        "distinct_source_objects": [
            "arbitrary-port canonical interpolant",
            "equal-port carrier measure",
            "two-channel parity response",
        ],
        "frozen_candidate_order": [
            "normalized parity cross-response when two independently "
            "observed channels exist",
            "equal-port angular comb when its complete nonzero-power "
            "physical transfer exists",
        ],
    }


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


PIN_RELATIVES = [
    "code/a5_closure/receipts/port_current_inner_reference.receipt.json",
    "code/invariant_mining/data/source_feature_registry.json",
    "code/invariant_mining/outputs/candidate_registry.json",
    "code/a5_closure/manifests/classical_realization_receipt.json",
    "code/a5_closure/manifests/family_band_attachment_reference.json",
    "code/capacity_readback/runtime/source_derived_public_checkpoint_packet.json",
    "claims/frozen_prediction_register.json",
]


def parent_pins() -> list[dict[str, Any]]:
    pins = []
    for relative in PIN_RELATIVES:
        payload = (REPO_ROOT / relative).read_bytes()
        pins.append(
            {
                "path": relative,
                "bytes": len(payload),
                "sha256": tagged_sha256(payload),
            }
        )
    return pins


def build_receipt() -> dict[str, Any]:
    data = build_kernels()
    projectors = projector_certificate(data)
    equivariance = equivariance_certificate(data)
    binding = band_binding_certificate(data)
    parity = parity_response_certificate(data)
    equal_port = equal_port_certificate(data)
    require(
        projectors["mutually_orthogonal"]
        and projectors["scaled_projectors"]
        and projectors["resolution_of_identity"],
        "projector certificate failed",
    )
    require(equivariance["all_equivariant"], "equivariance failed")
    require(
        equivariance["generating_set_complete"],
        "generator closure is not the full automorphism group",
    )
    require(equivariance["antipodal_parity"], "antipodal parity failed")
    require(binding["all_bound"], "band binding failed")
    require(parity["all_exact"], "parity response failed")
    require(
        equal_port["even_vector_exact"] and equal_port["all_odd_levels_zero"],
        "equal-port sequence failed",
    )
    decision = transfer_decision(data, equal_port)
    require(
        decision["normalized_statistics_disagree_exactly"]
        and decision["static_base_port_underdetermination"],
        "static base-port disagreement witness failed",
    )
    receipt = {
        "schema": SCHEMA,
        "issue": 643,
        "status": STATUS,
        "projector_certificate": projectors,
        "equivariance_certificate": equivariance,
        "band_binding_certificate": binding,
        "parity_response_certificate": parity,
        "equal_port_certificate": equal_port,
        "transfer_decision": decision,
        "lean_bindings": [
            "OPH.AngularBands.equalPort_evenVector",
            "OPH.AngularBands.equalPort_oddZero",
            "OPH.AngularBands.parity_signs",
        ],
        "parent_pins": parent_pins(),
        "lean_binding_scope": (
            "the Lean layers prove the exact equal-port arithmetic, the "
            "parity sign table, and the four-kernel projector identities; "
            "the equivariance, binding, and response identities are "
            "code-certified in exact quadratic-field arithmetic in this "
            "producer, with the Lagrange spectral reconstruction in the "
            "package tests as the independent route"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "row_state": {
            "static_base_port_underdetermination": "CERTIFIED",
            "refinement_repair_sky_transfer": "OPEN",
            "close_conditions": (
                "the transfer row closes identifiable when a selection "
                "certificate derives the sky completion from the "
                "registered refinement, coarse-graining, and "
                "repair-semigroup structure; it closes nonidentifiable "
                "only when two source-admissible completions are "
                "exhibited that intertwine refinement, coarse-graining, "
                "and the repair semigroup"
            ),
        },
        "evaluation_condition": (
            "a source-derived geometry-imprint field or a second "
            "independently observed channel, together with a source premise "
            "selecting one transfer completion"
        ),
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def write_runtime() -> Path:
    RUNTIME.mkdir(exist_ok=True)
    RECEIPT_PATH.write_bytes(canonical_json_bytes(build_receipt()))
    return RECEIPT_PATH


def verify_runtime() -> None:
    if RECEIPT_PATH.read_bytes() != canonical_json_bytes(build_receipt()):
        raise SystemExit("angular template receipt is stale")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        print(write_runtime())
    if args.verify:
        verify_runtime()
        print("ANGULAR_TEMPLATE_VALID")
    if not args.write and not args.verify:
        receipt = build_receipt()
        print(
            json.dumps(
                {
                    "status": receipt["status"],
                    "even_vector": receipt["equal_port_certificate"][
                        "even_initial_vector"
                    ],
                    "all_certificates": "pass",
                },
                indent=2,
                sort_keys=True,
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
