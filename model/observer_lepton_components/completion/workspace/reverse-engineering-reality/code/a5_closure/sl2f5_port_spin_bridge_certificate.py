#!/usr/bin/env python3
"""Constructive exact bridge from SL(2,F5) to the OPH PORT-SPIN-LIFT.

This certificate composes two independently existing exact constructions:

* ``coset_carrier_certificate.py`` supplies ``SL(2,F5)``, its central
  quotient, and an explicit quotient action which relabels onto the sixty
  committed twelve-port rotations;
* ``super_tannakian_matter_lift_certificate.py`` independently reconstructs
  the exact SU(2) PORT-SPIN-LIFT over Q(sqrt(5)) from the #314 current carrier.

The two constructions use isomorphic but not identical carrier coordinate
labellings.  The bridge therefore derives an exact orientation-preserving
carrier relabelling before comparing the sixty-row actions; no row equality is
assumed across coordinate systems.

The bridge is not inferred from order, element-order profile, or uniqueness of
a non-split extension.  It chooses one exact presentation triple

    a^2 = b^3 = c^5 = abc = -I

in ``SL(2,F5)``, matches its three quotient port rows to the independently
computed spin lifts, exhausts the eight possible signs, and keeps the unique
sign choice satisfying the same presentation in the spin group.  The resulting
generator map is extended by words to all 120 source elements and checked on
all 120^2 products.  Finally the map to the committed port permutations is
checked to commute elementwise through both double covers.

The output is therefore an explicit finite-group isomorphism certificate
between the canonical ``SL(2,F5)`` source and the executable
``PORT-SPIN-LIFT`` matrix group, together with a faithful exact two-dimensional
complex representation.  It does not invoke McKay/E8, select phi, identify a
physical rotation group, or make a mass claim.  The PORT-SPIN-LIFT remains
conditional on the declared #566 current fixture from which it is reconstructed.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent

import coset_carrier_certificate as c784  # noqa: E402
import super_tannakian_matter_lift_certificate as m314  # noqa: E402

SCHEMA = "oph.sl2f5_port_spin_bridge.v1"

Matrix = c784.Matrix
SpinMatrix = list[list[m314.C5]]
PortRow = tuple[int, ...]


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise c784.CertificateError(code, message)


def spin_equal(left: Sequence[Sequence[m314.C5]], right: Sequence[Sequence[m314.C5]]) -> bool:
    return m314.c_is_zero(m314.csub(left, right))


def spin_neg(matrix: Sequence[Sequence[m314.C5]]) -> SpinMatrix:
    return m314.cscale(matrix, m314.C5(-m314.ONE, m314.ZERO))


def spin_pow(matrix: Sequence[Sequence[m314.C5]], exponent: int) -> SpinMatrix:
    out = m314.cidentity(2)
    for _ in range(exponent):
        out = m314.cmul(out, matrix)
    return out


def permutation_compose(left: PortRow, right: PortRow) -> PortRow:
    """Composition ``left`` after ``right``, matching the committed convention."""

    return tuple(left[right[i]] for i in range(len(left)))


def inverse_permutation(p: Sequence[int]) -> PortRow:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def relabel_port_row(row: PortRow, phi: Sequence[int]) -> PortRow:
    """Conjugate a permutation from carrier labels to committed labels.

    ``phi`` maps old carrier vertex labels to committed port labels, so the
    transported row is ``phi ∘ row ∘ phi^{-1}``.
    """

    inv = inverse_permutation(phi)
    return tuple(phi[row[inv[j]]] for j in range(len(row)))


def source_power(source: c784.Source, value: Matrix, exponent: int) -> Matrix:
    out = source.identity
    for _ in range(exponent):
        out = source.mul(out, value)
    return out


def presentation_solutions(source: c784.Source) -> list[tuple[Matrix, Matrix, Matrix]]:
    """All exact solutions of a^2 = b^3 = c^5 = abc = -I."""

    minus = source.minus_identity
    out: list[tuple[Matrix, Matrix, Matrix]] = []
    for a in source.elements:
        if source_power(source, a, 2) != minus:
            continue
        for b in source.elements:
            if source_power(source, b, 3) != minus:
                continue
            c = source.mul(source.inv(source.mul(a, b)), minus)
            if source_power(source, c, 5) == minus:
                out.append((a, b, c))
    return out


def placement_of_solution(
    quotient: c784.Quotient,
    solution: tuple[Matrix, Matrix, Matrix],
) -> c784.Placement:
    """The (C5,C3,C2) quotient placement carried by one presentation triple."""

    source = quotient.source
    a, b, c = solution
    return tuple(
        frozenset(
            quotient.canon(g)
            for g in c784.cyclic_subgroup(source.mul, source.identity, generator)
        )
        for generator in (c, b, a)
    )  # type: ignore[return-value]


def port_action_for_placement(
    quotient: c784.Quotient,
    placement: c784.Placement,
    phi: Sequence[int],
) -> Callable[[Matrix], PortRow]:
    """Return the explicit source-to-committed-port action fixed by ``phi``."""

    i_vertex, representatives = c784.coset_index(quotient, placement[0])
    port_to_vertex = {port: vertex for vertex, port in enumerate(phi)}

    def source_to_port(source_element: Matrix) -> PortRow:
        g = quotient.canon(source_element)
        return tuple(
            phi[
                i_vertex[
                    quotient.mul(g, representatives[port_to_vertex[port]])
                ]
            ]
            for port in range(12)
        )

    return source_to_port


def current_carrier_relabel(
    algebra: m314.CurrentAlgebra,
    fixture: c784.Fixture,
) -> tuple[PortRow, dict[str, int]]:
    """Derive an orientation-preserving relabelling of the #314 carrier."""

    isos = c784.graph_isomorphisms(algebra.carrier.adjacency, fixture.adjacency)
    require(
        bool(isos),
        "SPIN_CARRIER_RELABEL",
        "#314 carrier is not graph-isomorphic to the committed carrier",
    )

    committed = c784.relabel_faces(fixture.oriented_faces, tuple(range(12)))
    matching = [
        tuple(phi)
        for phi in isos
        if c784.relabel_faces(algebra.carrier.faces, phi) == committed
    ]
    require(
        bool(matching),
        "SPIN_CARRIER_ORIENTATION",
        "no carrier relabelling preserves the committed face orientation",
    )

    phi = min(matching)
    transported = {relabel_port_row(tuple(row), phi) for row in algebra.plus}
    require(
        transported == set(fixture.port_action),
        "SPIN_PORT_ACTION",
        "relabelled PORT-SPIN-LIFT rows differ from the committed port action",
    )
    return phi, {
        "graph_isomorphisms": len(isos),
        "orientation_matching_relabellings": len(matching),
    }


def build_spin_group_index(
    lifts: Mapping[PortRow, SpinMatrix],
) -> tuple[dict[tuple, SpinMatrix], dict[tuple, PortRow]]:
    """Index the full {+U_g,-U_g} spin group and its quotient port row."""

    elements: dict[tuple, SpinMatrix] = {}
    to_port: dict[tuple, PortRow] = {}
    for row, lift in lifts.items():
        for matrix in (lift, spin_neg(lift)):
            key = m314.matrix_key(matrix)
            if key in to_port:
                require(
                    to_port[key] == row,
                    "SPIN_PORT_AMBIGUOUS",
                    "one spin element is attached to two distinct port rows",
                )
            elements[key] = matrix
            to_port[key] = row
    require(len(elements) == 120, "SPIN_GROUP_ORDER", f"spin group has {len(elements)} elements")
    require(len(to_port) == 120, "SPIN_PORT_MAP", "spin-to-port map is not defined on all 120 lifts")
    return elements, to_port


def signed_generator_choices(
    base_lifts: Sequence[SpinMatrix],
) -> list[tuple[tuple[int, int, int], tuple[SpinMatrix, SpinMatrix, SpinMatrix]]]:
    """Sign choices satisfying A^2=B^3=C^5=ABC=-I exactly."""

    minus_identity = spin_neg(m314.cidentity(2))
    passing = []
    for signs in itertools.product((1, -1), repeat=3):
        images = tuple(
            matrix if sign == 1 else spin_neg(matrix)
            for sign, matrix in zip(signs, base_lifts)
        )
        A, B, C = images
        abc = m314.cmul(m314.cmul(A, B), C)
        if (
            spin_equal(spin_pow(A, 2), minus_identity)
            and spin_equal(spin_pow(B, 3), minus_identity)
            and spin_equal(spin_pow(C, 5), minus_identity)
            and spin_equal(abc, minus_identity)
        ):
            passing.append((signs, images))
    return passing


def extend_generator_map(
    source: c784.Source,
    generators: Sequence[Matrix],
    images: Sequence[SpinMatrix],
) -> dict[Matrix, SpinMatrix]:
    """Extend a generator assignment along a BFS word tree, checking consistency."""

    mapping: dict[Matrix, SpinMatrix] = {source.identity: m314.cidentity(2)}
    queue = [source.identity]
    for g in queue:
        for source_generator, image_generator in zip(generators, images):
            h = source.mul(g, source_generator)
            image = m314.cmul(mapping[g], image_generator)
            if h in mapping:
                require(
                    spin_equal(mapping[h], image),
                    "WORD_MAP_INCONSISTENT",
                    "two source words give different spin matrices",
                )
            else:
                mapping[h] = image
                queue.append(h)
    require(
        len(mapping) == len(source.elements) == 120,
        "SOURCE_GENERATION",
        f"presentation generators reach {len(mapping)} source elements",
    )
    return mapping


def matrix_as_list(matrix: Matrix) -> list[int]:
    return list(matrix)


def row_as_list(row: PortRow) -> list[int]:
    return list(row)


def build_certificate() -> dict[str, Any]:
    # --- Canonical SL(2,F5) source and its exact committed-port quotient -----
    coset_manifest = c784.load_json(MODULE_DIR / "manifests" / "coset_carrier_reference.json")
    source = c784.build_source(coset_manifest)
    quotient, quotient_summary = c784.central_quotient(source)
    fixture = c784.load_fixture()

    solutions = presentation_solutions(source)
    require(len(solutions) == 120, "PRESENTATION", f"expected 120 presentation solutions, got {len(solutions)}")
    source_generators = min(solutions)
    placement = placement_of_solution(quotient, source_generators)
    incidence = c784.build_incidence(quotient, placement)
    require(
        c784.classify(incidence) == "CARRIER",
        "PRESENTATION",
        "chosen presentation solution is not a carrier placement",
    )
    phi, relabel_summary = c784.relabel(quotient, placement, incidence, fixture)
    source_to_port = port_action_for_placement(quotient, placement, phi)

    source_port_rows = {source_to_port(g) for g in source.elements}
    require(
        source_port_rows == set(fixture.port_action),
        "SOURCE_PORT_ACTION",
        "SL2F5 quotient action differs from the committed port rows",
    )
    source_kernel = [g for g in source.elements if source_to_port(g) == tuple(range(12))]
    require(
        set(source_kernel) == {source.identity, source.minus_identity},
        "SOURCE_PORT_KERNEL",
        "source-to-port kernel is not exactly {+I,-I}",
    )

    # --- Recompute PORT-SPIN-LIFT, then transport its carrier coordinates ----
    matter_manifest = m314.load_json(MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json")
    upstream = m314.load_upstream(matter_manifest, MODULE_DIR)
    algebra = m314.CurrentAlgebra(upstream["current_manifest"], MODULE_DIR)
    spin = m314.spin_lift_certificate(algebra)
    raw_lifts: Mapping[PortRow, SpinMatrix] = spin["lifts"]
    require(
        set(raw_lifts) == {tuple(row) for row in algebra.plus},
        "SPIN_INTERNAL_ACTION",
        "PORT-SPIN-LIFT keys differ from the current fixture's own rotation rows",
    )

    current_phi, current_relabel_summary = current_carrier_relabel(algebra, fixture)
    lifts: dict[PortRow, SpinMatrix] = {}
    for raw_row, lift in raw_lifts.items():
        committed_row = relabel_port_row(tuple(raw_row), current_phi)
        require(
            committed_row not in lifts,
            "SPIN_CARRIER_RELABEL",
            "carrier relabelling identifies two distinct rotation rows",
        )
        lifts[committed_row] = lift
    require(
        set(lifts) == set(fixture.port_action),
        "SPIN_PORT_ACTION",
        "transported PORT-SPIN-LIFT rows differ from the committed port action",
    )
    spin_elements, spin_to_port = build_spin_group_index(lifts)

    # --- Match one exact presentation triple through the eight lift signs ----
    a, b, c = source_generators
    generator_rows = (source_to_port(a), source_to_port(b), source_to_port(c))
    require(
        tuple(c784.element_order(source.mul, source.identity, x) for x in source_generators) == (4, 6, 10),
        "PRESENTATION_ORDERS",
        "chosen SL2F5 presentation generators do not have orders (4,6,10)",
    )
    base_lifts = tuple(lifts[row] for row in generator_rows)
    choices = signed_generator_choices(base_lifts)
    require(
        len(choices) == 1,
        "GENERATOR_SIGN_CHOICE",
        f"expected one exact sign choice for the spin presentation, got {len(choices)}",
    )
    signs, spin_generators = choices[0]

    # --- Extend by words, then verify all 120^2 products ----------------------
    source_to_spin = extend_generator_map(source, source_generators, spin_generators)
    image_keys = {m314.matrix_key(matrix) for matrix in source_to_spin.values()}
    require(
        image_keys == set(spin_elements),
        "ISOMORPHISM_SURJECTIVE",
        "generator map does not hit exactly the 120 PORT-SPIN-LIFT elements",
    )

    product_checks = 0
    port_homomorphism_checks = 0
    for g in source.elements:
        for h in source.elements:
            gh = source.mul(g, h)
            require(
                spin_equal(source_to_spin[gh], m314.cmul(source_to_spin[g], source_to_spin[h])),
                "ISOMORPHISM_HOMOMORPHISM",
                "SL2F5-to-spin map fails a multiplication check",
            )
            product_checks += 1
            require(
                source_to_port(gh)
                == permutation_compose(source_to_port(g), source_to_port(h)),
                "SOURCE_PORT_HOMOMORPHISM",
                "SL2F5-to-port action fails a multiplication check",
            )
            port_homomorphism_checks += 1

    identity2 = m314.cidentity(2)
    minus_identity2 = spin_neg(identity2)
    require(
        spin_equal(source_to_spin[source.identity], identity2),
        "CENTER_MAP",
        "+I does not map to +I2",
    )
    require(
        spin_equal(source_to_spin[source.minus_identity], minus_identity2),
        "CENTER_MAP",
        "-I does not map to -I2",
    )

    # --- Verify the cover square element by element ---------------------------
    square_checks = 0
    for g in source.elements:
        spin_key = m314.matrix_key(source_to_spin[g])
        require(
            spin_key in spin_to_port,
            "SPIN_PORT_MAP",
            "isomorphism image is not a registered spin lift",
        )
        require(
            spin_to_port[spin_key] == source_to_port(g),
            "COVER_SQUARE",
            "the SL2F5-to-spin isomorphism does not commute with the port quotient",
        )
        square_checks += 1

    spin_kernel_keys = {key for key, row in spin_to_port.items() if row == tuple(range(12))}
    require(
        spin_kernel_keys == {m314.matrix_key(identity2), m314.matrix_key(minus_identity2)},
        "SPIN_PORT_KERNEL",
        "spin-to-port kernel is not exactly {+I2,-I2}",
    )

    return {
        "schema": SCHEMA,
        "construction": "explicit_generator_correspondence_with_derived_carrier_relabelling_and_exhaustive_group_check",
        "source": {
            "group": "SL(2,F5)",
            "order": len(source.elements),
            "center": ["+I", "-I"],
            "central_quotient_order": quotient_summary["order"],
            "presentation": "a^2 = b^3 = c^5 = abc = -I",
            "presentation_solution_count": len(solutions),
            "chosen_generators": [matrix_as_list(x) for x in source_generators],
            "chosen_generator_orders": [4, 6, 10],
        },
        "carrier_bridge": {
            "committed_port_rows": len(source_port_rows),
            "source_to_port_kernel": ["+I", "-I"],
            "coset_to_committed_relabel": relabel_summary,
            "current_to_committed_relabel": {
                **current_relabel_summary,
                "map": list(current_phi),
            },
            "chosen_generator_port_rows": [row_as_list(row) for row in generator_rows],
            "source_port_homomorphism_checks": port_homomorphism_checks,
        },
        "port_spin_lift": {
            "order": len(spin_elements),
            "center": ["+I2", "-I2"],
            "base_generator_lift_signs": list(signs),
            "order_profile": spin["order_profile"],
            "exact_field": "Q(sqrt(5), i)",
            "faithful_two_dimensional_complex_representation": True,
        },
        "isomorphism": {
            "source_elements_mapped": len(source_to_spin),
            "distinct_spin_images": len(image_keys),
            "multiplication_checks": product_checks,
            "center_map": "+I -> +I2; -I -> -I2",
            "bijective": True,
        },
        "commuting_cover_square": {
            "checks": square_checks,
            "source_kernel": ["+I", "-I"],
            "spin_kernel": ["+I2", "-I2"],
            "commutes_on_all_source_elements": True,
        },
        "claim_boundary": {
            "proves": (
                "an explicit exact group isomorphism from the canonical SL(2,F5) source "
                "to the executable PORT-SPIN-LIFT, after deriving the coordinate change "
                "between the #314 carrier labels and the committed Lean port labels"
            ),
            "method_not_used": [
                "same-order inference",
                "element-order-profile classification",
                "uniqueness of the non-split central extension",
            ],
            "does_not_prove": [
                "McKay or E8 correspondence",
                "selection or derivation of phi",
                "mass relation",
                "physical identification of the finite spin action",
                "source selection of the #566 charged-double-triplet current fixture",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify", "print"), nargs="?", default="verify")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    certificate = build_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8", newline="\n")
    if args.command == "print" or args.output is None:
        print(encoded, end="")


if __name__ == "__main__":
    main()
