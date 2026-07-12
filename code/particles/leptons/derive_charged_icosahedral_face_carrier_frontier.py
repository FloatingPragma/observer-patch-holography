#!/usr/bin/env python3
"""Derive the icosahedral face-corner C3 carrier and its exact frontier.

The icosahedron is a declared OPH screen-microphysics object.  Its twenty
outward-oriented triangular faces form the A5/C3 orbit, and each face has a
canonical three-corner cyclic permutation representation.  This supplies a
genuine geometric C3 carrier without requiring a nonuniform scalar port load.

The theorem supplies local rank-three geometric fibers, not one canonically
selected global three-dimensional space.  It does not identify a fiber with
the physical charged-generation space, nor does C3 symmetry select the
amplitude, phase, determinant, or mass scheme of a charged response.  Those
are recorded as explicit promotion gates.
No charged reference datum is consumed.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from derive_charged_family_non_singlet_port_attachment import icosahedron_vertices


ROOT = Path(__file__).resolve().parents[2]
SCREEN = (
    ROOT / "particles" / "hierarchy" / "certificates"
    / "R_screen_sieve_icosahedral_certificate.json"
)
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_icosahedral_face_carrier_frontier.json"
)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def icosahedral_incidence(tolerance: float = 1.0e-10) -> dict[str, Any]:
    points = icosahedron_vertices()
    pair_distances = [
        float(np.linalg.norm(points[left] - points[right]))
        for left, right in itertools.combinations(range(len(points)), 2)
    ]
    edge_length = min(value for value in pair_distances if value > tolerance)
    edges = {
        (left, right)
        for left, right in itertools.combinations(range(len(points)), 2)
        if abs(float(np.linalg.norm(points[left] - points[right])) - edge_length) < tolerance
    }
    faces = [
        triple
        for triple in itertools.combinations(range(len(points)), 3)
        if all(tuple(sorted(pair)) in edges for pair in itertools.combinations(triple, 2))
    ]
    degrees = [sum(index in edge for edge in edges) for index in range(len(points))]
    edge_face_counts = {
        edge: sum(set(edge).issubset(face) for face in faces)
        for edge in edges
    }
    return {
        "vertices": len(points),
        "edges": len(edges),
        "faces": len(faces),
        "euler_characteristic": len(points) - len(edges) + len(faces),
        "vertex_degrees": degrees,
        "face_sizes": [len(face) for face in faces],
        "edge_face_counts": list(edge_face_counts.values()),
        "edge_length": edge_length,
        "face_vertex_indices": [list(face) for face in faces],
        "face_corner_flag_count": sum(len(face) for face in faces),
    }


def cyclic_shift() -> np.ndarray:
    return np.asarray(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )


def hermitian_circulant(a: float, modulus: float, phase: float) -> np.ndarray:
    shift = cyclic_shift().astype(complex)
    coefficient = modulus * np.exp(1j * phase)
    return a * np.eye(3, dtype=complex) + coefficient * shift + np.conjugate(coefficient) * (shift @ shift)


def circulant_spectrum(a: float, modulus: float, phase: float) -> list[float]:
    return sorted(float(value) for value in np.linalg.eigvalsh(hermitian_circulant(a, modulus, phase)))


def analytic_circulant_spectrum(a: float, modulus: float, phase: float) -> list[float]:
    return sorted(
        a + 2.0 * modulus * math.cos(phase + 2.0 * math.pi * index / 3.0)
        for index in range(3)
    )


def _countermodel(a: float, modulus: float, phase: float) -> dict[str, Any]:
    matrix = hermitian_circulant(a, modulus, phase)
    spectrum = circulant_spectrum(a, modulus, phase)
    shift = cyclic_shift().astype(complex)
    return {
        "a": a,
        "modulus_b": modulus,
        "phase_b": phase,
        "spectrum": spectrum,
        "positive_simple": min(spectrum) > 0.0 and len({round(value, 14) for value in spectrum}) == 3,
        "commutator_with_C3_shift_norm": float(np.linalg.norm(matrix @ shift - shift @ matrix)),
    }


def build_artifact(
    screen: dict[str, Any],
    screen_sha256: str | None = None,
) -> dict[str, Any]:
    incidence = icosahedral_incidence()
    orbit = screen.get("orbit_stabilizer", {})
    screen_premise_checks = {
        "screen_status_closed": (
            screen.get("status") == "closed_on_declared_triangulated_screen_branch"
        ),
        "screen_certificate_check_passes": (
            screen.get("checks", {}).get("icosahedral_orbit_has_twelve_vertices") is True
        ),
        "orbit_size_is_12": orbit.get("orbit_size") == 12,
        "rotation_group_order_is_60": orbit.get("group_order") == 60,
        "rotation_group_is_A5": "A5" in str(orbit.get("group", "")),
        "vertex_stabilizer_order_is_5": orbit.get("fivefold_stabilizer_order") == 5,
    }
    if not all(screen_premise_checks.values()):
        raise ValueError("face carrier requires the closed A5 twelve-vertex screen certificate")

    shift = cyclic_shift()
    c3_checks = {
        "R_cubed_is_identity": bool(np.allclose(np.linalg.matrix_power(shift, 3), np.eye(3))),
        "R_is_not_identity": not bool(np.allclose(shift, np.eye(3))),
        "A5_face_stabilizer_order": 60 // incidence["faces"],
        "face_stabilizer_matches_C3": 60 // incidence["faces"] == 3,
    }
    incidence_checks = {
        "V_E_F_are_12_30_20": (
            incidence["vertices"], incidence["edges"], incidence["faces"]
        ) == (12, 30, 20),
        "sphere_euler_characteristic_is_2": incidence["euler_characteristic"] == 2,
        "every_vertex_is_five_valent": set(incidence["vertex_degrees"]) == {5},
        "every_face_has_three_corners": set(incidence["face_sizes"]) == {3},
        "every_edge_belongs_to_two_faces": set(incidence["edge_face_counts"]) == {2},
    }
    model_a = _countermodel(1.0, 0.20, 0.10)
    model_b = _countermodel(1.0, 0.27, 0.23)
    analytic_check = bool(
        np.allclose(
            circulant_spectrum(1.0, 0.20, 0.10),
            analytic_circulant_spectrum(1.0, 0.20, 0.10),
        )
    )
    orientation_check = bool(
        np.allclose(
            circulant_spectrum(1.0, 0.20, 0.10),
            circulant_spectrum(1.0, 0.20, -0.10),
        )
    )
    selector_checks = {
        "countermodels_are_positive_simple": (
            model_a["positive_simple"] and model_b["positive_simple"]
        ),
        "countermodels_commute_with_C3": (
            model_a["commutator_with_C3_shift_norm"] < 1.0e-12
            and model_b["commutator_with_C3_shift_norm"] < 1.0e-12
        ),
        "numeric_spectrum_matches_analytic_cosine_law": analytic_check,
        "orientation_reversal_preserves_unordered_spectrum": orientation_check,
    }
    checks_pass = (
        all(screen_premise_checks.values())
        and all(incidence_checks.values())
        and all(c3_checks.values())
        and all(selector_checks.values())
    )

    return {
        "artifact": "oph_charged_icosahedral_face_carrier_frontier",
        "status": "CLOSED_GEOMETRIC_C3_FACE_CARRIER_PHYSICAL_ATTACHMENT_AND_VALUE_LAWS_OPEN",
        "source_only": True,
        "charged_reference_data_consumed": False,
        "public_charged_mass_promotion_allowed": False,
        "screen_source": "R_screen_sieve_icosahedral_certificate.json",
        "screen_source_sha256": screen_sha256,
        "screen_premise_checks": screen_premise_checks,
        "incidence": incidence,
        "incidence_checks": incidence_checks,
        "face_orbit_theorem": {
            "face_orbit": "A5/C3",
            "A5_order": 60,
            "face_count": incidence["faces"],
            "face_stabilizer_order": 3,
            "fiber": "the three corners of one outward-oriented triangular face",
            "statement": (
                "The declared icosahedral screen canonically carries a twenty-face A5/C3 "
                "orbit. The C3 stabilizer cyclically permutes the three corners of a face, "
                "giving an A5-associated bundle of local rank-three geometric fibers. Changing "
                "the face acts by A5 transport; no hidden face or port ID is required for "
                "the isomorphism class. For a fiber operator transported A5-equivariantly "
                "from one face, every face representative has the same conjugacy-invariant "
                "unordered spectrum."
            ),
            "global_flag_boundary": (
                "The 20*3=60 face-corner flags form a free transitive A5 set, hence their "
                "linearization is the 60-dimensional regular A5 representation. They do not "
                "canonically collapse to one global three-dimensional family space. Selecting "
                "or trivializing a face fiber requires an additional quotient-visible section, "
                "connection, or physical intertwiner."
            ),
            "orientation_boundary": (
                "Reversing the fiber orientation exchanges R and R^2. A Hermitian circulant "
                "is complex-conjugated and retains the same unordered real spectrum."
            ),
            "checks": c3_checks,
        },
        "circulant_commutant_theorem": {
            "statement": (
                "Every Hermitian operator on one face fiber that commutes with the regular "
                "C3 shift has the form C=aI+bR+conjugate(b)R^2. Its three eigenvalues are "
                "a+2|b|cos(delta+2*pi*k/3)."
            ),
            "real_parameter_count": 3,
            "shape_parameter_count_after_common_scale": 2,
            "failed_selector_countermodels": [model_a, model_b],
            "checks": selector_checks,
            "conclusion": (
                "The face carrier supplies the correct dimension and cyclic covariance, but "
                "C3 symmetry alone selects neither |b|/a nor delta. Distinct positive simple "
                "spectra obey exactly the same face symmetry."
            ),
        },
        "relation_to_twelve_port_W5_theorem": (
            "The face-corner bundle is a valid alternative incidence carrier. It does not "
            "contradict the W5 no-go, which applies to linear moments of scalar weights on the "
            "twelve vertices. Incidence geometry can supply a geometric C3 fiber without a "
            "nonuniform port scalar, but still requires a physical charged attachment and "
            "dynamical spectral value law."
        ),
        "conditional_physical_completion_theorem": {
            "if_supplied": [
                "an OPH quotient-visible attachment Xi from the face-corner bundle to the physical charged-generation multiplicity space",
                "a source law fixing the Hermitian-circulant amplitude and phase",
                "a normalized charged determinant-line character",
                "a declared running coordinate and RG/threshold/pole conversion",
            ],
            "then": (
                "the simple fiber spectrum gives an ID-independent unordered charged triple; "
                "a further attachment relative to neutrino or other family observables fixes "
                "the physical family lines"
            ),
        },
        "pro_completion_v1_boundary": {
            "validated_here": [
                "icosahedral V=12, E=30, F=20 incidence",
                "face stabilizer C3 and the three-corner cyclic carrier",
                "Hermitian-circulant spectral algebra",
                "face-ID independence of the unordered spectrum",
            ],
            "not_derived_by_face_incidence": [
                "delta=beta_EW/(2*N_c*N_g) or its endpoint correction",
                "the one-bit amplitude balance or its leakage correction",
                "det(Y_e)=sqrt(2)/6^14 or det(M_e)=v^3/(2*6^14)",
                "the kappa, chi_rho, and zeta_delta finite-incidence equalizers",
                "identification with the on-shell pole-mass coordinate",
            ],
        },
        "claim_boundary": (
            "The icosahedron is integral to OPH microphysics and the face-corner C3 carrier is "
            "a genuine new local geometric lemma. The numerical completion remains a conditional, "
            "retrospectively discovered conjecture until the listed source maps are proved."
        ),
        "checks_pass": checks_pass,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, default=SCREEN)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    screen_raw = args.screen.read_bytes()
    artifact = build_artifact(
        json.loads(screen_raw),
        screen_sha256=sha256(screen_raw),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": artifact["status"], "checks_pass": artifact["checks_pass"]}, indent=2))
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
