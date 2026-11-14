"""Exact native family anisotropy, before a physical Yukawa attachment.

The protected-core reference and the rank-three slow response band are actual
existing public constructions. Their one-body compression is not invariant
under all carrier rotations. This supplies a concrete covariant symmetry
breaker; assigning it to a physical mass interaction is a further choice.
No measured particle mass or coupling is used.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from .. import slow_geometry
from ...native_modular_reference import local_preparation
from ...observer import native


ROOT = Path(__file__).resolve().parents[3]
FAMILY = ROOT / "model/observer_particles/results/family_naturality_operators.json"
FAMILY_RECEIPT = FAMILY.with_name("family_naturality.json")
MANIFEST = ROOT / "code/a5_closure/manifests/echosahedral_federation_reference.json"


def _transpose(matrix):
    return list(map(list, zip(*matrix, strict=True)))


def _scale(c, coefficient, matrix):
    return [[c.q5_mul(coefficient, x) for x in row] for row in matrix]


@lru_cache(maxsize=16)
def _parse_projector_entry(expression):
    # Only trusted retained algebra is parsed; its digest is checked below.
    import sympy as s
    value = s.expand(s.sympify(expression))
    b = value.coeff(s.sqrt(5))
    a = s.expand(value - b*s.sqrt(5))
    if not a.is_Rational or not b.is_Rational:
        raise ValueError("the retained projector left Q(sqrt(5))")
    return F(str(a)), F(str(b))


def _chart_map(source_edges, source_faces, target):
    """An explicit oriented incidence isomorphism, without identifying labels."""
    edges = {tuple(sorted(int(p[1:]) for p in edge)) for edge in target["edges"]}
    adjacency = [set() for _ in range(12)]
    other = [set() for _ in range(12)]
    for i, j in source_edges:
        adjacency[i].add(j); adjacency[j].add(i)
    for i, j in edges:
        other[i].add(j); other[j].add(i)

    def cyclic(face):
        face = tuple(face)
        k = face.index(min(face))
        return face[k:] + face[:k]

    faces = {cyclic(tuple(int(p[1:]) for p in face)) for face in target["oriented_faces"]}

    def extend(assignment, used):
        position = len(assignment)
        if position == 12:
            if {cyclic(tuple(assignment[p] for p in face)) for face in source_faces} == faces:
                return tuple(assignment)
            return None
        for image in range(12):
            if image in used or any((q in adjacency[position]) != (assignment[q] in other[image])
                                    for q in range(position)):
                continue
            result = extend([*assignment, image], used | {image})
            if result is not None:
                return result
        return None

    result = extend([], set())
    if result is None:
        raise ArithmeticError("native and retained family carriers have no oriented incidence map")
    return result


@lru_cache(maxsize=1)
def _data():
    c = slow_geometry.source()  # Validates the existing carrier source hashes.
    probabilities = local_preparation.reference()
    states = native.states(binary=True)
    means = tuple(sum((p*x[i] for x, p in zip(states, probabilities, strict=True)), F())
                  for i in range(12))
    if sum(probabilities, F()) != 1 or sum(means, F()) != 2:
        raise ArithmeticError("native reference normalization or total charge changed")
    projector = c.slow_band_projector_exact()
    basis = [list(row) for row in slow_geometry.raw_frame()]
    gram = c.q5_matmul(_transpose(basis), basis)
    metric = gram[0][0]
    if gram != [[metric if i == j else c.q5(0) for j in range(3)] for i in range(3)]:
        raise ArithmeticError("source coordinate frame lost its isotropic band metric")
    left = _scale(c, c.q5_inv(metric), _transpose(basis))
    if c.q5_matmul(basis, left) != [list(row) for row in projector]:
        raise ArithmeticError("coordinate frame is not the actual slow band")

    receipt = json.loads(FAMILY_RECEIPT.read_text())
    if hashlib.sha256(FAMILY.read_bytes()).hexdigest() != receipt["artifacts"][FAMILY.name]:
        raise ArithmeticError("retained family projector changed")
    packet = json.loads(FAMILY.read_text())
    retained = [[_parse_projector_entry(x) for x in row] for row in packet["source_band_projector"]]
    carrier = json.loads(MANIFEST.read_text())["carrier"]
    chart = _chart_map(c.seams(), c.oriented_faces(), carrier)
    if any(projector[i][j] != retained[chart[i]][chart[j]] for i in range(12) for j in range(12)):
        raise ArithmeticError("incidence chart map does not intertwine the two rank-three bands")

    def compress(diagonal):
        return c.q5_matmul(left, [[c.q5_mul(c.q5(diagonal[i]), x) for x in basis[i]]
                                 for i in range(12)])

    feature = compress(means)
    spectra = tuple(sorted((feature[i][i] for i in range(3)), key=c.q5_float))
    if any(feature[i][j] != c.q5(0) for i in range(3) for j in range(3) if i != j):
        raise ArithmeticError("the declared core no longer diagonalizes in the source frame")
    expected = (c.q5(F(61, 450)), c.q5(F(41, 225), F(-7, 750)),
                c.q5(F(41, 225), F(7, 750)))
    if spectra != expected or any(c.q5_sign(x) <= 0 for x in spectra):
        raise ArithmeticError("source family anisotropy spectrum changed")
    if any(c.q5_sign(c.q5_sub(b, a)) <= 0 for a, b in zip(spectra, spectra[1:])):
        raise ArithmeticError("family feature eigenvalues are not strictly ordered")
    conjugates, invariant, state_invariant = [], 0, 0
    for permutation in c.rotations():
        inverse = tuple(permutation.index(i) for i in range(12))
        rotation = c.q5_matmul(left, [basis[inverse[i]] for i in range(12)])
        moved_means = tuple(means[inverse[i]] for i in range(12))
        moved = compress(moved_means)
        conjugate = c.q5_matmul(c.q5_matmul(rotation, feature), _transpose(rotation))
        if moved != conjugate:
            raise ArithmeticError("native family compression is not covariant")
        conjugates.append(moved)
        invariant += moved == feature
        state_invariant += moved_means == means
    twirl = [[c.q5(0) for _ in range(3)] for _ in range(3)]
    for moved in conjugates:
        twirl = [[c.q5_add(twirl[i][j], c.q5_mul(c.q5(F(1, 60)), moved[i][j]))
                  for j in range(3)] for i in range(3)]
    if twirl != [[c.q5(F(1, 6)) if i == j else c.q5(0) for j in range(3)] for i in range(3)]:
        raise ArithmeticError("the full rotation twirl is not the invariant scalar")
    return c, means, projector, feature, spectra, chart, invariant, state_invariant, twirl


def report():
    """JSON-safe exact source-feature report; no physical mass is inferred."""
    c, means, projector, feature, spectra, chart, invariant, state_invariant, twirl = _data()
    paths = {Path(__file__), FAMILY, FAMILY_RECEIPT, MANIFEST,
             ROOT / "model/common_source/slow_geometry.py", slow_geometry.SPEC,
             slow_geometry.ARCHIVE / "oph_exact/carrier.py",
             slow_geometry.ARCHIVE / "oph_fpe/dynamics/self_readback_repair_closure.py"}
    paths.update(ROOT / name for name in local_preparation.SOURCE_PATHS)
    return {
        "schema": "oph.particle_identification.native_family_anisotropy.v1",
        "status": "EXACT_NATIVE_COVARIANT_NONSCALAR_FAMILY_FEATURE",
        "source": {"protected_core": list(local_preparation.CORE), "binary_states": 66,
                   "ideal_reference_contract": True, "finite_preparation_errors_propagated": False,
                   "core_and_preparation_are_existing_operational_choices": True,
                   "occupations": list(map(str, means)), "empirical_mass_or_coupling_inputs": []},
        "family_band": {"rank": 3, "native_to_retained_family_chart": list(chart),
                        "oriented_faces_and_all144_projector_entries_intertwined": True,
                        "chart_isomorphism_is_a_physical_or_unique_identification": False,
                        "physical_family_attachment_is_contract": True},
        "compression": {"definition": "F = (B^T B)^(-1) B^T diag(E_sigma[n_i]) B",
                        "native_coordinate_basis_is_metric_orthogonal": True,
                        "matrix_Qsqrt5": [[c.q5_str(x) for x in row] for row in feature],
                        "eigenvalues_Qsqrt5": list(map(c.q5_str, spectra)),
                        "eigenvalue_floats_diagnostic": list(map(c.q5_float, spectra)),
                        "distinct_and_positive_exact": True},
        "symmetry": {"proper_rotations_checked": 60, "all_covariance_identities_exact": True,
                     "feature_stabilizer_size": invariant, "reference_stabilizer_size": state_invariant,
                     "full_A5_invariance": False,
                     "full_group_twirl_Qsqrt5": [[c.q5_str(x) for x in row] for row in twirl],
                     "unbroken_family_naturality_theorem_contradicted": False,
                     "no_native_symmetry_breaker_is_available": False},
        "if_Yukawa_endomorphism_is_declared_equal_to_F": {
            "relative_singular_values_Qsqrt5": [c.q5_str(c.q5_mul(x, c.q5_inv(spectra[0]))) for x in spectra],
            "canonical_common_kinetic_and_Higgs_attachment_required": True,
            "assigning_F_to_Yukawa_is_architecture_choice": True,
            "one_common_F_selects_distinct_up_down_lepton_textures": False},
        "interpretation": {"source_covariant_symmetry_breaker_computed": True,
                           "source_selected_physical_Yukawa_intertwiner": False,
                           "charged_propagating_pole_representation_identified": False,
                           "physical_chirality_selected": False,
                           "absolute_scale_or_measured_masses_derived": False,
                           "Standard_Model_particle_identification_completed": False},
        "source_bindings": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in sorted(paths)},
    }


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, allow_nan=False))
