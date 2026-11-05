"""Hostile finite/algebraic controls; these do not test all continuous maps."""

import copy
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import equivariant_degree_crosscheck as producer
import verify_equivariant_degree_crosscheck as independent


@pytest.fixture(scope="module")
def packet():
    return json.loads(producer.RECEIPT.read_text())


def test_independent_e8_character_and_basepoint_controls(packet):
    result = independent.verify(packet)
    assert result["minimum_homogeneous_degrees"] == {"Q": 1, "Q_prime": 7}
    assert result["checked_through"] == 94


def test_common_orientation_reversal_and_marked_generator_square(packet):
    for residues in packet["necessary_residues_mod_30"].values():
        assert sorted((-x) % 30 for x in residues) == residues
    for prime in (2, 3, 5):
        for a in range(1, prime):
            for b in range(1, prime):
                for k in range(1, prime):
                    assert (b*pow(a, -1, prime)) % prime == ((k*b)*pow(k*a, -1, prime)) % prime
    assert 2*pow(1, -1, 5) % 5 == 2
    assert 4*pow(2, -1, 5) % 5 == 2  # coherent generator square
    assert 4*pow(1, -1, 5) % 5 != 2  # target-only square changes the typing


@pytest.mark.parametrize("mutation", ["wrong_Qprime", "wrong_class", "omit_lift", "multiplicity",
    "wrong_numerator", "wrong_graph", "common_zero", "ambient_sign", "universal_promotion",
    "physical_promotion", "false_provenance", "implementation_hash"])
def test_rejects_algebraic_and_scope_forgery(packet, mutation):
    bad = copy.deepcopy(packet)
    if mutation == "wrong_Qprime":
        bad["characters"][7] = bad["characters"][1]
    elif mutation == "wrong_class":
        bad["class_Q_traces"][5] = [-1, 1]
    elif mutation == "omit_lift":
        bad["class_sizes"][1] = 0
    elif mutation == "multiplicity":
        bad["multiplicities"]["Q_prime"][3] = 1
    elif mutation == "wrong_numerator":
        bad["hilbert_numerator_exponents"]["Q_prime"] = [1, 11, 19, 29]
    elif mutation == "wrong_graph":
        bad["mckay_adjacency"][6][7] = 0
    elif mutation == "common_zero":
        bad["hms"]["denominator"] = bad["hms"]["numerator"]
    elif mutation == "ambient_sign":
        bad["minimum_absolute_degree_from_congruences"]["twin"] = 1
    elif mutation == "universal_promotion":
        bad["continuous_theorem"] = "CONTINUOUS_TWIN_DEGREE_FLOOR_PROVED"
    elif mutation == "physical_promotion":
        bad["physical_selection"] = True
    elif mutation == "false_provenance":
        bad["bernhard_original_proof_and_scripts"] = "CHECKED"
    else:
        bad["implementation_sha256"]["equivariant_degree_crosscheck.py"] = "0"*64
    with pytest.raises(ValueError):
        independent.verify(bad)


@pytest.fixture(scope="module")
def marked_geometry():
    # This verifier reconstructs coordinates in Q(phi) and uses Cramer's rule.
    # It imports neither the frame producer nor the source-audit producer.
    import ast
    import verify_galois_port_frame_independent as frame
    import verify_port_support_source_audit as incidence

    reference = frame.load(frame.REFERENCE)
    frame.verify(reference)
    body = (frame.ROOT / frame.SOURCE).read_text().split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    faces = [list(f) for f in ast.literal_eval(body[body.index("["):].strip())]
    group = incidence.search(faces)
    assert len(group) == 60
    vectors = {name: frame.matrix(reference["family"][name]["vectors"], 12, 3)
               for name in ("plus", "minus")}
    return frame, incidence, faces, group, vectors


def test_same_abstract_normalizer_swaps_realized_poles_in_both_frames(marked_geometry):
    frame, incidence, faces, group, vectors = marked_geometry
    identity = tuple(range(12))
    edges = sorted({tuple(sorted((a, b))) for f in faces for a in f for b in f if a != b})

    def order(g):
        value = identity
        for n in range(1, 61):
            value = incidence.cp(g, value)
            if value == identity:
                return n
        raise AssertionError("non-permutation action")

    for g in group:
        prime = order(g)
        if prime not in (2, 3, 5):
            continue
        powers, value = set(), identity
        for _ in range(prime):
            powers.add(value)
            value = incidence.cp(g, value)
        normalizer = [n for n in group if incidence.cp(incidence.cp(n, g), incidence.inv(n)) in powers]
        assert len(normalizer) == {2: 4, 3: 6, 5: 10}[prime]
        objects = [(i,) for i in range(12)] if prime == 5 else edges if prime == 2 else faces
        fixed = {tuple(sorted(o)) for o in objects if {g[i] for i in o} == set(o)}
        assert len(fixed) == 2
        pole = min(fixed)
        swaps = [n for n in normalizer if n not in powers]
        for n in swaps:
            assert order(n) == 2
            if prime == 2:
                assert incidence.cp(n, g) == incidence.cp(g, n)
            else:
                assert incidence.cp(incidence.cp(n, g), incidence.inv(n)) == incidence.inv(g)
            for v in vectors.values():
                axis = [frame.total(v[i][j] for i in pole) for j in range(3)]
                moved = [frame.total(v[n[i]][j] for i in pole) for j in range(3)]
                assert axis != [frame.ZERO]*3
                assert moved == [frame.neg(x) for x in axis]
                if prime == 2:
                    # The other half-turn's edge-axis is perpendicular, in
                    # both faithful representations; inversion alone is empty.
                    other_edge = next(e for e in edges if {n[i] for i in e} == set(e))
                    other_axis = [frame.total(v[i][j] for i in other_edge) for j in range(3)]
                    assert frame.dot(axis, other_axis) == frame.ZERO


def test_live_marked_rotation_traces_and_radial_orientation(marked_geometry):
    frame, incidence, faces, group, vectors = marked_geometry
    marked = (0, 2, 4, 1, 6, 8, 3, 5, 10, 7, 9, 11)
    assert marked in group
    results = {}
    for name, v in vectors.items():
        basis = v[:3]  # Stored as columns for Cramer's rule.
        denominator = frame.determinant(basis)

        def solve(columns, rhs):
            det = frame.determinant(columns)
            return [frame.div(frame.determinant([rhs if j == k else col
                    for j, col in enumerate(columns)]), det) for k in range(3)]

        unit_coefficients = [solve(basis, [frame.scalar(int(i == j)) for i in range(3)])
                             for j in range(3)]
        for g in group:
            columns = [[frame.total(frame.mul(c[k], v[g[k]][j]) for k in range(3))
                        for j in range(3)] for c in unit_coefficients]
            assert frame.determinant(columns) == frame.ONE
            assert all(frame.dot(columns[i], columns[j]) == frame.scalar(int(i == j))
                       for i in range(3) for j in range(3))
            assert all([frame.total(frame.mul(columns[k][j], v[i][k]) for k in range(3))
                        for j in range(3)] == v[g[i]] for i in range(12))
            if g == marked:
                trace = frame.total(columns[i][i] for i in range(3))
                assert trace == (frame.PHI if name == "plus" else frame.sigma(frame.PHI))

        assert all(frame.determinant([v[i] for i in f]) == denominator for f in faces)
        assert frame.sign(denominator) == -1
        if name == "plus":
            # Each listed plane supports the convex hull, with the origin on
            # its interior side. The source realization really is a sphere.
            for face in faces:
                def side(w):
                    return frame.determinant([[frame.sub(v[i][j], w[j]) for j in range(3)]
                                              for i in face])
                assert all(frame.sign(side(v[i])) == frame.sign(denominator)
                           for i in range(12) if i not in face)

        result = []
        for direction in ([1, 2, 4], [-3, 5, 2], [7, -2, 3]):
            hits, signs = [], []
            for index, face in enumerate(faces):
                columns = [v[i] for i in face]
                c = solve(columns, list(map(frame.scalar, direction)))
                assert all(x != frame.ZERO for x in c)
                if all(frame.sign(x) > 0 for x in c):
                    hits.append(index)
                    signs.append(-frame.sign(frame.determinant(columns)))
            assert signs == [1]*(1 if name == "plus" else 7)
            result.append(hits)
        results[name] = result
    assert results["plus"][0] == [2]
    assert results["minus"][0] == [0, 3, 6, 8, 16, 17, 18]
