"""Independent matrix and count checks for bound three-triplet records."""

from itertools import permutations, product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1.sectors import TripleSector
from model.executable_model import carrier


HERE = Path(__file__).parent


def receipt(name):
    return json.loads((HERE / name).read_text())


def antisymmetrizer():
    basis = list(product(range(3), repeat=3))
    a = np.zeros((27, 27))
    for p in permutations(range(3)):
        inv = sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
        for i, t in enumerate(basis):
            a[basis.index(tuple(t[k] for k in p)), i] += (-1) ** inv / 6
    np.testing.assert_allclose(a @ a, a, atol=1e-14)
    assert np.linalg.matrix_rank(a) == 1
    return a


def test_faces_on_full_determinant_and_exact_escape():
    data = receipt("results.json")
    original = data["strengths"]["8"]["original"]
    p = antisymmetrizer()
    groups = original["groups"]
    c = carrier.cell()
    triples = [
        tuple(sorted(groups[i][t[i]] for i in range(3)))
        for t in product(range(3), repeat=3)
    ]
    faces = set(c.faces())
    assert np.isclose(sum(p[i, i] for i, t in enumerate(triples) if t in faces), 2 / 3)
    assert sum(data["geometry_histogram"].values()) == 2217600
    assert (
        sum(int(k) * n for k, n in data["geometry_histogram"].items())
        == 2217600 * 6 * 20 // 220
    )
    assert data["geometry_histogram"]["5"] == data["geometry_histogram"]["6"] == 0
    spare = set(original["spare_sites"])
    cut = sum((a in spare) != (b in spare) for a, b in c.seams)
    assert np.isclose(cut / 12, original["exact_spare_hopping_norm_squared"])


def test_partition_certificates_and_native_integer_transfer():
    c = carrier.cell()
    sector = TripleSector(c)
    h = 2 * sector.hamiltonian(8).toarray() + 288 * np.eye(sector.dim)
    data = receipt("results.json")["strengths"]["8"]["modular_observability"]
    assert len(data) == 437 and sum(row["orbit"] for row in data) == 42735
    assert all(row["rank_b_ba"] == row["columns"] for row in data)
    exceptional = [r for r in data if r["rank_b"] < r["columns"]]
    assert len(exceptional) == 3
    for row in exceptional + data[::19]:
        membership = {s: i for i, g in enumerate(row["groups"]) for s in g}
        inside = np.array(
            [len({membership[int(v)] for v in t}) == 3 for t in sector.tri]
        )
        b = h[np.ix_(~inside, inside)]
        a = h[np.ix_(inside, inside)]
        assert np.linalg.matrix_rank(b) == row["rank_b"]
        assert np.linalg.matrix_rank(np.vstack([b, b @ a])) == inside.sum()


def test_finite_filter_and_same_copy_detector_counts():
    j = receipt("frames.json")
    c = carrier.cell()
    sector = TripleSector(c)
    en, v = eigh(sector.hamiltonian(8).toarray() + 144 * np.eye(sector.dim))
    launch = sector.rank(*j["face_for_energy_filter"])
    spec = j["finite_filter"]
    omega = (en - j["energy"]) * spec["step"]
    gains = np.ones_like(omega, dtype=complex)
    ok = abs(omega) > 1e-10
    gains[ok] = np.expm1(-1j * spec["timing_arms"] * omega[ok]) / (
        spec["timing_arms"] * np.expm1(-1j * omega[ok])
    )
    coeff = gains * v[launch]
    admission = np.linalg.norm(coeff) ** 2
    assert abs(admission - spec["source_acceptance"]) < 1e-8
    s = spec["records"]["32"]
    wave = v @ (np.exp(-32j * (en - j["energy"])) * coeff / np.sqrt(admission))
    tri = np.array([sector.rank(*sorted(t)) for t in product(*j["groups"])]).reshape(
        tuple(map(len, j["groups"]))
    )
    frames = [np.array(x) for x in j["sharp_energy_filtered"]["frames"]]
    face = sector.agreed == 3
    placed = np.einsum(
        "ijk,ia,jb,kc->abc", (wave * face)[tri], *frames, optimize=True
    ).ravel()
    pj = np.vdot(placed, antisymmetrizer() @ placed).real
    np.testing.assert_allclose(pj, s["face_and_singlet"], atol=1e-8)
    assert abs(s["counts"][0] / j["acquisitions"] - pj * admission) < 0.002
    exact = j["sharp_energy_filtered"]["face_then_colour"]
    d = spec["trace_distance_from_ideal_at_all_free_times"]
    assert (exact["face_and_singlet"] - d) / (exact["face_and_triplets"] + d) > 0.437
    assert s["dephased_joint_singlet"] < 0.05
    assert (
        j["sharp_energy_filtered"]["held_turns"]["half_turn_all"]["dark_given_three"]
        > 0.15
    )
