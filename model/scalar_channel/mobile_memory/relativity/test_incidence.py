"""Independent checks of the local incidence and its detector records."""

import numpy as np
from scipy import sparse

from ....contrast_memory_spherical_v1 import carrier
from .incidence import (
    CAGE_MASS,
    CHARGE_MASS,
    cell_incidence,
    dirac_operator,
    positive_embedding,
    positive_mode,
    run,
)


def test_signed_dirac_branches_and_positive_mode_against_literal_matrix():
    incidence = sparse.csr_matrix([[1.0, 0.0], [-1.0, 0.0]])
    mass = 1.25
    direct = np.array([[mass, 0, 1, 0], [0, mass, -1, 0],
                       [1, -1, -mass, 0], [0, 0, 0, -mass]])
    np.testing.assert_array_equal(dirac_operator(incidence, mass).toarray(), direct)
    energy = np.sqrt(mass ** 2 + 2)
    np.testing.assert_allclose(np.linalg.eigvalsh(direct), [-energy, -mass, mass, energy])
    for eigenvalue, seed in ((0.0, [1, 1]), (2.0, [1, -1])):
        value, state = positive_mode(incidence, mass, eigenvalue, np.asarray(seed) / np.sqrt(2))
        np.testing.assert_allclose(direct @ state, value * state, atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(state), 1, atol=1e-14)
    embedding = positive_embedding(incidence, mass)
    np.testing.assert_allclose(embedding.T @ embedding, np.eye(2), atol=1e-14)


def test_weighted_incidence_keeps_the_seam_column_and_its_orientation():
    car = carrier.federation(0)
    unweighted, seams = cell_incidence(car)
    seam = seams[0]
    blocked, blocked_seams = cell_incidence(car, {})
    weighted, weighted_seams = cell_incidence(car, {seam: 4})
    assert seams == blocked_seams == weighted_seams
    assert blocked.nnz == 0 and weighted.nnz == 2
    np.testing.assert_array_equal(weighted.toarray()[:, 0], 2 * unweighted.toarray()[:, 0])


def test_actual_carrier_shells_and_matched_observer_records():
    result = run()
    assert [row["cells"] for row in result["all_written"]] == [20, 80, 320]
    for row in result["all_written"]:
        assert row["charge"]["mass"] == CHARGE_MASS
        assert row["cage"]["mass"] == CAGE_MASS
        for key in ("charge", "cage"):
            shell = row[key]
            assert shell["components"] == 1
            assert shell["positive_band_size"] == row["cells"]
            assert shell["unpaired_negative_edge_modes"] == row["seams"] - row["cells"] + 1
            assert shell["square_residual"] < 1e-12
    for moving in result["addressed"].values():
        assert moving["written_receiver_gain"] > 0.05
        np.testing.assert_allclose(moving["written"][0], moving["blocked"][0], atol=1e-14)
        for row in moving["blocked"]:
            np.testing.assert_allclose(row, moving["blocked"][0], atol=1e-14)
    ring = result["ring"]
    assert len(ring["cell_ring"]) >= 3
    assert any(row["finite_spacing_defect"] > 0 for row in ring["records"])
    for row in ring["records"]:
        assert abs(row["group_speed"]) <= 1
