from fractions import Fraction
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import expm

from model.executable_model import carrier
from . import refinement


def test_refined_cage_detectors_keep_all_outcomes():
    replay = refinement.run()
    assert replay["cells_with_exact_cages_by_level"] == [20, 80, 320]
    for step, cells, seams in zip(replay["steps"], (20, 80), (30, 120), strict=True):
        assert step["coarse_cells"] == cells
        assert step["fine_cells"] == 4 * cells
        assert step["coarse_seams"] == seams
        assert step["fine_crossings"] == 2 * seams
        for fiber in step["seam_fibers"]:
            assert len(fiber["fine"]) == 2
            assert all(
                [step["parent_of_fine_cell"][child] for child in edge] == fiber["coarse"]
                for edge in fiber["fine"]
            )
        assert step["no_proper_parent_local_invariant_subspace"]
        for key in (
            "fine_cell_records", "fine_parent_records", "lifted_coarse_cell_records",
            "energy_matched_coarse_records", "zero_link_cell_records", "zero_link_parent_records",
        ):
            assert np.all(np.asarray(step[key]) >= 0)
            np.testing.assert_allclose(sum(step[key]), 1, atol=1e-12, rtol=0)
        assert step["zero_link_parent_records"][step["receiver_parent"]] == 0
        assert step["receiver_difference"] >= float(Fraction(step["receiver_difference_lower_bound"]))

    # Read the same finite record by dense exponentiation with a separate
    # cell builder and constructor, without using the certificate matrices.
    first = replay["steps"][0]
    fine = carrier.federation(1)
    matrix = np.zeros((fine.cells, fine.cells))
    for a, b in fine.external:
        left, right = a // 12, b // 12
        matrix[left, left] += 1
        matrix[right, right] += 1
        matrix[left, right] -= 1
        matrix[right, left] -= 1
    state = np.zeros(fine.cells)
    state[:4] = 0.5
    direct = expm(-1j * float(Fraction(first["wait"])) * matrix) @ state
    np.testing.assert_allclose(abs(direct) ** 2, first["fine_cell_records"], atol=1e-12, rtol=0)


def test_refinement_rejects_rewired_seams_even_at_fixed_degrees():
    coarse = carrier.federation(0)
    fine = carrier.federation(1)
    coarse_edges = refinement._cell_seams(coarse)
    fine_edges = refinement._cell_seams(fine)
    crossings = [edge for edge in fine_edges if edge[0] // 4 != edge[1] // 4]
    for (a, b), (c, d) in combinations(crossings, 2):
        parents = {a // 4, b // 4, c // 4, d // 4}
        switched = tuple(sorted((a, c))), tuple(sorted((b, d)))
        if (len(parents) == 4 and all(edge not in fine_edges for edge in switched)
                and tuple(sorted((a // 4, c // 4))) not in coarse_edges):
            mutated = [edge for edge in fine_edges if edge not in ((a, b), (c, d))]
            mutated.extend(switched)
            with pytest.raises(ValueError, match="exactly two actual fine crossings"):
                refinement.certify_step(coarse.cells, coarse_edges, fine.cells, mutated)
            return
    pytest.fail("the refinement fixture lacks a degree-preserving wrong crossing")
