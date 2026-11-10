"""Exact completeness and actual-current checks for all physical sectors."""
from collections import defaultdict

import sympy as s

from .basis import (block_labels, character_multiplicities, complete_metadata,
                    coordinate_rows, diagonal_basis, local_current_column,
                    physical_basis, raising_matrix)


def test_all_sector_exact_Gauss_dimensions_and_independent_characters():
    expected = {0: (4140, 1053, 153), -1: (2496, 718, 96), 1: (3358, 783, 103)}
    for boundary, (diagonal, physical, largest) in expected.items():
        data = complete_metadata(boundary)
        assert data['diagonal_dimension'] == diagonal
        assert data['physical_dimension'] == physical
        assert max(row['physical_dimension'] for row in data['blocks']) == largest
        assert not data['numerical_rank_tolerance_used']
        assert data['independent_SU2_character_multiplicities_match']
        counts = character_multiplicities(boundary)
        for ell, spin in block_labels(boundary):
            b, gram = physical_basis(boundary, ell, spin)
            m, _ = raising_matrix(boundary, ell, spin)
            assert b.cols == counts.get((ell, spin), 0)
            assert m * b == s.zeros(m.rows, b.cols)
            assert b.extract(coordinate_rows(boundary, ell, spin), range(b.cols)) == s.eye(b.cols)
            assert gram == b.H * b
            assert all(value.is_Rational for value in b.todok().values())


def test_actual_twelve_currents_and_retained_adjacent_weight_images():
    # The largest block of each boundary sector tests actual full Gauss
    # currents, including both lowering and raising terms and boundary sign.
    for boundary, ell, spin in ((0, 2, 0), (-1, 1, 1), (1, 3, 1)):
        keys = diagonal_basis(boundary, ell, spin)
        b, _ = physical_basis(boundary, ell, spin)
        m, row_keys = raising_matrix(boundary, ell, spin)
        assert m.rows > 0 and b.cols < len(keys)
        key_set = set(keys)
        assert all(target not in key_set for _, target in row_keys)
        for site in (0, 1):
            for current in range(12):
                residual = defaultdict(lambda: s.Integer(0))
                for (row, col), coefficient in b.todok().items():
                    for target, value in local_current_column(keys[row], site, current, boundary).items():
                        residual[target, col] += coefficient * value
                assert all(s.expand(value) == 0 for value in residual.values())
        if boundary:
            # Omitting the physical boundary charge fails the source-only
            # hypercharge constraint in every retained column.
            for key in keys:
                assert local_current_column(key, 1, 11, 0) == {key: -s.Integer(boundary)}


def test_reflected_boundary_domains_and_integer_shorthand():
    for q in (0, -1, 1):
        assert block_labels(q) == block_labels((0, q))
        for ell, spin in block_labels(q):
            assert diagonal_basis(q, ell, spin) == diagonal_basis((0, q), ell, spin)
    right = complete_metadata(-1)
    left = complete_metadata((1, 0))
    assert left['diagonal_dimension'] == 2496
    assert left['physical_dimension'] == 718
    reflected = {(4 - row['L'], -row['twoSz']): row['physical_dimension'] for row in right['blocks']}
    assert {(row['L'], row['twoSz']): row['physical_dimension'] for row in left['blocks']} == reflected
    assert left['boundary_hypercharges'] == [1, 0]
    keys = diagonal_basis((1, 0), 3, -1)
    for key in keys:
        assert local_current_column(key, 0, 11, (1, 0)) == {}
        assert local_current_column(key, 0, 11, 0) == {key: -s.Integer(1)}
