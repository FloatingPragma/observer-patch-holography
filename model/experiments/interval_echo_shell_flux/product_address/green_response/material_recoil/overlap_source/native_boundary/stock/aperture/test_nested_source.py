from itertools import permutations

import numpy as np
import pytest

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    nested_source as nested,
)


def test_selected_minor_determinant_keeps_the_reported_pivot_order():
    rows = np.array([[0, 0, 0], [0, 0, 11], [0, 5, 1], [13, 1, 2], [0, 1, 7]])
    ids, det = nested.minor(rows)
    a = rows[ids]
    exact = 0
    for p in permutations(range(3)):
        inversions = sum(p[i] > p[j] for i in range(3) for j in range(i + 1, 3))
        exact += (-1) ** inversions * np.prod(a[np.arange(3), p])
    assert int(exact) % nested.PRIME == det != 0
    with pytest.raises(ValueError):
        nested.minor(rows[:, [0, 0, 2]])
