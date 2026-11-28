from fractions import Fraction as F

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    clock_release_time as time,
)


def test_sparse_release_derivatives_against_nested_full_commutators():
    h = np.array([[9, -2, 1], [-2, 4, 3], [1, 3, 8]])
    observation = np.diag([1, -1, 2])

    def action(matrix):
        def multiply(vector):
            full = np.array([vector.get(k, 0) for k in range(3)])
            return dict(enumerate((matrix @ full).tolist()))

        return multiply

    for parity in (1, -1):
        initial = {0: 1, 1: parity}
        second, fourth = time.derivatives(initial, action(h), action(observation))
        vector = np.array([1, parity, 0])
        obj = observation
        for order in range(1, 5):
            obj = h @ obj - obj @ h
            if order % 2 == 0:
                expected = F(
                    (-1) ** (order // 2) * int(vector @ obj @ vector), 2 * 4**order
                )
                assert expected == (second if order == 2 else fourth)
