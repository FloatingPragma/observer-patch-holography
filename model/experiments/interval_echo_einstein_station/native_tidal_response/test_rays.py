"""Check the source and matched propagation on an unreduced native tensor product."""

import unittest

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.experiments.material_inertia.dark_collar import label
from .measure import ARMS, TAGS


def apply_native(state, local, axes):
    shape = state.shape
    moved = np.moveaxis(state, axes, (0, 1, 2))
    output = (local @ moved.reshape(local.shape[1], -1)).reshape(moved.shape)
    return np.moveaxis(output, (0, 1, 2), axes).reshape(shape)


class RestrictedRayTest(unittest.TestCase):
    def test_contact_without_projection(self):
        # A four-port chain parks both ends. Mobile sites each lose one native seam.
        # Keep every memory and load word in the matrix action to detect sector leakage.
        law = candidate.native_law
        exchange, contrast = law.local_term_parts()
        h = exchange + candidate.KAPPA * contrast
        empty = np.eye(1 + law.N_LAB)[0]
        degrees = np.array([1, 2, 2, 1])
        site_cost = np.zeros((3, 3, 3, 3))
        for loads in np.ndindex(site_cost.shape):
            site_cost[loads] = candidate.KAPPA * ((6 - degrees) @ np.square(loads))
        for arm, (magnitude, coin) in ARMS.items():
            memories = np.einsum("i,j,k->ijk", label(magnitude, 1, coin), empty, label(0, 1))
            for q in TAGS:
                with self.subTest(arm=arm, tag=q):
                    basis = []
                    for a, b in [(q, 0), (0, q)]:
                        loads = np.zeros((3, 3, 3, 3))
                        loads[1, a, b, 1] = 1
                        basis.append(memories[..., None, None, None, None]*loads)
                    result = site_cost*basis[0]
                    for axes in [(0, 3, 4), (1, 4, 5), (2, 5, 6)]:
                        result += apply_native(basis[0], h, axes)
                    source_gap = candidate.KAPPA * magnitude**2
                    loss = source_gap if arm == "coin_even" and q == 1 else 0.
                    scalar = candidate.KAPPA * (10 + 5*q*q) + source_gap
                    expected = (scalar - loss + .5)*basis[0] - .5*basis[1]
                    np.testing.assert_allclose(result, expected, atol=2e-13, rtol=0)


if __name__ == "__main__":
    unittest.main()
