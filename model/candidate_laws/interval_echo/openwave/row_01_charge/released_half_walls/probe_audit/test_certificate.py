"""Check all deposited word and operator witnesses, including an in-domain control."""

import json
from math import comb

import numpy as np

from . import certificate


def test_all_censuses_and_operators_replay():
    expected = json.loads((certificate.HERE / "census.json").read_text())
    assert json.loads(json.dumps(certificate.main())) == expected
    for check in expected:
        n = check["size"]
        assert check["two_wall_words"] == 1
        assert check["same_inventory_words"] == comb(n - 4, n // 2 - 3)
        assert check["opposite_inventory_words"] == comb(n - 4, n // 2 - 5)
        swaps = check["witnesses"][0]["swap_bonds_and_resulting_words"]
        assert len(swaps) == n - 6
        assert all(len(certificate.wall_pairs(w, n)) == 4 for _, _, w in swaps)


def test_probe_resolves_an_exchange_through_its_aperture():
    n, level = 14, 0
    ring, (left, right) = certificate.apertures(level, n)
    settings = (left[0], right[0])
    start = certificate.source_census(n, settings, 1)["witnesses"][0]["source_word"]
    # The endpoint's outer neighbour belongs to the alternating domain.
    j, k = left[0], (left[0] - 1) % n
    assert ((start >> j) ^ (start >> k)) & 1
    altered = start ^ (1 << j) ^ (1 << k)

    def read(word):
        loads = {p: (word >> i) & 1 for i, p in enumerate(ring)}
        return certificate.rail_matrix(
            certificate.federation(level), left[1:], loads, 1
        )

    assert start.bit_count() == altered.bit_count()
    np.testing.assert_allclose(
        read(altered) - read(start), np.diag([8, -8]), atol=2e-13
    )


def test_independent_probabilities_and_changed_counts():
    from .check_counts import check

    raw = json.loads(
        (certificate.HERE.parent / "material_probe/counts.json").read_text()
    )
    result = check(raw)
    assert result["covered"] == result["predictions"] == 7424
    raw["measurements"][0][1] = [raw["shots"], 0, 0, 0, 0, 0, 0, 0]
    mutated = check(raw)
    assert mutated["covered"] < mutated["predictions"]
