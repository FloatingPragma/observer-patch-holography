from fractions import Fraction as Q

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    closed_energy as energy,
    closed_pairing as pair,
)


def test_closed_marginals_determine_one_port_and_full_support_density():
    _, words, cells, _, _, _, _, _ = pair.setup()
    target = energy.least_change(cells)
    middle = words[:, 1]
    # Here the equality is necessary word by word, including before a free
    # evolution or any choice of weighting within the closed sharp support.
    np.testing.assert_array_equal(cells == 0, middle != 1)
    assert sum((target - 1) * (middle - 1) ** 2) / len(words) == Q(1, 14)
    assert np.all((target - 1) @ words[:, :4] == 0)
    board, homes, seams, terms, private, diagonal = energy.returned(words)
    np.testing.assert_array_equal(diagonal, terms.sum(0) + private.sum(1))
    prices = energy.all_supported_prices(words, cells, board, homes, seams, diagonal)
    actual = sum((target - 1) * diagonal) / (4 * len(words))
    assert Q(prices["coherent_lower"]) <= Q(prices["sharp_min"]) <= actual < 0
    assert 0 < Q(prices["sharp_max"]) <= Q(prices["coherent_upper"])


def test_same_spectator_predicts_double_cell_even_with_unknown_initial_filling():
    for words in pair.source.classes()[1]:
        mags = pair.stock.MAG[words[:, 4:]]
        both = (mags[:, :2] == 2).all(1)
        observed = (words[:, :4] == 1).all(1) & (mags[:, 2] != 2)
        np.testing.assert_array_equal(both, observed)
