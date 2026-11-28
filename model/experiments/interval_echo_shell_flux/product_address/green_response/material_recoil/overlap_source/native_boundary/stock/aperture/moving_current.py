"""Calibrate the three resolved material currents on stationary ordinary preparations."""

from fractions import Fraction as Q

import numpy as np

from . import moving_service as lab


def run(numbers, before, after, board):
    _, half = lab.gates()
    _, _, h, h0 = lab.energy.pair.stock.kernels()
    basis = ((0, 1), (0, 2), (1, 2))
    costs = np.zeros((3, 99), int)
    for old in range(99):
        a, b = divmod(old % 9, 3)
        if a != b:
            i = basis.index(tuple(sorted((a, b))))
            costs[i, old] = sum(
                int(half[new, old]) for new in range(99) if new % 9 != old % 9
            )
    delta = numbers @ ((after - before) @ costs.T) / 8
    assert delta[0] > 0 > delta[1] and delta[2] < 0
    # Use the *sharp uniformly distributed* stationary ordinary preparations
    # for both cost and movement. Includes the native private degree energies.
    rng = np.random.default_rng(42193639)
    units = {}
    for t in (1, 2):
        native = lab.energy.pair.stock.kernels()[2][:9, :9] / 4
        ix = [t, 3 * t]
        line, modes = np.linalg.eigh(native[np.ix_(ix, ix)])
        np.testing.assert_allclose(line, [8 * t * t, 8 * t * t + 1])
        np.testing.assert_allclose(modes[0] ** 2, 0.5)
        # One of six ports on the equipped site is sampled uniformly. For a
        # private slot the energy is sharp; live slots resolve the even/odd line.
        site = rng.multinomial(2**25, np.full(lab.V, 1 / lab.V))
        live = sum(
            int(rng.binomial(int(k), int(d) / 6))
            for k, d in zip(site, board.degree, strict=True)
        )
        odd = int(rng.binomial(live, 0.5))
        price = 48 * t * t + Q(sum(int(x) for x in board.degree), 2 * lab.V)
        estimate = 48 * t * t + 6 * odd / sum(site)
        error = 6 * np.sqrt(np.log(400) / (2 * sum(site)))
        assert abs(estimate - float(price)) <= error
        # One of this many exchanges occurs per empty-board stationary sweep.
        events = Q(int(sum(numbers)), lab.V)
        units[str(t)] = dict(
            energy=str(price),
            jump_count=str(events),
            price=str(price / events),
            odd_lines=odd,
            shots=int(sum(site)),
            measured_energy_99=[estimate - error, estimate + error],
        )
    p1, p2 = [Q(units[str(t)]["price"]) for t in (1, 2)]
    cap = p1 * delta[0] + p2 * delta[1]
    assert cap < 0
    # A different experiment reads every physical bond and background degree
    # term on the same stationary configurations. Material moves conserve
    # the multiset of loads, even on the written seams.
    bond = numbers @ (after - before) @ h.diagonal() / 4
    # The one-site marginal on either end of an ordinary edge is uniform over sites.
    n2 = np.tile(np.repeat(np.arange(3) ** 2, 3), 11)
    private = 8 * int(np.sum(6 - board.degree)) * ((after - before)[0] @ n2)
    exchange = numbers @ (after - before) @ h0.diagonal() / 4
    return dict(
        jump_changes=dict(zip(("0-1", "0-2", "1-2"), map(str, delta), strict=True)),
        sharp_ordinary_standards=units,
        required_double_to_single_price=str(p2 / p1),
        largest_ratio_permitting_positive_with_free_mixed_price=str(
            -delta[0] / delta[1]
        ),
        maximum_nonnegative_mixed_price_response=str(cap),
        squared_load_current_change=str(delta[0] + 4 * delta[1] + delta[2]),
        stationary_whole_energy_change=str(bond + private),
        stationary_exchange_energy_change=str(exchange),
    )
