import itertools
import json
import math
from pathlib import Path

import numpy as np
import pytest

from . import refinement as exp


def test_kernel_agrees_with_complete_small_alphabets():
    counts, ordinary = [0, 1, 3], 2
    expected = np.zeros((3, 3))
    for start, tiles in enumerate(counts):
        words = list(itertools.product(*(range(ordinary + (j < tiles)) for j in range(3))))
        for target, permitted in enumerate(counts):
            for word in words:
                end = start if ordinary in word[permitted:] else target
                expected[start, end] += 1 / (len(counts) * len(words))
    actual = exp.sector_kernel(counts, ordinary, 'reversible')
    np.testing.assert_allclose(actual, expected, atol=1e-15)
    weights = np.array([ordinary**(3-n) * (ordinary+1)**n for n in counts])
    microscopically_weighted = weights[:, None] * actual
    np.testing.assert_allclose(microscopically_weighted, microscopically_weighted.T)


@pytest.mark.parametrize('branch', exp.BRANCHES)
def test_small_ball_frame_response(branch):
    ell = 1 / 32
    flat_v, flat_a = exp.spatial_measures(0, ell)
    for u in exp.FRAME:
        stretch_squared = float(1 + u[1]**2 + u[2]**2)
        for curvature in exp.CURVATURES:
            a = exp.equal_volume_surface(branch, curvature, ell, u, flat_v)
            contraction = (3 if branch == 'maximally_symmetric' else stretch_squared) * curvature
            assert abs((flat_a - a) / (4 * math.pi / 15 * ell**4) - contraction) < 0.04
    fine = exp.product_measures(1, .25, 9, order=24)
    np.testing.assert_allclose(exp.product_measures(1, .25, 9), fine, rtol=2e-8)


def test_terminal_law_shares_one_ensemble_and_replays_the_observed_gates():
    receipt = json.loads(Path(__file__).with_name('refinement.json').read_text())
    for p in receipt['panels']:
        alphabet = p['ordinary_symbols']
        n = np.array(p['full_tiles'])
        for law in exp.LAWS:
            equilibrium = np.exp((n - n.max()) * math.log1p(1 / alphabet))
            if law == 'clear_on_contraction':
                equilibrium[:] = 1
            equilibrium /= equilibrium.sum()
            terminal = np.linalg.matrix_power(exp.sector_kernel(n, alphabet, law), 1024)
            np.testing.assert_allclose(terminal, np.broadcast_to(equilibrium, terminal.shape),
                                       rtol=1e-12)
            for obs in p['observers']:
                for launch in range(2):
                    reading = exp.read(p['ell'], p['geometry_budget'], obs['records'],
                                       p['terminal_sector_bit_counts'][law][launch])
                    assert reading == obs['interpretations'][law][launch]
        assert p['max_held_interior_z'] < 7
        if p['ell'] == exp.RADII[-1]:
            assert p['max_abs_area_over_ell4'] > 1
            assert p['max_abs_residual_over_ell4']['clear_on_contraction'] > 1
            if p['branch'] == 'maximally_symmetric':
                assert p['max_abs_residual_over_ell4']['reversible'] < 0.25
            else:
                for idx in (0, 5, 6):
                    r = p['observers'][idx]['interpretations']['reversible'][0]
                    assert np.max(np.abs(r['residuals_over_ell4'])) < 0.25
                r = p['observers'][7]['interpretations']['reversible'][0]
                assert r['residuals_over_ell4'][-1][0] == pytest.approx(-32*math.pi/15, abs=.09)
