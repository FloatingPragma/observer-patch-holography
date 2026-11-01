import math

import pytest

from .readout import infer_charge, infer_gap


@pytest.mark.parametrize('duration,gap', [(1., .1), (2., .5), (.1, 10.), (4., 0.)])
def test_inverse_is_independent_of_target_and_retains_error_interval(duration, gap):
    probability = math.cos(duration*gap/2)**2
    result = infer_gap(probability, duration, 1e-7)
    assert abs(result['central_gap_GeV']-gap) < 1e-12
    assert result['interval_GeV'][0] <= gap <= result['interval_GeV'][1]
    assert 'alias_assumption' in result


@pytest.mark.parametrize('p,t,e', [(1.1, 1, .1), (.2, 0, .1), (.2, 1, -1), (.2, 1, 1), (float('nan'), 1, .1)])
def test_invalid_source_inversion_rejected(p, t, e):
    with pytest.raises(ValueError):
        infer_gap(p, t, e)


def test_alias_is_a_real_premise_not_implicitly_resolved():
    p = math.cos(2.)**2  # true phase gap*t=4 lies beyond the declared first lobe.
    assert abs(infer_gap(p, 1., 0.)['central_gap_GeV']-4.) > 1


@pytest.mark.parametrize('charge', [-1., -2/3, -1/3, 0., 1/3, 2/3, 1.])
def test_charge_pointer_inverse_resolves_both_signs(charge):
    result = infer_charge((1-math.sin(charge))/2, 1e-7)
    assert abs(result['central_charge']-charge) < 1e-14
    assert result['charge_interval'][0] <= charge <= result['charge_interval'][1]
