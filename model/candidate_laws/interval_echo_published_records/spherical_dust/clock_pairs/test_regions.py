from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import json
import math

import numpy as np

from .. import survey
from . import census, modular_reader, regions


def test_clock_graph_locates_an_independently_launched_geodesic():
    with localcontext() as ctx:
        ctx.prec = 65
        h, s, gamma = D("0.45"), D(1)/128, D(3)
        boost = (gamma*gamma-1).sqrt()
        direction = (D(3)/5, D(0), -D(4)/5)
        time, position = survey.ray(s, h, gamma, boost, direction)
        args = (np.array([position], dtype=float), np.zeros(3),
                np.array([3., 0., 0., math.sqrt(8)]), 1., float(h))
        inside, clock, scale = regions.material_clock(*args, float(s)+1e-10)
        excluded = regions.material_clock(*args, float(s)-1e-10)[0]
        assert inside[0] and not excluded[0]
        assert abs(clock[0]-float(survey.conformal(time, h))) < 2e-13
        assert abs(scale[0]-float(survey.scale(time, h))) < 2e-13


def test_retained_pair_has_one_positive_finite_law_at_a_deep_cutoff():
    source = regions.Band(140).source(F(19, 33))
    assert source.weight > 0
    assert source.shift(0, 0, F(1, 2**140)) < 0
    assert regions.Band(140).source(F(0)).shift(0, 0, F(1, 2**140)) == 0
    assert modular_reader.verify() == math.log(3)
    period, arms = census.verify()
    assert all(sum(row["word_counts"]) == period for row in arms.values())


def test_deep_moving_aperture_preserves_endpoint_lags_at_nonzero_epochs():
    band = regions.Band(30)
    velocity = np.array([2., 1., 1., 1.])

    def aperture(epoch, timed):
        return regions._band_read(band, np.zeros(3), velocity, 1., 0., epoch,
                                  4*band.width, F(1, 3), 9, timed)

    reference, count = aperture(0, True)
    frozen, frozen_count = aperture(0, False)
    assert count == frozen_count and count > 0 and reference > 1.5*frozen
    # An epoch shift changes either rounded lag by at most one fine tick.
    # A huge epoch additionally tests that it is never added to an offset in float.
    for epoch in (F(1, 32), -F(1, 32), D(2)**80):
        with np.errstate(over="raise", invalid="raise"):
            signal, retained = aperture(epoch, True)
        assert retained == count
        assert math.isclose(signal, reference, rel_tol=1/band.source(F(1, 3)).grid.bins)


def test_pair_batch_accepts_signed_lags_wider_than_machine_integers():
    pair = regions.Band(68).source(F(1, 3))
    half = pair.grid.bins*pair.grid.subphases//2
    assert half > 2**64
    batch = pair.entropy_from_tick_differences([-half, 0, half])
    for value, delay in zip(batch, (F(-1, 2), 0, F(1, 2)), strict=True):
        assert math.isclose(value, pair.entropy_change(0, delay, 1), rel_tol=1e-13)


def test_held_diagonal_separates_independent_controls():
    data = json.loads((regions.HERE/"regions.json").read_text())
    rows = data["held_out"]
    assert len(data["ten_frames"]) == 10
    assert max(abs(o["entropy_over_r4"]/o["bulk_over_r4"]-1) for o in data["ten_frames"]) < 0.025
    assert abs(rows["held"]["balance_over_r4"]) < 0.05
    assert rows["frozen_pair_clock"]["balance_over_r4"] < -1.1
    assert rows["blind_collar_work"]["balance_over_r4"] < -1.1
    assert rows["half_collar_force"]["balance_over_r4"] > 0.5
    assert rows["blank"]["entropy_over_r4"] == rows["blank"]["area_over_r4"] == 0
