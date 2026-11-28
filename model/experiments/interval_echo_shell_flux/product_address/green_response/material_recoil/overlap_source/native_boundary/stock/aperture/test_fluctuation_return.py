from fractions import Fraction as F

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    fluctuation_return as meter,
)


def test_homogeneous_return_and_moving_source_with_all_sites_active():
    for n in (3, 4, 9):
        for h in (F(0), F(1, 4), F(1)):
            assert meter.exact_first_step(n, h) == meter.response(n, h)
        assert (
            abs(meter.moving_response(np.ones(n)) - float(meter.response(n, F(1, 4))))
            < 1e-12
        )


def test_static_coverage_bound_on_every_small_layout():
    n = 7
    for mask in range(1 << n):
        sites = np.array([(mask >> i) & 1 for i in range(n)])
        assert (
            meter.local_response(sites / 4)
            <= float(meter.local_cap(n, int(sites.sum()))) + 1e-12
        )
