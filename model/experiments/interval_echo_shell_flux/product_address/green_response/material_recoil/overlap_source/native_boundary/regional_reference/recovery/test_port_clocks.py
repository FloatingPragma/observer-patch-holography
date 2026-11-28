"""Independent checks of integer rounding, signs and the experimental decoder."""

import json

import numpy as np
from scipy.linalg import expm

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.regional_reference.recovery import (
    port_enclosure as e,
)


def test_integer_products_and_phases():
    rng = np.random.default_rng(124)
    re, im = rng.integers(-e.GRID, e.GRID, size=(2, 5, 5))
    r, i = re.astype(object), im.astype(object)
    for a, b in zip(e.square(re, im), (r @ r - i @ i, r @ i + i @ r), strict=True):
        np.testing.assert_array_equal(a, (b + e.GRID // 2) // e.GRID)
    h = np.array([[3, 2, -1], [2, 7, 0], [-1, 0, 4]], dtype=np.int64)
    r, i, err = e.block_exponential(h, 4)
    expected = expm(-3j * (h - 4 * np.eye(3)))
    assert np.linalg.norm((r + 1j * i) / e.GRID - expected, ord=2) < float(err)


def test_endpoint_density_and_every_address():
    data = json.loads((e.clocks.HERE / "port_enclosure.json").read_bytes())
    law, ports, _ = e.clocks.clock_laws(8)
    weights = np.array([8357, 0, -91643, 0]) / 100000
    for k in [*ports, "matched", "unmatched"]:
        score = (
            weights @ (law[k][-1] - law["baseline"][-1]) * data["interval_denominator"]
        )
        bounds = (
            data["ports"][k.removeprefix("port_")]
            if k in ports
            else data["endpoints"][{"matched": "f", "unmatched": "g"}[k]]
        )
        assert bounds[0] <= score <= bounds[1]
    assert len(ports) == len(data["ports"]) == 24
    assert data["maximum_port_upper"] <= -24800 and data["source_lower"] >= 245300000
