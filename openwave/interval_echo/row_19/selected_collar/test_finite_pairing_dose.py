"""Reconstruct the registered finite change directly from its native lottery."""

import json

import numpy as np
from matched_einstein_limit.source_implication.full_algebra import station
from scipy.optimize import brentq

from openwave.interval_echo.row_19.selected_collar import finite_pairing_dose as dose


def test_exact_volume_root_and_finite_modular_identity():
    data = json.loads((dose.HERE / "finite_pairing_dose.json").read_bytes())
    axes = station.phase.AXES[1:]
    heads = [(np.eye(2) + s * a) / 2 for a in axes for s in (-1, 1)]
    atoms = np.array([np.kron(left, right) for left in heads for right in heads])

    def load_recipe(name):
        recipe = data["source_recipes"][name]
        kept, mixed = (np.asarray(recipe[key]) for key in ("kept", "dephased"))
        assert min(kept.min(), mixed.min()) > 0
        even, odd = (
            np.tensordot(kept + mixed / 2, atoms, axes=1),
            np.tensordot(mixed / 2, atoms, axes=1),
        )
        return np.stack((even, odd))

    x, p, w = map(load_recipe, ("reference", "independent", "compensator"))
    effects = station.radar_menus()[6]

    def read(y):
        return station.geometry(station.pair(effects, y))

    t = data["doses"][-1]["dose"]

    def mix(u):
        return (1 - t - u) * x + t * p + u * w

    a0, v0 = read(x)
    u = brentq(lambda z: read(mix(z))[1] - v0, *data["root_fraction_interval"])
    a, v = read(mix(u))
    assert abs(v - v0) < 1e-12
    lo, hi = data["radar_area_change_same_prepared_volume"]
    assert lo <= a - a0 <= hi < -data["minimum_area_contraction_same_prepared_volume"]
    held = mix(data["selected_fraction"])
    np.testing.assert_allclose(
        dose.heads.measure.marginals(held), dose.heads.measure.marginals(x), atol=1e-14
    )
    charge = station.pair(station.log_blocks(p) - station.log_blocks(x), held - x)
    cost = station.pair(station.log_blocks(held) - station.log_blocks(x), held)
    entropy = station.entropy(held) - station.entropy(x)
    assert abs(charge - cost - entropy) < 1e-13
    for key, value in (
        ("held_connected_modular_gain", charge),
        ("held_relative_entropy_to_original", cost),
        ("held_entropy_gain", entropy),
    ):
        lo, hi = data[key]
        assert 0 < lo <= value <= hi
