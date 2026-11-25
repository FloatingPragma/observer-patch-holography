import json
from pathlib import Path

import numpy as np
import pytest
from sympy import Poly, symbols

from . import certificate, measure


def test_divided_difference_and_rejection_of_real_spin():
    h = np.diag([1, 1, 3]).astype(object)
    x = symbols("x")
    mats, mod = certificate.residue_matrices(
        h, Poly((x - 1) * (x - 3), x), Poly(x - 1, x)
    )
    assert mod == [-1, 1]
    np.testing.assert_array_equal(mats[0], np.diag([-2, -2, 0]))
    # Positive control: a genuine spin channel with spin-insensitive hopping
    # transmits both input amplitudes at its shared positive energy.
    genuine = np.kron(np.array([[2, -1], [-1, 2]]), np.eye(2, dtype=int))
    with pytest.raises(AssertionError):
        certificate.certify_poles(genuine, [0, 1], [[2, 3]])


def test_all_physical_settings_and_every_sharp_line():
    result = measure.run()
    old = json.loads(Path(measure.__file__).with_name("results.json").read_text())
    assert (
        len(result["preparations"])
        == len({r["parities"] for r in result["preparations"]})
        == 64
    )
    assert result["positive_roots_by_multiplicity"] == {1: 464, 2: 24}
    assert result["broadband_negative_at_16"] == 32
    assert result["benches_with_repeated_positive"] == 12
    assert set(r["full_port_commutant"] for r in result["preparations"]) == {1}
    assert sum(r["zero_minor_identities"] for r in result["preparations"]) == 2016
    obs = result["native_readout"]
    assert (
        obs["partitions"] == 966
        and 0.00297 < obs["minimum_positive_crossing"] < 0.00298
    )
    for row, stored in zip(obs["banks"], old["native_readout"]["banks"], strict=True):
        assert row["joint_counts"] == stored["joint_counts"]
        assert sum(map(sum, row["joint_counts"])) == obs["starts_per_bank"]


def test_known_spin_and_flavor_controls():
    out = measure.known_controls()
    assert out["scalar_commutant"] == 1 and out["unsplit_spin_commutant"] == 4
    assert out["split_spin_commutant"] == 2 and out["three_flavor_commutant"] == 12
    assert out["spin_order_sign"] == -1 and out["split_energy_unique"]
    assert abs(out["unsplit_energy_read"]["pt_min"] + 1 / 9) < 1e-12
    assert abs(out["split_energy_read"]["pt_min"]) < 1e-12
    for key in ["split_energy_read", "unsplit_energy_read"]:
        assert abs(out[key]["mean_arrival"] - 2 / 9) < 1e-12
