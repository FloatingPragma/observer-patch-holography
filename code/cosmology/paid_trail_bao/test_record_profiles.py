from fractions import Fraction as Q
import json

import numpy as np

from . import neighbor_witness as neighbor
from . import witness_access as access

rp = neighbor.rp


def test_physical_recovery_matches_exact_counts_and_dilute_sign():
    bank = json.loads((rp.source.HERE / "record_density.json").read_bytes())
    x = np.array([0, 1 / 16, 1 / 4, 1 / 2, 3 / 4, 1])
    laws = rp.source_laws(x, bank)
    for menu, curve in laws.items():
        assert curve[0] == curve[-1] == 0
        for i, label in enumerate(("1/16", "1/4", "1/2"), 1):
            np.testing.assert_allclose(
                curve[i],
                bank["cases"][label][menu]["receiver_nats_per_attempt"],
                rtol=2e-13,
            )
        np.testing.assert_allclose(curve[2], curve[4], atol=1e-15)
    assert laws["none"][1] > 0 > laws["third_head"][1]
    threshold = rp.source.Q(neighbor.dilute_limit(bank)["threshold"])
    assert threshold == Q(570303, 900970)
    p = np.array([1e-7])
    sparse = rp.source_laws(p, bank)
    flags = sum(
        w * h / i["words"]
        for w, h, i in zip(
            rp.source.input_weights(p),
            bank["third_magnitude_two_words"],
            bank["classes"],
            strict=True,
        )
    )
    homogeneous = sparse["none"] + flags * (sparse["third_head"] - sparse["none"])
    assert sparse["none"][0] > 0 > homogeneous[0]


def test_joint_orbit_feasibility_requires_one_gain_and_one_mass_path():
    radius = [1.0, 2.0, 3.0]
    data = [Q(i) for i in (2, 3, 6)]
    same = [Q(i) for i in (0, 1, 1)]
    assert not access.constant_access(radius, data, data, same, same)["possible"]
    a = [Q(i) for i in (0, 1, 1)]
    b = [Q(i) for i in (0, 2, 4)]
    result = access.variable_access(radius, data, data, a, b)
    assert result["possible"] and list(map(Q, result["feasible_msun"])) == data
    assert not access.variable_access(radius, data, data, a, a)["possible"]
    # The only additional shell is negative. A decreasing mass remains admissible.
    result = access.variable_access(
        [1, 2],
        list(map(Q, [4, 2])),
        list(map(Q, [4, 2])),
        list(map(Q, [0, -1])),
        list(map(Q, [0, -1])),
    )
    assert result["possible"] and Q(result["lower"]) == Q(result["upper"]) == 2


def test_held_velocity_cannot_change_own_gain_or_repair_an_undefined_orbit():
    row = dict(observed_kms=[2.0, 3.0], baryons_kms=[1.0, 1.0], error_kms=[0.1, 0.1])
    inputs = dict(held=row, other=row)
    fields = dict(held=np.array([-3.0, 1.0]), other=np.array([3.0, 8.0]))

    def held_result():
        return rp.scale_scores(inputs, fields)["galaxies"]["held"]

    first = held_result()
    inputs["held"] = dict(row, observed_kms=[20.0, 30.0])
    second = held_result()
    assert first["gain_msun_per_nat_pc2"] == second["gain_msun_per_nat_pc2"] == 1.0
    assert first["undefined"] == second["undefined"] == 1
