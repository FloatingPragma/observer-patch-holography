"""Exercise the phase signatures using deposited observer records alone."""

from copy import deepcopy
import json

import pytest

from .measure import (
    HERE,
    read_hop,
    read_longitudinal,
    read_spatial,
    secants,
    synthetic_convex_control,
)


@pytest.fixture(scope="module")
def receipt():
    return json.loads((HERE / "records.json").read_text())


def test_bridge_read_uses_only_tapes(receipt):
    line = receipt["archived_arrival_line"]
    inputs = {name: line[name] for name in ("spatial", "temporal")}
    read = read_longitudinal(inputs, read_hop(receipt["blank_hop_counts"]))
    assert (
        read["spatial"]["radians_per_step"]
        == line["measured"]["spatial"]["radians_per_step"]
    )
    assert read["energy"] == line["measured"]["energy"]
    assert 16.4309 < read["energy"][0] < 16.43120649 < read["energy"][1] < 16.4315


def test_signed_parity_resolves_the_intensity_alias(receipt):
    spatial = deepcopy(receipt["archived_arrival_line"]["spatial"])
    first = read_spatial(spatial)
    for bank in spatial["contact_parity_counts"]:
        bank[0], bank[1] = bank[1], bank[0]
    flipped = read_spatial(spatial)
    assert first["cosine"] == pytest.approx([-x for x in flipped["cosine"][::-1]])
    assert first["radians_per_step"][1] < flipped["radians_per_step"][0]
    spatial["middle_count"][-1] -= 1
    with pytest.raises(ValueError, match="Incomplete spatial count bank"):
        read_spatial(spatial)


def test_secants_discriminate_the_positive_massive_control(receipt):
    hop = read_hop(receipt["blank_hop_counts"])
    for panel in receipt["native"]:
        held = [read_longitudinal(row, hop) for row in panel["readings"][-3:]]
        result = secants(held)
        assert result["positive_massive_convexity_rejected"]
        assert result["lower_slope_difference"] > 0.167
    assert not synthetic_convex_control()["positive_massive_convexity_rejected"]
