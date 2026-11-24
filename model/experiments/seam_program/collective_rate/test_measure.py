import json
from pathlib import Path

import numpy as np
import pytest

from model.experiments.seam_program.collective_rate.measure import (
    clock_records, gibbs_read, modular_transfer, nested_regions, spectrum,
)


@pytest.fixture(scope="module")
def active():
    return spectrum(8)


def test_population_records_and_short_time_laws(active):
    fillings = (.25, .5, .75)
    t = 1/512
    observations = clock_records(active, fillings, (t,))
    recorded = json.loads(Path(__file__).with_name("results.json").read_text())["clocks"]
    for f, row, archived in zip(fillings, observations, recorded[:3], strict=True):
        for key, value in row.items():
            assert value == pytest.approx(archived[key], abs=1e-12)
        assert row["exchange_quadrature"] / t == pytest.approx(1, abs=.0002)
        assert row["symmetric_quadrature"] / t == pytest.approx(50 - 68*f, abs=.04)
        assert row["adjacent_density_difference"] / t**2 == pytest.approx(.25, abs=.00004)


def test_pair_modular_transfer(active):
    transfer = modular_transfer()
    assert transfer["native_energy_slope"] == pytest.approx(1)
    assert all(value == 0 for value in transfer["modular_slopes"].values())
    control, warm = (gibbs_read(active, b) for b in (0, 1/16))
    assert control["swap_slope"] == pytest.approx(0, abs=1e-14)
    assert warm["filling"] == pytest.approx(.5)
    assert warm["minimum_regional_probability"] > .11
    # End exchange and charge conservation diagonalize the reduced Gibbs
    # density in the two parity words. Recover both slopes from populations.
    odd = warm["minimum_regional_probability"]
    even = odd + 2 * warm["exchange_coherence"]
    aligned = (1 - even - odd) / 2
    assert warm["swap_slope"] == pytest.approx(np.log(even / odd))
    assert warm["aligned_to_symmetric_slope"] == pytest.approx(np.log(aligned / even))
    assert warm["swap_slope"] == pytest.approx(.0593563906835)
    assert warm["aligned_to_symmetric_slope"] / 8 == pytest.approx(.148251700121)


def test_common_neighbors_change_the_two_modular_slopes(active):
    beta = 1/16384
    row = gibbs_read(active, beta)
    difference = row["aligned_to_symmetric_slope"] / 8 - row["swap_slope"]
    assert difference / beta**2 == pytest.approx(8.5, abs=.007)


def test_thermal_edge_inside_nested_graph_regions(active):
    beta = 1/16
    regions = nested_regions(active, beta)
    pair = gibbs_read(active, beta)
    assert [len(region["ports"]) for region in regions] == [2, 4, 8]
    assert regions[0]["swap_per_energy"] == pytest.approx(pair["swap_slope"], abs=1e-12)
    assert regions[0]["aligned_per_energy"] == pytest.approx(pair["aligned_to_symmetric_slope"] / 8,
                                                            abs=1e-12)
    gaps = [region["aligned_per_energy"] - region["swap_per_energy"] for region in regions]
    assert gaps == pytest.approx([.0888953094, .0126003337, -.00003281282], abs=1e-10)
    assert max(abs(regions[-1][key] - beta) / beta
               for key in ("swap_per_energy", "aligned_per_energy")) < .00054
    for region in regions:
        # Global charge conservation and reality leave u(XX+YY)+v ZZ after
        # removing one-port and identity terms. Compare its analytic spectrum.
        u = -(region["swap_per_energy"] - beta) / 4
        v = u - 4 * (region["aligned_per_energy"] - beta)
        norm = max(abs(v), abs(-v - 2*u), abs(-v + 2*u))
        assert region["fixed_marginal_operator_norm"] == pytest.approx(norm, abs=1e-12)
    assert regions[-1]["fixed_marginal_operator_norm"] == pytest.approx(.00002667039894, abs=1e-12)
