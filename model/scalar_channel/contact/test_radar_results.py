"""Independent checks of published physical statements against their archived inputs."""

from fractions import Fraction

import numpy as np

from . import certified_ranges, contact_response, intervention, material_radar, metric_certificate, radar
from .__main__ import HERE


def test_one_pulse_changes_both_public_marginals_and_not_only_the_fuel():
    spec = material_radar.protocol()
    with np.load(HERE / "radar_counts.npz", allow_pickle=False) as archive:
        result = intervention.read(spec, archive["intervention_counts"])
    assert result["menus_with_remote_and_material_response"] == 73
    assert result["seed_seam"] != result["receiver_seam"]
    first, other = result["simultaneous_absolute_lower_bounds"]
    assert first > .1859 and other > .0542


def test_exact_contact_visibility_and_early_propagation_orders():
    spec = material_radar.protocol()
    result = contact_response.run(spec)
    assert result["mark_cyclic_dimension"] == 37
    assert result["unaffected_at_every_time_dimension"] == 113
    chart, _ = contact_response.exact_chart(spec)
    base = chart.astype(object)
    marked = base.copy()
    marked[spec["mark"], spec["mark"]] += 2
    sites = spec["endpoints"]
    p = np.eye(len(base), dtype=object)[:, sites]
    q = p.copy()
    from math import factorial
    for order in range(14):
        difference = (q-p)[sites]
        for a, row in enumerate(result["station_propagation_amplitudes"]):
            for b, entry in enumerate(row):
                if order < entry["order"]:
                    assert difference[b, a] == 0
                if order == entry["order"]:
                    assert Fraction(difference[b, a], 2**order * factorial(order)) == Fraction(
                        entry["coefficient_of_minus_it_power"])
        p, q = base @ p, marked @ q


def test_certified_and_observed_geometry_both_exclude_five_percent():
    spec = material_radar.protocol()
    exact = certified_ranges.run(spec)
    with np.load(HERE / "tagged_counts.npz", allow_pickle=False) as archive:
        counts = archive["precision_counts"]
    for mark in (0, 1):
        record = exact["settings"][mark]
        scale = Fraction(1, exact["scale"])
        lo = np.array(record["distance_lower"], dtype=object) * scale
        hi = np.array(record["distance_upper"], dtype=object) * scale
        assert metric_certificate.certify((lo+hi)/2, lo, hi, numerical_guard=0, thresholds=[Fraction(1,20)]
                                          )["thresholds"][0]["excluded_subsets"] == (6, 3)[mark]
        means, hits = radar.distance_statistics(counts[mark], spec["dt"])
        # Independent of source probabilities: bound lattice conditioning using a
        # simultaneous lower confidence bound for each ordered success probability.
        n = 2**34
        successes = counts[mark, :, :, :-1, :-1].sum(axis=(-1, -2))
        assert (counts[mark].sum(axis=(-2, -1))[~np.eye(6, dtype=bool)] == n).all()
        # One-sided multiplicative Chernoff: if empirical p >= 2p then
        # Pr[X >= x] <= (e/4)^(x/2) for p <= x/(2n). Even x=1000
        # makes the union bound over 60 dyads less than 1e-80.
        assert np.min(successes[~np.eye(6,dtype=bool)]) >= 1000
        smallest = np.min(successes[~np.eye(6,dtype=bool)]) / (2*n)
        bias = 2*(spec["ticks"]-1)*spec["dt"]*(169/(2*radar.LATTICE)+1e-10)/smallest
        _, low, high = radar.confidence_envelopes(means, hits, spec["dt"], spec["ticks"],
                                                   np.full((6, 6), bias))
        verdict = metric_certificate.certify(means, low, high, thresholds=[Fraction(1,20)])
        assert verdict["thresholds"][0]["global_3d_fit_excluded"]
