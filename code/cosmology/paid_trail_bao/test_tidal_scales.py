from fractions import Fraction as Q
from itertools import product

import numpy as np

from . import tidal_scale_kernel as kernel
from . import tidal_neighbors as peers
from . import tidal_mixture as mixtures
from . import tidal_rotation as direct
from . import tidal_scales as scales


def test_all_peer_bins_and_no_peer_match_exhaustive_independent_boards():
    w = np.array([[0.0, 0.125, 0.5], [0.375, 0.0, 0.25]])
    rings = np.array([[0, 0, 1], [0, 1, 1]])
    bins = peers.search_bins(w.shape)
    tally = np.zeros((len(bins) + 1, 2))
    coords = np.argwhere(w)
    for bits in product((False, True), repeat=len(coords)):
        bit = np.array(bits)
        p = w[tuple(coords.T)]
        probability = np.prod(np.where(bit, p, 1 - p))
        positions = coords[bit]
        for y, x in positions:
            d2 = ((positions - [y, x]) ** 2).sum(1)
            first = min(d2[d2 > 0], default=bins[-1] + 1)
            tally[np.searchsorted(bins, first), rings[y, x]] += probability
    np.testing.assert_allclose(peers.peer_mass(w, rings, 2, 1, bins), tally, atol=2e-16)
    isolated = peers.peer_mass(
        np.array([[0.3]]), np.array([[0]]), 1, 2, peers.search_bins((1, 1))
    )
    np.testing.assert_allclose(isolated, [[0.6]])


def test_observer_polar_integral_and_tangent_limit():
    r = np.array([0.4, 0.9, 1.2, 2.8])
    k = [2, 3, 4, 5]
    K = [0, 0.6, 0.25, 0]
    for hop in (0.16, 0.45, 1.2):
        area = direct.annular_transfer(r, hop, k, K, 128, 2048) / hop**3
        result = kernel.annular_force(r, hop, k, K, 18)
        np.testing.assert_allclose(area, result, rtol=1e-4, atol=6e-5)
    # At the boundary of a filled disk the surviving half-plane
    # gives angular moment 2 in the small-scale limit.
    moment = kernel.integral(k, list(map(lambda x: Q(str(x)), K)), 2, Q(2), Q(5))
    result = kernel.annular_force(np.array([1.0]), 1e-8, k, K, 18)[0, 0]
    np.testing.assert_allclose(result, 2 * float(moment) / np.pi, rtol=1e-12)


def test_signed_enclosure_for_all_scales_on_nonmonotone_disk():
    r = np.array([1.0, 2.0, 4.0])
    cum = [[4.0, 4.5, 45.0]]
    dlo, dhi, unit = kernel.density_ticks(r, cum)
    L = kernel.length_edges(4)
    k, K = [1, 2, 3, 4], [Q(0), Q(1, 8), Q(1), Q(0)]
    bins = kernel.kernel_bins(k, K)
    weights = np.diff(np.r_[0.0, cum[0]])
    box = [kernel.disk_boxes(int(x), r.astype(int), dlo, dhi, L, 1, bins) for x in r]
    length = L / kernel.LENGTH

    def conv(value):
        return scales.physical_bounds(value, bins[-1], unit)

    for hop in [1e-8, 1e-4, 0.3334, 0.875, 1, 1.125, 1.37, 2.0, 2.76, 7.99]:
        result = kernel.annular_force(r, hop, k, K, 18) @ weights
        i = np.clip(np.searchsorted(length, hop) - 1, 0, len(L) - 2)
        for j, (lo, hi, flo, fhi) in enumerate(box):
            assert conv(lo[0, i])[0] - 1e-11 <= result[j] <= conv(hi[0, i])[1] + 1e-11
            assert conv(flo[0])[0] - 1e-11 <= result[j] <= conv(fhi[0])[1] + 1e-11
    # An exterior annulus alone can point outward.
    assert kernel.annular_force(r, 0.4, k, K)[0, -1] < 0


def test_same_scale_mass_cannot_be_spent_twice_at_different_radii():
    # One annulus can do one unit at either orbit, and may split its mass.
    # Each threshold alone is feasible, but three quarters at both is impossible.
    a = kernel.ANGLE
    boxes = np.array([[[a, 0]], [[0, a]]], np.int64)
    gain = [Q("3.141592653589793")]
    required = [Q(3, 4)] * 2
    upper = [Q(2)] * 2
    report = mixtures.support_witness(boxes, boxes, gain, required, upper)
    assert report["certified_exclusion"] and len(report["lower_weights"]) == 2
    assert Q(report["required_sum"]) > Q(report["outward_capacity"])
    feasible = mixtures.support_witness(boxes, boxes, gain, [Q(2, 5)] * 2, upper)
    assert not feasible["certified_exclusion"]


def test_hierarchy_moments_reproduce_direct_scale_integrals():
    profiles, mu = scales.source_data()
    knots, values, cal = profiles["source"]
    info = scales.hierarchy(profiles, mu)["source"]
    u_min, u_max = info["support"]
    # The point-source field is r^3 * K(r/L)/L^3 in enclosed mass units.
    leg, ws = np.polynomial.legendre.leggauss(24)
    for r in (0.73, 3.0, 11.0):
        log_mass, lin_mass = 0.0, 0.0
        for n in range(u_min, u_max):
            lo, hi = r / (n + 1), r / n
            L = lo + (leg + 1) * (hi - lo) / 2
            K = np.interp(r / L, knots, list(map(float, values)))
            a = (hi - lo) / 2 * ws * r**3 * K / L**3
            log_mass += sum(a / L)
            lin_mass += sum(a)
        D = cal[6]
        np.testing.assert_allclose(
            log_mass * float(mu / D),
            float(Q(info["log_density_enclosed_per_unit_log_mass"]["6"])),
            rtol=3e-14,
        )
        np.testing.assert_allclose(
            lin_mass * float(mu / D) / r,
            float(
                Q(
                    info["linear_density_enclosed_per_radius_per_unit_mass_per_length"][
                        "6"
                    ]
                )
            ),
            rtol=3e-14,
        )
