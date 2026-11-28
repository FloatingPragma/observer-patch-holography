"""Certify the native circular current-meter ceiling for every free time."""

from fractions import Fraction
from itertools import pairwise
import json
from math import inf, nextafter

import numpy as np
from scipy.linalg import eigh

from ..apparatus import fixture
from ..certificate import isolate
from ..measure import BASE, CY, ONE, SWAP
from ..recoil_search.guides import detour
from ..recoil_search.threshold import DELAYS, DT


HEAD_WORDS = 64
ENCLOSURE_PAD = Fraction(1, 20_000_000_000)
PI_UPPER = Fraction(22, 7)


def verify_layout(h, actions, ports, routes, marks, seams):
    """Check the cyclic preparation and cut identities on the admitted chart."""
    word = np.arange(h.shape[0])
    rotation = actions[CY]
    reflection = actions[SWAP]
    np.testing.assert_array_equal(rotation[rotation[rotation]], word)
    np.testing.assert_array_equal(reflection[reflection], word)
    np.testing.assert_array_equal(reflection[rotation[reflection]], rotation[rotation])
    for action in (rotation, reflection):
        assert (h[action][:, action] != h).nnz == 0
    assert (h != h.T.conj()).nnz == 0

    source = np.zeros(h.shape[0], dtype=np.int8)
    source[ONE] = [0, 1, -1]
    quadrature = np.zeros_like(source)
    quadrature[ONE] = [2, -1, -1]
    np.testing.assert_array_equal(
        source + source[rotation] + source[rotation[rotation]], 0
    )
    np.testing.assert_array_equal(source[reflection], -source)
    np.testing.assert_array_equal(
        source[rotation] - source[rotation[rotation]], quadrature
    )

    port_of_word = np.repeat(ports, HEAD_WORDS)
    mapped_ports = port_of_word[rotation].reshape(-1, HEAD_WORDS)
    np.testing.assert_array_equal(
        mapped_ports, np.broadcast_to(mapped_ports[:, :1], mapped_ports.shape)
    )
    port_rotation = dict(zip(ports, mapped_ports[:, 0], strict=True))
    rotated_routes = [[port_rotation[p] for p in route] for route in routes]
    assert len(routes) == 3
    assert rotated_routes in (routes[1:] + routes[:1], routes[2:] + routes[:2])
    assert all(len(route) == len(routes[0]) for route in routes)
    assert set(ports) == {port for route in routes for port in route}
    interiors = [port for route in routes for port in route[1:-1]]
    assert len(interiors) == len(set(interiors)) == len(ports) - 2
    blank_by_arm = [list(pairwise(route)) for route in routes]
    blank_edges = {edge for arm in blank_by_arm for edge in arm}
    marks_by_arm = [
        [edge for edge in marks if set(edge) <= set(route)] for route in routes
    ]
    assert len(marks) == sum(map(len, marks_by_arm)) == 6
    assert blank_edges.isdisjoint(marks)
    assert len(seams) == len(blank_edges) + len(marks)
    assert {edge for edge, *_ in seams} == blank_edges.union(marks)
    # A bit projector's only boundary is its own marked seam. Fixed orientation
    # makes its derivative one signed sum, without cancellation between rails.
    for edge, origin, target, magnitude in seams:
        np.testing.assert_array_equal(port_of_word[origin], edge[0])
        np.testing.assert_array_equal(port_of_word[target], edge[1])
        head_change = (origin % HEAD_WORDS) ^ (target % HEAD_WORDS)
        if edge in blank_edges:
            np.testing.assert_array_equal(head_change, 0)
            assert magnitude == 0.5
        else:
            bit = 1 << marks.index(edge)
            np.testing.assert_array_equal(head_change, bit)
            assert np.unique(origin & bit).size == 1
            assert magnitude == 4
    assert 2 * sum(len(origin) for _, origin, _, _ in seams) == (
        h.nnz - np.count_nonzero(h.diagonal())
    )

    levels = {port: j for route in routes for j, port in enumerate(route)}
    peak_marks_per_arm = 0
    cuts = len(routes[0]) - 1
    for cut in range(cuts):
        for blanks, written in zip(blank_by_arm, marks_by_arm, strict=True):
            assert sum(levels[u] <= cut < levels[v] for u, v in blanks) == 1
            mark_count = sum(
                min(levels[u], levels[v]) <= cut < max(levels[u], levels[v])
                for u, v in written
            )
            assert mark_count <= 1
            peak_marks_per_arm = max(peak_marks_per_arm, mark_count)
    return cuts, peak_marks_per_arm


def certify(
    h, actions, ports, routes, marks, seams, enclosure, center, target_amplitude
):
    """Bound detector record differences using the inherited spectral enclosure."""
    cuts, mark_crossings = verify_layout(h, actions, ports, routes, marks, seams)
    center = Fraction(center)
    tick = Fraction(DT)
    # Native H has half-integer entries and the row sums are exact in binary64.
    twice_entries = 2 * h.data
    np.testing.assert_array_equal(twice_entries, np.rint(twice_entries))
    assert np.max(abs(twice_entries)) * np.diff(h.indptr).max() < 2**53
    row_sum = Fraction(float(np.max(abs(h).sum(axis=1))))
    spectral_angle = tick * (row_sum + abs(center))
    assert (
        spectral_angle < 3
    )  # 3 < pi, so all energies are inside the no-alias interval.
    detuning = max(abs(Fraction(e) - center) for e in enclosure["interval"])
    detuning += ENCLOSURE_PAD
    source_floor = abs(Fraction(target_amplitude)) - Fraction(
        enclosure["eigenvector_error"]
    )
    source_floor -= ENCLOSURE_PAD
    phase = tick * (DELAYS - 1) * detuning / 2
    quadrature_floor = 1 - phase**2 / 2
    assert source_floor > 0 and quadrature_floor > 0
    target_floor = source_floor * quadrature_floor
    residual_ceiling = PI_UPPER / (DELAYS * tick * target_floor)
    blank_ceiling = 2 * (Fraction(1, 3) + mark_crossings) * residual_ceiling
    marked_ceiling = residual_ceiling / 4
    return {
        "cyclic_ordinal_cuts": cuts,
        "maximum_marked_crossings_per_arm": mark_crossings,
        "spectral_angle_upper": nextafter(float(spectral_angle), inf),
        "target_amplitude_lower": nextafter(float(target_floor), -inf),
        "normalized_residual_upper": nextafter(float(residual_ceiling), inf),
        "conditional_blank_meter_score_upper": nextafter(float(blank_ceiling), inf),
        "conditional_marked_meter_score_upper": nextafter(float(marked_ceiling), inf),
    }


def run():
    """Replay the local decoder and three inherited clock certificates."""
    # Import lazily so study.py can use certify without a module import cycle.
    from .study import THRESHOLD, bonds, verify_local_meter

    verify_local_meter()
    original = json.loads(THRESHOLD.read_text(encoding="utf-8"))
    assert original["clock_tick"] == DT and original["clock_delays"] == DELAYS
    stock = fixture()
    cases = []
    for i, threshold_case in enumerate(original["lengths"]):
        _, ports, legs, h, actions, lift, small = detour(i)
        routes = [[stock["hub"], *leg, stock["tip"]] for leg in legs]
        seams, _ = bonds(h, actions, ports, routes, stock["marked"])
        energies, vectors = eigh(small, driver="evd")
        line = int(energies.searchsorted(16.0))
        center = threshold_case["lowest_and_next_harmonic"][0]["energy"]
        assert abs(float(energies[line]) - center) < 1e-11
        # A finer dyadic grid stabilizes the strict ceiling across solvers.
        # The enclosure still checks that all sparse dot products fit exactly.
        enclosure = isolate(h, actions, lift, vectors, energies, line, dyadic_bits=45)
        beta = BASE[:, 0] @ (lift[ONE] @ vectors[:, line])
        certificate = certify(
            h,
            actions,
            ports,
            routes,
            stock["marked"],
            seams,
            enclosure,
            center,
            float(beta),
        )
        assert certificate["cyclic_ordinal_cuts"] == threshold_case["length"]
        cases.append(certificate)
    return cases


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
