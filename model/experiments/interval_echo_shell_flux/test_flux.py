from fractions import Fraction
import json
from math import log, sqrt

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from model.executable_model import federation
from .measure import GRID, HERE, SHOTS, stationary


def panels():
    return json.loads((HERE / "records.json").read_text())["panels"]


def error(hits):
    # Hoeffding bounds p above by p_hat + h. Bernstein then uses that bound
    # in the per-event variance p(1-p) <= p. The paired union has size 4K.
    x = log(4 * 10000 / 0.01)
    h = sqrt(x / (2 * SHOTS))
    return sqrt(2 * (hits / SHOTS + h) * x / SHOTS) + 2 * x / (3 * SHOTS)


def signed(out, inside):
    return (out - inside) / SHOTS, error(out) + error(inside)


def directed_counts(panel, arm, capture=None):
    ns = panel["neighbor_slots"]
    flow = {}
    total = 0
    for i, (bins, slots) in enumerate(zip(arm["edge_bank"], ns, strict=True)):
        total += sum(bins)
        if slots and i != capture:
            assert len(bins) == 2 * len(slots) + 1
            for j, hits in zip(slots, bins[2::2], strict=True):
                flow[i, j] = hits
        else:
            assert len(bins) == 1
    assert total == SHOTS
    return flow


def cut(flow, region):
    outward = sum(k for (a, b), k in flow.items() if a in region and b not in region)
    inward = sum(k for (a, b), k in flow.items() if b in region and a not in region)
    return signed(outward, inward)


def test_full_native_pointer_and_both_addressed_exposures():
    exchange, contrast = native_law.local_term_parts()
    gate = exchange + KAPPA * contrast
    ix = [native_law.SRC.index(q) for q in ((2, 0), (0, 2))]
    for generator, time, expected in (
        (gate, np.pi, [[0, 1], [1, 0]]),
        (gate, np.pi / 2, [[0.5, 0.5], [0.5, 0.5]]),
        (KAPPA * contrast, np.pi / 2, np.eye(2)),
    ):
        cols = expm(-1j * time * generator)[:, ix]
        np.testing.assert_allclose(abs(cols[ix]) ** 2, expected, atol=3e-13, rtol=0)
        np.testing.assert_allclose(np.delete(cols, ix, axis=0), 0, atol=3e-13, rtol=0)
    echo = native_law.echo()
    consensus = np.tile(
        [int(x == (1, 1)) for x in native_law.SRC], native_law.DIM // native_law.N_SRC
    )
    for load in (0, 2):
        source = np.eye(native_law.DIM)[:, native_law.SRC.index((load, 0))]
        published = echo @ source
        np.testing.assert_allclose(
            np.vdot(published, consensus * published),
            int(load == 2),
            rtol=0,
            atol=3e-13,
        )
        np.testing.assert_allclose(echo @ published, source, atol=3e-13, rtol=0)


def test_unrestricted_steady_source_and_independently_counted_shells():
    # Grounded path 0--1--drain gives two resistors in series. The center
    # conducts 1/49 doubles per tick, including the separate reservoir tick.
    small = stationary([[1], [0, 2], []])
    assert Fraction(small["central"], small["denominator"]) == Fraction(1, 49)
    for row in panels():
        cert = row["certificate"]
        pi, z, g, denominator = (
            cert[k]
            for k in ("stationary_numerators", "potential", "central", "denominator")
        )
        ns = row["neighbor_slots"]
        level = row["layers"]
        arm = row["arms"]["native"]
        q = Fraction(g, denominator)
        tv = float(Fraction(arm["burn_tv_upper"]))
        assert tv < 5e-12 and arm["burn_roundoff_units"] / GRID < 1e-35
        assert min(pi) > 0 and sum(pi) == denominator
        flow = directed_counts(row, arm)
        source, drain = (arm[k][0] for k in ("source_bank", "drain_bank"))
        for observed in (source, drain):
            assert abs(observed / SHOTS - float(q)) < error(observed) + tv
        # Check the complete rational source equation independently, including
        # arbitrary positive integer levels and each unconstrained drain mass.
        for i, adj in enumerate(ns):
            if adj:
                assert sum(z[i] - z[j] for j in adj) == g * (i == 0)
                assert pi[i] == 16 * z[i]
                value, eps = cut(flow, {i})
                assert abs(value - float(q) * (i == 0)) < eps + 2 * tv
            else:
                assert z[i] == 0
                assert pi[i] == sum(z[j] for j, adj in enumerate(ns) if i in adj)
        for i, adj in enumerate(ns):
            for j in adj:
                k = flow[i, j]
                assert abs(k / SHOTS - z[i] / denominator) <= error(k) + tv
        for n in range(1, 7):
            region = {i for i, d in enumerate(level) if d <= n}
            net, eps = cut(flow, region)
            assert abs(net - source / SHOTS) < eps + error(source) + 4 * tv


def test_distance_and_flux_shape_against_independent_controls():
    for row in panels():
        board = federation(row["level"])
        linkread = row["full_contact_survey"]
        assert linkread["seams_offered"] == len(board.seams)
        assert linkread["per_edge_forward_reverse_hits"] == [32, 32]
        assert row["parked_survey"]["per_edge_forward_reverse_hits"] == [0, 0]
        seen = row["all_origins_first_two_cuts"]
        assert sum(k for _, _, k in seen) == board.n
        assert all(second != 4 * first and first > 0 for first, second, _ in seen)
        root = row["ports"].index(0)
        depth = row["layers"]
        ns = row["neighbor_slots"]
        flow = directed_counts(row, row["arms"]["native"])
        inward_witness = False
        for n in range(1, 7):
            intervals = []
            for i in range(len(ns)):
                if depth[i] != n:
                    continue
                for j in ns[i]:
                    if depth[j] == n + 1:
                        net, eps = signed(flow[i, j], flow.get((j, i), 0))
                        intervals.append((net - eps, net + eps))
            assert max(lo for lo, _ in intervals) > min(hi for _, hi in intervals)
            inward_witness |= any(hi < 0 for _, hi in intervals)
        assert inward_witness
        parked = row["arms"]["parked"]
        assert parked["source_bank"][0] == parked["drain_bank"][0] == 0
        assert not any(directed_counts(row, parked).values())
        stopped = row["arms"]["source_off"]
        assert stopped["source_bank"][0] == 0
        assert stopped["drain_bank"][0] / SHOTS < 5e-7
        assert sum(stopped["edge_bank"][root]) < SHOTS / 1000


def test_hidden_inner_collector_is_detected_by_independent_emission_meter():
    for row in panels():
        data = row["interior_collector"]
        capture = data["capture_slot"]
        assert row["layers"][capture] == 3
        flow = directed_counts(row, data, capture)
        source, inner, outer = [
            data[k][0] for k in ("source_bank", "inner_bank", "drain_bank")
        ]
        assert source > outer > 0 and inner > 0
        assert abs((source - inner - outer) / SHOTS) < sum(
            map(error, (source, inner, outer))
        )
        for shell in range(3, 7):
            net, eps = cut(flow, {i for i, d in enumerate(row["layers"]) if d <= shell})
            assert net + eps < source / SHOTS - error(source)
            assert abs(net - (source - inner) / SHOTS) < eps + error(source) + error(
                inner
            )


def test_same_source_has_steady_divergence_and_equal_outward_edges_in_the_cell():
    case = json.loads((HERE / "records.json").read_text())["antipodal"]
    assert case["survey"]["seams_offered"] == 30
    for p in case["banks"]:
        ns, d = p["neighbor_slots"], p["layers"]
        cert, arm = p["certificate"], p["records"]
        pi, z, g, denom = (
            cert[k]
            for k in ("stationary_numerators", "potential", "central", "denominator")
        )
        tv = float(Fraction(p["burn_tv_upper"]))
        assert tv < 1e-14
        assert sum(pi) == denom
        source = arm["source_bank"][0]
        assert abs(source / SHOTS - g / denom) < error(source) + tv
        flow = directed_counts(p, arm)
        for i, ne in enumerate(ns):
            if ne:
                assert sum(z[i] - z[j] for j in ne) == (g if i == 0 else 0)
        for layer, cardinality in ((1, 10), (2, 5)):
            edges = [
                (i, j)
                for i, ne in enumerate(ns)
                if d[i] == layer
                for j in ne
                if d[j] == layer + 1
            ]
            assert len(edges) == cardinality
            realized = set()
            intervals = []
            for i, j in edges:
                exact = Fraction(z[i] - z[j], denom)
                mean, eps = signed(flow[i, j], flow.get((j, i), 0))
                assert abs(mean - float(exact)) < eps + 2 * tv and mean > eps
                intervals.append((mean - eps, mean + eps))
                realized.add(exact)
            if not p["cut"]:
                assert realized == {Fraction(1, 49 * cardinality)}
                assert max(lo for lo, _ in intervals) < min(hi for _, hi in intervals)
            else:
                assert len(realized) > 1
                assert max(lo for lo, _ in intervals) > min(hi for _, hi in intervals)
