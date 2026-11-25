import json
from math import log, sqrt
from itertools import product
import numpy as np
from scipy.linalg import expm
from ..measure import gates
from .measure import HERE, apply_offer, menu, native_law, KAPPA


def read_summary():
    j = json.loads((HERE / "records.json").read_text())
    thermometers = {}
    for row in j["calibrations"]:
        count = sum(n for _, n in row["words"])
        value = sum(word.bit_count() * n for word, n in row["words"]) / (
            count * len(row["owners"])
        )
        assert count + sum(row["first_refusals"]) == j["calibration_shots"]
        thermometers[row["body"], row["exposure"]] = value, count
    stats = []
    for row in j["contacts"]:
        size = len(row["owners"])
        before = np.zeros(2**size, dtype=int)
        after = np.zeros_like(before)
        pairs = np.zeros((2**size, 2**size), int)
        for tag, old, new, count in row["records"]:
            before[old] += count
            after[new] += count
            pairs[old, new] += count
            assert old.bit_count() == new.bit_count()
        n = before.sum()
        assert n + sum(row["first_refusals"]) == j["shots"]
        delta = (after - before) / n
        flow = delta[1::2].sum()
        j_upper = None
        if row["second"] == "B":
            j_upper = sum(delta[i] for i in range(len(delta)) if i & 6 == 6)
        stats.append((row, n, max(abs(delta)), flow, j_upper))
    return j, thermometers, stats


def test_complete_histories_and_independent_heat_prediction():
    j, t, stats = read_summary()

    def allowance(n, width=1):
        return width * sqrt(log(2 * 2048 / 0.01) / (2 * n))

    for (body, b), (p, n) in t.items():
        assert abs(p - 1 / (1 + 2**b)) < allowance(n)
    for r, n, largest, heat, both in stats:
        p, n1 = t[r["first"], 1]
        q, n2 = t[r["second"], r["exposure"] if r["second"] == "C" else 1]
        tags = r["tags"]
        contact = sum(tag[0] == "swap" and tag[1] == 0 for tag in tags) / len(tags)
        if r["control"] == "exchange_off":
            contact = 0
        assert abs(heat - contact * (q - p)) < allowance(n, 2) + contact * (
            allowance(n1) + allowance(n2)
        )
        if r["copy"] and r["control"] == "native":
            assert abs(both + 1 / 36) < allowance(n, 2)
            assert abs(both) > allowance(n, 2)
            assert abs(heat) < allowance(n, 2)
        elif r["exposure"] == 1 or r["control"] != "native":
            assert largest < allowance(n, 2)
        else:
            assert abs(heat) > allowance(n, 2)
    for r in j["coherence"]:
        assert (
            sum(untouched + switched for untouched, _, switched in r["banks"])
            == j["shots"]
        )
        changed = sum(switched for _, _, switched in r["banks"]) / j["shots"]
        assert abs(changed - (0.5 if r["mode"] == "native" else 0)) < allowance(
            j["shots"]
        )


def test_contact_as_an_operator_in_the_material_energy_basis():
    exchange, cl = native_law.local_term_parts()
    raw = expm(-1j * np.pi * (exchange + KAPPA * cl))[:9, :9]
    # Two transversal native contacts swap the complete two-dimer state.
    op = (
        np.kron(raw, raw)
        .reshape((3,) * 8)
        .transpose(0, 2, 1, 3, 4, 6, 5, 7)
        .reshape(81, 81)
    )
    assert tuple(native_law.SRC) == tuple(product(range(3), repeat=2))
    loc = np.zeros((9, 2))
    loc[native_law.SRC.index((2, 0))] = [1, 1]
    loc[native_law.SRC.index((0, 2))] = [1, -1]
    basis = np.kron(loc, loc) / 2
    expected = basis[:, [0, 2, 1, 3]]
    np.testing.assert_allclose(op @ basis, expected, atol=3e-13)
    h = gates("native")[0]
    np.testing.assert_allclose(np.linalg.eigvalsh(h[:2, :2]), [32, 33], atol=3e-13)


def test_background_dependent_odds_expose_a_correlated_body():
    # Multiplying these exact integer weights by a common factor is immaterial.
    # Body B is half independent bits, half the copy of a single sign.
    weights = np.array([[10, 2], [2, 4]])
    assert weights[0].sum() == 2 * weights[1].sum()
    assert weights[0, 1] * weights[1, 0] != weights[0, 0] * weights[1, 1]
    tags = menu(1, 2)
    assert len(tags) == 8 and len(set(tags)) == 6
    assert apply_offer((0, 1, 1), ("swap", 0, 1), "native") == (1, 0, 1)


def test_joint_fixed_diagonal_orbits_and_the_local_phase_read():
    import networkx as nx

    for left, right in ((1, 1), (1, 2), (2, 3), (3, 3)):
        tags = menu(left, right)
        words = list(product(range(2), repeat=left + right))
        graph = nx.Graph()
        graph.add_nodes_from(words)
        for w in words:
            for tag in tags:
                graph.add_edge(w, apply_offer(w, tag, "native"))
        sectors = list(nx.connected_components(graph))
        assert len(sectors) == left + right + 1
        assert all(len({sum(w) for w in sector}) == 1 for sector in sectors)
        # Each off-diagonal matrix unit anticommutes with at least one offered local Z.
        phase_sites = [tag[1] for tag in tags if tag[0] == "phase"]
        assert all(
            any(v[i] != w[i] for i in phase_sites)
            for n, v in enumerate(words)
            for w in words[:n]
        )
