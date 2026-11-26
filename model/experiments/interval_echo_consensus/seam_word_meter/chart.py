"""Subsequent preparation check for the common affine chart in the frozen protocol.

On each independently prepared data/stock pair, the protected native pulse is
a common phase times exp(i t SWAP/2) on the span of 02 and 20. At t = pi it
permutes every ordered ternary pair exactly. Conjugating one local pulse by
port permutations thus gives U_data,stock tensor identity on every other
material factor, even when those spectators are entangled. Distinct stock
ports leave the twelve coordinate preparations independent.
"""

from fractions import Fraction
from math import asin, isfinite, pi, sqrt

import networkx as nx
import numpy as np

from . import apparatus as app


def exchange_word(path):
    """Native swaps realizing an endpoint transposition with the interior restored."""
    edges = tuple(zip(path, path[1:]))
    return edges + edges[-2::-1]


def stock_layout():
    """Route disjoint data/stock inputs to one shuttered writing seam and back."""
    shared = app.layout()
    board = app.federation(0)
    head = (60, 61)
    stocks = tuple(range(48, 60))
    reserved = set(shared["calendar"])
    for station in shared["stations"]:
        reserved.update((station["head"], *station["standards"]))
    graph = nx.Graph(board.seams)
    graph.remove_nodes_from(reserved)
    if not graph.has_edge(*head):
        raise ValueError("The writing head needs a native blank contact")
    routes = []
    for data, stock in enumerate(stocks):
        paths = [
            nx.shortest_path(graph, data, head[0]),
            nx.shortest_path(graph, stock, head[1]),
        ]
        labels = list(range(board.n))
        for path in paths:
            for a, b in exchange_word(path):
                labels[a], labels[b] = labels[b], labels[a]
        expected = list(range(board.n))
        for first, second in ((data, head[0]), (stock, head[1])):
            expected[first], expected[second] = expected[second], expected[first]
        if labels != expected:
            raise ArithmeticError("Routed writing pulse would use a spectator")
        routes.append(paths)
    return {"head": list(head), "stock_ports": list(stocks), "paths_to_head": routes}


def prepare_chart(values, lower, upper):
    """Return the local instrument readings for twelve sharp 02 stocks in one chart."""
    values = tuple(values)
    if len(values) != app.N or not all(isfinite(x) for x in (*values, lower, upper)):
        raise ValueError(
            "Twelve finite coordinates and two finite chart endpoints are required"
        )
    if lower >= upper or any(x < lower or x > upper for x in values):
        raise ValueError(
            "Every coordinate must belong to the declared increasing chart"
        )
    # Exact ratios of the supplied finite values avoid overflow in wide charts.
    lo, hi = Fraction(lower), Fraction(upper)
    sharp = app.LAW.SRC.index((0, 2))
    crossed = app.LAW.SRC.index((2, 0))
    masks = [np.kron(p, np.eye(3)) for p in app.projectors()]
    ideal_support = np.zeros(9, dtype=bool)
    ideal_support[[sharp, crossed]] = True
    rows = []
    word_count = 1
    for x in values:
        fraction = (Fraction(x) - lo) / (hi - lo)
        # Near the upper endpoint use the small untransferred amplitude,
        # avoiding loss from rounding the larger probability up to one.
        t = (
            2 * asin(sqrt(float(fraction)))
            if fraction <= Fraction(1, 2)
            else pi - 2 * asin(sqrt(float(1 - fraction)))
        )
        state = app.pulse(t)[:, sharp]
        if np.linalg.norm(state[~ideal_support]) > app.TOL:
            raise ArithmeticError("Physical pulse left the sharp 02/20 pair span")
        probabilities = [float(np.linalg.norm(mask @ state) ** 2) for mask in masks]
        p = float(fraction)
        np.testing.assert_allclose(probabilities, [1 - p, 0, p], rtol=0, atol=app.TOL)
        # Normalizing the three measured effects only removes matrix roundoff.
        pn = np.asarray(probabilities) / sum(probabilities)
        stock_fraction = Fraction(float((pn[1] + 2 * pn[2]) / 2))
        measured = float((1 - stock_fraction) * lo + stock_fraction * hi)
        supp = [0] if fraction == 0 else [2] if fraction == 1 else [0, 2]
        word_count *= len(supp)
        rows.append(
            {
                "requested": float(x),
                "contact_time": float(t),
                "sharp_data_support": supp,
                "ternary_read": probabilities,
                "decoded_scalar": measured,
            }
        )
    return {
        "endpoints": [float(lower), float(upper)],
        "coordinates": rows,
        "joint_sharp_work_words": word_count,
    }


def scalar_chart_report():
    """Immutable examples of the source and two common calibrations, without sampling."""
    values = (-3, -2.75, -sqrt(2), -1, -0.5, 0, 0.375, 1, sqrt(2), 2, 3, 5)
    return {
        "provenance": "Subsequent implementation check of the preregistered common affine chart. Entries are native operator probabilities, without stochastic draws.",
        "initial_data_stock_word": [0, 2],
        "routing": stock_layout(),
        "charts": [prepare_chart(values, -3, 5), prepare_chart(values, -5, 7)],
    }
