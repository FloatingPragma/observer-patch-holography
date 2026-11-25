"""Round-trip threshold lengths on the retained joint source, with blank-path clocks."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.optimize import brentq

from . import flight, tetrahedron
from .. import live_modular, live_radar, live_source, live_survey
from ..source import certificate, measure as source
from ....source_implication.common_front import measure as front

HERE = Path(__file__).parent
ORDERS = (1, 2, 3, 4)
IDENTITY = np.eye(8, dtype=np.int64)
Z = np.r_[np.ones(4), -np.ones(4)]
POWERS = (Fraction(4, 5) ** 4, Fraction(5, 4) ** 4)


def hop(board, route, charge):
    addresses = live_radar.station_rows(board, range(board.n)).ravel()
    h = source.system(board, route, charge)[0][addresses][:, addresses].tocoo()
    off = (h.row // 8 != h.col // 8) & (h.data != 0)
    matrix = sparse.csr_matrix((-h.data[off], (h.row[off], h.col[off])), shape=h.shape)
    if matrix.data.min() <= 0:
        raise ValueError("shortest-path positivity needs nonnegative hopping blocks")
    return matrix


def blocks(matrix, a, b, order):
    maximum = int(abs(matrix).sum(axis=1).max())
    if maximum**order >= 2**63:
        raise OverflowError("integer power exceeds the checked sparse-product bound")
    power = matrix**order
    i = 8 * a[:, None, None] + np.arange(8)[None, :, None]
    j = 8 * b[:, None, None] + np.arange(8)[None, None, :]
    return np.asarray(
        power[
            np.broadcast_to(i, (len(a), 8, 8)).ravel(),
            np.broadcast_to(j, (len(a), 8, 8)).ravel(),
        ]
    ).reshape(-1, 8, 8)


def effect(block, order):
    # The second Gram product may exceed int64 even when the path powers fit.
    m = block.astype(object)
    u, v = m.T @ m, m @ m.T
    return u @ u + v @ v, 2 * 20 ** (4 * order)


def outer_triangle(d):
    slack, ports = max(
        (float(d[i, k] - d[i, j] - d[j, k]), [flight.STATIONS[t] for t in (i, j, k)])
        for j in range(len(d))
        for i, k in combinations((t for t in range(len(d)) if t != j), 2)
    )
    return dict(excess=slack, ports=ports, uniform_error_lower=max(0.0, slack) / 3)


def scalar_types(effects):
    """Identify local operators by exact rational equality, with a fixed denominator."""
    expected = {
        "blank": IDENTITY,
        "marked": 2048 * IDENTITY,
        "strain": np.diag([POWERS[0]] * 4 + [POWERS[1]] * 4),
    }
    result = []
    for op, den in effects:
        labels = [
            k for k, v in expected.items() if np.array_equal(op, den * v.astype(object))
        ]
        if len(labels) != 1:
            raise ValueError("unexpected reciprocal single-seam effect")
        result.append(labels[0])
    return result


def blank(board, route, order, a, b, effects, index, hops):
    # Require a unique shortest path and enough separation to miss both heads and the bridge.
    eligible = np.flatnonzero(hops[list(route)][:, a].min(axis=0) > order)
    multiplicity = np.asarray(
        (board.adjacency().astype(np.int64) ** order)[a[eligible], b[eligible]]
    ).ravel()
    k = eligible[np.flatnonzero(multiplicity == 1)[0]]
    op, den = effects[index[k]]
    np.testing.assert_array_equal(op, den * IDENTITY.astype(object))
    return [int(a[k]), int(b[k])]


def survey():
    board = source.formation.federation(1)
    route = source.formation.region(board)
    hops = shortest_path(board.adjacency(), unweighted=True, directed=False)
    negative, positive = (hop(board, route, q) for q in (-1, 1))
    ref, changes, _ = live_survey.retained(True)
    stations = np.array(flight.STATIONS)
    direct = np.zeros((8, 8))
    raw = {k: np.zeros((8, 8)) for k in ("reference", *changes)}
    panels = []
    triangle_operators = {}
    for order in ORDERS:
        a, b = np.where(np.triu(hops == order, 1))
        m = blocks(negative, a, b, order)
        np.testing.assert_array_equal(
            m.transpose(0, 2, 1), blocks(positive, a, b, order)
        )
        if not np.all(np.any(m != 0, axis=(1, 2))):
            raise ValueError("the declared hop order exceeds a leading contact order")
        menu, index = np.unique(m.reshape(-1, 64), axis=0, return_inverse=True)
        operators = [effect(x.reshape(8, 8), order) for x in menu]
        obs = np.stack([np.array(op, dtype=float) / den for op, den in operators])
        mean = np.einsum("kij,ji->k", obs, ref).real
        if np.min(mean) <= 0:
            raise ValueError(
                "the common blank-clock threshold needs a positive leading coefficient"
            )
        weights = order * mean ** (-1 / (4 * order + 2))
        d, _, matrix = front.graph(board.n, stations, weights[index], (a, b))
        if connected_components(matrix)[0] != 1:
            raise ValueError("fixed-order relay graph is disconnected")
        panel = dict(
            order=order,
            links=len(a),
            blank_standard=blank(board, route, order, a, b, operators, index, hops),
            effect_types=len(np.unique(obs.reshape(-1, 64), axis=0)),
            minimum_coefficient=float(min(mean)),
            relay_distances=d.tolist(),
            four_gram_eigenvalues=np.linalg.eigvalsh(front.gram(d[:4, :4])).tolist(),
        )
        if order == 4:
            panel["tetrahedron"] = tetrahedron.match(
                board, stations, obs, index, (a, b), ref, changes
            )
            panel["eight_gram_eigenvalues"] = np.linalg.eigvalsh(front.gram(d)).tolist()
        if order == 1:
            types = scalar_types(operators)
            classes = [types[k] for k in index]
            panel["link_classes"] = dict(Counter(classes))
            panel["exceptional_links"] = {
                tag: [
                    [int(i), int(j)]
                    for i, j, k in zip(a, b, classes, strict=True)
                    if k == tag
                ]
                for tag in ("marked", "strain")
            }
        ids = {(int(p), int(q)): k for k, (p, q) in enumerate(zip(a, b, strict=True))}
        for i, j in combinations(range(8), 2):
            key = tuple(sorted((stations[i], stations[j])))
            if key not in ids:
                continue
            k = index[ids[key]]
            if key in ((7, 42), (7, 39), (39, 42)):
                triangle_operators[key] = (order, *operators[k])
            direct[i, j] = direct[j, i] = weights[k]
            for name, rho in {"reference": ref, **changes}.items():
                raw[name][i, j] = raw[name][j, i] = float(np.trace(obs[k] @ rho).real)
        panels.append(panel)
    order_matrix = hops[stations][:, stations].astype(int)
    if np.any(direct[np.triu_indices(8, 1)] <= 0):
        raise ValueError("some survey contacts lack a blank-clock length")
    return dict(
        orders=panels,
        stations=stations.tolist(),
        hop_order=order_matrix.tolist(),
        direct_lengths=direct.tolist(),
        triangle=outer_triangle(direct),
        certified_triangle=certified_triangle(triangle_operators),
        direct_coefficients={k: m.tolist() for k, m in raw.items()},
    )


def certified_triangle(operators):
    """Rational checks on the three threshold lengths, before rounding any root."""
    import sympy as sp

    (formed, _), ferr, _ = live_source.formation_center()
    _, _, err, grams = live_source.pump_center(True)
    ref = live_source.rational_read(grams, formed, (2, 1))
    trace_error = ferr + 4 * err + 4 * err**2
    for pair, minimum, bound in (
        ((7, 42), True, sp.Rational(298, 100)),
        ((7, 39), False, sp.Rational(174, 100)),
        ((39, 42), False, sp.Integer(1)),
    ):
        order, op, den = operators[pair]
        matrix = sp.Matrix(op) / den
        if matrix == sp.eye(8):
            assert order == bound
            continue
        mean = sp.re(sp.trace(matrix * ref))
        radius = (
            max(sum(abs(val) for val in matrix.row(i)) for i in range(8)) * trace_error
        )
        target = (order / bound) ** (4 * order + 2)
        if minimum:
            assert mean + radius < target
        else:
            assert mean - radius > target
    return dict(ports=[7, 39, 42], excess_lower=0.24, uniform_error_lower=0.08)


def z_expectation(state):
    return np.dot(Z, state.diagonal().real)


@lru_cache(maxsize=1)
def z_intervals():
    """Independent rational enclosures for the implicit-function denominator and numerator."""
    import sympy as sp

    (formed, product), ferr, _ = live_source.formation_center()
    _, _, colerr, grams = live_source.pump_center(True)
    native = 4 * colerr + 4 * colerr**2
    rotation, rot_error = live_modular.rotation_centers(formed, ferr)["right"]
    directions = dict(
        decorrelated=(product - formed, 3 * ferr + 2 * native),
        right=(rotation, rot_error + (1 + rot_error) * native),
    )
    out = {}
    for name, (inp, err) in directions.items():
        state = live_source.rational_read(grams, inp, (2, 1))
        nominal = sum((1 if i < 4 else -1) * sp.re(state[i, i]) for i in range(8))
        out[name] = certificate.enclosed(nominal, err)
    return out


def matched():
    """A native one-sided partial decorrelation with the entire one-hop metric fixed."""
    intervals = z_intervals()
    dz, rot = (
        [Fraction(x) for x in intervals[tag]] for tag in ("decorrelated", "right")
    )
    if not max(rot) < 0:
        raise ValueError("right marginal rotation cannot locally keep strain fixed")
    alpha = [min(-x / y for x in dz for y in rot), max(-x / y for x in dz for y in rot)]
    mod = json.loads((HERE.parent / "live_modular.json").read_text())["modular_slopes"]
    km, kr = ([Fraction(x) for x in mod[tag]] for tag in ("decorrelated", "right"))
    slopes = [x + a * y for x in km for a in alpha for y in kr]
    if min(slopes) <= Fraction(59, 1000):
        raise ValueError("retained entropy obstruction lost its certified margin")
    board = source.formation.federation(1)
    rho, product = source.form(board, source.formation.region(board), -1)
    g = live_survey.rotations(rho)["right"]
    ref, _, _ = live_survey.retained(True)
    baseline = z_expectation(ref)

    def entropy(state):
        eigenvalues = np.linalg.eigvalsh(state)
        return -np.dot(eigenvalues, np.log(eigenvalues))

    base_s = entropy(ref)
    checks = []
    for t in (1e-4, 1e-5, 1e-6):
        trial = (1 - t) * rho + t * product

        def prepare(theta):
            # The traceless marginal-axis generator has eigenvalues ±1/2.
            unitary = np.cos(theta / 2) * np.eye(4) - 2j * np.sin(theta / 2) * g
            return live_source.propagate_heads(unitary @ trial @ unitary.conj().T)

        theta = brentq(
            lambda x: z_expectation(prepare(x)) - baseline,
            2 * t * float(alpha[0]),
            0.0,
            xtol=1e-13,
        )
        state = prepare(theta)
        checks.append(
            dict(
                partial_decorrelation=t,
                rotation=theta,
                strain_error=z_expectation(state) - baseline,
                entropy_secant=(entropy(state) - base_s) / t,
            )
        )
    return dict(
        strain_slopes=intervals,
        implicit_rotation_slope=[float(x) for x in alpha],
        retained_entropy_slope_lower=0.059,
        finite_checks=checks,
    )


def run():
    return dict(**survey(), fixed_one_hop_metric=matched())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE / "threshold.json").write_text(result)
    else:
        print(result, end="")
