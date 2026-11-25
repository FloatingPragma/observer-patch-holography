"""Angular image of a two-component launch bias on a five-port material screen."""

from pathlib import Path
import argparse
import json

import numpy as np
from scipy.sparse.linalg import expm_multiply

from ..common_front import measure as front
from ..shared_radar.aperture import measure as aperture

HERE = Path(__file__).resolve().parent
UNIT = 10**9
# Rational reading weights for the regular pentagon's first angular harmonic.
PIXELS = np.array(
    [
        (UNIT, 0),
        (309016994, 951056516),
        (-809016994, 587785252),
        (-809016994, -587785252),
        (309016994, -951056516),
    ],
    dtype=np.int64,
)
BEARINGS = PIXELS / UNIT


def ring(board, gateway, marked):
    pixels = set(p for p in board.neighbours[gateway] if p // 12 == gateway // 12)
    cycle = [marked, min(pixels.intersection(board.neighbours[marked]))]
    while len(cycle) < 5:
        following = (pixels - set(cycle)).intersection(board.neighbours[cycle[-1]])
        assert len(following) == 1
        cycle.extend(following)
    assert cycle[-1] in board.neighbours[marked]
    return cycle


def screens(board, route):
    x, y = ring(board, route[1], route[0]), ring(board, route[2], route[3])
    return x, y


def indices(ports):
    return (4 * np.array(ports)[:, None] + np.arange(4)).ravel()


def effects(board, route, rings, charge, time, strength=8):
    """Forward only. Every launch has five exclusive hits or one outside outcome."""
    sx, sy = map(indices, rings)
    grams = []
    for odd in (False, True):
        if strength == 8:
            h = front.native.generator(board, route, charge, odd) / 2
        else:
            h = front.optics.ranging.generator(board, route, charge, strength, odd)
        w = np.zeros((h.shape[0], len(sx)))
        w[sx] = np.eye(len(sx))
        if time is None:
            assert int(np.max(abs(2 * h).sum(axis=1))) ** 6 < 2**63
            h2 = (2 * h).astype(np.int64)
            assert (2 * h - h2).count_nonzero() == 0
            w = w.astype(np.int64)
            for _ in range(3):
                assert not np.any(w[sy])
                w = h2 @ w
            assert np.max(np.abs(w[sy])) < 1000
        else:
            ht = -1j * time * h
            w = expm_multiply(ht, w, traceA=ht.diagonal().sum())
        blocks = w[sy].reshape(5, 4, 5, 4).transpose(0, 2, 1, 3)
        if time is None:
            for dest in range(5):
                for src in range(5):
                    np.testing.assert_array_equal(
                        blocks[dest, src], blocks[min(dest, 1), min(src, 1)]
                    )
        grams.append(blocks.conj().transpose(0, 1, 3, 2) @ blocks)
    return np.stack(grams)


def probabilities(table, state):
    return np.einsum("byxij,bji->yx", table, state).real


def image(p):
    throughput = p.sum()
    x = p.sum(axis=0) @ BEARINGS / throughput
    y = p.sum(axis=1) @ BEARINGS / throughput
    jacobian = 2 * (BEARINGS.T @ p @ BEARINGS / throughput - np.outer(y, x))
    return jacobian, throughput / 5


def area_derivative(p, d):
    m = p.sum()
    x, y = p.sum(axis=0) @ BEARINGS / m, p.sum(axis=1) @ BEARINGS / m
    dm = d.sum() / m
    dx, dy = (
        d.sum(axis=0) @ BEARINGS / m - dm * x,
        d.sum(axis=1) @ BEARINGS / m - dm * y,
    )
    j = image(p)[0]
    dj = 2 * (
        BEARINGS.T @ (d - dm * p) @ BEARINGS / m - np.outer(dy, x) - np.outer(y, dx)
    )
    return float(
        j[1, 1] * dj[0, 0]
        + j[0, 0] * dj[1, 1]
        - j[0, 1] * dj[1, 0]
        - j[1, 0] * dj[0, 1]
    )


def read(table, state, tangents, charges):
    p = probabilities(table, state)
    j, throughput = image(p)
    slopes = {
        name: area_derivative(p, probabilities(table, t))
        for name, t in tangents.items()
    }
    return dict(
        hit_or_leading_weight=throughput,
        jacobian=j.tolist(),
        signed_area_gain=float(np.linalg.det(j)),
        singular_values=np.linalg.svd(j, compute_uv=False).tolist(),
        area_slopes=slopes,
        observed_to_predicted={
            key: slopes[key]
            * charges["ordinary"]["full"]
            / (slopes["ordinary"] * charges[key]["full"])
            for key in ("unpairing", "phase")
        }
        if abs(slopes["ordinary"]) > 1e-16
        else {},
    )


def check_leading(table, product_state):
    for src in range(5):
        for dst in range(5):
            np.testing.assert_array_equal(
                table[:, dst, src], table[:, min(1, dst), min(1, src)]
            )
    # For each launch and parity the blank/blank effect is I. The two separate
    # effects act on opposite coin factors, and their joint equals their product.
    for raw in table:
        np.testing.assert_array_equal(raw[1, 1], np.eye(4))
        np.testing.assert_array_equal(raw[0, 0], raw[1, 0] @ raw[0, 1])
        assert np.array_equal(raw[0, 1], np.kron(np.eye(2), raw[0, 1, :2, :2]))
        assert np.array_equal(raw[1, 0], np.kron(raw[1, 0, ::2, ::2], np.eye(2)))
    np.testing.assert_array_equal(table[0, 0, 1], table[1, 0, 1])
    j = image(probabilities(table, product_state))[0]
    assert np.linalg.norm(j) < 1e-12


def gateway_census():
    """Test that every first path between the two pixel rings uses their bridge."""
    result = []
    for level in (0, 1, 2):
        board = front.optics.pairing.formation.federation(level)
        g = board.adjacency().astype(np.int64)
        early = g + g @ g
        three = g**3
        count = 0
        for src, dst in board.external:
            a = [p for p in board.neighbours[src] if p // 12 == src // 12]
            b = [p for p in board.neighbours[dst] if p // 12 == dst // 12]
            assert len(a) == len(b) == 5
            assert not early[np.ix_(a, b)].nnz
            np.testing.assert_array_equal(
                three[np.ix_(a, b)].toarray(), np.ones((5, 5))
            )
            count += 1
        result.append(
            dict(cells=board.cells, bridges=count, unique_three_hop_pairs=25 * count)
        )
    return result


def run():
    board = front.optics.pairing.formation.federation(1)
    route = front.optics.pairing.formation.region(board)
    rings = screens(board, route)
    states = aperture.survey.preparations()
    chi, tangents, charges, _ = aperture.source_and_changes()
    data = {}
    for charge in (-1, 1):
        table = effects(board, route, rings, charge, None)
        check_leading(table, states["pairing"])
        r = read(table, chi, tangents, charges)
        # Here the zero determinant and its derivative follow algebraically.
        r["leading_distinct_input_columns"] = 2
        weights = probabilities(table, chi)
        left, right, joint = map(float, (weights[1, 0], weights[0, 1], weights[0, 0]))
        covariance = joint - left * right
        gain = 50 * covariance / weights.sum() ** 2
        np.testing.assert_allclose(
            r["jacobian"], np.diag([gain, 0]), atol=2e-15, rtol=0
        )
        r["covariance"] = covariance
        r["angular_gain_from_covariance"] = gain
        r["exact_signed_area"] = 0
        r["exact_area_slopes"] = {key: 0 for key in tangents}
        del r["observed_to_predicted"]
        values = {}
        for time in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
            values[str(time)] = read(
                effects(board, route, rings, charge, time), chi, tangents, charges
            )
        control = read(
            effects(board, route, rings, charge, 4.0, strength=0),
            chi,
            tangents,
            charges,
        )
        data[str(charge)] = dict(
            leading=r, finite=values, zero_contrast_at_four=control
        )
    return dict(
        gateway_census=gateway_census(),
        source_screen=rings[0],
        receiver_screen=rings[1],
        bearing_numerators=PIXELS.tolist(),
        bearing_denominator=UNIT,
        full_modular_slopes={k: v["full"] for k, v in charges.items()},
        charges=data,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(report)
    else:
        print(report, end="")
