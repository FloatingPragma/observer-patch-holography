"""Two charges and five gateway views of the same encounter-generated head source."""

import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ....beam_response import measure as camera
from model.experiments.interval_echo_einstein_station.shared_geometry_probe.flight import (
    prepared_states,
)

HERE = Path(__file__).resolve().parent
SHOTS = 2**48
DWELL = 4
STOCKS = ("independent", "half_formed", "formed", "ordinary")
ALPHA = 0.01 / (2 * 5 * 2 * 4 * 2 * 25)
BEARINGS = camera.BEARINGS


def mounts(board, route):
    """All different gateways of the two source-owner cells, in native address order."""
    hosts = {x // 12 for x in route}
    primary = tuple(sorted(route[1:3]))
    edges = [e for e in board.external if any(x // 12 in hosts for x in e)]
    assert len(edges) == 5
    masks = {}
    for edge in sorted(edges, key=lambda e: (e != primary, e)):
        rings = []
        for i, g in enumerate(edge):
            mark = (
                route[0 if i == 0 else 3]
                if edge == primary
                else min(h for h in board.neighbours[g] if h // 12 == g // 12)
            )
            rings.append(camera.ring(board, g, mark))
        masks[edge] = rings
    return masks


def kernels(board, route, stations, charge):
    indices = camera.indices(stations)
    out = []
    for odd in (False, True):
        h = camera.front.native.generator(board, route, charge, odd)
        assert (h - h.T).nnz == 0
        initial = np.zeros((h.shape[0], len(indices)))
        initial[indices] = np.eye(len(indices))
        ht = -0.5j * DWELL * h
        w = expm_multiply(ht, initial, traceA=ht.diagonal().sum())[indices]
        out.append(w.reshape(len(stations), 4, len(stations), 4).transpose(0, 2, 1, 3))
    return np.array(out)


def responses(waves, rho, tau):
    result = []
    for coin, f in ((tau, 0), ((rho + tau) / 2, 0), (rho, 0), (tau, 0.5)):
        tx = waves @ coin @ waves.conj().swapaxes(-1, -2)
        # Compare the two full coin pairs on every flight, irrespective of
        # the independent reversal outcomes. In these stocks the conditional
        # coin state on both reversal sectors is identically the same coin.
        p = np.einsum("ij,pyxji->pyx", coin.T, tx).real / 4
        result.append((1 - f) * p[0] + f * p[1])
    return np.array(result)


def jacobian_derivative(base, change):
    w = base.sum()
    dw = change.sum() / w
    x = base.sum(0) @ BEARINGS / w
    y = base.sum(1) @ BEARINGS / w
    dx = change.sum(0) @ BEARINGS / w - dw * x
    dy = change.sum(1) @ BEARINGS / w - dw * y
    return 2 * (
        BEARINGS.T @ (change - dw * base) @ BEARINGS / w
        - np.outer(dy, x)
        - np.outer(y, dx)
    )


def source_entropies(rho):
    """Modular first derivatives for the two physically fixed references, in nats."""
    tau = np.kron(*camera.front.optics.pairing.formation.marginals(rho))

    def modular(state):
        eig, basis = np.linalg.eigh(state)
        if eig.min() <= 0:
            raise ValueError("finite modular tangent needs faithful support")
        return (basis * -np.log(eig)) @ basis.conj().T

    def mean(k, state):
        return float(np.trace(k @ state).real)

    formed_k = modular(rho)
    product_k = modular(tau)
    sf = mean(formed_k, rho)
    sp = mean(product_k, tau)
    return dict(
        formed_entropy=sf,
        product_entropy=sp,
        mutual_information=sp - sf,
        at_product=mean(product_k, rho - tau),
        at_formed=mean(formed_k, rho - tau),
    )


def pencil(tables):
    p, half, full, ordinary = tables
    # Native probabilities are quadratic in the shared prepared coin,
    # and linear in its left odd-reversal fraction on the registered family.
    f = jacobian_derivative(p, 4 * half - full - 3 * p)
    o = jacobian_derivative(p, 2 * (ordinary - p))
    g = camera.image(p)[0]

    def area_slope(delta):
        return float(
            delta[0, 0] * g[1, 1]
            + delta[1, 1] * g[0, 0]
            - delta[0, 1] * g[1, 0]
            - delta[1, 0] * g[0, 1]
        )

    prop = np.linalg.solve(o, f)
    return dict(
        formed_jacobian=f.tolist(),
        ordinary_jacobian=o.tolist(),
        base_jacobian=g.tolist(),
        energy_from_area=area_slope(f) / area_slope(o),
        formed_area_slope=area_slope(f),
        ordinary_area_slope=area_slope(o),
        relative_trace=float(np.trace(prop)),
        relative_determinant=float(np.linalg.det(prop)),
        formed_det=float(np.linalg.det(f)),
        ordinary_det=float(np.linalg.det(o)),
        formed_to_independent_gain=float(
            np.linalg.det(camera.image(full)[0]) / np.linalg.det(g)
        ),
    )


def count(tables, rng):
    bank = []
    lower = []
    upper = []
    for table in tables:
        rows = np.array(
            [rng.multinomial(SHOTS, np.r_[col, 1 - sum(col)]) for col in table.T]
        )
        hits = rows[:, :-1].T
        lo = np.nan_to_num(beta.ppf(ALPHA, hits, SHOTS - hits + 1))
        hi = np.nan_to_num(beta.ppf(1 - ALPHA, hits + 1, SHOTS - hits), nan=1.0)
        assert np.all(lo <= table) and np.all(table <= hi)
        bank.append(rows.tolist())
        lower.append(lo)
        upper.append(hi)
    return bank, (np.asarray(lower), np.asarray(upper))


def run():
    board = camera.front.optics.pairing.formation.federation(1)
    route = camera.front.optics.pairing.formation.region(board)
    rho, tau, _ = prepared_states(board, route, 1)
    frame = mounts(board, route)
    stations = sorted({p for rings in frame.values() for ring in rings for p in ring})
    rng = np.random.default_rng(41826116)
    records = {}
    for charge in (-1, 1):
        panels = responses(kernels(board, route, stations, charge), rho, tau)
        series = {}
        for edge, rings in frame.items():
            i, j = ([stations.index(p) for p in s] for s in rings)
            rows = {}
            vals = []
            for name, dst, src in (("forward", j, i), ("reverse", i, j)):
                t = panels[:, dst, :][:, :, src]
                vals.append(t)
                bank, _ = count(t, rng)
                rows[name] = dict(counts=bank, population=pencil(t))
            np.testing.assert_allclose(
                vals[0], vals[1].transpose(0, 2, 1), atol=3e-14, rtol=0
            )
            series[str(edge)] = rows
        records[str(charge)] = series
    return dict(
        rings={str(k): v for k, v in frame.items()},
        stocks=list(STOCKS),
        dwell=DWELL,
        starts_per_setting=SHOTS,
        spectrum=np.linalg.eigvalsh(rho).tolist(),
        state_real=rho.real.tolist(),
        state_imag=rho.imag.tolist(),
        charges=records,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
