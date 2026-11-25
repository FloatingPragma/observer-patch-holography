"""Regional modular energy and its boundary flow on the central native bridge."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv
from scipy import sparse

from ...full_algebra import certificate as c
from ...full_algebra.reuse import invisible
from .. import bridge
from ..connection import equal_radar
from ..traction.certificate import four_tick_columns
from ..translation.certificate import compress

HERE = Path(__file__).resolve().parent
NAMES = ("occupation", "coherence", "quadrature", "imbalance")


def source():
    info = equal_radar.tangent_source()[0]
    words = info["fixed_state_centers"]
    frame = [c.prior.sym_matrix(x) for x in c.iso.basis(False)]
    centers = [
        sp.eye(4) / 4
        + sum(
            (sp.Rational(n, words["denominator"]) * a for n, a in zip(row, frame)),
            sp.zeros(4),
        )
        for row in words["even_then_odd"]
    ]
    return (
        centers,
        sp.Rational(info["conditional_HS_error_upper"]),
        iv.mpf(info["matched_odd_weight"]),
    )


def mean(blocks, ref, error):
    centers, radius, p = ref
    base = []
    dev = []
    for h, m in zip(blocks, centers):
        cc = h - h.trace() * sp.eye(4) / 4
        base.append(c.ri(c.inner(h, m)))
        dev.append(iv.sqrt(c.ri(c.inner(cc, cc))) * c.ri(radius))
    return (
        base[0]
        + p * (base[1] - base[0])
        + iv.mpf([-1, 1]) * (max(e.b for e in dev) + c.ri(error))
    )


def modular(data):
    d, x, y, z = data
    u = iv.sqrt(x * x + y * y + z * z)
    high, low = (d + u) / 2, (d - u) / 2
    assert low.a > 0 and high.b < 1 and d.b < 1 and u.a > 0
    beta = -iv.ln(high / low) / (2 * u)
    a = iv.ln(1 - d) - iv.ln(high * low) / 2
    # K shifted to vanish on the outside outcome, in canonical load Pauli basis.
    return [a, beta * x, beta * y, beta * z], [1 - d, high, low]


def meters(h, endpoints):
    group = np.arange(8)
    outer = np.setdiff1d(np.arange(h.shape[0]), group)
    local = sparse.lil_matrix(h.shape, dtype=np.int64)
    local[:8, :8] = np.kron([[17, -1], [-1, 17]], np.eye(4, dtype=int))
    local = local.tocsr()
    onsite = sparse.lil_matrix(h.shape, dtype=np.int64)
    onsite[:8, :8] = h[:8, :8] - local[:8, :8]
    onsite = onsite.tocsr()
    assert onsite[:4, 4:8].count_nonzero() == 0 and onsite[4:8, :4].count_nonzero() == 0
    hop = h.copy().tolil()
    hop[:8, :8] = 0
    hop[np.ix_(outer, outer)] = 0
    hop = hop.tocsr()
    assert len(endpoints) == 12 and h.shape == (48, 48)
    return dict(hopping=hop, head_onsite=onsite, inside=local)


def local_effects(board, route, rings, charge, weights, *, include_flow=True):
    sites = list(route[1:3])
    outer = sorted(
        {
            int(y)
            for p, t in board.seams
            for x, y in ((p, t), (t, p))
            if x in sites and y not in sites
        }
    )
    points = sites + outer
    idx = bridge.beam.indices(points)
    launches = bridge.beam.indices(sum(rings, []))
    matrices = (
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1], [1, 0]]),
        np.diag([1, -1]),
    )
    labels = ("reference", "hopping", "head_onsite", "inside", "total_rate")
    readings = {name: [[] for _ in matrices] for name in labels}
    for odd in (False, True):
        h = bridge.native.generator(board, route, charge, odd)
        w = four_tick_columns(h, idx, launches, weights)
        cut = h[idx, :][:, idx]
        terms = meters(cut, points)
        for j, small in enumerate(matrices):
            mat = sparse.lil_matrix(cut.shape, dtype=int)
            mat[:8, :8] = np.kron(small, np.eye(4, dtype=int)).astype(int)
            mat = mat.tocsr()
            imaginary = j == 2
            readings["reference"][j].append(compress(w, mat, -sp.I if imaginary else 1))
            if not include_flow:
                continue
            full = cut @ mat - mat @ cut
            pieces = {key: obj @ mat - mat @ obj for key, obj in terms.items()}
            assert (sum(pieces.values()) - full).count_nonzero() == 0
            for k, result in [*pieces.items(), ("total_rate", full)]:
                assert max(abs(result).sum(axis=1).A1) < 256
                # h=2H. For O=i n, i[H,O]=-[h,n]/2. Otherwise multiply [h,n]/2 by i.
                readings[k][j].append(
                    compress(w, result, -2 if imaginary else -2 * sp.I)
                )
    return readings, points


def run():
    iv.dps = 32
    board = bridge.native.iso.pairing.formation.federation(1)
    route = bridge.native.iso.pairing.formation.region(board)
    rings = bridge.beam.screens(board, route)
    weights = bridge.native.prior.weights()[7]
    ref = source()
    template = equal_radar.json.loads(
        (equal_radar.HERE / "force_entropy.json").read_text()
    )
    witness = sp.Matrix(
        [sp.Rational(x) for x in template["direction_center_even_then_odd"]]
    )
    corr = [c.prior.sym_matrix(x) for x in c.iso.basis(True)]
    # The previous certified correction plus a binary64 coordinate conversion bound.
    witness_error = sp.Rational(".000000568")
    results = {}
    for q in (-1, 1):
        table, points = local_effects(board, route, rings, q, weights)
        values = {
            item: [mean(blocks, ref, sp.Rational("3e-15")) for blocks in byread]
            for item, byread in table.items()
        }
        k, spectrum = modular(values["reference"])
        stress_weights = [k[0] + 17 * k[1], -2 * k[1], 2 * k[2], k[3]]
        rates = {
            label: sum(a * b for a, b in zip(k, v))
            for label, v in values.items()
            if label != "reference"
        }
        assert (
            0 in rates["inside"]
            and 0 in rates["total_rate"] - rates["hopping"] - rates["head_onsite"]
        )
        m = sp.Matrix(
            [
                [c.inner(block, b) for block in blocks for b in corr]
                for blocks in table["reference"]
            ]
        )
        tangents = []
        for i in range(4):
            nom = c.ri((m[i, :] * witness)[0])
            err = invisible.frobenius(m[i, :]) * c.ri(witness_error) + c.ri("1e-14")
            tangents.append(nom + iv.mpf([-1, 1]) * err)
        entropy_slope = sum(v * w for v, w in zip(tangents, k))
        results[str(q)] = dict(
            populations=dict(zip(NAMES, map(str, values["reference"]))),
            eigenvalues=list(map(str, spectrum)),
            modular_coefficients_canonical=list(map(str, k)),
            work_basis_coefficients=dict(
                zip(
                    ("occupation", "strain_work", "current_1_to_38", "imbalance"),
                    map(str, stress_weights),
                )
            ),
            flux_rates_nats_per_tick={n: str(s) for n, s in rates.items()},
            earlier_equal_ranges_and_forces_witness=dict(
                dimer_tangents=list(map(str, tangents)),
                regional_entropy_slope=str(entropy_slope),
            ),
            endpoint_neighborhood=points,
        )
    return dict(
        ticks=4, central_bridge=route[1:3], ten_equal_launches=rings, charges=results
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    output = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "regional.json").write_text(output)
    else:
        print(output, end="")
