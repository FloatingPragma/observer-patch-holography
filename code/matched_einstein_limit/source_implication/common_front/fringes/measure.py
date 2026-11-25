"""Remote spatial interference with a finite clock-resolved Fourier receiver."""

import argparse
import json
from itertools import combinations, permutations
from math import comb
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from ...action_response import measure as action
from ...action_response.four_ticks import raw_effect, volume_ratio

HERE = Path(__file__).resolve().parent
PORTS = (0, 5, 38, 39)
# Equal bins 1/2 to 6; unweighted clock points and three registered phase advances.
PHASORS = {"zero": 1, "quarter_turn": sp.I, "rational_turn": (3 + 4 * sp.I) / 5}


def matrices(board, route, q, native, co):
    positions = np.array([4 * p + j for p in PORTS for j in range(4)])
    sectors = []
    for odd in (False, True):
        h = native.generator(board, route, q, odd)
        u = native.kernel_panels(h, positions, co)
        ys = []
        for y in (2, 3):
            series = []
            for n in range(12):
                k = [
                    sp.Matrix(u[0, n, y, x]) + sp.I * sp.Matrix(u[1, n, y, x])
                    for x in (0, 1)
                ]
                series.append(
                    (k[1].H * k[0]).applyfunc(sp.expand) / native.prior.FRAC**2
                )
            ys.append(series)
        sectors.append(ys)
    return sectors


def records(blocks, source, reference_only=False):
    def quadrature(ms):
        # A whole 4x4 block has norm error <1e-19. Multiplying two
        # contractions costs <3e-19. Clock weights are bounded by 36,
        # and the matched tangent trace norms are below 999.
        error = sp.Rational("2e-17")
        tangent_allowance = iv.mpf(["-2e-14", "2e-14"])
        values = {"reference": raw_effect(ms, source[0], sp.Rational("6e-7"), error)}
        if reference_only:
            return values
        values.update(
            {
                k: j + tangent_allowance
                for k, j in action.source_jets(ms, source).items()
            }
        )
        for key, tgt, err in (
            ("raw_o", "ordinary", "12e-7"),
            ("raw_p", "unpair", "18e-7"),
        ):
            values[key] = raw_effect(ms, source[1][tgt], sp.Rational(err), 2 * error)
        return values

    re = quadrature([(m + m.H) / 2 for m in blocks])
    im = quadrature([(m - m.H) / (2 * sp.I) for m in blocks])
    return {k: iv.mpc(re[k], im[k]) for k in re}


def determinant(matrix):
    n = len(matrix)
    value = iv.mpf(0)
    for perm in permutations(range(n)):
        product = iv.mpf(
            (-1) ** sum(a > b for i, a in enumerate(perm) for b in perm[i + 1 :])
        )
        for i, j in enumerate(perm):
            product *= matrix[i][j]
        value += product
    return value


def frame_witnesses(rows):
    out = {}
    for label, keys in (
        ("vacancy", [k for k in rows if k[0] == "-1"]),
        ("excess", [k for k in rows if k[0] == "1"]),
        ("both", list(rows)),
    ):
        witnesses = []
        for group in combinations(keys, 5):
            # The first column is identically 1, so subtract the first row and
            # expand along that column. This transformation preserves rank.
            matrix = [
                [rows[k][j] - rows[group[0]][j] for j in range(1, 5)] for k in group[1:]
            ]
            det = determinant(matrix)
            if 0 not in det:
                witnesses.append((float((abs(det).a / abs(det).b).a), group, det))
        if witnesses:
            _, witness, cert = max(witnesses, key=lambda x: x[0])
            out[label] = dict(
                nonsingular_subsets=len(witnesses),
                total_subsets=comb(len(keys), 5),
                rows=[list(x) for x in witness],
                determinant=str(cert),
            )
        else:
            out[label] = dict(nonsingular_subsets=0, total_subsets=comb(len(keys), 5))
    return out


def run():
    iv.dps = 32
    native = action.native
    board = native.iso.pairing.formation.federation(1)
    route = native.iso.pairing.formation.region(board)
    assert PORTS[1] in board.neighbours[PORTS[0]]
    co = native.prior.weights()
    source = action.front.inputs()
    rho = volume_ratio(source)
    out = {}
    frame_rows = {}
    for q in (-1, 1):
        raw = matrices(board, route, q, native, co)
        row = {}
        for y, port in enumerate(PORTS[2:]):
            channel = {}
            for name, z in PHASORS.items():
                moments = []
                for order in (0, 1, 2):
                    blocks = [
                        sum(
                            (
                                sp.Rational(n, 2) ** order * z**n * m
                                for n, m in enumerate(sector[y], 1)
                            ),
                            sp.zeros(4),
                        ).applyfunc(sp.expand)
                        / 12
                        for sector in raw
                    ]
                    moments.append(records(blocks, source, reference_only=order == 2))
                f, g, g2 = moments
                base = f["reference"]
                assert abs(base).a > 0
                clock = g["reference"] / base
                slopes = {
                    k: ((g[k] - clock * f[k]) / base).real
                    for k in ("ordinary", "pairing", "raw_o", "raw_p")
                }
                omega = 2 * iv.atan2(
                    action.front.full.ri(sp.im(z)), action.front.full.ri(sp.re(z))
                )
                omega_dispersion = -omega * (g2["reference"] / base - clock**2).imag
                # From t_phys = b+a*t at fixed proper frequency, D tau acquires
                # D b + D a (tau + omega d tau/d omega). Keep both terms
                # independently to cover phase-indexed and frequency-indexed fits.
                without_rotation = slopes["raw_p"] - slopes["raw_o"] / rho
                frame_rows[(str(q), str(port), name)] = [
                    iv.mpf(1),
                    clock.real,
                    omega_dispersion,
                    slopes["ordinary"],
                    without_rotation,
                ]
                entry = dict(
                    omega_delay_dispersion=str(omega_dispersion),
                    pairing_minus_ordinary_over_volume_ratio=str(without_rotation),
                    visibility=str(abs(base)),
                    delay=str(clock.real),
                    delay_slopes={k: str(slopes[k]) for k in ("ordinary", "pairing")},
                )
                if 0 not in slopes["ordinary"]:
                    direct = slopes["pairing"] / slopes["ordinary"]
                    shared = (
                        1 / rho
                        + (slopes["raw_p"] - slopes["raw_o"] / rho) / slopes["ordinary"]
                    )
                    ratio = iv.mpf([max(direct.a, shared.a), min(direct.b, shared.b)])
                    entry["collar_over_ordinary"] = str(ratio)
                channel[name] = entry
            row[str(port)] = channel
        out[str(q)] = row
    witnesses = frame_witnesses(frame_rows)
    return dict(
        rotation_ratio_ordinary_over_pairing=str(rho),
        affine_clock_witnesses=witnesses,
        launch_ports=list(PORTS[:2]),
        receive_ports=list(PORTS[2:]),
        clock_ticks=[n / 2 for n in range(1, 13)],
        phase_advance_per_half_tick={k: str(v) for k, v in PHASORS.items()},
        charges=out,
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
