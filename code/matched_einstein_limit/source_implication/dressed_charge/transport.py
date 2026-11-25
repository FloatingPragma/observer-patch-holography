"""Conserved real polynomial extension on the propagating native 80-cell carrier."""

from pathlib import Path
import argparse
import json

import numpy as np
import sympy as sp
from scipy import sparse

from . import measure as compact


ORDERS = (0, *range(2, 21))
UPPER = np.triu_indices(4)


def powers(h, port):
    """Exact integer columns, using sparse graph edges and arbitrary-precision accumulators."""
    h = h.tocoo()
    wave = np.zeros((h.shape[0], 4), dtype=object)
    wave[4 * port : 4 * port + 4] = np.eye(4, dtype=int)
    for j in range(max(ORDERS) + 1):
        if j in ORDERS:
            yield j, wave
        if j < max(ORDERS):
            nxt = np.zeros_like(wave)
            np.add.at(nxt, h.row, h.data[:, None] * wave[h.col])
            wave = nxt


def sector(h, k):
    """Verify independent symmetric self and cross coordinates and replay both target signs."""
    p, s, r = slice(20, 24), slice(4, 8), slice(152, 156)
    np.testing.assert_array_equal(h[s, p].toarray(), -np.eye(4, dtype=int))
    np.testing.assert_array_equal(h[r, s].toarray(), -np.eye(4, dtype=int))
    centered = h - h[20, 20] * sparse.eye(h.shape[0], dtype=np.int64, format="csr")
    columns, control = [], []
    for _, w in powers(centered, 5):
        local, crossing = w[p], w[r]
        np.testing.assert_array_equal(local, local.T)
        columns.append([*local[UPPER], *(crossing + crossing.T)[UPPER]])
        control.append(w[s].ravel().tolist())
    matrix = sp.Matrix(columns).T
    minor = int(matrix.det()) % 1009
    assert minor
    # Exact arithmetic on the floating center's binary-rational input. The invertibility
    # proof above independently covers arbitrary real targets, including the true optimum.
    real = sp.Matrix([sp.Rational(float(a)) for a in k.real[UPPER]])
    rhs = sp.Matrix.hstack(*(real.col_join(2 * sign * real) for sign in (-1, 0, 1)))
    fitted = matrix.inv() * rhs
    assert matrix * fitted == rhs
    # A single off-diagonal beta-independent matrix element witnesses the scale
    # of these particular degree-20 lifts. Its value is a lower bound on operator norm.
    other = sp.Matrix(control).T * fitted
    lower = [float(max(abs(x) for x in other[:, j])) for j in range(3)]
    return dict(
        minor_mod_1009=minor, nominal_norm_lower_for_cross_minus_zero_plus=lower
    )


def real_pairing_certificate(x):
    """Enclose the true raw unpairing response of the real part of the frozen logarithm."""
    from ..full_algebra import certificate as cert
    from ..full_algebra.native_precision import prior
    from ..full_algebra.enclosure import iv

    iv.dps = 28
    x = [prior.sym_matrix(h) for h in x]
    x = [(h + h.H) / 2 for h in x]
    tr = sum(h.trace() for h in x)
    x = [h / tr for h in x]
    total = x[0] + x[1]
    right = sp.Matrix(
        2, 2, lambda i, j: sum(total[2 * k + i, 2 * k + j] for k in range(2))
    )
    left = [
        sp.Matrix(2, 2, lambda i, j: sum(z[2 * i + k, 2 * j + k] for k in range(2)))
        for z in x
    ]
    central = iv.mpf(0)
    for z, block in zip(x, left):
        d = (sp.kronecker_product(block, right) - z).applyfunc(sp.re)
        vals, err = prior.entropy_coordinates(z, 0, sp.Rational(".02"), [d])
        central += cert.ri_band(vals[0], 2 * err)
    radius = cert.ri("2e-7")
    error = (2 / (cert.ri(".02") - radius) + 3 * iv.sqrt(8) * 4) * radius
    enclosure = central + iv.mpf([-1, 1]) * error
    assert enclosure.a > cert.ri(".23105") and enclosure.b < cert.ri(".23112")
    return str(enclosure)


def calibrated_bridge_control(board, route, state, ordinary, unpair):
    """Keep the reference and ordinary bridge responses fixed while changing the held one."""
    from ..common_front import measure as front

    ruler = front.Ruler(board, route, 3)
    link = ruler.ops[:, ruler.ids[(0, 39)]] / 2
    pair = compact.station.pair
    coefficient = float(pair(link, state))
    # At the reference the direct shortest length is coefficient ** (-1/6).
    length_gradient = -link / (6 * coefficient ** (7 / 6))
    identity = np.stack([np.eye(4), np.eye(4)])
    energy = np.stack([np.zeros((4, 4)), np.eye(4)])
    ordinary_energy = float(pair(energy, ordinary))
    assert ordinary_energy > 0.2975
    gain = float(pair(length_gradient, ordinary)) / ordinary_energy
    control = length_gradient - gain * energy
    control -= float(pair(control, state)) * identity
    assert abs(pair(energy, unpair)) < 1e-12
    contrasts = [
        -float(pair(control, direction)) / 4 for direction in (state, ordinary, unpair)
    ]
    assert max(abs(x) for x in contrasts[:2]) < 1e-12
    # The common-front proof encloses the true length derivative on the same
    # survey-volume-matched unpairing. Centering and subtracting energy vanish there.
    source = json.loads((front.HERE / "certificate.json").read_text())
    low, high = source["3"]["matched_slopes"]["unpair"]
    outward = [-0.009052, -0.008940]
    assert outward[0] < -high / 4 < contrasts[2] < -low / 4 < outward[1]
    return dict(
        calibration_link=[0, 39],
        order=3,
        length_per_ordinary_energy=gain,
        nominal_reference_and_ordinary_and_unpairing_changes=contrasts,
        exact_source_unpairing_current_change_enclosure=outward,
    )


def run():
    form = compact.station.iso.pairing.formation
    board = form.federation(1)
    route = form.region(board)
    # A unique two-hop path reaches a different owner from launch 5.
    paths = [
        (v, w) for v in board.neighbours[5] for w in board.neighbours[v] if w // 12 != 0
    ]
    assert paths == [(1, 38)]
    chi, k = compact.read_source()
    left = np.trace(chi.reshape(2, 2, 2, 2, 2), axis1=2, axis2=4)
    right = np.trace(chi.sum(axis=0).reshape(2, 2, 2, 2), axis1=0, axis2=2)
    tangent = np.stack([np.kron(block, right) for block in left]) - chi
    from .. import signed_transfers as signed

    source = json.loads(Path(signed.__file__).with_suffix(".json").read_text())
    angle = next(
        row["angle_per_unit"] for row in source["tests"] if row["source"] == "collar"
    )
    bloch = np.array([np.trace(right @ a).real for a in signed.phase.AXES[1:]])
    g = np.kron(
        np.eye(2),
        sum(v * a for v, a in zip(bloch, signed.phase.AXES[1:]))
        / np.linalg.norm(bloch)
        / 2,
    )
    rotation = -1j * (g @ chi - chi @ g)
    matched = tangent + angle * rotation
    ordinary_angle = next(
        row["angle_per_unit"] for row in source["tests"] if row["source"] == "ordinary"
    )
    ordinary = chi.sum(axis=0) / 2 - chi + ordinary_angle * rotation
    results = {}
    for q in (-1, 1):
        results[str(q)] = [
            sector(compact.native_precision.generator(board, route, q, o), k[int(o)])
            for o in (False, True)
        ]
    return dict(
        orders=list(ORDERS),
        sector_order=["even", "odd"],
        sectors=results,
        receiving_cell=3,
        unique_two_hop_path=[5, 1, 38],
        cross_target="Re K",
        ordinary_calibrated_control=calibrated_bridge_control(
            board, route, chi, ordinary, matched
        ),
        true_raw_real_pairing_nats=real_pairing_certificate(chi),
        nominal_initial_current_rate=-float(compact.station.pair(k.real, chi)) / 4,
        nominal_raw_pairing_current_rate=-float(compact.station.pair(k.real, tangent))
        / 4,
        nominal_fixed_volume_pairing_current_rate=-float(
            compact.station.pair(k.real, matched)
        )
        / 4,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        Path(__file__).with_name("transport.json").write_text(record)
    else:
        print(record, end="")
