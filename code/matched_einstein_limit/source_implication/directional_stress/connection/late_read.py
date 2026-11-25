"""Recover the hidden force contrast by a later single-port hit and two phase quadratures."""

import argparse
import json
from math import isqrt
from itertools import product

import mpmath as mp
import numpy as np
import sympy as sp

from .. import bridge
from ..traction.certificate import four_tick_columns
from . import measure


def mm(x):
    return mp.matrix(
        [
            [mp.mpc(str(z.real), str(z.imag)) for z in row]
            for row in np.asarray(x, dtype=complex)
        ]
    )


def mp_wave(a, b):
    # Rational fixed-point real/imaginary numerators. Conversion retains far
    # more digits than the native kernel bound.
    return mp.matrix(
        [
            [mp.mpc(int(x), int(y)) / bridge.native.prior.FRAC for x, y in zip(rx, ry)]
            for rx, ry in zip(a, b)
        ]
    )


def invertibility(a, b):
    fixed = (sp.Matrix(a) + sp.I * sp.Matrix(b)) / bridge.native.prior.FRAC
    rec = fixed.inv()
    square = sp.expand((rec.H * rec).trace())
    n = isqrt(int(sp.ceiling(square))) + 1
    assert n < 10**18
    return f"1/{2 * n}"


def trace(m):
    return sum(m[i, i] for i in range(m.rows))


def calibrate(a, b, basis, keep):
    """Product-input process tomography fixes both amplitude matrices up to one common phase."""
    columns = mp.matrix(4, 8)
    for i in range(4):
        for j in range(4):
            columns[i, j] = a[i, j]
            columns[i, j + 4] = b[i, j]
    single = (
        np.array([1, 0]),
        np.array([0, 1]),
        np.array([1, 1]) / np.sqrt(2),
        np.array([1, 1j]) / np.sqrt(2),
    )
    records = {}
    archive = {}
    for label in product(range(4), repeat=3):
        # Independent preparations of port rail, left coin and right coin.
        vec = np.kron(np.kron(single[label[0]], single[label[1]]), single[label[2]])
        state = mm(vec[:, None])
        state /= mp.sqrt(trace(state.H * state))
        output = columns * state
        density = output * output.H
        paulis = [mp.re(trace(p * density)) for p in basis]
        records[label] = sum((v * p / 4 for v, p in zip(paulis, basis)), mp.zeros(4))
        if keep:
            archive["".join(map(str, label))] = [mp.nstr(x, 34) for x in paulis]
    first = records[0, 0, 0]
    k = max(range(4), key=lambda z: mp.re(first[z, z]))
    anchor = mp.sqrt(mp.re(first[k, k]))
    assert anchor > mp.mpf(".3")
    # |0><0| and |1><0| in the four input projector basis.
    dual = ((1, 0, 0, 0), ((-1 + mp.j) / 2, (-1 + mp.j) / 2, 1, -mp.j))
    fitted = mp.zeros(4, 8)
    for j in range(8):
        bits = [(j >> 2) & 1, (j >> 1) & 1, j & 1]
        for label, rec in records.items():
            weight = mp.fprod(dual[v][z] for v, z in zip(bits, label))
            if weight:
                for i in range(4):
                    fitted[i, j] += weight * rec[i, k] / anchor
    phase = columns[k, 0] / anchor
    assert mp.norm(fitted * phase - columns) < mp.mpf("1e-40")
    return fitted[:, :4], fitted[:, 4:], archive


def decode_local(a, b, x, y, coin, dual, basis, keep):
    """Compute positive conditional hit records first. Reconstruct from their two quadratures."""
    incident = [[x, y] for x, y in zip(x, y)]
    base = (
        sum(
            (a * v * coin * v.H * a.H + b * t * coin * t.H * b.H for v, t in incident),
            mp.zeros(4),
        )
        / 10
    )
    transported = sum((a * v * coin * t.H * b.H for v, t in incident), mp.zeros(4)) / 10
    detected = {}
    archived = {}
    for k, phase in (("cos", 1), ("sin", mp.j)):
        interference = mp.conj(phase) * transported + phase * transported.H
        hands = []
        for h in (base + interference, base - interference):
            hands.append([mp.re(trace(pauli * h)) for pauli in basis])
        difference = [r - s for r, s in zip(*hands)]
        detected[k] = sum((r * p / 4 for r, p in zip(difference, basis)), mp.zeros(4))
        if keep:
            archived[k] = {
                "correlated": [mp.nstr(v, 36) for v in hands[0]],
                "anticorrelated": [mp.nstr(v, 36) for v in hands[1]],
            }
    response = mp.re(trace(dual * (detected["cos"] + mp.j * detected["sin"])))
    return response, archived


def run():
    from scipy import sparse

    mp.mp.dps = 52
    atlas = measure.translation.Atlas(1)
    board = atlas.board
    route = bridge.native.iso.pairing.formation.region(board)
    patch = measure.translation.chosen(atlas)
    shifted = measure.translation.map_matrix(board, patch)
    partners = [int(shifted[:, 4 * p].nonzero()[0][0]) // 4 for p in route]
    pairs = [
        (p, t)
        for i, p in enumerate(route)
        for t in partners[2 * (i // 2) : 2 * (i // 2) + 2]
    ]
    starts = sum((list(z) for z in bridge.beam.screens(board, route)), [])
    coeff = bridge.native.prior.weights()[7]
    states = bridge.action.preparations()
    sigma = [
        np.eye(2),
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    basis = [mm(np.kron(a, b)) for a in sigma for b in sigma]
    angles = {
        ("ordinary" if row["source"] == "ordinary" else "pairing"): row[
            "angle_per_unit"
        ]
        for row in json.loads(
            (bridge.HERE.parent / "signed_transfers.json").read_text()
        )["tests"]
    }

    def slopes(sums):
        return {
            k: float(
                sums[k] - sums["equilibrium"] + t * (sums["plus"] - sums["minus"]) / 2
            )
            for k, t in angles.items()
        }

    out = {}
    for q in (-1, 1):
        grand = {k: mp.mpf(0) for k in states}
        selected = {k: mp.mpf(0) for k in states}
        raw = {}
        for odd in (False, True):
            h, _, _, _, force4 = measure.controls(board, route, patch, q, odd)
            receivers = sorted(set(starts + list(route) + partners))
            rr = bridge.beam.indices(receivers)
            w = four_tick_columns(h, rr, rr, coeff)

            def block(p, t):
                dst, src = receivers.index(p), receivers.index(t)
                return w[:, 4 * dst : 4 * dst + 4, 4 * src : 4 * src + 4]

            tables = {}
            rebuild = sparse.lil_matrix(h.shape, dtype=int)
            for p, t in pairs:
                aa = block(p, p)
                bb = block(p, t)
                a, b = mp_wave(*aa), mp_wave(*bb)
                early_p = [mp_wave(*block(p, launch)) for launch in starts]
                early_t = [mp_wave(*block(t, launch)) for launch in starts]
                f = force4[4 * t : 4 * t + 4, 4 * p : 4 * p + 4].toarray()
                force = mm(f) / 4
                rebuild[4 * t : 4 * t + 4, 4 * p : 4 * p + 4] = f
                rebuild[4 * p : 4 * p + 4, 4 * t : 4 * t + 4] = f.T
                ca, cb, survey = calibrate(a, b, basis, (p, t) == (0, 229))
                dual = (cb.H**-1) * force * (ca**-1)
                recorded = {}
                for name, st in states.items():
                    keep = (p, t) == (0, 229) or name == "equilibrium"
                    response, history = decode_local(
                        a, b, early_p, early_t, mm(st[int(odd)]), dual, basis, keep
                    )
                    # Opposite lotteries flip the entire pair contribution, so
                    # their difference is twice the original expectation.
                    grand[name] += response / 2
                    if (p, t) == (0, 229):
                        selected[name] += response
                    recorded[name] = {"force_contrast": mp.nstr(response, 30)}
                    if keep:
                        recorded[name]["port_pauli"] = history
                row = dict(
                    rigorous_minimum_singular_value_lower_bounds=[
                        invertibility(*aa),
                        invertibility(*bb),
                    ],
                    numerical_work_dual_norm=float(mp.norm(dual)),
                    local_readouts=recorded,
                )
                if (p, t) == (0, 229):
                    walks = np.zeros((h.shape[0], 4), dtype=np.int64)
                    walks[4 * t : 4 * t + 4] = np.eye(4, dtype=np.int64)
                    assert int(abs(h).sum(axis=1).max()) ** 9 < 2**63
                    for _ in range(9):
                        assert not np.any(walks[4 * p : 4 * p + 4])
                        walks = h @ walks
                    first = sp.Matrix(walks[4 * p : 4 * p + 4])
                    assert first.det() == 1048576
                    row["calibration_product_inputs"] = survey
                    row["earliest_block"] = first.tolist()
                    row["integer_determinant"] = int(first.det())
                tables[f"{p}_{t}"] = row
            assert (rebuild.tocsr() - force4).count_nonzero() == 0
            raw[str(int(odd))] = tables
        independent = json.loads((measure.HERE / "results.json").read_text())[
            "charges"
        ][str(q)]["after_original_flight"]
        original = json.loads(
            (measure.HERE.parent / "translation/certificate.json").read_text()
        )[str(q)]["marked_force"]
        for sums, truth in ((selected, independent), (grand, original)):
            for name, val in {
                "reference": float(sums["equilibrium"]),
                **slopes(sums),
            }.items():
                lo, hi = json.loads(truth[name])
                assert lo < val < hi
        out[str(q)] = dict(
            reconstructed_witness_reference=float(selected["equilibrium"]),
            reconstructed_witness_slopes=slopes(selected),
            reconstructed_full_force=float(grand["equilibrium"]),
            reconstructed_full_slopes=slopes(grand),
            parities=raw,
        )
    return dict(
        appointment=4,
        extra_native_flight=4,
        final_ports=list(route),
        pair_count=len(pairs),
        calibration_product_settings_per_channel=64,
        phase_quadratures=[0, "pi/2"],
        pauli_order=["I", "X", "Y", "Z"],
        readout_type="simulation population expectations",
        charges=out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False, default=int) + "\n"
    if args.write:
        (measure.HERE / "late_read.json").write_text(result)
    else:
        print(result, end="")
