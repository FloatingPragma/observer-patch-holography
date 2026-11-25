"""Recover the sharp native optical image acceleration from terminal local material observations."""

import argparse
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from ...jacobi_screen import measure as optical
from . import measure

PAULI = (
    np.eye(2),
    np.array([[0, 1], [1, 0]]),
    np.array([[0, -1j], [1j, 0]]),
    np.diag([1, -1]),
)
HEADS = np.array([np.kron(a, b) for a in PAULI for b in PAULI])


def material_read(op, w):
    """Each real/imaginary port-interference reading is conditioned on a local Pauli result."""
    co = op.tocoo()
    pairs = sorted(
        {(i // 4, j // 4) for i, j in zip(co.row, co.col) if i // 4 <= j // 4}
    )
    effect = np.zeros((10, 4, 4), dtype=complex)
    count = 0
    for p, r in pairs:
        block = op[4 * p : 4 * p + 4, 4 * r : 4 * r + 4].toarray()
        coeff = np.einsum("aij,ji->a", HEADS, block) / 4
        np.testing.assert_array_equal(np.einsum("a,aij->ij", coeff, HEADS), block)
        for weight, pol in zip(coeff, HEADS):
            if weight == 0:
                continue
            mx = w[p].conj().swapaxes(-1, -2) @ pol @ w[r]
            if p == r:
                assert weight.imag == 0
                effect += weight.real * mx
            else:
                # X measurement has upper-right matrix element +1; Y has -i.
                effect += weight.real * (mx + mx.conj().swapaxes(-1, -2))
                effect -= weight.imag * (-1j * mx + 1j * mx.conj().swapaxes(-1, -2))
            count += int(weight.real != 0) + int(weight.imag != 0)
    return effect, count


def receiver(h, w, y):
    diag = np.zeros(h.shape[0], dtype=int)
    diag[4 * y : 4 * y + 4] = 1
    pop = sparse.diags(diag, dtype=np.int64, format="csr")
    vel = h @ pop - pop @ h
    acc = vel @ h - h @ vel
    # Check operator-by-operator sums of all phased incident physical bonds.
    one = h @ pop - pop @ h @ pop
    zports = sorted(set(one.nonzero()[0] // 4))
    currents = []
    for z in zports:
        keep = np.zeros_like(diag)
        keep[4 * z : 4 * z + 4] = 1
        incoming = sparse.diags(keep, dtype=int) @ one
        currents.append(incoming - incoming.T)
    assert (sum(currents) - vel).count_nonzero() == 0
    assert (sum(j @ h - h @ j for j in currents) - acc).count_nonzero() == 0
    effects, counts = [], []
    for op, scale in ((pop, 1), (1j * vel, 0.5), (acc, 0.25)):
        v, num = material_read(op, w)
        effects.append(scale * v)
        counts.append(num)
    return effects, counts


def contractions(jets, state, variation):
    coords = np.c_[np.ones(5), measure.bridge.beam.BEARINGS]
    probabilities = np.einsum("kbyxij,bji->kyx", jets, state).real
    slopes = np.einsum("kbyxij,bji->kyx", jets, variation).real
    m = coords.T @ probabilities @ coords
    dm = coords.T @ slopes @ coords
    im = optical.matrix_jets(m)
    _, b, c = optical.optical_matrices(im)
    db, dc = optical.deformation_response(im, optical.response_jets(m, dm))
    return np.array([b.trace(), -c.trace(), db.trace(), -dc.trace()])


def decode_images(raw, time_reads, rho, variation):
    vals = {}
    for tag, rs, cs in (
        ("forward", slice(5, 10), slice(5)),
        ("reverse", slice(5), slice(5, 10)),
    ):
        exact = contractions(raw[:, :, rs, cs], rho, variation)
        stencil = {}
        for dt, blocks in time_reads.items():
            jets = blocks[:, :, rs, cs]
            vals_st = contractions(jets, rho, variation)
            stencil[str(dt)] = vals_st.tolist()
        finest = np.array(stencil[str(min(time_reads))])
        coarse = np.array(stencil[str(2 * min(time_reads))])
        extrapolated = (4 * finest - coarse) / 3
        np.testing.assert_allclose(extrapolated, exact, rtol=3e-6, atol=2e-6)
        vals[tag] = dict(
            material_terminal_ward_read=exact.tolist(),
            arrival_time_stencils=stencil,
            extrapolated_arrivals=extrapolated.tolist(),
        )
    return vals


def run(arm=None):
    board = measure.bridge.native.iso.pairing.formation.federation(1)
    route = measure.bridge.native.iso.pairing.formation.region(board)
    rings = measure.bridge.beam.screens(board, route)
    nodes = sum(rings, [])
    idx = measure.bridge.beam.indices(nodes)
    central, _, _ = measure.source()
    p = 0.1343
    rho = np.array(
        [np.array(s, dtype=complex) * z for s, z in zip(central, (1 - p, p))]
    )
    if arm is None:
        arm = json.loads((measure.HERE / "control.json").read_text())
    co = np.array(arm["center_coordinates"]).reshape(2, 15)
    variation = np.einsum("bk,kij->bij", co, measure.c.iso.basis(False))
    for t in [-0.01, 0.01]:
        assert min(np.linalg.eigvalsh(z).min() for z in rho + t * variation) > 0
    output = {}
    banks = []
    for q in (-1, 1):
        halves = []
        time_reads = {s: [] for s in (1 / 256, 1 / 512, 1 / 1024)}
        counts = []
        for odd in (False, True):
            h = measure.bridge.native.generator(board, route, q, odd)
            start = np.zeros((h.shape[0], 40))
            start[idx] = np.eye(40)

            def flight(t, w):
                z = (-0.5j * t) * h
                return expm_multiply(z, w, traceA=z.diagonal().sum())

            central_w = flight(4, start)
            w = central_w.reshape(board.n, 4, 10, 4).transpose(0, 2, 1, 3)
            sector = []
            for y in nodes:
                reads, cs = receiver(h, w, y)
                sector.append(reads)
                counts.append(cs)
            halves.append(np.moveaxis(sector, 1, 0))
            for dt in time_reads:
                pjet = []
                for s in (-1, 0, 1):
                    a = flight(s * dt, central_w) if s else central_w
                    a = a[idx].reshape(10, 4, 10, 4).transpose(0, 2, 1, 3)
                    pjet.append(a.conj().swapaxes(-1, -2) @ a)
                minus, base, plus = pjet
                time_reads[dt].append(
                    [base, (plus - minus) / (2 * dt), (plus - 2 * base + minus) / dt**2]
                )
        raw = np.array(halves).swapaxes(0, 1)
        time_reads = {dt: np.array(v).swapaxes(0, 1) for dt, v in time_reads.items()}
        banks.append((raw, time_reads))
        vals = decode_images(raw, time_reads, rho, variation)
        output[str(q)] = dict(
            phase_and_head_settings_per_receiver_max=np.max(counts, axis=0).tolist(),
            directional_images=vals,
        )
    return dict(
        dwell=4,
        launches=rings,
        nominal_reused_odd_weight=p,
        order=[
            "expansion",
            "focusing",
            "source_expansion_slope",
            "source_focusing_slope",
        ],
        charges=output,
        fair_untagged_images=decode_images(
            (banks[0][0] + banks[1][0]) / 2,
            {dt: (banks[0][1][dt] + banks[1][1][dt]) / 2 for dt in time_reads},
            rho,
            variation,
        ),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "screen_ward.json").write_text(report)
    else:
        print(report, end="")
