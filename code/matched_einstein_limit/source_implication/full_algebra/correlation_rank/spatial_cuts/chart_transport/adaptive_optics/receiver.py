"""Synthesize a native terminal Pauli lottery using only independent material inputs."""

from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import iv, mp

from .. import held, transfer

HERE = Path(__file__).parent
TOTAL = 2**48
SHIFT = TOTAL // 2048
cuts = transfer.cuts
station = held.station


def emit_atoms():
    pauli = [(np.eye(2) + s * p) / 2 for p in station.phase.AXES[1:] for s in (-1, 1)]
    return np.array([np.kron(a, b) for a in pauli for b in pauli])


def dyadic_lottery(values, offset):
    # Largest-remainder rounding conserves the sum exactly. These are native
    # preparation probabilities, each specified by 48 fair bits.
    ratios = [Fraction.from_float(float(x)) for x in values]
    scaled = [v * TOTAL / sum(ratios) for v in ratios]
    counts = np.array([v.__floor__() for v in scaled], dtype="int64")
    order = sorted(
        range(len(values)), key=lambda i: scaled[i] - int(counts[i]), reverse=True
    )
    counts[order[: TOTAL - int(counts.sum())]] += 1
    # Holding the second atom fixed, both X signs and both Z signs sum to I.
    # The endpoints choose opposite decompositions of the same kept operator.
    counts[[0, 6]] += offset
    counts[[24, 30]] -= offset
    if min(counts) <= 0 or counts.sum() != TOTAL:
        raise ValueError("the inferred preparation is outside the positive lottery")
    return counts


def source_stock(counts):
    atoms = emit_atoms()
    # All projector entries are Gaussian integers divided by four. Numerators
    # stay below 2^53 since each group of positive counts has mass <= TOTAL.
    n = np.array([np.rint(4 * atoms.real), np.rint(4 * atoms.imag)], dtype="int64")
    kept = n.transpose(0, 2, 3, 1) @ counts[:36]
    mixed = n.transpose(0, 2, 3, 1) @ counts[36:]
    both = np.stack((2 * kept + mixed, mixed))
    return (both[:, 0] + 1j * both[:, 1]) / (8 * TOTAL)


def calibrate(keys, rows, observations, offset):
    b = np.array([rows[a][j] for a, j in keys])
    inv = transfer.center(b) ** -1
    raw = mp.matrix([mp.mpf(v.numerator) / v.denominator for v in observations])
    learned = np.array((inv * raw).tolist(), dtype=float).ravel()
    fair = np.tile(np.eye(4), (2, 1, 1)) / 8
    density = fair + held.blocks(learned)
    prep = station.preparation(density)
    counts = dyadic_lottery(prep["kept"] + prep["dephased"], offset)
    stock = source_stock(counts)
    # Verify the entire solver and finite switch approximation in native
    # Hilbert-Schmidt coordinates. No input density is available to this reader.
    j = station.iso.basis(True)
    loc = np.asarray(transfer.full.local_basis(), dtype=complex)
    basis = np.concatenate((j, loc))
    axes = np.rint(basis.real * 2).astype("int64") + 1j * np.rint(
        basis.imag * 2
    ).astype("int64")
    # Dyadic entries of 2*stock have <= 53 bits, but summing Pauli products
    # may need more. Take the Pauli traces in exact Fraction arithmetic.
    a = stock * 8 * TOTAL
    value = []
    per = []
    for block in a:
        comp = []
        for op in axes:
            numerator = sum(
                int(x.real) * int(y.real) - int(x.imag) * int(y.imag)
                for x, y in zip(block.flat, op.T.flat)
            )
            comp.append(iv.mpf(numerator) / (16 * TOTAL))
        per.append(comp)
    value.extend(per[0][:9] + per[1][:9] + per[0][9:] + per[1][9:])
    contrast = round(np.trace(a[0]).real) - round(np.trace(a[1]).real)
    value.append(iv.mpf(int(contrast)) / (8 * TOTAL) / iv.sqrt(8))
    raw_iv = np.array([iv.mpf(v.numerator) / v.denominator for v in observations])
    R = transfer.intervals(transfer.publish(inv))
    defect = cuts.norm(np.eye(31, dtype=int) - R @ b)
    residual = cuts.norm(np.array(value, dtype=object) - R @ raw_iv) + defect
    assert residual.b < cuts.ri("1e-10")
    return (
        counts,
        stock.transpose(0, 2, 1),
        dict(
            min_count=int(min(counts)),
            max_count=int(max(counts)),
            counts=counts.tolist(),
            implementation_and_inverse_error_upper=transfer.upper(residual),
            source_gain=transfer.upper(cuts.norm(R)),
        ),
    )
