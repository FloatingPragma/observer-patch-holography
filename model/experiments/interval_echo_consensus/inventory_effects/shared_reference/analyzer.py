"""Independent forecasts and sector-resolved analyzer settings from Fourier minors."""

import numpy as np


def fourier(m):
    return np.exp(2j * np.pi * np.outer(np.arange(m), np.arange(m)) / m) / np.sqrt(m)


def occupied(mask, m):
    return [i for i in range(m) if mask & (1 << i)]


def minor(u, y, x):
    a, b = occupied(y, len(u)), occupied(x, len(u))
    if len(a) != len(b):
        return 0.0
    return np.linalg.det(u[np.ix_(a, b)]) if a else 1.0


def receiver(m, alice):
    count = m - alice.bit_count()
    if count in (0, m):
        return np.eye(m), {mask: int(count == m) for mask in range(1 << m)}
    # Columns correspond to local one-particle addresses or to their missing
    # addresses for the single-hole description of two loads on three rails.
    full = (1 << m) - 1
    masks = [1 << j if count == 1 else full ^ (1 << j) for j in range(m)]
    transform = fourier(m)
    components = np.zeros((2, m), complex)
    for j, mask in enumerate(masks):
        src = full ^ mask
        components[src & 1, j] = minor(transform, alice, src)
    components /= np.linalg.norm(components, axis=1)[:, None]
    pos, neg = (
        (components[0] + components[1]) / np.sqrt(2),
        (components[0] - components[1]) / np.sqrt(2),
    )
    rows = [pos.conj(), neg.conj()]
    if m == 3:
        rows.append(np.cross(pos, neg))
    chart = np.array(rows)
    # wedge^(m-1) U in hole order is det(U) D conj(U) D.
    if count != 1:
        signs = (-1.0) ** np.arange(m)
        chart = signs[:, None] * chart.conj() * signs
    # 1 denotes a positive-phase guess. All complement outputs are kept.
    return chart, {mask: int(mask != masks[1]) for mask in range(1 << m)}
