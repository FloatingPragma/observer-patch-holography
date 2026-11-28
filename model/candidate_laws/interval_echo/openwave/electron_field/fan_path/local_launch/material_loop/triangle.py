"""Build the shuttered material triangle from its actual microscopic contacts."""

from fractions import Fraction
from math import factorial

import numpy as np
from scipy.linalg import expm


MATERIAL = 4
METER_SIZE = 65
HUB = [3, 0]  # Ordered reported words a,b; the oriented material chain runs b to a.
I4 = np.eye(4)


def geometry(path):
    face = path.field.rings[0]
    ports = [path.cells[-1][i] for i in face]
    assert ports == [1045, 1049, 1053]
    mask = path.field.masks[0]
    cw = path.field.clockwise[0]
    all_bits = sum(path.field.bits.values())
    a, b = (334, 365) if path.sign == 1 else (334 ^ all_bits, 365 ^ all_bits)
    assert b == a ^ mask
    vertices = [face[0], face[1], face[2], face[0]]
    words = [b]
    for u, v in zip(vertices, vertices[1:]):
        words.append(words[-1] ^ path.field.bits[tuple(sorted((u, v)))])
    assert words[-1] == a
    basis = list(zip(vertices, words, strict=True))
    spectators = sorted({w[:-1] for p, w in path.basis if p == ports[0] and w[-1] == a})
    # Rebuild every allowed triangle term, including its Laplacian diagonal.
    sub = np.zeros((4, 4))
    for i, (x, w) in enumerate(basis):
        for (u, v), bit in path.field.bits.items():
            if x not in (u, v) or u not in face or v not in face:
                continue
            if ((x == u) == bool(w & bit)) == (path.sign == 1):
                j = basis.index((v if x == u else u, w ^ bit))
                sub[j, i] -= MATERIAL
                sub[i, i] += MATERIAL
        if w & mask in (cw, mask ^ cw):
            j = basis.index((x, w ^ mask))
            sub[j, i] -= 0.5
            sub[i, i] += 0.5
        sub[i, i] += 8 * (6 - int(path.board.degree[path.cells[-1][x]]))
    np.testing.assert_array_equal(sub, block(0.5, 0))
    off = ~np.eye(4, dtype=bool)
    for upstream in spectators:
        rows = [path.ids[path.cells[-1][x], (*upstream, w)] for x, w in basis]
        np.testing.assert_array_equal(path.h[rows][:, rows].toarray()[off], sub[off])
    return {
        "port_word_in_oriented_chain_order": [[path.cells[-1][x], w] for x, w in basis],
        "upstream_word_combinations": len(spectators),
        "degrees": [int(path.board.degree[p]) for p in ports],
    }


def block(g, phase):
    """Flux is on the physical 1049 -> 1053 transfer, including its meter."""
    h = np.diag([4 + g, 16, 16, 4 + g]).astype(complex)
    for i, j, rate in [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, g)]:
        h[j, i] = -rate * (np.exp(1j * phase) if i == 1 else 1)
        h[i, j] = h[j, i].conjugate()
    return h


def counted_blocks(g, wait, size):
    """Fourier inversion yields final (register, material/word, input word)."""
    eig, v = np.linalg.eigh([block(g, 2 * np.pi * k / size) for k in range(size)])
    unitary = (v * np.exp(-1j * wait * eig)[:, None, :]) @ v.conj().transpose(0, 2, 1)
    out = np.fft.fft(unitary[:, :, HUB], axis=0) / size
    np.testing.assert_allclose(
        out.reshape(-1, 2).conj().T @ out.reshape(-1, 2), np.eye(2), atol=2e-14, rtol=0
    )
    np.testing.assert_allclose(
        out.sum(axis=0), expm(-1j * wait * block(g, 0))[:, HUB], atol=2e-13, rtol=0
    )
    return out


def shifted_apply(a, g):
    """Direct meter recurrence for (H - 12I), independent of Fourier signs."""
    b = np.empty_like(a)
    b[:, 0] = (g - 8) * a[:, 0] - 4 * a[:, 1] - g * a[:, 3]
    b[:, 1] = 4 * a[:, 1] - 4 * a[:, 0] - 4 * np.roll(a[:, 2], -1, axis=0)
    b[:, 2] = 4 * a[:, 2] - 4 * np.roll(a[:, 1], 1, axis=0) - 4 * a[:, 3]
    b[:, 3] = (g - 8) * a[:, 3] - 4 * a[:, 2] - g * a[:, 0]
    return b


def taylor_check(g, wait, size):
    """32 unitary substeps. The scalar global shift is restored at the end."""
    wave = np.zeros((size, 4, 2), complex)
    wave[0, HUB] = np.eye(2)
    for _ in range(32):
        term = wave
        partial = wave.copy()
        for n in range(1, 25):
            term = (-1j * wait / (32 * n)) * shifted_apply(term, g)
            partial += term
        wave = partial
    return wave * np.exp(-12j * wait)


def exact_leading_determinant():
    # Integer native 2H with meter initially zero. Each endpoint entry first
    # appears before wraparound. (-i)^n converts the H^n/n! coefficients to U.
    cols = np.zeros((9, 4, 2), dtype=np.int64)
    cols[0, [0, 3]] = np.eye(2, dtype=np.int64)
    hits = {}
    for n in range(6):
        for m in (-1, 1):
            for i, row in enumerate([0, 3]):
                for j in (0, 1):
                    entry = int(cols[m, row, j])
                    if entry:
                        hits.setdefault(
                            (m, i, j), (n, Fraction(entry, 2**n * factorial(n)))
                        )
        # At native rate 1/2, 2*(H-12I) and 2H have exact integer entries.
        nxt = np.empty_like(cols)
        nxt[:, 0] = 9 * cols[:, 0] - 8 * cols[:, 1] - cols[:, 3]
        nxt[:, 1] = (
            32 * cols[:, 1] - 8 * cols[:, 0] - 8 * np.roll(cols[:, 2], -1, axis=0)
        )
        nxt[:, 2] = (
            32 * cols[:, 2] - 8 * np.roll(cols[:, 1], 1, axis=0) - 8 * cols[:, 3]
        )
        nxt[:, 3] = 9 * cols[:, 3] - 8 * cols[:, 2] - cols[:, 0]
        cols = nxt
    for m in (-1, 1):
        a, b, c, d = (hits[m, i, j] for i in (0, 1) for j in (0, 1))
        assert a[0] + d[0] == b[0] + c[0] == 8
        assert a[1] * d[1] - b[1] * c[1] == Fraction(4**6, 4 * 2880)
    return str(Fraction(4**6, 4 * 2880))


def get_counters():
    records = {}
    for g in (0.5, 0.25, 0.0):
        for eighths in (0, 1, 3, 8):
            wait = eighths * np.pi / 8
            k = counted_blocks(g, wait, METER_SIZE)
            direct = taylor_check(g, wait, METER_SIZE)
            assert np.linalg.norm((k - direct).reshape(-1, 2), ord=2) < 2e-13
            larger = counted_blocks(g, wait, 2 * METER_SIZE - 1)
            np.testing.assert_allclose(
                k[[-1, 0, 1]], larger[[-1, 0, 1]], atol=2e-13, rtol=0
            )
            records[g, eighths] = k
    return records


def spatial_current(i, j):
    jmat = np.zeros((4, 4), complex)
    jmat[j, i] = 1j * MATERIAL
    jmat[i, j] = -1j * MATERIAL
    return jmat


def side_analyzers():
    """A quarter swap on each isolated contact followed by its two sharp port reads."""
    banks = []
    for i in range(3):
        j = i + 1
        pair = I4[:, [i, j]]
        differential = np.outer(pair[:, 1], pair[:, 1]) - np.outer(
            pair[:, 0], pair[:, 0]
        )
        swap = np.outer(pair[:, 0], pair[:, 1]) + np.outer(pair[:, 1], pair[:, 0])
        # Equalize the two on-site voltages for the supplied analyzer only.
        turn = expm(1j * np.pi / 4 * swap)
        np.testing.assert_allclose(
            MATERIAL * turn.conj().T @ differential @ turn,
            spatial_current(i, j),
            atol=3e-15,
        )
        banks.append((i, j, turn))
    return banks


def stationary_meter():
    x = 1 + np.sqrt(2)
    g = MATERIAL * x**2 / (2 * x - 1)
    target = 8 + 4 * np.sqrt(2)
    h = block(g, 0)
    e, v = np.linalg.eigh(h)
    np.testing.assert_allclose(e[1:3], target, atol=1e-13, rtol=0)
    project = v[:, 1:3] @ v[:, 1:3].conj().T
    filters = [
        (I4 + expm(-1j * np.pi * (h - target * I4) / abs(outside - target))) / 2
        for outside in (e[0], e[3])
    ]
    pass_both = filters[1] @ filters[0]
    np.testing.assert_allclose(pass_both, project, atol=1e-13, rtol=0)
    currents = [spatial_current(i, i + 1) for i in range(3)]
    for side in currents[1:]:
        np.testing.assert_allclose(
            project @ side @ project,
            project @ currents[0] @ project,
            atol=5e-14,
            rtol=0,
        )
    j = currents[1]
    slopes = np.linalg.eigvalsh(v[:, 1:3].conj().T @ j @ v[:, 1:3])
    np.testing.assert_allclose(slopes.sum(), 0, atol=1e-13)
    native = np.linalg.eigvalsh(block(0.5, 0))
    np.testing.assert_allclose(
        native,
        sorted([8 - 4 * np.sqrt(2), 8 + 4 * np.sqrt(2), 4, 21]),
        atol=1e-13,
        rtol=0,
    )
    probes = []
    for phase in (1 / 32, 1 / 64):
        tuned = np.linalg.eigvalsh(block(g, phase))
        original = np.linalg.eigvalsh(block(0.5, phase))
        for rate, values in ((g, tuned), (0.5, original)):
            np.testing.assert_allclose(
                values, np.linalg.eigvalsh(block(rate, -phase)), atol=1e-13, rtol=0
            )
        splitting = float((tuned[2] - tuned[1]) / (2 * phase))
        assert abs(splitting / slopes[1] - 1) < 0.01
        probes.append(
            {
                "phase_magnitude": phase,
                "native_four_lines_either_sign": original.tolist(),
                "tuned_four_lines_either_sign": tuned.tolist(),
                "split_over_twice_phase": splitting,
            }
        )
    native_e, native_v = np.linalg.eigh(block(0.5, 0))
    np.testing.assert_allclose(native_e, native, atol=2e-14)
    native_lines = [np.outer(vv, vv.conj())[:, HUB] for vv in native_v.T]
    return {
        "ring_rate": g,
        "target_energy": target,
        "native_energies": native.tolist(),
        "tuned_energies": e.tolist(),
        "persistent_material_current_extrema": slopes.tolist(),
        "external_phase_comparison": probes,
    }, (filters[0][:, HUB], pass_both[:, HUB], j, expm(-1j * np.pi * h), native_lines)
