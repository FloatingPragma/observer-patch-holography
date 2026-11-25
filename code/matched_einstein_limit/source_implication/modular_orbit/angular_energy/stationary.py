"""Exact bound on waiting-time dephasing from a rational commuting approximation."""

import numpy as np
import sympy as sp
from mpmath import iv

from model.experiments.interval_echo_angular import modular
from model.experiments.interval_echo_einstein_station import cap


def commuting_distance_squared(words, state, ham):
    """Weighted least squares on exact orbital indicators of symmetric matrices."""
    lookup = {word: i for i, word in enumerate(words)}
    permutations = [
        [lookup[tuple(sorted(g[k] for k in w))] for w in words]
        for g in modular.rotations()
    ]
    size = len(words)
    seen = set()
    basis = []
    entries = []
    for i in range(size):
        for j in range(i, size):
            if (i, j) in seen:
                continue
            orbit = {(p[i], p[j]) for p in permutations} | {
                (p[j], p[i]) for p in permutations
            }
            seen.update(orbit)
            indicator = np.zeros((size, size), dtype=np.int64)
            for a, b in orbit:
                assert state[a, b] == state[i, j]
                indicator[a, b] = 1
            basis.append(indicator)
            entries.append(int(state[i, j]))
    basis = np.array(basis)
    comm = (ham @ basis - basis @ ham).reshape(len(basis), -1)
    equations = sp.Matrix(np.unique(comm.T, axis=0).tolist())
    kernel = sp.Matrix.hstack(*equations.nullspace())
    weight = sp.diag(*[int(m.sum()) for m in basis])
    source = sp.Matrix(entries)
    reduced = kernel.T * weight
    error = source - kernel * (reduced * kernel).inv() * (reduced * source)
    return (error.T * weight * error)[0]


def dephasing_bound(low, high, den):
    # Native mixed-endpoint energy is 8I + (I-SWAP)/2.
    ham = np.zeros(low.shape, dtype=np.int64)
    for edge in cap.cell().seams:
        a, b = cap.seam_indices(edge)
        ham[a, a] += 17
        ham[b, b] += 17
        ham[a, b] -= 1
        ham[b, a] -= 1
    error2 = sp.Integer(0)
    for sl in (slice(0, 1), slice(1, 13), slice(13, 79)):
        error2 += commuting_distance_squared(cap.WORDS[sl], low[sl, sl], ham[sl, sl])
    trace_norm = (
        iv.sqrt(79)
        * (
            iv.sqrt(iv.mpf(int(error2.p)) / int(error2.q))
            + iv.sqrt(iv.mpf(int(((high - low) ** 2).sum())))
        )
        / den
    )
    return trace_norm


def waiting_protocol(reference, hamiltonian):
    """Construct finite fair waits from isolated roots of the native integer Hamiltonian."""
    twice = (2 * hamiltonian).astype(np.int64)
    np.testing.assert_array_equal(twice, 2 * hamiltonian)
    result = {}
    final = np.zeros_like(reference)
    for name, sl in (
        ("vacuum", slice(0, 1)),
        ("one", slice(1, 13)),
        ("two", slice(13, 79)),
    ):
        block = sp.Matrix(twice[sl, sl].tolist())
        roots = (
            block.charpoly().as_poly().sqf_part().intervals(eps=sp.Rational(1, 10**18))
        )
        bands = []
        for (low, high), multiplicity in roots:
            assert multiplicity == 1
            bands.append(
                iv.mpf([int(low.p), int(low.p)]) / int(low.q)
                + iv.mpf([0, 1]) * iv.mpf(int((high - low).p)) / int((high - low).q)
            )
        eig, u = np.linalg.eigh(hamiltonian[sl, sl])
        gaps = eig[:, None] - eig
        waiting = np.ones(gaps.shape, dtype=complex)
        max_time = iv.mpf(0)
        appointments = []
        for j, high in enumerate(bands):
            for low in bands[:j]:
                assert high.a > low.b
                duration = 2 * iv.pi / (high - low)
                mid = float(duration.mid)
                appointments.append(mid)
                max_time += duration
                waiting *= (1 + np.exp(-1j * gaps * mid)) / 2
        active = u.T @ reference[sl, sl] @ u
        stopped = u @ (waiting * active) @ u.T
        assert np.max(abs(stopped.imag)) < 1e-13
        final[sl, sl] = stopped.real
        result[name] = dict(
            levels=len(roots), fair_waits=len(appointments), duration_sum=str(max_time)
        )
    return final, result
