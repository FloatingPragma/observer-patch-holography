"""Compare independent local preparations, signed reads and recapture of the angular code."""

from math import factorial
from pathlib import Path
import json

import numpy as np
from scipy.linalg import expm

from . import measure as m

RIM = (7, 8, 9, 10, 11)
SHOTS = 2**20
TIMES = (0, 128, 512, 2048, 8192, 16384)


def local_emit(target, sites):
    """Prepare on the chosen owner's contiguous rim. No pulse uses the bridge."""
    u = np.eye(12, dtype=complex)
    station = 1 if sites[0] < 5 else 8
    source = sites[0]
    for j, (a, b) in enumerate(zip(sites, sites[1:])):
        tail = np.linalg.norm(target[list(sites[j + 1 :])])
        if tail > 1e-13 or abs(target[a]) > 1e-13:
            u = m.contact(a, b, 2 * np.arctan2(tail, abs(target[a]))) @ u
        if abs(target[a]) > 1e-13:
            u = m.phase(a, np.angle(target[a]) - np.angle(u[a, source]), station) @ u
    a = sites[-1]
    if abs(target[a]) > 1e-13:
        u = m.phase(a, np.angle(target[a]) - np.angle(u[a, source]), station) @ u
    np.testing.assert_allclose(u[:, source], target, atol=1e-12)
    return u


def pair_meter(plus, minus, sites):
    # First local inverse sends plus to the first port. A second inverse
    # on just the tail sends the orthogonal state to the second port.
    first = local_emit(plus, sites)
    v = first.conj().T @ minus
    assert abs(v[sites[0]]) < 1e-12
    second = local_emit(v, sites[1:])
    unit = second.conj().T @ first.conj().T
    np.testing.assert_allclose(unit @ plus, np.eye(12)[:, sites[0]], atol=1e-12)
    np.testing.assert_allclose(unit @ minus, np.eye(12)[:, sites[1]], atol=1e-12)
    return unit


def bases(side):
    z, y = m.mode(1, side), m.mode(-1, side)
    return [
        (z, y),
        ((z + y) / np.sqrt(2), (z - y) / np.sqrt(2)),
        ((z + 1j * y) / np.sqrt(2), (z - 1j * y) / np.sqrt(2)),
    ]


def hamiltonian():
    h = m.native()[0].copy()
    block = m.native()[1]
    for a, b in m.CUT:
        h[np.ix_([a, b], [a, b])] -= block
    return h


def exact_channel(h):
    """Nonzero determinant coefficient over Q(zeta_5), witnessed modulo 31."""
    prime = 31
    zeta = 2
    assert pow(zeta, 5, prime) == 1 and zeta != 1
    assert sum(zeta**k for k in range(5)) % prime == 0
    left = np.zeros((12, 2), int)
    right = np.zeros((2, 12), int)
    for j in range(5):
        for n, s in enumerate((1, -1)):
            left[j, n] = pow(zeta, (s * j) % 5, prime)
            right[n, j + 7] = pow(zeta, (-s * j) % 5, prime) * pow(5, -1, prime) % prime
    integer = np.rint(2 * h).astype(int)
    moments = []
    power = left
    for k in range(25):
        moments.append((right @ power) * pow(factorial(k), -1, prime) % prime)
        coeff = (
            sum(
                int(moments[j][0, 0]) * int(moments[k - j][1, 1])
                - int(moments[j][0, 1]) * int(moments[k - j][1, 0])
                for j in range(k + 1)
            )
            % prime
        )
        if coeff:
            return dict(prime=prime, fifth_root=zeta, order=k, nonzero_residue=coeff)
        power = integer @ power % prime
    raise AssertionError("No nonzero coefficient in the declared expansion")


def sampled_witness(reads):
    # Original-starts maps reconstructed from six local input/output axes.
    z = np.diag([1, -1])
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    density = []
    for row in reads[:6]:
        banks = np.asarray(row["recapture_counts"])
        trace = np.mean(banks[:, :2].sum(axis=1)) / SHOTS
        diff = (banks[:, 0] - banks[:, 1]) / SHOTS
        density.append(
            (trace * np.eye(2) + sum(v * a for v, a in zip(diff, (z, x, y)))) / 2
        )
    cross = (density[2] - density[3] + 1j * (density[4] - density[5])) / 2
    J = np.block([[density[0], cross], [cross.conj().T, density[1]]]) / 2
    pt = J.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)
    # Simultaneous two-sided Hoeffding bounds over 18 parities and 18
    # success indicators. The Frobenius norm is invariant under PT.
    ep = np.sqrt(2 * np.log(72 / 0.01) / SHOTS)
    local = np.sqrt(((ep / 2) ** 2 + 3 * ep**2) / 2)
    error = np.sqrt(2 * local**2 + 2 * (2 * local) ** 2) / 2
    return dict(
        min_pt=float(np.linalg.eigvalsh(pt)[0]), simultaneous99_radius=float(error)
    )


def rotation_echo(h, t, rng):
    # Recombine by a SECOND forward evolution through the same open bridge.
    u = expm(-1j * t * h)
    turn = np.eye(12, dtype=complex)
    live = {tuple(sorted(e)) for e in m.EDGES if e not in m.CUT}
    for j in reversed(range(7, 11)):
        assert (j, j + 1) in live
        turn = m.transpose(j, j + 1) @ turn
    for j in range(5):
        np.testing.assert_allclose(
            turn[:, 7 + j], np.eye(12)[:, 7 + (j + 1) % 5], atol=1e-12
        )
    # The four pulses act identically on the whole source owner and far hub.
    # Each addressed physical head stays blank on every participating load.
    np.testing.assert_allclose(turn[:, :7], np.eye(12)[:, :7], atol=1e-12)
    np.testing.assert_allclose(np.linalg.matrix_power(turn, 5), np.eye(12), atol=1e-12)
    source = local_emit(m.mode(1), tuple(range(5)))
    launched = source[:, 0]
    m.shutter(launched)
    rows = []
    for n in (0, 1, 5, 10, "supplied_remote_minus"):
        action = (
            np.linalg.matrix_power(turn, n)
            if isinstance(n, int)
            else np.diag([1] * 6 + [-1] * 6)
        )
        terminal = u @ action @ u @ launched
        terminal[5] = 0  # Negative source-hub read is retained as a miss.
        m.shutter(terminal)
        p = float(abs((source.conj().T @ terminal)[0]) ** 2)
        rows.append(
            dict(
                turns=n,
                source_return_probability=p,
                positive=int(rng.binomial(SHOTS, np.clip(p, 0, 1))),
            )
        )
    return rows


def recapture():
    rng = np.random.default_rng(63231)
    h = hamiltonian()
    originals = bases(0)
    target = bases(1)
    meters = [pair_meter(*pair, RIM) for pair in target]
    v = np.column_stack(originals[0])
    w = np.column_stack(target[0])
    uniform = local_emit(m.mode(0), tuple(range(5)))[:, 0]
    phase_sources = []
    for sign in (1, -1):
        u = uniform.copy()
        for j in range(5):
            u = m.phase(j, sign * 2 * np.pi * j / 5) @ u
        np.testing.assert_allclose(u, m.mode(sign), atol=2e-12)
        phase_sources.append(u)
    rows = []
    for t in TIMES:
        u = expm(-1j * h * t)
        code = w.conj().T @ u @ v
        singular = np.linalg.svd(code, compute_uv=False)
        trials = []
        for name, src in [
            (f"direct_{axis}_{i}", local_emit(state, tuple(range(5))))
            for axis, pair in enumerate(originals)
            for i, state in enumerate(pair)
        ]:
            trials.append((name, src[:, 0]))
        trials.extend((f"phase_grating_{i}", s) for i, s in enumerate(phase_sources))
        reads = []
        for name, state in trials:
            m.shutter(state)
            y = u @ state
            # A sharp test for vacancy at the far hub. The 'filled' branch keeps
            # all angular coherence and permits the reversible tip shutter.
            ready = y.copy()
            ready[6] = 0
            m.shutter(ready, 8, 6)
            counts = []
            for decoder in meters:
                z = decoder @ ready
                p = abs(z[list(RIM[:2])]) ** 2
                law = np.maximum([*p, 1 - sum(p)], 0)
                counts.append(rng.multinomial(SHOTS, law / law.sum()).tolist())
            # Independently record signed vacancy arrival and the current on
            # the intact first receiving rim bond from a native half exchange.
            arrival = sum(abs(y[6:]) ** 2)
            current = m.contact(7, 8, np.pi / 2) @ y
            p = abs(current[[7, 8]]) ** 2
            law = np.maximum([*p, 1 - sum(p)], 0)
            reads.append(
                dict(
                    source=name,
                    recapture_counts=counts,
                    signed_arrivals=int(rng.binomial(SHOTS, np.clip(arrival, 0, 1))),
                    current_bond=rng.multinomial(SHOTS, law / law.sum()).tolist(),
                )
            )
        rows.append(
            dict(
                time=t,
                squared_singular_values=(singular**2).tolist(),
                det_abs=float(abs(np.linalg.det(code))),
                sampled_witness=sampled_witness(reads),
                current_energy_variance=float(
                    np.linalg.norm(
                        h @ v[:, 0] - v[:, 0] * (v[:, 0].conj() @ h @ v[:, 0])
                    )
                    ** 2
                ),
                reads=reads,
            )
        )
    return dict(
        geometry=m.fixture(True),
        source_axes=["winding", "real", "imaginary"],
        starts=SHOTS,
        exact_determinant=exact_channel(h),
        records=rows,
        after_flight_turns={str(t): rotation_echo(h, t, rng) for t in (2048, 16384)},
    )


if __name__ == "__main__":
    Path(__file__).with_name("release.json").write_text(
        json.dumps(recapture(), indent=2) + "\n"
    )
