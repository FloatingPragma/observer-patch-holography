"""Collective dispersion and physical domain-wall clicks in locally paired material."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm, svd

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from model.executable_model.carrier import seam_digest

HERE = Path(__file__).resolve().parent
K = float(candidate.KAPPA)
# Magnitude of the blank contact transition between the two dimer orientations.
_blank_law = candidate.native_law
H = -float(_blank_law.seam_term(0)[_blank_law.SRC.index((0, 1)), _blank_law.SRC.index((1, 0))])
SEEDS = (312141, 312142, 312143)
DT, SAMPLES = 0.08, 16384
# Importance sampling just one physical product-Pauli setting from each
# normalized expansion at each endpoint. One complete trial has product +/-1.
SHOTS = 2**32


def native_checks():
    law = candidate.native_law
    exchange, contrast = law.local_term_parts()
    basis = np.eye(law.DIM)
    words = [(0, 0), (0, 1), (1, 0), (1, 1)]
    empty = law.echo() @ basis[:, law.SRC.index((0, 0))]
    source = sum(basis[:, law.SRC.index(w)] for w in [(2, 0), (0, 2)]) / np.sqrt(2)
    drawn = law.echo() @ source
    idle = empty.reshape(-1, law.N_SRC)[:, law.SRC.index((0, 0))]
    collar = drawn.reshape(-1, law.N_SRC)[:, law.SRC.index((1, 1))]
    assert np.isclose(np.dot(collar, collar), 1)
    block = np.array([[0.5, -0.5], [-0.5, 0.5]])
    dimer = basis[:, [law.SRC.index(w) for w in words[1:3]]]
    residuals = [
        np.max(abs(exchange @ dimer - dimer @ block)),
        np.max(abs(contrast @ dimer - dimer)),
    ]
    for word in words:
        st = np.eye(law.N_SRC)[:, law.SRC.index(word)]
        for memory, cost in ((idle, 0), (collar, 4 * (1 - np.prod(word)))):
            col = np.kron(memory, st)
            residuals.extend(
                [np.max(abs(exchange @ col)), np.max(abs(contrast @ col - cost * col))]
            )
    assert max(residuals) < 1e-12
    # Unknown-neighbour spin echo gives a diagonal rotation on either dimer.
    # Ordering here is occupation of the exposed link endpoint (0 or 1).
    p = np.diag([0, 1])
    x = np.array([[0, 1], [1, 0]])
    ident = np.eye(4)
    two = 4 * K * (ident - np.kron(p, p))
    gate = expm(-1j * np.pi / (8 * K) * two)
    flip = np.kron(np.eye(2), x)
    effect = gate @ flip @ gate @ flip
    ideal = np.kron(np.diag([-1, -1j]), np.eye(2))
    residuals.append(np.max(abs(effect - ideal)))
    assert max(residuals) < 1e-12
    return dict(columns_per_seam=law.DIM, max_amplitude_residual=float(max(residuals)))


def grow(board, seed):
    edges = np.asarray(board.seams)
    rng = np.random.default_rng(seed)
    partners = []
    for stage in range(2):
        a, b = edges.T
        if stage == 0:
            ok = board.degree[a] == board.degree[b]
        else:
            ok = (partners[0][a] >= 0) & (partners[0][b] >= 0) & (partners[0][a] != b)
        mate = np.full(board.n, -1)
        for idx in rng.permutation(np.flatnonzero(ok)):
            a, b = edges[idx]
            if mate[a] == -1 and mate[b] == -1:
                mate[a], mate[b] = b, a
        partners.append(mate)
    dimers, links = partners
    seen = set()
    out = []
    for a in sorted(np.flatnonzero(dimers >= 0), key=lambda a: links[a] >= 0):
        if a in seen:
            continue
        port, route = a, []
        while port >= 0 and port not in seen:
            other = dimers[port]
            route.append((int(port), int(other)))
            seen.update([port, other])
            port = links[other]
        out.append((port == a, route))
    return out


def dual(n, j):
    # Sites in the dual open chain are the n+1 bonds including the fixed ends.
    # Staggering or the sign of the off-diagonal coupling fixes the k origin.
    m = n + 1
    c = 2 * (j * np.eye(m) - H * np.eye(m, k=-1))
    u, energies, vh = svd(c)
    return u, energies, vh.T


def exact_check(n):
    """Build the physical 2^n orientation matrix with *its* two end fields."""
    states = np.arange(2**n)
    z = 1 - 2 * ((states[:, None] >> np.arange(n)) & 1)
    mat = np.diag(K * (np.sum(z[:, :-1] * z[:, 1:], axis=1) + z[:, 0] - z[:, -1]))
    for i in range(n):
        mat[states, states ^ (1 << i)] -= H
    ev = np.linalg.eigvalsh(mat)
    _, energies, _ = dual(n, K)
    occ = (np.arange(2 ** (n + 1))[:, None] >> np.arange(n + 1)) & 1
    selected = occ[np.sum(occ, axis=1) % 2 == n % 2] - 0.5
    error = float(max(abs(ev - np.sort(selected @ energies))))
    assert error < 1e-11
    # An even Majorana bilinear is one physical product of local dimer Paulis.
    j0, k0 = 1, n - 1
    left, right = j0, k0
    extended = np.column_stack([np.ones(len(states)), z, -np.ones(len(states))])
    tau = -extended[:, :-1] * extended[:, 1:]
    W = np.zeros_like(mat, dtype=complex)
    mask = sum(1 << i for i in range(left, right))
    W[states ^ mask, states] = 1j * np.prod(tau[:, left:right], axis=1)
    np.testing.assert_allclose(W, W.conj().T)
    T = 2.35
    propagator = expm(-1j * T * mat)
    direct = float(
        np.trace(propagator.conj().T @ W @ propagator @ W).real / len(states)
    )
    u, energies, _ = dual(n, K)
    aa = (u * np.cos(T * energies)) @ u.T
    jw = float(np.linalg.det(aa[np.ix_([left, right], [left, right])]))
    assert abs(direct - jw) < 1e-11

    # Check physical single-bond load-product probabilities as well as energies.
    # The parity permits a common extra wall at the far end on even n.
    def sharp(q):
        marks = np.zeros(n + 1, dtype=int)
        marks[q] = 1
        marks[-1] = 1 - n % 2
        spins = np.cumprod(2 * marks - 1)[:-1]
        return int(((spins < 0) * (1 << np.arange(n))).sum())

    full_z = np.column_stack([np.ones(2**n), z, -np.ones(2**n)])
    indicators = (1 + full_z[:, :-1] * full_z[:, 1:]) / 2
    t = 2.35
    prop = expm(-1j * t * mat)
    density = (abs(prop[:, sharp(1)]) ** 2 - abs(prop[:, sharp(3)]) ** 2) @ indicators
    kernel = walls(n, K, t)
    residual = float(max(abs(density - (kernel[:, 1] - kernel[:, 3]))))
    assert residual < 1e-11
    return dict(
        dimers=n,
        physical_states=2**n,
        max_spectrum_error=error,
        density_error=residual,
        even_string_trace_error=abs(direct - jw),
    )


def walls(n, j, t):
    u, e, v = dual(n, j)
    co, si = np.cos(t * e), np.sin(t * e)
    a, b = (u * co) @ u.T, (v * co) @ v.T
    ab, ba = (u * si) @ v.T, -(v * si) @ u.T
    return a * b - ab * ba


def wedge_norm(f, g):
    return np.abs(np.triu(np.outer(f, g) - np.outer(g, f), 1)).sum()


class CountStream:
    def __init__(self, sample_output):
        self.sample_output = sample_output
        self.hashes = {}
        self.sample = {}

    def keep(self, key, clicks):
        k = np.asarray(clicks, dtype="<u8", order="C")
        self.hashes[key] = hashlib.sha256(k.tobytes()).hexdigest()
        if self.sample_output and (
            key.startswith("clock") or key.startswith(("p1_2_", "p1_3_", "p2_3_"))
        ):
            self.sample[key] = k

    def finish(self, n):
        if self.sample_output:
            np.savez_compressed(
                self.sample_output, n=n, shots=SHOTS, clock_shots=2**20, **self.sample
            )


def counts(signal, normalization, rng, stream, key):
    prob = np.clip((1 + signal / normalization) / 2, 0, 1)
    clicks = rng.binomial(SHOTS, prob)
    stream.keep(key, clicks)
    return normalization * (2 * clicks / SHOTS - 1)


def read_clock(times, rng, stream):
    # Isolated blank dimer. The second independent copy gets a pi/2 native
    # exchange exposure at the terminal. Keep every load-left / load-right read.
    shots = 2**20
    cr = rng.binomial(shots, (1 + np.cos(times)) / 2)
    sr = rng.binomial(shots, (1 + np.sin(times)) / 2)
    stream.keep("clock_cos", cr)
    stream.keep("clock_sin", sr)
    cos, sin = 2 * cr / shots - 1, 2 * sr / shots - 1
    phase = np.unwrap(np.arctan2(sin, cos))
    slope, intercept = np.polyfit(np.arange(len(times)), phase, 1)
    return dict(
        dt_observed=float(slope),
        max_step=float(np.max(abs(np.diff(phase)))),
        residual_rms=float(
            np.sqrt(np.mean((phase - intercept - slope * np.arange(len(times))) ** 2))
        ),
        cos_sin_trials_per_hold=shots,
        cycles=float((phase[-1] - phase[0]) / (2 * np.pi)),
    )


def projection(trace, idx, dt):
    # Sub-bin peak interpolated from the same three Hann bins on every trial.
    powers = abs(np.fft.rfft(trace * np.hanning(len(trace))))
    if idx is None:
        start = len(powers) // 2
        idx = start + np.argmax(powers[start:-1])
    a, b, c = np.log(np.maximum(powers[idx - 1 : idx + 2], 1e-300))
    offset = (a - c) / (2 * (a - 2 * b + c))
    return float(2 * np.pi * (idx + np.clip(offset, -0.5, 0.5)) / (len(trace) * dt))


def spectroscopy(n, rng, sample_path=None):
    stream = CountStream(sample_path)
    m = n + 1
    u, e, _ = dual(n, K)
    end = min(8, max(4, (m + 1) // 4))
    channels = list(range(1, end + 1))
    pairs = [(1, b) for b in channels[1:]] + [(2, 3)] + [(2, b) for b in channels[3:]]
    spatial = np.sqrt(2 / (m + 1)) * np.sin(
        np.pi * np.outer(np.arange(1, m + 1), channels) / (m + 1)
    )
    times = DT * np.arange(SAMPLES)
    clock = read_clock(times, rng, stream)
    observed_times = clock["dt_observed"] * np.arange(SAMPLES)
    u_off, e_off, _ = dual(n, 0)
    co_off = np.cos(times[:, None] * e_off)
    co = np.cos(times[:, None] * e)
    # Physical two-time transfer arrays (t, source sine profile, position).
    propagated = np.stack([(co * (u.T @ f)) @ u.T for f in spatial.T], axis=1)
    off = np.stack([(co_off * (u_off.T @ f)) @ u_off.T for f in spatial.T], axis=1)
    controls = {}
    pair_sums = {}
    planes = {}
    for a, b in pairs:
        fa, fb = spatial[:, a - 1], spatial[:, b - 1]
        pa, pb = propagated[:, a - 1], propagated[:, b - 1]
        src_norm = wedge_norm(fa, fb)
        aa, ab, bb = pa @ fa, pa @ fb, pb @ fb
        key = f"p{a}_{b}_"
        trace = counts(aa * bb - ab**2, src_norm**2, rng, stream, key + "sum")
        freq = projection(trace, None, clock["dt_observed"])
        ca, cb = off[:, a - 1], off[:, b - 1]
        changed = counts(
            (ca @ fa) * (cb @ fb) - (ca @ fb) ** 2,
            src_norm**2,
            rng,
            stream,
            key + "changed",
        )
        win = np.hanning(SAMPLES)

        def amplitude(data, w):
            return float(
                2
                * abs(np.dot(win * data, np.exp(-1j * w * observed_times)))
                / win.sum()
            )

        controls[f"{a},{b}"] = dict(
            native_sum_amplitude=amplitude(trace, freq),
            changed_at_same_sum=amplitude(changed, freq),
            changed_at_two=amplitude(changed, 2),
        )
        pair_sums[a, b] = freq
        # Two receiver profiles, scanned over every wall position. These are
        # signed product-Pauli expectations, with the same all-outcome sampling.
        demod = np.hanning(SAMPLES) * np.exp(-1j * freq * observed_times)
        plane = []
        for ref, profile in [
            (fa, aa[:, None] * pb - ab[:, None] * pa),
            (fb, ab[:, None] * pb - bb[:, None] * pa),
        ]:
            norm = src_norm * (abs(ref).sum() - abs(ref))
            observed = counts(
                profile, norm, rng, stream, key + "scan" + str(len(plane))
            )
            phase = 2 * np.sum(demod[:, None] * observed, axis=0) / demod.real.size
            gphase = np.angle(phase[np.argmax(abs(phase))])
            plane.append(np.real(phase * np.exp(-1j * gphase)))
        planes[a, b] = np.linalg.qr(np.asarray(plane).T)[0]
    e1 = (pair_sums[1, 2] + pair_sums[1, 3] - pair_sums[2, 3]) / 2
    rows = []
    for i in channels:
        ps = [basis for pair, basis in planes.items() if i in pair]
        # The intersection of the two independently read two-mode planes singles
        # out one standing wave. Fit its interior three-point spatial recurrence.
        projector = sum(p @ p.T for p in ps)
        _, vec = eigh(projector, subset_by_index=[m - 1, m - 1])
        w = vec[:, 0]
        mid = w[1:-1]
        cos_k = np.dot(mid, w[:-2] + w[2:]) / (2 * np.dot(mid, mid))
        k = float(np.arccos(np.clip(cos_k, -1, 1)))
        energy = e1 if i == 1 else pair_sums[1, i] - e1
        lattice = float(np.sqrt(2 - 2 * np.cos(k)))
        rest = 2 * (K - H)
        speed2 = 4 * K * H
        omega_lattice = float(np.sqrt(rest**2 + speed2 * lattice**2))
        omega_bloch = float(np.sqrt(rest**2 + speed2 * k**2))
        tangent = rest + speed2 * lattice**2 / (2 * rest)
        rows.append(
            dict(
                spatial_launch=i,
                k=k,
                graph_gradient_wave_number=lattice,
                omega=energy,
                exact_lattice_forecast=omega_lattice,
                continuum_forecast=omega_bloch,
                matched_lattice_parabola=float(tangent),
            )
        )
    closures = [
        abs(pair_sums[2, b] - (pair_sums[1, b] + pair_sums[1, 2] - 2 * e1))
        for b in channels[3:]
    ]
    stream.finish(n)
    return dict(
        count_sha256=stream.hashes,
        pair_sum_lines={f"{a},{b}": v for (a, b), v in pair_sums.items()},
        withheld_pair_sum_max_error=float(max(closures)),
        clock=clock,
        same_spatial_source_changed_contrast=controls,
        modes=rows,
    )


def transport(route, j, seed):
    n = len(route)
    source, ref = n // 3, 1
    rng = np.random.default_rng(seed)
    # Wall j sees the right endpoint of dimer j-1 and the left of dimer j.
    cross = {i for i in range(1, n) if route[i - 1][1] // 12 != route[i][0] // 12}
    # Cells touched by an actual bond reader. Exclude source and ref cells.
    cells = {i: {route[i - 1][1] // 12, route[i][0] // 12} for i in range(1, n)}
    far = [
        i
        for i in range(1, n)
        if cells[i].isdisjoint(cells[ref] | cells[source])
        and min(abs(i - source), abs(i - ref)) >= 3
    ]
    times = (0, 2, 4, 8)
    data = []
    for t in times:
        ker = walls(n, j, t)
        # Expected wall density for any initial product of t=0 wall indicators.
        # Both fixed-parity preparations share their spectator when required.
        initial = np.zeros(n + 1)
        initial[-1] = 1 - n % 2
        first, second = initial.copy(), initial.copy()
        first[source], second[ref] = 1, 1
        p = 0.5 + ker @ (first - 0.5)
        q = 0.5 + ker @ (second - 0.5)
        clicks = rng.binomial(2**20, np.clip(p, 0, 1)) - rng.binomial(
            2**20, np.clip(q, 0, 1)
        )
        far_signal = max((abs((p - q)[i]) for i in far), default=0)
        data.append(
            dict(
                t=t,
                distant_max_density_difference=float(far_signal),
                distant_max_count_difference=int(
                    max((abs(clicks[i]) for i in far), default=0)
                ),
                differences=clicks.tolist(),
            )
        )
    return dict(
        source=source,
        reference=ref,
        intercell_links=sorted(cross),
        far_indices=far,
        trials_per_arm=2**20,
        flights=data,
    )


def run(sample_path=None):
    out = dict(
        native=native_checks(),
        physical_spin_check=[exact_check(i) for i in (5, 6)],
        dt=DT,
        samples=SAMPLES,
        shots_per_time_setting=SHOTS,
        boards=[],
    )
    for level in (0, 1, 2):
        board = federation(level)
        for seed in SEEDS:
            components = grow(board, seed)
            route = max((path for closed, path in components if not closed), key=len)
            obs = spectroscopy(
                len(route),
                np.random.default_rng([seed, level, 47]),
                sample_path if level == 0 and seed == SEEDS[0] else None,
            )
            out["boards"].append(
                dict(
                    level=level,
                    seed=seed,
                    seam_digest=seam_digest(board),
                    dimers=len(route),
                    cells=len(set(np.ravel(route) // 12)),
                    route=route,
                    rings=int(sum(x for x, _ in components)),
                    components=len(components),
                    observer=obs,
                    native_transport=transport(route, K, seed),
                    disabled_contrast_transport=transport(route, 0, seed),
                )
            )
    return out


def verify(actual, saved):
    """Check discrete words exactly and floating reductions to numeric precision."""
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for key in actual:
            verify(actual[key], saved[key])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for x, y in zip(actual, saved, strict=True):
            verify(x, y)
    elif isinstance(saved, float):
        np.testing.assert_allclose(actual, saved, atol=2e-9, rtol=2e-9)
    else:
        assert actual == saved, (actual, saved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--sample-counts",
        type=Path,
        help="Save raw representative clocks, pair lines and spatial scans as compressed NPZ",
    )
    args = parser.parse_args()
    result = run(args.sample_counts)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    else:
        verify(result, json.loads((HERE / "results.json").read_text()))
    for b in result["boards"]:
        print(b["level"], b["seed"], b["dimers"], b["cells"], b["observer"]["modes"])
