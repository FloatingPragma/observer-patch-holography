"""Local charge density on the winding guide, with no reference arm or energy subtraction."""

import argparse
import hashlib
import io
import json
import numpy as np
from scipy.optimize import minimize_scalar

from . import geometry, measure, phase

HERE, N, TIMES = measure.HERE, measure.N, phase.TIMES
PROTOCOL = "6de8823cc131eef7855c79cfbac880a28da1444137f6793ba37e5345b8b5ba10"
SHOTS = 2**20
FIRST = (69, 64, 57, 49, 42, 39)
TRAIN, HELD = np.arange(4, 35), np.arange(35, 68)


def states(j):
    v = phase.marked_source(j)
    mirror = v.reshape(-1, 2)[:, ::-1].ravel()
    s, a = v + mirror, v - mirror
    return np.column_stack([s / np.linalg.norm(s), a / np.linalg.norm(a)])


def samples(p, shots, rng):
    flat = p.reshape(-1, p.shape[-1])
    result = np.array([rng.multinomial(shots, row / row.sum()) for row in flat])
    return result.reshape(p.shape).astype(np.int32)


def regression(h, y):
    c = np.linalg.lstsq(h, y, rcond=None)[0]
    return c, np.mean((y - h @ c) ** 2)


def bounded(search, objective):
    i = np.argmin([objective(x) for x in search])
    return minimize_scalar(
        objective,
        bounds=(search[max(0, i - 1)], search[min(i + 1, len(search) - 1)]),
        method="bounded",
        options=dict(xatol=1e-13),
    ).x


def time_basis(t, omega):
    return np.column_stack((np.ones(len(t)), np.cos(omega * t)))


def spatial_basis(x, k):
    return np.column_stack((np.ones(len(x)), np.cos(k * x), np.sin(k * x)))


def decode(prob):
    y = prob[:, 2 * TRAIN]
    t = TIMES[::2]

    def loss(omega):
        return regression(time_basis(t, omega), y[::2])[1]

    base = bounded(np.linspace(0.0001, 2 * np.pi - 0.0001, 513), loss)
    coeff = regression(time_basis(t, base), y[::2])[0]
    aliases = [base, 4 * np.pi - base]

    # Quarter-step odd holds on training ports resolve the even-hold cosine alias.
    def alias_residual(f):
        return np.sqrt(np.mean((time_basis(TIMES[1::2], f) @ coeff - y[1::2]) ** 2))

    aliases.sort(key=alias_residual)
    freq, other = aliases

    def spatial_loss(k):
        return regression(spatial_basis(TRAIN, k), coeff[1])[1]

    k = bounded(np.linspace(0.0001, np.pi - 0.0001, 513), spatial_loss)
    profiles = regression(spatial_basis(TRAIN, k), coeff.T)[0]
    forecast = time_basis(TIMES[1::2], freq) @ (spatial_basis(HELD, k) @ profiles).T
    error = np.sqrt(np.mean((prob[1::2, 2 * HELD] - forecast) ** 2))
    amplitude = np.sqrt(np.mean(coeff[1] ** 2))
    return dict(
        frequency=freq,
        wave_number=k,
        other_temporal_alias=other,
        alias_rms=[alias_residual(freq), alias_residual(other)],
        held_rms=error,
        modulation_rms=amplitude,
        resolved=bool(amplitude > 0.0005),
        spatial_pass=bool(error < 0.00022),
        lattice_error=freq - np.sqrt(49 + 32 * np.sin(k / 4) ** 2),
        continuum_error=freq - np.sqrt(49 + 2 * k * k),
        folded_upper_error=freq - np.sqrt(49 + 32 * np.cos(k / 4) ** 2),
        critical_error=freq - 2 * abs(np.sin(k / 4)),
    )


def native(route, kappa):
    h, ports = geometry.microscopic(route, kappa)
    energy, basis = np.linalg.eigh(h / 8)
    return ports, energy, basis


def scans(raw, n8, n1, rng):
    results = []
    for kappa, data, indices in (
        (8, n8, (*FIRST, *(73 - j for j in FIRST))),
        (1, n1, FIRST),
    ):
        for j in indices:
            state = states(j).sum(axis=1) / np.sqrt(2)
            p = measure.release(*data, state, TIMES)
            y = samples(p, SHOTS, rng)
            raw[f"line_{kappa}_{j}"] = y
            row = decode(y / SHOTS)
            theory = decode(p)
            results.append(dict(kappa=kappa, source=j, sampled=row, exact=theory))
    subspace = states(57)
    h = (n8[2] * n8[1]) @ n8[2].T
    state = subspace @ np.linalg.eigh(subspace.T @ h @ subspace)[1][:, 0]
    prob = measure.release(*n8, state, TIMES)
    shots = samples(prob, SHOTS, rng)
    raw["one_band"] = shots
    control = dict(
        maximum_exact_temporal_range=float(np.ptp(prob, axis=0).max()),
        sampled=decode(shots / SHOTS),
    )
    return results, control


def packets(raw, data, rng):
    ports = data[0]
    cells = np.arange(N) // 2
    t = np.array([0, 6, 12, 18])
    out = []
    for q in (0.4, -0.4, 0.8, -0.8, 1.2, -1.2):
        a = measure.packet(ports, q, 1)
        b = measure.packet(ports, -q, -1)
        # The finite packets can overlap. Normalize both complementary coherent sources.
        waves = [(a + s * b) / np.linalg.norm(a + s * b) for s in (1, -1)]
        p = [measure.release(*data, state, t) for state in waves]
        incoherent = (measure.release(*data, a, t) + measure.release(*data, b, t)) / 2
        y = samples(np.array([*p, incoherent, incoherent]), SHOTS, rng)
        raw[f"packet_{q}"] = y
        mix = y[2:] / SHOTS
        coherent = y[:2] / SHOTS
        # Weight sign-preparations by their physical norms to implement the exact phase average.
        weights = [np.linalg.norm(a + s * b) ** 2 / 4 for s in (1, -1)]
        mean = np.einsum("i,its->ts", weights, coherent) @ cells
        velocity = np.polyfit(t, mean, 1)[0]
        forecast = np.mean(
            [measure.group_velocity(ports, state)[0] for state in (a, b)]
        )
        fourier = np.exp(-1j * 2 * abs(q) * cells[::2])
        fringe = (coherent[0, :, ::2] - coherent[1, :, ::2]) @ fourier
        null = (mix[0, :, ::2] - mix[1, :, ::2]) @ fourier
        out.append(
            dict(
                q=q,
                physical_signed_speed=velocity,
                predicted_speed=forecast,
                total_sign_weights=weights,
                fringe_modulus=abs(fringe).tolist(),
                fringe_phase=np.angle(fringe).tolist(),
                mixture_modulus=abs(null).tolist(),
                passed=bool(
                    np.sign(velocity) == np.sign(forecast)
                    and abs(velocity - forecast) < 0.03
                    and abs(fringe[0]) > 0.15
                    and max(abs(null)) < 0.01
                ),
            )
        )
    return out


def sharp(raw, data, klist):
    ports, energy, basis = data
    a = np.flatnonzero(ports % 2 == 0)
    start = ports[a] // 2
    rng = np.random.default_rng(306144)
    counts, exact = [], []
    for t in TIMES:
        amp = (basis[a] * np.exp(-1j * t * energy)) @ basis.T
        p = geometry.physical(amp, ports, N)
        counts.append(samples(p, 8192, rng))
        exact.append(p)
    y = np.array(counts)
    raw["sharp"] = y
    true = np.array(exact)[:, :, ::2]
    release = y[:, :, ::2] / 8192
    gap = energy[:, None] - energy[None, :]
    interband = (energy[:, None] > 4.5) != (energy[None, :] > 4.5)
    out = []
    for k in klist:
        obs = np.exp(1j * k * (np.arange(N // 2)[None, :] - start[:, None]))

        def corr(dist):
            return np.einsum("tia,ia->t", dist, obs).real / len(a)

        # The real correlation uses frequency differences with positive quadratic weights.
        dk = np.zeros(len(ports), dtype=complex)
        dk[a] = np.exp(1j * k * start)
        weights = abs((basis.T * dk) @ basis) ** 2 / len(a)
        omega = np.sqrt(49 + 32 * np.sin(k / 4) ** 2)
        kept = weights[interband].sum()
        away = weights[interband & (abs(abs(gap) - omega) > 0.025)].sum() / kept
        normal = np.cos(np.outer(TIMES, gap.ravel())) @ weights.ravel()
        np.testing.assert_allclose(normal, corr(true), atol=5e-13)
        collapsed = np.cos(np.outer(TIMES, gap[~interband])) @ weights[
            ~interband
        ] + kept * np.cos(omega * TIMES)
        out.append(
            dict(
                k=k,
                massive_frequency=omega,
                total_interband_weight=kept,
                interband_fraction_outside_point025=away,
                exact_single_line_error=max(abs(normal - collapsed)),
                sampled_single_line_error=max(abs(corr(release) - collapsed)),
                sampled_full_law_error=max(abs(corr(release) - normal)),
            )
        )
    return out


def run():
    assert hashlib.sha256((HERE / "DENSITY.md").read_bytes()).hexdigest() == PROTOCOL
    route, _ = measure.source.material.geometry()
    n8, n1 = native(route, 8), native(route, 1)
    raw = {}
    rng = np.random.default_rng(305144)
    lines, control = scans(raw, n8, n1, rng)
    movement = packets(raw, n8, rng)
    correlations = sharp(
        raw, n8, [x["sampled"]["wave_number"] for x in lines[: len(FIRST)]]
    )
    archive = io.BytesIO()
    np.savez_compressed(archive, **raw)
    data = archive.getvalue()
    out = dict(
        protocol_sha256=PROTOCOL,
        counts_sha256=hashlib.sha256(data).hexdigest(),
        laws=[8, 1],
        single_band=control,
        lines=lines,
        packets=movement,
        sharp=correlations,
    )
    return json.dumps(out, indent=2, sort_keys=True) + "\n", data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    write = parser.parse_args().write
    txt, data = run()
    for fn, blob in (("density.json", txt.encode()), ("density_counts.npz", data)):
        path = HERE / fn
        if write:
            path.write_bytes(blob)
        else:
            assert path.read_bytes() == blob, f"replay mismatch: {fn}"
    print(
        "Winding density acquisition written."
        if write
        else "Winding density acquisition verified."
    )
