"""Fourier content of same-charge single-band recoil on the independent 82-port material."""

import argparse
from hashlib import sha256
import io
import json
import numpy as np
from scipy.linalg import eigh, expm

from . import geometry, measure, preparation

HERE = measure.HERE
PROTOCOL = "da6a3358528dd67c5beba8675efcb01ae475962955f7aa3b06f5d41495baad25"
TIMES = 2.0 ** np.arange(-5, 16)
N_CLOCK, N_SPACE, N_METER = 2**20, 2**34, 2**38
EVENTS = 30000


def tolerance(n):
    return np.sqrt(np.log(2 * EVENTS / 0.01) / (2 * n))


def meter(raw):
    rng = np.random.default_rng(310144)
    terms = measure.source.contacts()
    estimates = {}
    for name in ("blank", "even"):
        d, h = np.array(terms[name]) / 8
        local = np.array([[d, h], [h, d]])
        times = (1 / 8, 3 / 16)
        probs = [abs(expm(-1j * local * t)[1, 0]) ** 2 for t in times]
        hit = np.array([rng.binomial(N_METER, p) for p in probs])
        raw["meter_" + name] = np.column_stack((hit, N_METER - hit))

        def interval(p, t):
            return (
                np.arcsin(
                    np.sqrt(np.clip(p + np.array([-1, 1]) * tolerance(N_METER), 0, 1))
                )
                / t
            )

        bounds = interval(hit[0] / N_METER, times[0])
        held = np.sin(bounds * times[1]) ** 2
        estimates[name] = dict(
            exact_swap=probs,
            band=bounds.tolist(),
            held_frequency_band=interval(hit[1] / N_METER, times[1]).tolist(),
            held_pass=bool(
                held[0] - tolerance(N_METER)
                < hit[1] / N_METER
                < held[1] + tolerance(N_METER)
            ),
        )
    return estimates


def harmonic_interval(n, a, b, L):
    # c_l are the power series coefficients of (1-z)^(1/2).
    rl, rh = a[0] / b[1], a[1] / b[0]
    c = [1.0]
    for i in range(1, 38):
        c.append(c[-1] * (i - 1.5) / i)
    lo, hi = 0.0, 0.0
    for order in range(30):
        scalar = 2 * c[n + order] * c[order]
        band = (
            np.array([b[0] * rl ** (2 * order + n), b[1] * rh ** (2 * order + n)])
            * scalar
        )
        lo += min(band)
        hi += max(band)
    tail = 2 * b[1] * rh ** (60 + n) / (1 - rh**2)
    alias = (
        4
        * b[1]
        * (rh ** (2 * L - n) + rh ** (2 * L + n))
        / ((1 - rh**2) * (1 - rh ** (2 * L)))
    )
    return [lo - tail - alias - 1e-13, hi + tail + alias + 1e-13]


def frequency(prob, rng):
    shots = rng.binomial(N_CLOCK, np.clip(prob, 0, 1))
    z = (2 * shots[0] / N_CLOCK - 1) + 1j * (2 * shots[1] / N_CLOCK - 1)
    angles = np.angle(z)
    eps = np.arcsin(2 * np.sqrt(2) * tolerance(N_CLOCK))
    phase = angles[0]
    assert abs(phase) + eps < np.pi
    margin = np.pi
    for value in angles[1:]:
        predicted = 2 * phase
        phase = value + 2 * np.pi * round((predicted - value) / (2 * np.pi))
        margin = min(margin, np.pi - abs(phase - predicted) - 3 * eps)
    assert margin > 0
    return phase / TIMES[-1], eps / TIMES[-1], margin, shots


def spatial(ref, wave, rng, L):
    counts = []
    for sign in (1, -1):
        state = (ref + sign * wave) / np.sqrt(2)
        p = abs(state) ** 2
        physical = p[:-2].copy()
        physical[:2] += p[-2:]
        counts.append(rng.multinomial(N_SPACE, physical / physical.sum()))
    cell = np.arange(1, L - 1)
    diff = (counts[0] - counts[1])[2 * cell] / N_SPACE
    errors = []
    for i in range(1, L):
        k = np.pi * i / L
        mat = np.column_stack((np.sin(k * cell), np.cos(k * cell)))
        fit = mat @ np.linalg.lstsq(mat, diff, rcond=None)[0]
        errors.append(np.sqrt(np.mean((diff - fit) ** 2)))
    allowed = np.flatnonzero(np.array(errors) <= 2 * tolerance(N_SPACE)) + 1
    return dict(
        best_index=int(np.argmin(errors)) + 1,
        accepted_indices=allowed.tolist(),
        best_rms=min(errors),
        runner_up_rms=sorted(errors)[1],
    ), counts


def run():
    assert sha256((HERE / "HARMONICS.md").read_bytes()).hexdigest() == PROTOCOL
    survey = json.loads((HERE / "geometry.json").read_text())
    route = max((r["route"] for r in survey["flights"] if r["level"] == 2), key=len)
    assert len(route) == 82
    h, _ = geometry.microscopic(route)
    h = h / 8
    L = len(h) // 2
    e, v = eigh(h)
    x = np.arange(len(h))
    rest = (-1.0) ** (x // 2) / np.sqrt(len(x))
    top = (-1.0) ** x / np.sqrt(len(x))
    np.testing.assert_allclose(h @ rest, 8 * rest, atol=1e-14)
    hs = []
    raw = {}
    local = meter(raw)
    rng = np.random.default_rng(309144)
    for i in range(1, L + 1):
        q = np.pi * i / L
        if i == L:
            val, wave = 9.0, top
        else:
            marked = np.zeros(2 * L)
            nodes = np.arange(1, L)
            aperture = np.sin(min(i, L - i) * np.pi * nodes / L)
            if i < L - i:
                aperture *= (-1.0) ** (nodes + 1)
            marked[1:-1] = np.repeat(aperture, 2)
            mirror = marked.reshape(-1, 2)[:, ::-1].ravel()
            s, a = marked + mirror, marked - mirror
            frame = np.column_stack([s / np.linalg.norm(s), a / np.linalg.norm(a)])
            val, vec = eigh(frame.T @ h @ frame)
            val, wave = val[-1], frame @ vec[:, -1]
        res = np.linalg.norm(h @ wave - val * wave)
        assert res < 3e-13
        if i < L:
            shape, raw[f"space_{i}"] = spatial(rest, wave, rng, L)
        else:
            shape = dict(
                best_index=L, accepted_indices=[L], known_stationary_endpoint=True
            )
        start = (rest + wave) / np.sqrt(2)
        analysis = np.stack(
            [(rest + wave) / np.sqrt(2), (rest - 1j * wave) / np.sqrt(2)]
        )
        for target in (start, *analysis):
            np.testing.assert_allclose(
                preparation.compile_state(target, measure.source.contacts()),
                target,
                atol=1e-12,
            )
        # Full evolution and one-projector read, independent of the frequency used by the decoder.
        signal = v @ (np.exp(-1j * np.outer(e, TIMES)) * (v.T @ start)[:, None])
        p = abs(analysis.conj() @ signal) ** 2
        measured, df, margin, raw[f"phase_{i}"] = frequency(p, rng)
        # A scalar band with a free offset and nonzero first/second hoppings, same reader.
        alt = 0.61 * (1 - np.cos(q)) - 0.09 * (1 - np.cos(2 * q))
        palt = (1 + np.array([np.cos(alt * TIMES), np.sin(alt * TIMES)])) / 2
        foil, _, _, raw[f"foil_{i}"] = frequency(palt, rng)
        hs.append(
            dict(
                index=i,
                q=q,
                shape=shape,
                native_recoil=val - 8,
                sampled_recoil=measured,
                phase_bound_99=df,
                unwrap_margin=margin,
                source_residual=res,
                foil_recoil=foil,
            )
        )
    a, b = local["blank"]["band"], local["even"]["band"]
    harm = []
    recoils = np.r_[0, [z["sampled_recoil"] for z in hs]]
    foil = np.r_[0, [z["foil_recoil"] for z in hs]]
    exact = np.r_[0, [z["native_recoil"] for z in hs]]
    for n in (3, 4, 5):
        coeff = np.cos(n * np.pi * np.arange(L + 1) / L) * 2 / L
        coeff[[0, -1]] /= 2
        dx = abs(coeff[1:]).sum() * df + 2e-10
        value = coeff @ recoils
        other = coeff @ foil
        predicted = harmonic_interval(n, a, b, L)
        harm.append(
            dict(
                harmonic=n,
                exact=coeff @ exact,
                sampled=value,
                sampling_interval=[value - dx, value + dx],
                max_numerical_allowance=2e-10,
                contact_forecast=predicted,
                foil=other,
                foil_interval=[other - dx, other + dx],
            )
        )
    out = dict(
        protocol_sha256=PROTOCOL,
        route=route,
        clock_exposures=TIMES.tolist(),
        simultaneous_events_upper=EVENTS,
        shots_per_phase=N_CLOCK,
        shots_per_spatial_histogram=N_SPACE,
        shots_per_meter=N_METER,
        meter=local,
        mode_rows=hs,
        harmonics=harm,
        all_spatial_certified=all(
            z["shape"]["accepted_indices"] == [z["index"]] for z in hs
        ),
        max_recoil_error=max(abs(z["sampled_recoil"] - z["native_recoil"]) for z in hs),
    )
    archive = io.BytesIO()
    np.savez_compressed(archive, **raw)
    data = archive.getvalue()
    out["counts_sha256"] = sha256(data).hexdigest()
    return json.dumps(out, indent=2, sort_keys=True) + "\n", data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    write = p.parse_args().write
    text, blob = run()
    for name, data in (
        ("harmonics.json", text.encode()),
        ("harmonics_counts.npz", blob),
    ):
        path = HERE / name
        if write:
            path.write_bytes(data)
        else:
            assert path.read_bytes() == data, f"Recoil harmonics mismatch: {name}"
    print(
        "Charge recoil harmonics written."
        if write
        else "Charge recoil harmonics verified."
    )
