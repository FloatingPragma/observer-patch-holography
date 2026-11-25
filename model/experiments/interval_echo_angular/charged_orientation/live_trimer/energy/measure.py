"""Energy-resolved amplitudes and an independently counted intercell beat in the live trimer."""

from argparse import ArgumentParser
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from scipy.sparse.csgraph import connected_components
from sympy import Poly, symbols
from model.candidate_laws.interval_echo.openwave.row_11_spectrum.compact_parents import (
    metrology,
)

from .. import measure as trimer
from . import certificate

HERE = Path(__file__).resolve().parent
FLIGHT = 131072
SHOTS = 2**18
HALF = 512
STEP = 4096
COINS = tuple(product((0, 1), repeat=3))
IN = ((0, 2, 0, 1, 1, 1), (0, 0, 2, 1, 1, 1))
OUT = ((1, 1, 1, 0, 2, 0), (1, 1, 1, 0, 0, 2))


def dynamics(sign):
    h, ix = trimer.sector(sign, True, True)
    xs = [ix[x, (0, 0, 0)] for x in IN]
    _, labels = connected_components(h - np.diag(h.diagonal()), directed=False)
    comp = np.flatnonzero(labels == labels[xs[0]])
    assert labels[xs[0]] == labels[xs[1]]
    xs = np.searchsorted(comp, xs)
    h = h[np.ix_(comp, comp)]
    cyc = certificate.cyclic(h, int(xs[0]))
    other = certificate.cyclic(h, int(xs[1]))
    assert cyc["monic_charpoly_ascending"] == other["monic_charpoly_ascending"]
    cyc["second_rail_index"] = int(xs[1])
    cyc["second_rail_det_mod_prime"] = other["det_mod_prime"]
    if sign == 1:
        # Reflection interchanges the two written inside heads and reverses
        # the bridge coin. It removes the near-degenerate numerical mixing.
        rev = np.array([ix[w[3:] + w[:3], (c[2], 1 - c[1], c[0])] for w, c in ix])
        rev = np.searchsorted(comp, rev[comp])
        np.testing.assert_array_equal(h, h[np.ix_(rev, rev)])
        first = np.flatnonzero(np.arange(len(comp)) < rev)
        assert 2 * len(first) == len(comp)
        es = []
        vs = []
        par = []
        for s in (-1, 1):
            j = (
                np.eye(len(h))[:, first] + s * np.eye(len(h))[:, rev[first]]
            ) / np.sqrt(2)
            en, vec = eigh(j.T @ h @ j / 8)
            es.extend(en)
            vs.extend((j @ vec).T)
            par.extend([s] * len(en))
        order = np.argsort(es)
        e = np.array(es)[order]
        v = np.array(vs)[order].T
        mirror = np.array(par)[order]
    else:
        e, v = eigh(h / 8)
        mirror = np.zeros(len(e), int)
    source = v[xs, :].T
    dest = np.zeros((len(COINS), 2, len(e)))
    for a, c in enumerate(COINS):
        for i, w in enumerate(OUT):
            row = ix[w, c]
            j = np.searchsorted(comp, row)
            if j < len(comp) and comp[j] == row:
                dest[a, i, :] = v[j]
    amps = np.einsum("cnk,ki->kcni", dest, source)
    charges = np.array([sum(w[3:]) - 3 for w, c in ix])[comp]
    return h, cyc, e, v, mirror, source, amps, charges


def process(mats):
    bell = mats.transpose(0, 2, 1).reshape(8, 4) / np.sqrt(2)
    rho = bell.T @ bell.conj()
    pt = rho.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    return float(np.trace(rho).real), float(np.linalg.eigvalsh(pt)[0])


def bands(e, source, amp):
    nominal = {}
    for a, b in [(66, 67), (72, 73), (81, 84), (101, 102), (66, 84)]:
        mask = (e >= a) & (e < b)
        src = source[mask]
        moved = np.einsum("k,kcij->cij", np.exp(-1j * FLIGHT * e[mask]), amp[mask])
        mean, pt = process(moved)
        nominal[f"{a}-{b}"] = dict(
            lines=int(mask.sum()),
            input_weight_extrema=np.linalg.eigvalsh(src.T @ src).tolist(),
            ideal_coherent_mixed_arrival=mean,
            ideal_choi_pt_min=pt,
        )
    return nominal


def gap_read(e, amp, pair):
    # Independent, complete end-of-flight load reads at positively jittered
    # appointments; the two jitter offsets are uniform in [0,STEP].
    g = e[:, None] - e[None, :]
    window = np.sinc(g * STEP / (2 * np.pi)) ** 2
    times = STEP * np.arange(1, 2 * HALF + 1)
    prop = np.exp(-1j * e[:, None] * times)
    rng = np.random.default_rng(672619)
    tip = []
    for rail in (0, 1):
        src = amp[:, :, rail, 0]
        weights = (src @ src.T) * window
        tip.append(np.sum(prop.conj() * (weights @ prop), axis=0).real)
    gram = np.einsum("ica,jca->ij", amp[:, :, :, 0], amp[:, :, :, 0]) * window
    i, j = pair
    baseline = float(np.trace(gram))
    tone = 2 * float(gram[i, j])
    others = np.triu(gram, 1)
    others[i, j] = 0
    remaining = 2 * float(np.sum(abs(others)))
    clicks = []
    for p, q in zip(*tip, strict=True):
        probs = np.round([p, q, 1 - p - q], 11)
        clicks.append(rng.multinomial(SHOTS, probs).tolist())
    clicks = np.array(clicks)
    reads = clicks[:, :2].sum(axis=1) / SHOTS

    def fit(t, y):
        power = abs(np.fft.rfft(y - y.mean()))
        k = int(np.argmax(power[1:]) + 1)
        spacing = 2 * np.pi / (len(t) * STEP)

        def slope(f):
            design = np.stack([np.ones(len(t)), np.cos(f * t), np.sin(f * t)], axis=1)
            return np.linalg.norm(
                design @ np.linalg.lstsq(design, y, rcond=None)[0] - y
            )

        return minimize_scalar(
            slope,
            bounds=((k - 0.5) * spacing, (k + 0.5) * spacing),
            method="bounded",
            options={"xatol": 1e-15},
        ).x

    shot = np.sqrt(np.log(2 * len(times) / 0.01) / (2 * SHOTS))
    dominant = baseline + tone * np.cos((e[j] - e[i]) * times)
    observed_residual = float(max(abs(reads - dominant)))
    assert observed_residual <= remaining + shot
    train = fit(times[:HALF], reads[:HALF])
    held = fit(times[HALF:], reads[HALF:])
    return dict(
        shots_per_wait=SHOTS,
        step=STEP,
        stationary_mean=baseline,
        smoothed_main_cosine=tone,
        other_tone_absolute_bound=remaining,
        max_observed_residual=observed_residual,
        simultaneous_99_shot_bound=float(shot),
        independent_uniform_jitters=[0, STEP],
        first_and_second_tip_else=clicks.tolist(),
        first_half_frequency=float(train),
        held_half_frequency=float(held),
    )


def charge_read(e, v, source, charge):
    rng = np.random.default_rng(672620)
    out = []
    for t in (0, 65536, 131072):
        psi = v @ (np.exp(-1j * e * t) * source[:, 0])
        levels = np.arange(-2, 2)
        p = np.array([np.sum(abs(psi[charge == q]) ** 2) for q in levels])
        counts = rng.multinomial(SHOTS, np.round(p, 12))
        out.append(
            dict(
                wait=t,
                signed_load_bins=counts.tolist(),
                observed_net_charge=float(counts @ levels / SHOTS),
            )
        )
    return out


def run():
    cases = {}
    arrays = {}
    for par in (1, -1):
        h, cyc, e, v, mirror, source, amp, charges = dynamics(par)
        bounds, d, err = certificate.ordered_bound(h, e, v)
        cases[str(par)] = dict(
            cyclicity=cyc, ordered_certificate=bounds, ideal_bands=bands(e, source, amp)
        )
        # A simple line is represented in a unique reflection sector in
        # the even case. No frequency tolerance is used to merge levels.
        arrays[f"energies_{par}"] = e
        arrays[f"source_weights_{par}"] = source[:, :, None] * source[:, None, :]
        arrays[f"amplitudes_{par}"] = amp
        arrays[f"reflection_{par}"] = mirror
        if par == 1:
            ticks = bounds["energy_ticks"]
            idx = np.argmax(np.sum(source**2, axis=1) * ((e > 66) & (e < 67)))
            i = idx if abs(e[idx + 1] - e[idx]) < abs(e[idx] - e[idx - 1]) else idx - 1
            gap = int(d[i + 1] - d[i]) / ticks
            cases[str(par)].update(
                main_pair=[int(i), int(i + 1)],
                energies=e[i : i + 2].tolist(),
                main_gap_inclusion=[gap - 2 * float(err), gap + 2 * float(err)],
                main_pair_mirror=mirror[i : i + 2].tolist(),
                pair_cosine_coefficient=2
                * float(np.sum(amp[i, :, :, 0] * amp[i + 1, :, :, 0])),
                count_spectroscopy=gap_read(e, amp, (i, i + 1)),
                integrated_boundary_current=charge_read(e, v, source, charges),
            )
    x = symbols("x")
    polynomials = [
        Poly.from_list(
            c["cyclicity"]["monic_charpoly_ascending"][::-1], x, modulus=65521
        )
        for c in cases.values()
    ]
    common = polynomials[0].gcd(polynomials[1])
    assert common.degree() == 0
    clock = metrology.clock(np.random.default_rng(672621))
    unit = metrology.read(clock)
    return dict(
        cases=cases,
        flight=FLIGHT,
        cross_parity_mod_p_gcd_ascending=[
            int(z) % 65521 for z in common.all_coeffs()[::-1]
        ],
        dimer_clock=dict(
            angular_frequency=unit,
            starts=metrology.CLOCK_SHOTS,
            waits=metrology.WAITS.tolist(),
            counts=clock.tolist(),
        ),
    ), arrays


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opts = parser.parse_args()
    data, arrays = run()
    if opts.write:
        (HERE / "results.json").write_text(json.dumps(data, indent=1) + "\n")
        np.savez_compressed(HERE / "lines.npz", **arrays)
    else:
        print(json.dumps(data, indent=1))


if __name__ == "__main__":
    main()
