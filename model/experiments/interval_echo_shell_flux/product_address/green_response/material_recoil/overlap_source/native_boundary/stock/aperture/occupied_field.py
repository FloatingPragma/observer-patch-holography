"""Measure whether the extra source-to-current response is a scalar gradient."""

import argparse
from hashlib import sha256
import json
from math import factorial

import numpy as np
from mpmath import iv

from . import occupied_transfer as exposure
from . import occupied_release as flight

HERE = flight.HERE
T = 64
READS = 2**50
EPS = np.sqrt(flight.RISK / (2 * READS))
HSTEP = 1 / 512
COEFF = np.array(
    [1 / 280, -4 / 105, 1 / 5, -4 / 5, 0, 4 / 5, -1 / 5, 4 / 105, -1 / 280]
)
VERIFY = 2**58


def ratio(num, den):
    return flight.standard.quotient(num, den)


def shell_bands(df, signal, obs):
    # No direct written potential acts across an exterior edge. Recover its
    # common drop if the summed coherence supports division on both arms.
    if not len(df):
        return None
    norm = len(df)
    z = signal[:, df].sum(1)
    out = []
    for w in range(2):
        den = z[w] + np.array([-1, 1]) * norm * EPS
        m = obs[w, df].sum(0)
        out.append(ratio(m, den))
    return dict(
        coherences=z.tolist(), arm_drops=out, exact_zero=bool(np.all(obs[:, df] == 0))
    )


def probe(before, marks, mesh, r, old, rng):
    h0 = before / 2
    u, v = mesh
    n = len(r)
    sectors = np.array(flight.SECTORS)
    coeff = sectors @ (marks[:, u] - marks[:, v])
    cal = old["decoded"]["ordinary_energy_bands"]
    energy_error = max(
        max(abs(x - round(sum(pair) / 2)) for x in pair) for pair in cal.values()
    )
    assert energy_error < 0.001
    source = np.array(flight.sector_prior(old)) @ flight.CONDITIONAL
    sectors_waves = []
    densities = []
    specs = []
    for a, b in sectors:
        h = h0 + np.diag(a * marks[0] + b * marks[1])
        e, basis = np.linalg.eigh(h)
        amps = basis.T @ (np.ones(n) / np.sqrt(n))
        wave = basis @ (amps * np.exp(-1j * T * e))
        rho = np.outer(wave, wave.conj())
        densities.append(rho)
        sectors_waves.append(wave)
        specs.append((e, basis, amps))
        direct = ((h @ rho - rho @ h) - (h0 @ rho - rho @ h0))[u, v].real
        np.testing.assert_allclose(
            direct,
            (a * (marks[0, u] - marks[0, v]) + b * (marks[1, u] - marks[1, v]))
            * rho[u, v].real,
            atol=2e-14,
        )
    waves = np.array(sectors_waves)
    z = (waves[:, u] * waves[:, v].conj()).real
    avail = abs(waves[:, u]) ** 2 + abs(waves[:, v]) ** 2
    probs = np.stack([(avail + 2 * z) / 2, (avail - 2 * z) / 2, 1 - avail], axis=-1)
    assert probs.min() > -1e-13
    signal = np.zeros((2, len(u)))
    num = np.zeros((2, len(u), 2))
    raw = []
    for arm, w in enumerate(source):
        readings = []
        for k in range(len(u)):
            law = (w[:, None] * probs[:, k]).ravel()
            count = rng.multinomial(
                READS, np.maximum(law, 0) / np.maximum(law, 0).sum()
            ).reshape(-1, 3)
            freq = (count[:, 0] - count[:, 1]) / (2 * READS)
            signal[arm, k] = sum(freq)
            m = coeff[:, k] @ freq
            # The uncertain energy is weighted by the probability that the
            # defect was in either exposed endpoint arm of the fringe.
            met = (count[:, :2].sum() / READS + EPS) / 2
            err = (
                max(abs(coeff[:, k])) * EPS
                + energy_error * sum(abs(marks[:, u[k]] - marks[:, v[k]])) * met
            )
            num[arm, k] = m + np.array([-err, err])
            readings.append(count)
        raw.append(np.array(readings))
    drops = []
    for k in range(len(u)):
        arms = [
            ratio(num[b, k], signal[b, k] + np.array([-EPS, EPS])) for b in range(2)
        ]
        drops.append(
            None
            if any(q is None for q in arms)
            else flight.lab.interval_output(iv.mpf(arms[1]) - iv.mpf(arms[0]))
        )
    edges = {
        tuple(sorted((int(x), int(y)))): k
        for k, (x, y) in enumerate(zip(u, v, strict=True))
    }
    cycles = []
    for i, j in edges:
        for k in np.flatnonzero((h0[i] < 0) & (h0[j] < 0)):
            if k <= j:
                continue
            loop = [edges[i, j], edges[j, k], edges[i, k]]
            if any(drops[e] is None for e in loop):
                continue
            total = (
                iv.mpf(drops[loop[0]]) + iv.mpf(drops[loop[1]]) - iv.mpf(drops[loop[2]])
            )
            lo, hi = flight.lab.interval_output(total)
            cycles.append((max(lo, -hi), [i, j, int(k)], loop, [lo, hi]))
    best = max(cycles, key=lambda v: v[0])
    # Physically reassign the two conserved heads jointly to another material
    # copy in the same arm. Independent sampling of the measured marginals
    # leaves the actual ordinary fringe intact and removes its association.
    control = []
    for k in best[2]:
        bands = []
        counts = []
        for w in source:
            num_head = coeff[:, k]
            flag = w @ probs[:, k]
            dist = (w[:, None] * flag).ravel()
            sample = rng.multinomial(READS, dist / dist.sum()).reshape(-1, 3)
            counts.append(sample.tolist())
            fringe = (sample[:, 0] - sample[:, 1]) / (2 * READS)
            z = sum(fringe)
            zz = z + np.array([-EPS, EPS])
            dwell = (
                max(abs(num_head)) * EPS
                + energy_error
                * sum(abs(marks[:, u[k]] - marks[:, v[k]]))
                * (sample[:, :2].sum() / READS + EPS)
                / 2
            )
            band = ratio(num_head @ fringe + np.array([-dwell, dwell]), zz)
            if band is not None:
                assert band[0] <= w @ num_head <= band[1]
            bands.append(band)
        control.append(
            dict(edge=[int(u[k]), int(v[k])], counts=counts, arm_drops=bands)
        )
    assert np.max(abs((source[1] - source[0]) @ coeff)) < 2e-14
    shell = {}
    for radius in flight.CUTS:
        if radius == 0:
            continue
        ids = np.flatnonzero((np.minimum(r[u], r[v]) == radius) & (r[u] != r[v]))
        assert not coeff[:, ids].any()
        shell[str(radius)] = shell_bands(ids, signal, num)
    # Validate the selected loop on new material-only copies. The time slope
    # uses eight independent imaginary-coherence reads and the kinetic
    # subtraction uses local real interference, with no source labels.
    checks = []
    for k in best[2]:
        checks.append(
            independent(k, u[k], v[k], h0, specs, source, densities, signal[:, k], rng)
        )
    if all(x["increment"] is not None for x in checks):
        loop = (
            iv.mpf(checks[0]["increment"])
            + iv.mpf(checks[1]["increment"])
            - iv.mpf(checks[2]["increment"])
        )
        verified = flight.lab.interval_output(loop)
    else:
        verified = None
    assert verified is not None and max(verified[0], best[3][0]) <= min(
        verified[1], best[3][1]
    )
    return dict(
        reassigned_controls=control,
        intact_triangles=sum(
            all(tuple(sorted(e)) in edges for e in ((x, y), (y, z), (z, x)))
            for x, y, z in flight.lab.motion.cold.federation(0).faces()
        ),
        resolved_triangles=len(cycles),
        triangles_excluding_zero=sum(m > 0 for m, *_ in cycles),
        greatest_margin_triangle=best[1],
        predicted_circulation=best[3],
        material_only_circulation=verified,
        material_checks=checks,
        loop_head_and_fringe_counts=[raw[b][best[2]].tolist() for b in range(2)],
        complete_head_and_fringe_sha256=[
            sha256(row.astype("<i8").tobytes()).hexdigest() for row in raw
        ],
        exterior=shell,
    )


def read(rng, prob):
    # Each 2^50-copy block is exactly representable by the generator.
    samples = rng.multinomial(2**50, prob / prob.sum(), size=VERIFY // 2**50)
    return samples.sum(0).tolist()


def independent(index, u, v, h0, specs, source, densities, coherence, rng):
    d = np.zeros_like(h0)
    d[:, v] += h0[u, :]
    d[u, :] -= h0[:, v]
    operator = (d + d.T) / 2
    p = []
    width = 0.0
    for e, basis, amps in specs:
        width = max(width, float(e[-1] - e[0]) + 1e-8)
        times = T + HSTEP * np.arange(-4, 5)
        waves = basis[[u, v]] @ (amps[:, None] * np.exp(-1j * np.outer(e, times)))
        avail = (abs(waves) ** 2).sum(0)
        z = (waves[0] * waves[1].conj()).imag
        p.append(
            np.stack([(avail + 2 * z) / 2, (avail - 2 * z) / 2, 1 - avail], axis=-1)
        )
    p = np.einsum("bs,stq->btq", source, np.array(p))
    eps = np.sqrt(flight.RISK / (2 * VERIFY))
    trunc = (
        sum(abs(COEFF) * abs(np.arange(-4, 5)) ** 9)
        * width**9
        * HSTEP**8
        / (2 * factorial(9))
    )
    values = []
    banks = []
    for b, w in enumerate(source):
        raw = [read(rng, z) for z in p[b]]
        samples = np.array(raw) / VERIFY
        timed = COEFF @ (samples[:, 0] - samples[:, 1]) / (2 * HSTEP)
        rho = sum(v * mat.real for v, mat in zip(w, densities, strict=True))
        norm, prob = exposure.flags(operator, rho)
        kinetic = read(rng, prob)
        moment = norm * (np.array(kinetic) / VERIFY) @ [1, -1, 1, -1, 0]
        residual = -timed - moment
        err = eps * (sum(abs(COEFF)) / HSTEP + 2 * norm) + trunc + 2e-10
        values.append(
            ratio(
                residual + np.array([-err, err]), coherence[b] + np.array([-EPS, EPS])
            )
        )
        banks.append(dict(imaginary_flags=raw, kinetic_flags=kinetic))
    inc = (
        None
        if any(q is None for q in values)
        else flight.lab.interval_output(iv.mpf(values[1]) - iv.mpf(values[0]))
    )
    return dict(
        edge=[int(u), int(v)],
        observed_real_coherence=coherence.tolist(),
        independent_banks=banks,
        derivative_error_bound=trunc,
        arm_drops=values,
        increment=inc,
    )


def run():
    board, old, edges, rad, _ = flight.chart()
    rng = np.random.default_rng(631222)
    result = {}
    for tag, load in (("vacancy", 0), ("excess", 2)):
        h, marks = flight.apparatus(board, edges, load, -1, 1)
        result[tag] = probe(h, marks, np.where(np.triu(h, 1) < 0), rad, old, rng)
    return dict(
        sha256={
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("OCCUPIED_FIELD.md", "occupied_release.json", "pair_contact.json")
        },
        time=T,
        head_trials_per_edge=READS,
        independent_trials_per_meter=VERIFY,
        step=HSTEP,
        results=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    out = run()
    path = HERE / "occupied_field.json"
    if opt.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        flight.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(path.read_bytes()), out
        )
    for tag, row in out["results"].items():
        print(
            tag,
            row["greatest_margin_triangle"],
            row["predicted_circulation"],
            row["material_only_circulation"],
        )
