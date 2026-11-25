"""Recoil, a frozen-grating read, and two physical length holdouts on one doublet."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import triu
from scipy.stats import beta as beta_law

from .. import measure as source
from .guides import detour
from .impulse import certify

HERE = Path(__file__).resolve().parent
SEED = 8051223
STARTS = 2**50
KICKS = (-0.25, -0.125, 0, 0.125, 0.25)
CUTOFF = 0.01 / 200  # simultaneous Hoeffding/Clopper-Pearson failures across all banks


def clock(ev, spec):
    z = (ev - spec["center"]) * spec["tick"] / (2 * np.pi)
    n = spec["number_delays"]
    return np.exp(-1j * np.pi * (n - 1) * z) * np.sinc(n * z) / np.sinc(z)


def binomial_interval(raw, total):
    a = CUTOFF / 2
    low = beta_law.ppf(a, raw, total - raw + 1) if raw else 0.0
    high = beta_law.isf(a, raw + 1, total - raw) if raw < total else 1.0
    return [float(low), float(high)]


def length_read(count, spec, rng):
    replace, _, legs, _, _, q, small = detour(count)
    ev, v = eigh(small, driver="evd")
    beta = source.BASE[:, 0] @ (q[source.ONE] @ v)
    response = clock(ev, spec)
    nearest = int(np.argmin(abs(ev - spec["center"])))
    # Assemble a conservative spectral exclusion bound from the complete
    # numerical diagonalization of the accessible standard row. For invertible
    # V, a real eigenvalue of H lies within ||V^{-1}|| ||HV-VD|| of D.
    ortho = np.linalg.norm(v.T @ v - np.eye(len(ev)), ord="fro")
    resid = np.linalg.norm(small @ v - v * ev, ord="fro") / np.sqrt(1 - ortho)
    # Loose outward allowance also covers the dense product's summation error
    # via gamma_n |A| |V| with the Frobenius/operator product inequality.
    eps = np.finfo(float).eps
    gamma = (len(ev) + 2) * eps / (1 - (len(ev) + 2) * eps)
    resid += (
        gamma * np.linalg.norm(abs(small), ord="fro") * np.linalg.norm(v, ord="fro") * 4
    )
    assert ortho < 1e-8 and resid < 1e-6 and np.ptp(ev) * spec["tick"] < 2 * np.pi
    admission = np.linalg.norm(response * beta) ** 2 / 3
    encoded = int(rng.binomial(STARTS, 1 / 3))
    clock_events = int(rng.binomial(encoded, min(1.0, 3 * admission)))
    return dict(
        arm_steps=len(legs[0]) + 1,
        added_ports=list(replace.values()),
        nearest_row_energy=float(ev[nearest]),
        distance=abs(ev[nearest] - spec["center"]),
        exclusion_allowance=float(resid),
        source_probability=admission,
        starts=STARTS,
        source_refused=STARTS - encoded,
        clock_refused=encoded - clock_events,
        accepted=clock_events,
        probability_interval=binomial_interval(clock_events, STARTS),
    )


def energy_matchings(h):
    f = source.fixture()
    upper = triu(h, k=1).tocoo()
    groups = []
    for a, b in f["blank"] + f["marked"]:
        i, j = sorted([f["ports"].index(a), f["ports"].index(b)])
        mask = (upper.row // 64 == i) & (upper.col // 64 == j)
        rows, cols, hs = upper.row[mask], upper.col[mask], upper.data[mask]
        assert len(hs) in (32, 64) and np.unique(hs).size == 1
        occupied = {a, b}
        for used, edges in groups:
            if not used & occupied:
                break
        else:
            used, edges = set(), []
            groups.append((used, edges))
        used.update(occupied)
        edges.append((rows, cols, float(hs[0])))
    return [edges for _, edges in groups]


def energy_read(psi, groups, p_src, rng):
    # A -pi/2 relative phase on the second endpoint followed by the native
    # half contact maps spatial ± superpositions to the two load counters.
    # The marked active space is resolved by a full return contact after the
    # first read. On blank contacts ± cost 8/9, on even active marks −8/0.
    results = []
    support = np.array([-8, 8, 9, 0])
    for group in groups:
        weights = np.zeros(4)
        for a, b, hop in group:
            pplus = 0.5 * np.sum(abs(psi[a] + psi[b]) ** 2)
            if hop == -0.5:
                weights[1] += pplus
                weights[2] += 0.5 * np.sum(abs(psi[a] - psi[b]) ** 2)
            else:
                assert hop == -4
                weights[0] += pplus
        weights[-1] = 1 - sum(weights[:-1])
        assert min(weights) > -2e-13
        weights = np.maximum(weights, 0)
        n = int(rng.binomial(STARTS, p_src))
        bins = rng.multinomial(n, weights / sum(weights))
        results.append(
            dict(
                source_misses=STARTS - n,
                energy_bins=bins.tolist(),
                exact=float(support @ weights),
                estimate=float(support @ bins / n),
                radius=17 * np.sqrt(np.log(2 / CUTOFF) / (2 * n)),
            )
        )
    return dict(
        banks=results,
        estimate=sum(r["estimate"] for r in results),
        exact=sum(r["exact"] for r in results),
        radius=sum(r["radius"] for r in results),
    )


def current_read(psi, x, groups, p_src, rng):
    # A bare native half-contact (the energy analyzer with its relative dwell
    # omitted) makes the receiving-minus-sending load counts equal to
    # 2 Im(a* b). A full reverse contact distinguishes blocked marked rails.
    banks = []
    support = np.array([-8.0, -0.5, 0.0, 0.5, 8.0])
    for group in groups:
        weights = np.zeros(5)
        for a, b, hop in group:
            speed = -hop * (x[b] - x[a])
            assert np.all(speed == speed[0])
            v = float(speed[0])
            for value, amplitudes in [
                (v, 1j * psi[a] + psi[b]),
                (-v, psi[a] + 1j * psi[b]),
            ]:
                weights[support.tolist().index(value)] += 0.5 * np.sum(
                    abs(amplitudes) ** 2
                )
        weights[2] = 1 - sum(weights[[0, 1, 3, 4]])
        assert min(weights) > -2e-13
        weights = np.maximum(weights, 0)
        n = int(rng.binomial(STARTS, p_src))
        bins = rng.multinomial(n, weights / sum(weights))
        banks.append(
            dict(
                source_misses=STARTS - n,
                current_bins=bins.tolist(),
                exact=float(support @ weights),
                estimate=float(support @ bins / n),
                radius=16 * np.sqrt(np.log(2 / CUTOFF) / (2 * n)),
            )
        )
    return dict(
        banks=banks,
        estimate=sum(row["estimate"] for row in banks),
        exact=sum(row["exact"] for row in banks),
        radius=sum(row["radius"] for row in banks),
    )


def position_read(psi, x, p_src, rng):
    p = np.bincount(x, weights=abs(psi) ** 2)
    n = int(rng.binomial(STARTS, p_src))
    counts = rng.multinomial(n, p / p.sum())
    return dict(
        misses=STARTS - n,
        port_ordinal_counts=counts.tolist(),
        estimate=float(np.arange(len(counts)) @ counts / n),
        radius=float(x.max() * np.sqrt(np.log(2 / CUTOFF) / (2 * n))),
    )


def run():
    rng = np.random.default_rng(SEED)
    previous = json.loads((HERE.parent / "results.json").read_text())
    spec, idx = previous["clock"], previous["selected_standard_index"]
    ev, _, beta, q, v, h, actions = source.build()
    F = clock(ev, spec)
    initial = q @ (v @ (beta * F))
    adm = np.vdot(initial, initial).real
    psi = initial / np.sqrt(adm)
    z = q @ v[:, idx]
    cycle = actions[source.CY]
    other = (z[cycle] - z[cycle[cycle]]) / np.sqrt(3)
    Q = np.stack([z, other], axis=1)
    np.testing.assert_allclose(Q.T @ Q, np.eye(2), atol=1e-12)
    np.testing.assert_allclose(h @ Q, ev[idx] * Q, atol=4e-12)
    f = source.fixture()
    length = len(f["legs"][0])
    x = np.repeat([0, length + 1] + list(range(1, length + 1)) * 3, 64)
    vx = 1j * (h @ (x[:, None] * Q) - x[:, None] * (h @ Q))
    projected_current = float(np.max(abs(Q.T @ vx)))
    assert projected_current < 2e-12
    groups = energy_matchings(h)
    terms = triu(h, k=1).tocoo()
    delta = abs(x[terms.row] - x[terms.col])
    assert set(delta) == {1, 2}
    c = [
        -2
        * sum(
            terms.data[delta == d] * z[terms.row[delta == d]] * z[terms.col[delta == d]]
        )
        for d in [1, 2]
    ]
    kicks = []
    for k in KICKS:
        phase = np.exp(1j * k * x)
        columns = Q.T @ (phase[:, None] * Q)
        np.testing.assert_allclose(columns, columns[0, 0] * np.eye(2), atol=3e-12)
        state = phase * psi
        g = v.T @ (q.T @ state)
        np.testing.assert_allclose(q @ (v @ g), state, atol=1e-12)
        energies = energy_read(state, groups, adm / 3, rng)
        mean = np.vdot(state, h @ state).real
        assert abs(energies["exact"] - mean) < 1e-10
        velocity = 2 * np.vdot(x * state, h @ state).imag
        current = current_read(state, x, groups, adm / 3, rng)
        assert abs(current["exact"] - velocity) < 1e-10
        waits = [0, 0.125, 0.25, 1.0, 4.0, 16.0]
        positions = []
        for t in waits:
            wave = q @ (v @ (g * np.exp(-1j * (ev - ev[idx]) * t)))
            positions.append(dict(time=t, **position_read(wave, x, adm / 3, rng)))
        second = np.linalg.norm(F * g) ** 2
        admitted = int(rng.binomial(STARTS, adm / 3))
        passed = int(rng.binomial(admitted, min(1.0, second)))
        kicks.append(
            dict(
                k=k,
                energy=energies,
                velocity=float(velocity),
                position=positions,
                conditional_repeat_filter=float(second),
                repeat_accept=passed,
                repeat_reject=admitted - passed,
                unprepared=STARTS - admitted,
                measured_current=current,
                ideal_line_survival=float(abs(columns[0, 0]) ** 2),
            )
        )
    return dict(
        seed=SEED,
        original_starts_per_bank=STARTS,
        source_probability=adm / 3,
        original_pole=float(ev[idx]),
        phase_gate=certify(0.125),
        ideal_recoil_coefficients=c,
        ideal_inverse_longitudinal_curvature=1 / (c[0] + 4 * c[1]),
        projected_current=projected_current,
        pinned_control=abs(2 * np.vdot(x * state, h.diagonal() * state).imag),
        kicks=kicks,
        lengths=[length_read(count, spec, rng) for count in (0, 1, 2)],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = run()
    if args.write:
        (HERE / "readings.json").write_text(json.dumps(data, indent=2) + "\n")
    else:
        print(
            json.dumps(
                {
                    k: data[k]
                    for k in ("original_pole", "ideal_recoil_coefficients", "lengths")
                },
                indent=2,
            )
        )
