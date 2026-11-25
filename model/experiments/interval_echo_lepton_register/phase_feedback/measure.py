"""Resolve three successful-tick phase rates by native interferometry on one module."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta as beta_dist

from . import circuit
from .. import measure as marked

HERE = Path(__file__).resolve().parent
SHOTS = 2**42
SEED = 704315
POWERS = range(25)
DWELLS = (1 / 128, 1 / 16)
# Both signed quadratures have range [-1,1]. Count all 300 quadrature banks.
RAD = np.sqrt(2 * np.log(2 * len(POWERS) * 2 * 3 * len(DWELLS) / 0.01) / SHOTS)


def instrument():
    sites, _, _, _, _, F, _, _ = circuit.build()
    centers = sites[:, 0]
    tables = marked.terminal_labels()
    checks = sorted(set(range(12)) - set(centers))
    assert len(checks) == 9 and all(
        [v for v in table if v not in checks] == sorted(centers) for table in tables
    )
    meters = marked.routes(0)
    links = [
        (a, b) if a < 12 else (b, a)
        for a, b in marked.federation(0).external
        if min(a, b) < 12
    ]
    gate, ref = next(
        (a, b) for a, b in links if b not in [r["reference"] for r in meters]
    )
    # A blank external half swap splits the physical left/right vacancy; matched
    # reversal reads (1 ± e^{i phi} z^n)/2. The phase shifts of the actual
    # splitter cancel through the reversal. Internal gates have vacuum as idle.
    half = expm(-0.5j * np.pi * marked.seam_block(1))
    np.testing.assert_allclose(abs(half[:, 1]) ** 2, 0.5, atol=1e-13)
    z = 0.7 + 0.1j
    for angle in (0, np.pi / 2):
        w = half.conj().T @ np.diag([np.exp(1j * angle) * z, 1]) @ half[:, 1]
        np.testing.assert_allclose(
            abs(w[::-1]) ** 2,
            [
                abs(1 + np.exp(1j * angle) * z) ** 2 / 4,
                abs(1 - np.exp(1j * angle) * z) ** 2 / 4,
            ],
            atol=1e-13,
        )
    # Encode each character with the same real ternary split; the input
    # label chooses only its Fourier phases. A routed full swap connects the
    # splitter gateway and the first central port.
    for k in range(3):
        ops = circuit.preparation(centers, [1 / 3, 1 / 3, 1 / 3])
        ops += [
            ("phase", int(a), 0, 2 * np.pi * i * k / 3) for i, a in enumerate(centers)
        ]
        u = circuit.compile_word(ops)
        np.testing.assert_allclose(u[centers, centers[0]], F[:, k], atol=2e-13)
    return dict(
        reference_port=ref,
        signal_gate=gate,
        route=list(marked.route(gate, int(centers[0]))),
        nine_rejected_addresses=checks,
        pointer_references=[r["reference"] for r in meters],
    )


def infer(rows, bound, amplitude_error=0.0):
    # Use only original-start counts, the supplied common kick and the
    # positive-effect prior 0<=omega<=epsilon. A low-visibility bank can
    # contribute no unique phase inequality and remains in the raw record.
    lo, hi = 0.0, bound
    for row in rows:
        x, y = [(v[0] - v[1]) / SHOTS for v in row["counts"]]
        r = np.hypot(x, y)
        delta = np.sqrt(2) * RAD + amplitude_error
        if r <= delta:
            continue
        half = np.arcsin(delta / r)
        theta = np.arctan2(y, x)
        n = row["ticks"]
        lower = int(np.ceil((n * lo - theta - half) / (2 * np.pi)))
        upper = int(np.floor((n * hi - theta + half) / (2 * np.pi)))
        if lower > upper:
            raise ArithmeticError("Inconsistent measured phases.")
        if lower == upper:
            center = theta + 2 * np.pi * lower
            lo = max(lo, (center - half) / n)
            hi = min(hi, (center + half) / n)
    return [float(lo), float(hi)]


def bank(z, eps, rng):
    rows = []
    for exponent in POWERS:
        n = 2**exponent
        response = z**n
        outcomes = []
        for angle in (0, np.pi / 2):
            r = np.exp(1j * angle) * response
            p = np.array(
                [
                    abs(1 + r) ** 2 / 4,
                    abs(1 - r) ** 2 / 4,
                    max(0.0, (1 - abs(r) ** 2) / 2),
                ]
            )
            outcomes.append(rng.multinomial(SHOTS, p / p.sum()).tolist())
        rows.append(dict(ticks=n, counts=outcomes))
    rates = infer(rows, eps)
    theta = np.array(rates)
    recovered = np.sin(theta) / (np.sin(theta) + np.sin(eps - theta))
    return dict(
        phase_rate99=rates, recovered_effect99=recovered.tolist(), observations=rows
    )


def ratios(bounds):
    light = np.array(bounds[2])
    return [[1.0, 1.0], *(list(np.asarray(bounds[k]) / light[::-1]) for k in [1, 0])]


def run():
    rng = np.random.default_rng(SEED)
    sites, rot, theta, U, back, F, p, word = circuit.build()
    groups = []
    for eps in DWELLS:
        z = circuit.tick(-eps)
        readings = [bank(val, eps, rng) for val in z]
        groups.append(
            dict(
                common_exposure=eps,
                signal_dwell_phase=-eps,
                phase_rates=(-np.angle(z)).tolist(),
                exact_finite_phase_ratios=(np.angle(z)[::-1] / np.angle(z[2])).tolist(),
                recovered_effect_ratio99=ratios(
                    [r["recovered_effect99"] for r in readings]
                ),
                measured_phase_ratio99=ratios([r["phase_rate99"] for r in readings]),
                characters=readings,
            )
        )
    # Undoing the herald changes which physical operator is iterated. Without
    # it, unitarily conjugating the phase gate gives just the two imposed phases.
    full = (
        back
        @ np.diag(
            [np.exp(1j * DWELLS[0]) if i in sites[:, 0] else 1 for i in range(12)]
        )
        @ U
    )
    vals = np.linalg.eigvals(full)
    assert (
        sum(abs(vals - 1) < 3e-12) == 9
        and sum(abs(vals - np.exp(1j * DWELLS[0])) < 3e-12) == 3
    )
    flip = circuit.tick(np.pi)
    assert np.max(abs(np.angle(flip))) < 1e-14 and min(abs(flip)) > 0
    # Repeat the same drive with the existing finite exclusive-flag calibration
    # instead of the analytic 2/9 probability. The empirical pulse setting is
    # a calibration input and has its own shot uncertainty.
    archived = json.loads((HERE.parent / "results.json").read_text())
    mu = archived["measured_xor"]
    n = archived["shots_per_bank"]
    k = round(n * mu)
    bracket = beta_dist.ppf([0.005, 0.995], [k, k + 1], [n - k + 1, n - k])

    def ratio(angle):
        pp = np.sort(
            (np.sqrt(2) + 2 * np.cos(angle + 2 * np.pi * np.arange(3) / 3)) ** 2 / 36
        )
        return (pp / pp[0]).tolist()

    return dict(
        shots_per_quadrature=SHOTS,
        seed=SEED,
        simultaneous_hoeffding_radius=RAD,
        material=instrument(),
        physical_module=sites.tolist(),
        proper_rotation=rot,
        native_compiled_operations=len(word),
        exclusive_probability=theta,
        static_probabilities=p.tolist(),
        weak_limit_ratios=(p[::-1] / p[2]).tolist(),
        shared_dwell_banks=groups,
        phase_pi_control_per_tick=np.angle(flip).tolist(),
        unheralded_distinct_phases=[0, DWELLS[0]],
        measured_xor=mu,
        archived_xor99=bracket.tolist(),
        measured_xor_ratios=ratio(mu),
        ratios_at_xor_interval_endpoints=[ratio(a) for a in bracket],
        unforced_release=circuit.release(),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    d = run()
    if ap.parse_args().write:
        (HERE / "readings.json").write_text(json.dumps(d, indent=2) + "\n")
    else:
        print(
            d["weak_limit_ratios"],
            [g["measured_phase_ratio99"] for g in d["shared_dwell_banks"]],
        )
