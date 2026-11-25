"""Compare counted physical pulses with modular forecasts from two joint clock references."""

import argparse
import json
from itertools import product
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation
from model.experiments.interval_echo_clock.modular_reference import reference_reader
from model.experiments.interval_echo_clock.reference_exchange import measure as exchange

HERE = Path(__file__).resolve().parent
N = 2**29
PULSES = 16
ANGLE = np.pi / 16
X_BASIS = np.array([[1, 1], [1, -1]]) / np.sqrt(2)


def sample(rng, mass):
    mass = np.asarray(mass)
    np.testing.assert_allclose(mass.sum(), 1, rtol=0, atol=3e-11)
    if min(mass) < -1e-12:
        raise ValueError("Negative terminal probability")
    mass = np.round(np.maximum(mass, 0), 13)
    return rng.multinomial(N, mass / mass.sum())


def outcomes(rho, first, second):
    return np.array(
        [
            np.trace(rho @ np.kron(a, b)).real
            for a, b in product(
                [np.eye(2) - first, first], [np.eye(2) - second, second]
            )
        ]
    )


def acquire():
    board = federation(0)
    launches, settings, pulse = reference_reader.instruments(board)
    references = []
    for _, state in launches:
        orbit = [np.linalg.matrix_power(pulse, n) @ state for n in range(32)]
        references.append(sum(np.outer(v[:2], v[:2].conj()) for v in orbit) / 32)
    rho_a, rho_b = references
    herald = launches[1][0]
    first, shifted, third = [(u.conj().T @ e @ u)[:2, :2] for e, u in settings]
    # Exchange.gate uses the opposite ordering of the two localized words.
    full = exchange.gate(np.pi)
    rails = exchange.contacts(board, [0, 1])
    exchange.certify_route(rails, np.pi, full)
    u = full[np.ix_(exchange.SECTOR, exchange.SECTOR)][::-1, ::-1]
    source = np.kron(rho_a, rho_b)
    joined = (source + u @ source @ u.conj().T) / 2
    marginal = (rho_a + rho_b) / 2
    sources = {
        "shared_exchange": (herald, joined),
        "independent_exchange": (herald**2, np.kron(marginal, marginal)),
    }
    sampler = np.random.default_rng(912304)
    calibration = {}
    for key, (ready, rho) in sources.items():
        calibration[key] = sample(
            sampler, [1 - ready, *(ready * outcomes(rho, third, third))]
        ).tolist()
    e = np.diag(X_BASIS.T @ third @ X_BASIS).real
    response = np.vstack([1 - e, e])
    # Flight uses the native published fair tick and invariant blank hand.
    _, rates, _ = clock._primitives()
    np.testing.assert_allclose(list(rates.values()), 0.5, atol=1e-13)
    ids = [clock.law.SRC.index(word) for word in [(2, 0), (0, 2)]]
    h0, h1 = clock.law.local_term_parts()
    columns = np.eye(clock.law.DIM)[:, ids]
    for h in [h0, h1]:
        np.testing.assert_allclose(
            h @ columns, columns @ (columns.T @ h @ columns), atol=1e-13
        )
    probe = np.kron(np.outer(launches[0][1][:2], launches[0][1][:2].conj()), rho_a)
    counts = {}
    for law, h in {
        "native": h0 + clock.KAPPA * h1,
        "contrast_only": clock.KAPPA * h1,
    }.items():
        impulse = expm(-1j * ANGLE * columns.T @ h @ columns)
        for analyzer, effect in enumerate((first, shifted)):
            stack = []
            for n in range(PULSES + 1):
                hand = np.linalg.matrix_power(impulse, n)
                flight = np.kron(hand, hand)
                out = outcomes(flight @ probe @ flight.conj().T, effect, third)
                stack.extend(comb(PULSES, n) * (0.5**PULSES) * out)
            counts[f"{law}_{analyzer}"] = sample(sampler, stack).reshape(-1, 4).tolist()
    return dict(
        shots=N,
        response=response.tolist(),
        calibration=calibration,
        clock_counts=counts,
    )


def margin(n):
    # Fewer than 1,000 distinct real calibration and terminal bin probabilities.
    # Conditional on the observed count, shots are i.i.d. in the selected cohort.
    return np.sqrt(np.log(2000 / 0.01) / (2 * n)) + 3e-10


def reference(bank, response):
    raw = np.array(bank[1:]).reshape(2, 2)
    inverse = np.linalg.inv(response)
    weights = inverse @ (raw / raw.sum()) @ inverse.T
    error = margin(raw.sum()) * np.outer(abs(inverse).sum(1), abs(inverse).sum(1))
    low, high = weights - error, weights + error
    if (low <= 0).any():
        raise ValueError("Unresolved reference population")

    def gaps(w):
        return np.log(w[0] / w[1])

    def limits(lower, upper):
        return np.array([np.log(lower[0] / upper[1]), np.log(upper[0] / lower[1])])

    local = gaps(weights.sum(1))
    ll, lh = limits(low.sum(1), high.sum(1))
    if ll <= 0:
        raise ValueError("Unresolved marginal gap")
    gl, gh = limits(low, high)
    if (gl <= 0).any():
        raise ValueError("Unresolved joint gap")
    return weights, gaps(weights) / local, np.array([gl / lh, gh / ll])


def forecast(rates, bounds, n, analyzer, response):
    scale = n * ANGLE
    hit = (1 + np.cos(scale * rates + analyzer * np.pi / 2)) / 2
    outcomes = np.vstack([1 - hit, hit]) @ response.T / 2
    angle_unc = np.max(abs(bounds - rates), axis=0) * scale
    # Each cosine hit probability is 1/2-Lipschitz in its angle.
    radii = np.tile((angle_unc @ response.T) / 4, (2, 1))
    cov = (hit[0] - hit[1]) * (response[1, 0] - response[1, 1]) / 4
    cov_error = angle_unc.sum() * abs(response[1, 0] - response[1, 1]) / 8
    return outcomes, radii, cov, cov_error


def summarize(raw):
    response = np.array(raw["response"])
    answer = {}
    for ref, bank in raw["calibration"].items():
        weights, rates, bounds = reference(bank, response)
        rows = []
        for n, analyzer in product([4, 8, 12], [0, 1]):
            pred, err, cov, cov_err = forecast(rates, bounds, n, analyzer, response)
            quads = np.array(raw["clock_counts"][f"native_{analyzer}"][n]).reshape(2, 2)
            observed = quads / quads.sum()
            eps = margin(quads.sum())
            # Product error <= 4 eps for two [0,1] marginals obtained from two bins each.
            obs_cov = observed[1, 1] - observed[1].sum() * observed[:, 1].sum()
            rows.append(
                dict(
                    accepted=n,
                    analyzer=analyzer,
                    clock_frequencies=observed.tolist(),
                    forecast=pred.tolist(),
                    joint_radius=(err + eps).tolist(),
                    disjoint_entries=int((abs(pred - observed) > err + eps).sum()),
                    measured_covariance=float(obs_cov),
                    measured_covariance_error=float(5 * eps),
                    forecast_covariance=float(cov),
                    forecast_covariance_error=float(cov_err),
                )
            )
        answer[ref] = dict(
            x_basis_weights=weights.tolist(),
            relative_rates=rates.tolist(),
            relative_rate_limits=bounds.tolist(),
            held=rows,
        )
    return answer


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    output = p.parse_args().output
    raw = acquire()
    report = {"records": raw, "inference": summarize(raw)}
    if output:
        output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    else:
        from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
            compare,
        )

        compare(report, json.loads((HERE / "records.json").read_bytes()), "joint clock")
    for name, data in report["inference"].items():
        print(
            name,
            data["relative_rates"],
            [row["disjoint_entries"] for row in data["held"]],
        )


if __name__ == "__main__":
    main()
