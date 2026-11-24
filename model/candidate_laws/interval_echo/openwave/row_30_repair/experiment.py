"""Single-port return of a spatially uniform probe after native sharp repairs."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import connected_components, laplacian
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_28_acceleration.experiment import (
    history,
    layout,
)
from model.candidate_laws.interval_echo.openwave.row_30_formed.experiment import (
    synthesize,
)
from model.executable_model import federation

HERE = Path(__file__).resolve().parent
K = candidate.KAPPA
SHOTS = 2**22
TIMES = (1 / 16, 1 / 4)
# Two marginal intervals for each of two controls and two exposures on 27 boards.
TAIL = 0.01 / (2 * 2 * 2 * 2 * 27)


def check_signed_collars():
    """Native generator factors into a clock on the memory and identity on the tag."""
    law = candidate.native_law
    columns = [
        (1 + law.LABELS.index((coin, sign))) * law.N_SRC + law.SRC.index(load)
        for coin in (0, 1)
        for sign in (-2, 2)
        for load in ((0, 0), (1, 0), (0, 1))
    ]
    basis = np.eye(law.DIM)[:, columns]
    kinetic, contrast = law.local_term_parts()
    flip_term = np.array([[1.0, -1.0], [-1.0, 1.0]]) / 2
    mem = np.kron(np.eye(2), flip_term)
    np.testing.assert_allclose(
        kinetic @ basis, basis @ np.kron(mem, np.eye(3)), atol=1e-13
    )
    np.testing.assert_allclose(contrast @ basis, 4 * basis, atol=1e-13)
    for loads in ((2, 0), (0, 2)):
        written = law.echo()[:, law.SRC.index(loads)].reshape(-1, law.N_SRC)
        expected = np.zeros_like(written)
        for coin in (0, 1):
            expected[
                1 + law.LABELS.index((coin, loads[0] - loads[1])), law.SRC.index((1, 1))
            ] = 2**-0.5
        np.testing.assert_allclose(written, expected, atol=1e-14)


def interval(n):
    lo = 0 if n == 0 else beta.ppf(TAIL, n, SHOTS - n + 1)
    hi = 1 if n == SHOTS else beta.ppf(1 - TAIL, n + 1, SHOTS - n)
    return float(lo), float(hi)


def source(level, seed):
    n, edges, ranks, _, _ = layout(level)
    rng = np.random.default_rng(seed)
    word = 2 * rng.integers(0, 2, n, dtype=np.int8)
    epoch, signs, last = history(n, edges, ranks, word, 0.25, rng)
    assert np.array_equal(epoch != 0, abs(signs) == 2)
    written = np.bincount(edges[epoch > 0].ravel(), minlength=n)
    assert written.max() <= 1
    a, b = edges[epoch == 0].T
    graph = csr_matrix((np.ones(2 * len(a)), (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    count, label = connected_components(graph, directed=False)
    largest = np.argmax(np.bincount(label))
    uniform = []
    for c in range(count):
        sites = np.flatnonzero(label == c)
        if len(np.unique(written[sites])) == 1:
            uniform.append(
                [
                    int(len(sites)),
                    int(len(np.unique(sites // 12))),
                    int(written[sites[0]]),
                ]
            )
    sites = np.flatnonzero(label == largest)
    sub = graph[sites][:, sites]
    owners = np.flatnonzero(federation(level).degree[sites] == 5)
    launch = int(owners[0])
    psi, err = synthesize(
        np.ones(len(sites)) / np.sqrt(len(sites)),
        sub,
        federation(level).degree[sites],
        launch,
    )
    return (
        sub,
        written[sites],
        psi,
        dict(
            seed=seed,
            cells=federation(level).cells,
            home=int(sites[launch]),
            ports=len(sites),
            total_components=count,
            uniform_components=uniform,
            last_commit=last,
            zeros=int((written[sites] == 0).sum()),
            ones=int(written[sites].sum()),
            compilation_error_bound=bool(err < 1e-11),
        ),
    )


def measure(level, seed):
    graph, written, wave, panel = source(level, seed)
    spatial = laplacian(graph).tocsr()
    assert np.linalg.norm(spatial @ wave) < 2e-11
    rng = np.random.default_rng([30601, level, seed])
    records = []
    for law, coupling in (("registered", K), ("exchange", 0)):
        # All scalar energy offsets cancel in the return probability.
        generator = spatial / 2 - diags(coupling * written, dtype=float)
        variance = float(coupling**2 * np.var(written))
        calculated = generator @ wave
        calculated -= np.vdot(wave, calculated) * wave
        np.testing.assert_allclose(
            np.vdot(calculated, calculated).real, variance, atol=2e-12
        )
        for t in TIMES:
            evolved = expm_multiply(-1j * t * generator, wave)
            fidelity = abs(np.vdot(wave, evolved)) ** 2
            if not -2e-12 <= 1 - fidelity <= 1:
                raise ArithmeticError("invalid return probability")
            # The changed-law uniform wave is exactly stationary.
            miss = float(1 - fidelity) if coupling else 0.0
            hand = np.sin(t / 2) ** 2
            probabilities = np.outer([1 - hand, hand], [1 - miss, miss]).ravel()
            outcomes = rng.multinomial(SHOTS, probabilities.round(13))
            lower, upper = interval(int(outcomes[2:].sum()))
            tlo, thi = 2 * np.arcsin(np.sqrt([lower, upper]))
            mlo, mhi = interval(int(outcomes[1] + outcomes[3]))
            # For a spectral distribution of width W, the Taylor remainder is bounded
            # by E[(E-E')^4] t^4/24 <= W² Var(E) t^4/12.
            width = 6 + coupling
            forecast = [
                float(max(0, variance * tlo**2 - width**2 * variance * thi**4 / 12)),
                float(min(1, variance * thi**2)),
            ]
            compatible = mlo <= forecast[1] and mhi >= forecast[0]
            records.append(
                dict(
                    law=law,
                    nominal_exposure=t,
                    hand_return_counts=outcomes.tolist(),
                    variance=round(variance, 10),
                    measured_time=[tlo, thi],
                    measured_miss=[mlo, mhi],
                    predicted_miss=forecast,
                    compatible=bool(compatible),
                )
            )
    panel["measurements"] = records
    return panel


def run():
    check_signed_collars()
    return dict(
        shots=SHOTS,
        outcomes=["hand0_hit", "hand0_miss", "hand1_hit", "hand1_miss"],
        panels=[
            measure(level, seed)
            for level in range(3)
            for seed in [*range(280210 + level * 10, 280218 + level * 10), 301121]
        ],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = run()
    file = HERE / "results.json"
    serialized = json.dumps(report, indent=1, allow_nan=False) + "\n"
    if args.write:
        file.write_text(serialized)
    elif json.loads(file.read_text()) != report:
        raise ValueError("Archived record differs")
    if not all(m["compatible"] for p in report["panels"] for m in p["measurements"]):
        raise ValueError(
            "Short-time prediction contradicted by returned clock and detector records"
        )
    print("Verified 27 repair histories and 108 joint clock/home-port measurements.")


if __name__ == "__main__":
    main()
