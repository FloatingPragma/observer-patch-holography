"""A deterministic, unmeasured return cycle on all twelve material rails."""

import argparse
import json

import numpy as np
from scipy.linalg import schur

from . import circuit
from .measure import HERE, POWERS, SHOTS, infer, ratios

EPS = 1 / 128
SEED = 704319


def run():
    rng = np.random.default_rng(SEED)
    sites, _, _, U, back, F, p, _ = circuit.build()
    central = sites[:, 0]
    reflect = circuit.compile_word(
        [("phase", i, 0, np.pi) for i in range(12) if i not in central]
    )
    kick = circuit.compile_word([("phase", int(i), 0, -EPS) for i in central])
    step = reflect @ back @ kick @ U
    np.testing.assert_allclose(step.conj().T @ step, np.eye(12), atol=2e-12)
    # Decode the full *physical* unitary before comparing with the principal
    # angle formula. Normalized Schur phases avoid accumulating norm roundoff
    # at the largest power. The counted return is at the final read only.
    d, basis = schur(step, output="complex")
    np.testing.assert_allclose(d, np.diag(np.diag(d)), atol=1e-12)
    phases = np.angle(np.diag(d))
    inputs = np.eye(12)[:, central] @ F
    mix = abs(basis.conj().T @ inputs) ** 2
    low = abs(phases) < EPS * 1.01
    assert low.sum() == 3 and sum(abs(np.exp(1j * phases) + 1) < 2e-12) == 6
    theory = EPS / 2 - np.arcsin((1 - 2 * p) * np.sin(EPS / 2))
    np.testing.assert_allclose(np.sort(-phases[low]), np.sort(theory), atol=5e-13)
    panels = []
    for j in range(3):
        wt = mix[:, j]
        own = np.argmin(abs(phases + theory[j]))
        assert 1 - wt[own] < np.tan(EPS / 2) ** 2 + 1e-13
        # The character input can populate its conjugate high-phase partner.
        # Its weight is bounded without learning the unknown low rate first.
        rows = []
        for exp in POWERS:
            n = 2**exp
            z = wt @ np.exp(1j * n * phases)
            banks = []
            for angle in (0, np.pi / 2):
                y = z * np.exp(1j * angle)
                prob = np.maximum(
                    [abs(1 + y) ** 2 / 4, abs(1 - y) ** 2 / 4, (1 - abs(y) ** 2) / 2],
                    0.0,
                )
                banks.append(rng.multinomial(SHOTS, prob / prob.sum()).tolist())
            rows.append(dict(ticks=n, counts=banks))
        observed = infer(rows, EPS, 2 * np.tan(EPS / 2) ** 2)
        assert observed[0] < theory[j] < observed[1]
        panels.append(
            dict(
                low_mode_weight=float(wt[own]),
                low_radians_per_tick=theory[j],
                observed_rate99=observed,
                observations=rows,
            )
        )
    return dict(
        shots_per_quadrature=SHOTS,
        seed=SEED,
        common_exposure=EPS,
        all_physical_phases=np.sort(phases).tolist(),
        characters=panels,
        measured_low_ratios=ratios([s["observed_rate99"] for s in panels]),
        exact_low_ratios=(theory[::-1] / theory[2]).tolist(),
        maximum_nuisance_amplitude=2 * np.tan(EPS / 2) ** 2,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    data = run()
    if ap.parse_args().write:
        (HERE / "coherent.json").write_text(json.dumps(data, indent=2) + "\n")
    else:
        print(data["measured_low_ratios"])
