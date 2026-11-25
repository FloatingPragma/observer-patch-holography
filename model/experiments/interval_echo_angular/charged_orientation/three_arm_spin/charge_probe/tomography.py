"""Compose two observer calibrations and predict new input and receiver bases."""

import json

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..apparatus import native_graph
from ..measure import BASE, ONE, X, Y, Z
from ..reciprocal import channel
from .measure import FLIGHT, HERE, SHOTS, local_assay

PAULI = np.array([X, Y, Z])
INPUTS = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
ANALYZERS = np.array([[2, 1, 2], [-1, 2, 2], [2, -2, 1]]) / 3
CALIBRATION_SHOTS = 2**22


def bounds(count, n):
    # All training, calibration and holdout marginals fit in 1,024 tails.
    tail = 0.01 / 1024
    return np.array(
        [
            beta.ppf(tail, count, n - count + 1) if count else 0,
            beta.isf(tail, count + 1, n - count) if count < n else 1,
        ]
    )


def trained_terms():
    train = json.loads((HERE.parent / "reciprocal.json").read_text())["forward"]
    terms = []
    for page in train["observations"]:
        for axis, setting in enumerate(page["settings"]):
            for sign, count in zip((1, -1), setting[:2], strict=True):
                terms.append(
                    (
                        page["input_axis"],
                        page["input_sign"],
                        axis,
                        sign,
                        bounds(count, train["shots_per_setting"]),
                    )
                )
    return terms


def predict(terms, source, analyzer):
    """Interval of an original-start probability. Inputs are counts and settings."""
    total = np.zeros(2)
    for axis, sign, target_axis, target_sign, interval in terms:
        w = (
            (1 / 6 + sign * source[axis] / 2)
            * (1 / 3 + target_sign * analyzer[target_axis])
            / 2
        )
        total += w * (interval if w >= 0 else interval[::-1])
    return np.clip(total, 0, 1)


def density(b):
    return (np.eye(2) + np.einsum("i,ijk->jk", b, PAULI)) / 2


def overlaps(a, b):
    return bool(a[0] <= b[1] and b[0] <= a[1])


def run():
    terms = trained_terms()
    rng = np.random.default_rng(269011)
    rows = []
    readings = []
    for charge in (-1, 1):
        h = native_graph(charge)
        columns = np.zeros((h.shape[0], 2))
        columns[ONE] = BASE
        K, *_ = channel(expm_multiply(-1j * FLIGHT * h, columns), 0)
        for probe in (-1, 1):
            # Quadrature -pi/2 gives an especially different probe response
            # for opposite source signs. Both output probabilities are kept.
            probs = local_assay(charge, probe, 8, -np.pi / 2)
            calibration = rng.multinomial(CALIBRATION_SHOTS, probs)
            calibrated = bounds(calibration[0], CALIBRATION_SHOTS)
            readings.append(
                dict(charge=charge, probe=probe, counts=calibration.tolist())
            )
            for i, b in enumerate(INPUTS):
                evolved = sum(k @ density(b) @ k.conj().T for k in K) / 9
                for j, m in enumerate(ANALYZERS):
                    # The native probabilities are used only to generate a
                    # fresh finite record. The predictions above use no H.
                    bright = float(np.trace(density(m) @ evolved).real)
                    joint = bright * probs[0]
                    code = float(np.trace(evolved).real)
                    count = rng.multinomial(
                        SHOTS,
                        [
                            joint,
                            (code - bright) * probs[0],
                            code * probs[1],
                            1 - code * sum(probs[:2]),
                        ],
                    )
                    forecast = predict(terms, b, m) * calibrated
                    memoryless = predict(terms, b, np.zeros(3)) * calibrated
                    acquired = bounds(count[0], SHOTS)
                    rows.append(
                        dict(
                            charge=charge,
                            probe=probe,
                            input=i,
                            analyzer=j,
                            counts=count.tolist(),
                            forecast99=forecast.tolist(),
                            observed99=acquired.tolist(),
                            unpolarized99=memoryless.tolist(),
                            compatible=overlaps(forecast, acquired),
                            unpolarized_compatible=overlaps(memoryless, acquired),
                        )
                    )
    return dict(
        calibration_shots=CALIBRATION_SHOTS,
        shots=SHOTS,
        seed=269011,
        training="../reciprocal.json (forward, time 32, negative charge only)",
        inputs=INPUTS.tolist(),
        analyzers=ANALYZERS.tolist(),
        stocked_probe_calibration=readings,
        count_columns=[
            "probe_bright_code_plus",
            "probe_bright_code_minus",
            "probe_dark_code",
            "all_other_starts",
        ],
        tested=len(rows),
        covered=sum(r["compatible"] for r in rows),
        unpolarized_covered=sum(r["unpolarized_compatible"] for r in rows),
        observations=rows,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = run()
    if parser.parse_args().write:
        (HERE / "tomography.json").write_text(json.dumps(result, indent=2) + "\n")
    print({k: result[k] for k in ("tested", "covered", "unpolarized_covered")})
