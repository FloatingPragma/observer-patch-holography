import json
import numpy as np
from scipy.linalg import expm
from . import orbital


def test_complete_native_spectrum_and_torus_curvature():
    d = 2 * np.pi / 10
    f = orbital.modes()
    for b, twist in ((0, 0), (1 / 32, 0), (0, 1 / 16)):
        h = orbital.generator(b, "registered", twist)
        r2 = orbital.fixture()[2]
        e = np.array(
            [
                -sum(
                    np.cos(
                        d * m * s
                        - orbital.CHARGE * b * r2 * np.sin(s * d) / 2
                        - d * twist * s
                    )
                    for s in (1, 2)
                )
                for m in range(10)
            ]
        )
        transformed = f.conj().T @ h @ f
        np.testing.assert_allclose(
            transformed, np.diag(e + np.mean(np.diag(transformed) - e)), atol=3e-12
        )
    report = json.loads((orbital.HERE / "orbital.json").read_text())
    for c in report["curvature"]:
        assert c["limit_lower"] < 5 * d**2 < c["limit_upper"]
        for row in c["normalized"]:
            flux = row["field"] * r2
            splitting = 2 * sum(
                np.sin(s * d) * np.sin(flux * np.sin(s * d) / 2) for s in (1, 2)
            )
            expected = abs(splitting / (flux * 5 * d**2))
            assert 0.5 < row["lower"] < expected < row["upper"] < 0.8


def test_counts_retain_missing_outcomes_and_source_sign():
    rng = np.random.default_rng(7)
    # Loss into a third port leaves an explicit complement in both quadratures.
    a, b, c = np.eye(3)
    h = np.outer(b, c) + np.outer(c, b)
    r = orbital.pair_counts(expm(-1j * np.pi * h / 2), a, b, rng)
    for counts in r:
        assert sum(counts) == orbital.SHOTS
        assert 0.499 < counts[2] / orbital.SHOTS < 0.501
    report = json.loads((orbital.HERE / "orbital.json").read_text())
    p = report["cases"][0]
    data = p["flux"]
    assert data[0]["splitting"] > 0 > data[-1]["splitting"]
    assert 2.8 < p["gaps"][1]["gap"] / p["gaps"][0]["gap"] < 2.9
    assert data[-1]["upper"] < 0.72 < report["cases"][-1]["flux"][-1]["lower"]
