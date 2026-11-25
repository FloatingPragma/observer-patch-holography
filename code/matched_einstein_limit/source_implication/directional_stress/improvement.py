"""A bounded head record and two conserved effective stresses with identical flat rest readings."""

import argparse
import json
from fractions import Fraction

import numpy as np
import sympy as sp
from mpmath import iv

from ..full_algebra import certificate as full
from ..null_optics.certificate import bounded_sources
from . import bridge


SIGMA = [
    sp.Matrix([[0, 1], [1, 0]]),
    sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    sp.diag(1, -1),
]


def witness():
    """Contrast of one coincidence analyzer aligned with the two reference marginals."""
    iv.dps = 30
    states, radius = bounded_sources(sp.eye(4) / 4)["state"]
    trace_state = sum(states, sp.zeros(4))

    def mean(a, b):
        return full.ri_band(full.inner(trace_state, sp.kronecker_product(a, b)), radius)

    left = [mean(p, sp.eye(2)) for p in SIGMA]
    right = [mean(sp.eye(2), p) for p in SIGMA]
    ln = iv.sqrt(sum(v**2 for v in left))
    rn = iv.sqrt(sum(v**2 for v in right))
    assert ln.a > 0 and rn.a > 0
    reference = sum(
        left[i] * right[j] * mean(a, b)
        for i, a in enumerate(SIGMA)
        for j, b in enumerate(SIGMA)
    ) / (ln * rn)
    product = ln * rn
    contrast = product - reference
    assert contrast.a > 0.19
    states = bridge.tangents()
    raw = states["reference"].sum(axis=0)
    sigma = [np.asarray(p, dtype=complex) for p in SIGMA]
    ell = np.array([np.trace(raw @ np.kron(p, np.eye(2))).real for p in sigma])
    r = np.array([np.trace(raw @ np.kron(np.eye(2), p)).real for p in sigma])
    a = sum(v * p for v, p in zip(ell / np.linalg.norm(ell), sigma))
    b = sum(v * p for v, p in zip(r / np.linalg.norm(r), sigma))
    w = np.kron(a, b)
    effect = np.array([w, w])
    reads = bridge.evaluate(effect, states)
    assert abs(reads["ordinary"]) < 1e-12 and reads["pairing"] in contrast
    return dict(
        left_length=str(ln),
        right_length=str(rn),
        reference=str(reference),
        unpaired=str(product),
        matched_pairing_response=str(contrast),
        matched_ordinary_response_exact=0,
        reference_analyzer_eigenvalues=[-1, 1],
        nominal_source=reads,
    )


def counts():
    """Estimate calibration axes and both joint moments in independent copies."""
    seed, n = 1462738, 2**42
    rng = np.random.default_rng(seed)
    states = bridge.action.preparations()
    eye = np.eye(2)
    sigma = [np.asarray(p, dtype=complex) for p in SIGMA]
    margin = iv.sqrt(2 * iv.ln(2 * 24 / iv.mpf(".01")) / n)
    observations = {}

    def read(label, a, b, state):
        mean = np.trace(state @ np.kron(a, b)).real
        plus = int(rng.binomial(n, (mean + 1) / 2))
        observations[label] = plus
        return (iv.mpf(2 * plus) / n - 1) + iv.mpf([-1, 1]) * margin

    ref = states["equilibrium"].sum(axis=0)
    left = [read(f"left_{i}", p, eye, ref) for i, p in enumerate(sigma)]
    right = [read(f"right_{i}", eye, p, ref) for i, p in enumerate(sigma)]
    norm = iv.sqrt(sum(v**2 for v in left) * sum(v**2 for v in right))
    moments = {}
    for name in ("equilibrium", "pairing"):
        rho = states[name].sum(axis=0)
        moments[name] = (
            sum(
                left[i] * right[j] * read(f"{name}_{i}_{j}", a, b, rho)
                for i, a in enumerate(sigma)
                for j, b in enumerate(sigma)
            )
            / norm
        )
    response = moments["pairing"] - moments["equilibrium"]
    assert response.a > 0.19
    return dict(
        seed=seed,
        trials_per_setting=n,
        plus_counts=observations,
        simultaneous_coverage=0.99,
        contrast=str(response),
    )


def twice_d(f, axis):
    return np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)


def grid(example, rest_projected=False):
    """Integer numerators for 8 I, for F = polynomial/2 in the local stencil."""
    n = 13
    t, x, y, z = np.indices((n,) * 4) - n // 2
    cutoff = np.maximum.reduce([abs(t), abs(x), abs(y), abs(z)]) <= 4
    polynomial = x * x - y * y if example == "quadrupole" else 2 * t * x
    twice_f = cutoff * polynomial
    metric = (-1, 1, 1, 1)
    hess = [[twice_d(twice_d(twice_f, i), j) for j in range(4)] for i in range(4)]
    box = sum(metric[i] * hess[i][i] for i in range(4))
    tensor = [
        [hess[i][j] - (metric[i] * box if i == j else 0) for j in range(4)]
        for i in range(4)
    ]
    if rest_projected:
        laplacian = sum(hess[i][i] for i in (1, 2, 3))
        for i in range(4):
            tensor[0][i] = tensor[i][0] = np.zeros_like(twice_f)
        for i in (1, 2, 3):
            tensor[i][i] = hess[i][i] - laplacian
    for j in range(4):
        assert not np.any(sum(metric[i] * twice_d(tensor[i][j], i) for i in range(4)))
        assert not np.any(tensor[0][j].sum(axis=(1, 2, 3)))
    null_sum = 0
    for direction in (-1, 1):
        kk = tensor[0][0] + 2 * direction * tensor[0][1] + tensor[1][1]
        for offset in range(n):
            ray = kk[np.arange(n), (offset + direction * np.arange(n)) % n]
            null_sum = max(null_sum, int(np.max(np.abs(ray.sum(axis=0)))))
    if not rest_projected:
        assert null_sum == 0
    local = np.array([[a[(6,) * 4] for a in row] for row in tensor])
    assert np.all(local % 8 == 0)
    local //= 8
    # Every point in the inner five-site-wide space-time patch has the same
    # values because the entire two-step differentiation stencil lies in the cutoff.
    assert all(
        np.all(a[4:9, 4:9, 4:9, 4:9] == 8 * local[i, j])
        for i, row in enumerate(tensor)
        for j, a in enumerate(row)
    )
    v = Fraction(1, 2)
    moving = (int(local[0, 0]) + 2 * v * int(local[0, 1]) + v**2 * int(local[1, 1])) / (
        1 - v**2
    )
    return dict(
        local_orthonormal_tensor=local.tolist(),
        local_energy_at_boost_half=str(moving),
        max_divergence_numerator=0,
        max_slice_charge_numerator=0,
        max_periodic_null_sum_numerator=null_sum,
        max_rest_energy_or_flux_numerator=int(np.max(np.abs(tensor[0]))),
    )


def run():
    quadrupole, flux = grid("quadrupole"), grid("mixed")
    rest = grid("quadrupole", rest_projected=True)
    assert rest["max_rest_energy_or_flux_numerator"] == 0
    assert rest["local_orthonormal_tensor"] == quadrupole["local_orthonormal_tensor"]
    assert quadrupole["local_orthonormal_tensor"] == np.diag([0, 1, -1, 0]).tolist()
    assert quadrupole["local_energy_at_boost_half"] == "1/3"
    assert flux["local_energy_at_boost_half"] == "4/3"
    return dict(
        source=witness(),
        simulated_counts=counts(),
        finite_chart={
            "shape": [13] * 4,
            "signature": [-1, 1, 1, 1],
            "quadrupole": quadrupole,
            "time_space_mixed": flux,
            "rest_projected_quadrupole": rest,
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (bridge.HERE / "improvement.json").write_text(result)
    else:
        print(result, end="")
