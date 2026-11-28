"""Read an ordinary thermal spectral weight with finite real-time Ramsey appointments."""

from fractions import Fraction as Q
import math

import numpy as np
from scipy.linalg import expm
from sympy import Matrix

from . import sector_thermal as thermal


def read(seed=722451):
    _, _, _, blocks, ops = thermal.apparatus()
    H = ops[0].sum(0) / 4
    beta = thermal.BETA
    rng = np.random.default_rng(seed)
    n = 2**32
    maximum = 1024
    k = np.arange(1, maximum + 1)
    risk = 0.002
    width = math.sqrt(2 * math.log(2 * len(blocks) / risk) / n)
    rows = []
    record_rng = np.random.default_rng(164307)
    for group in blocks:
        h = H[np.ix_(group, group)]
        diag = h.diagonal()
        off = abs(h).sum(1) - abs(diag)
        # Use integer row bounds and a margin beyond each physical energy.
        lo = math.floor(min(diag - off)) - 2
        hi = math.ceil(max(diag + off)) + 2
        span = hi - lo
        a = beta * span / math.pi
        end = math.exp(-beta * span)
        coeff = (2 * a / math.pi) * (1 - (-1.0) ** k * end) / (a * a + k * k)
        probs = np.r_[(1 - end) / (beta * span), coeff]
        # The Fourier cosine expansion of exp(-a*x) on [0,pi] gives
        # positive clock-time weights summing to one. The infinite tail
        # away from 0 and pi is bounded by Abel summation in each series.
        tail = (
            (2 * a / math.pi)
            * (1 + end)
            / ((maximum + 1) ** 2 + a * a)
            / math.sin(math.pi / span)
        )
        appointments = rng.multinomial(n, np.r_[probs, 1 - probs.sum()])
        base = expm((-1j * math.pi / span) * (h - lo * np.eye(len(group))))
        U = np.eye(len(group), dtype=complex)
        starts0 = record_rng.multinomial(appointments[0], np.full(len(h), 1 / len(h)))
        missed = record_rng.multinomial(appointments[-1], np.full(len(h), 1 / len(h)))
        all_starts = starts0 + missed
        signed_words = starts0.copy()
        signed = int(appointments[0])
        plus_total = signed
        for length, count in enumerate(appointments[1:-1], 1):
            U = U @ base
            starts = rng.multinomial(count, np.full(len(group), 1 / len(group)))
            # Interfere the ordinary clock with its parked reference arm.
            # Each sharp start has its own diagonal return amplitude.
            prob = (1 + U.diagonal().real) / 2
            assert min(prob) > -1e-10 and max(prob) < 1 + 1e-10
            bits = rng.binomial(starts, np.clip(prob, 0, 1))
            all_starts += starts
            signed_words += 2 * bits - starts
            plus = int(bits.sum())
            plus_total += plus
            signed += 2 * plus - int(count)
        assert np.max(abs(U.conj().T @ U - np.eye(len(group)))) < 1e-9
        z = signed / n
        err = width + tail + 1e-9
        assert z > err
        # Center each spectrum about its directly measurable ordinary
        # first moment. This cancels reversible energy offsets identically.
        offset = beta * (float(diag.mean()) - lo)
        rows.append(
            dict(
                trials=n,
                positive=plus_total,
                negative=n - int(appointments[-1]) - plus_total,
                unserved=int(appointments[-1]),
                energy_window=[lo, hi],
                clock_tick=math.pi / span,
                omitted_visibility_bound=tail,
                visibility=z,
                statistical_width=width,
                sharp_starts=all_starts.tolist(),
                sharp_signed_counts=signed_words.tolist(),
                centered_log_band=[
                    math.log(z - err) + offset,
                    math.log(z + err) + offset,
                ],
            )
        )
    return rows


def variance_meter(seed=445631):
    # Two independent sharp words from the same exchange component occupy
    # two clock arms. Their relative Ramsey visibility is |tr(U)/d|².
    # Quadrature weights cancel the next six even Taylor powers exactly.
    q = 7
    moments = Matrix([[(j + 1) ** (2 * i + 2) for j in range(q)] for i in range(q)])
    coeff = [Q(x) for x in moments.inv()[:, 0]]
    t = Q(1, 128)
    _, _, _, blocks, ops = thermal.apparatus()
    H = ops[0].sum(0) / 4
    n = 2**38
    rng = np.random.default_rng(seed)
    risk = 0.002
    err = math.sqrt(2 * math.log(2 * q * len(blocks) / risk) / n)
    result = []
    for group in blocks:
        h = H[np.ix_(group, group)]
        d = h.diagonal()
        rad = abs(h).sum(1) - abs(d)
        span = math.ceil(max(d + rad) - min(d - rad))
        step = expm(-1j * float(t) * h)
        U = np.eye(len(h), dtype=complex)
        estimate = 0.0
        error = 0.0
        counts = []
        for j, c in enumerate(coeff, 1):
            U = U @ step
            diag = U.diagonal()
            vis = np.real(diag[:, None] * diag.conj()[None, :]).ravel()
            starts = rng.multinomial(n, np.full(len(vis), 1 / len(vis)))
            plus = int(rng.binomial(starts, np.clip((1 + vis) / 2, 0, 1)).sum())
            # 1-cos(s*x) = (s*x)^2/2 + ..., with E[(E1-E2)^2]/2 = Var(E).
            average = 2 * plus / n - 1
            estimate += float(c / t**2) * (1 - average)
            trunc = Q((j * span) ** (2 * q + 2), math.factorial(2 * q + 2)) * t ** (
                2 * q + 2
            )
            error += float(abs(c) / t**2) * (err + float(trunc) + 2e-12)
            counts.append(plus)
        exact = float(((h * h).sum() / len(h)) - d.mean() ** 2)
        assert estimate - error <= exact <= estimate + error
        result.append(
            dict(positive=counts, bounds=[estimate - error, estimate + error])
        )
    return dict(
        trials_per_time=n,
        tick=str(t),
        coefficients=list(map(str, coeff)),
        joint_coverage=1 - risk,
        groups=result,
    )
