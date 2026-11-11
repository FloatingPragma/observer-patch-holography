"""A staged measurement in plain numpy: two Villain faces on one seam with pinned outer links.
The one free angle eps carries the weight w(eps) w(-eps); the twist replaces w(eps) by w(eps + 2 pi).
Exact value: the theta ratio with c = 1/2 (exact.theta_ratio). This exercises the estimator of the Monte Carlo."""
import math

import numpy as np

from .fits import estimate


def ln_villain(Fl, g2):
    k = np.arange(-3, 4)[:, None]
    return np.log(np.exp(-(Fl[None, :] - 4 * math.pi * k) ** 2 / (2 * g2)).sum(0))


def staged_pinned_pair(g2, N=12, chains=4000, therm=60, sweeps=120, blocks=20, seed=11):
    rng = np.random.default_rng(seed)
    eps = np.zeros(chains)
    a = 0.5 / N
    nodes = np.zeros((N + 1, blocks, 5))
    step = 1.2 * math.sqrt(g2)
    per = chains // blocks

    def lw(e, lam):
        return (1 - lam) * ln_villain(e, g2) + lam * ln_villain(e + 2 * math.pi, g2) + ln_villain(-e, g2)
    for k in range(N + 1):
        lam = k / N
        acc = np.zeros((blocks, 5))
        for it in range(therm + sweeps):
            new = eps + step * rng.uniform(-1, 1, chains)
            take = np.log(rng.random(chains)) < lw(new, lam) - lw(eps, lam)
            eps = np.where(take, new, eps)
            if it >= therm:
                d = ln_villain(eps + 2 * math.pi, g2) - ln_villain(eps, g2)
                cols = np.stack([d, np.exp(a * d), np.exp(-a * d), np.exp(2 * a * d), np.exp(-2 * a * d)], 1)
                acc += cols[:per * blocks].reshape(blocks, per, 5).mean(1)
        nodes[k] = acc / sweeps
    return estimate(nodes)
