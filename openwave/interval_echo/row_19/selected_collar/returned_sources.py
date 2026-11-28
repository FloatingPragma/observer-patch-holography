"""Infer regional connected modular responses from incoming and returned Pauli records."""

import numpy as np

from . import reference_read as heads

# Union bound also covers the timing decoder. The family contains fewer than
# 100,000 bounded estimates, including parity flags and signed Pauli scores.
FAMILY = 100_000
CONTRASTS = (
    ("ordinary", "reference", 1.0),
    ("turn_plus", "turn_minus", 0.5),
    ("source_plus", "source_minus", 1.0),
)


def probability_error(p, n):
    a = np.log(2 * FAMILY / 0.01) / n
    x = np.minimum(p, 1 - p)
    return np.sqrt(2 * a * (2 * x + 4 * a)) + 2 * a / 3 + 3e-14


def observed(tape, n, returned):
    C = np.asarray(tape)
    assert C.shape == (3, 3, 9 if returned else 8) and C.dtype.kind in "iu"
    assert np.all(C >= 0) and np.all(C.sum(axis=2) == n)
    p = C[:, :, :8].reshape(3, 3, 2, 2, 2) / n
    parity = p[0, 0].sum(axis=(1, 2))
    moments = np.empty((2, 4, 4))
    moments[:, 0, 0] = parity
    moments[:, 1:, 0] = np.einsum("ibst,s->bi", p[:, 0], heads.SIGNS)
    moments[:, 0, 1:] = np.einsum("jbst,t->bj", p[0], heads.SIGNS)
    moments[:, 1:, 1:] = np.einsum("ijbst,s,t->bij", p, heads.SIGNS, heads.SIGNS)
    a = np.log(2 * FAMILY / 0.01) / n
    # Every Pauli score is zero outside the indicated parity and success.
    # On the common event its variance is bounded by this parity upper limit.
    delta = np.sqrt(2 * a * (parity + probability_error(parity, n))) + 4 * a / 3 + 3e-14
    blocks = np.einsum("bij,ijkl->bkl", moments, heads.PRODUCTS) / 4
    return blocks, delta


def score(tapes, n, returned, dominance=None):
    stock = {key: observed(t, n, returned) for key, t in tapes.items()}
    ref, delta = stock["reference"]
    result = np.zeros((3, 2))
    norms = np.zeros((3, 2))
    unresolved = 0
    for b, X in enumerate(ref):
        err = 2 * delta[b]  # HS and operator error by orthogonality of the Pauli basis.
        head_pair = X.reshape(2, 2, 2, 2)
        p = X.trace().real
        Xl = np.trace(head_pair, axis1=1, axis2=3)
        Xr = np.trace(head_pair, axis1=0, axis2=2)
        mats = (X, Xl, Xr)
        radii = (err, 2 * err, 2 * err)
        ev = [np.linalg.eigh(M) for M in mats]
        good = p > delta[b] and all(
            vals[0] > e for (vals, _), e in zip(ev, radii, strict=True)
        )
        if good:
            logs = [(V * np.log(vals)) @ V.conj().T for vals, V in ev]
            K = (
                np.kron(logs[1], np.eye(2))
                + np.kron(np.eye(2), logs[2])
                - np.eye(4) * np.log(p)
                - logs[0]
            )
            dk = sum(
                e / (vals[0] - e) for (vals, _), e in zip(ev, radii, strict=True)
            ) + delta[b] / (p - delta[b])
            # Uncertainty in expectation of the nominal operator directly in
            # the measured signed Pauli coordinates.
            coefficients = abs(np.einsum("ijkl,lk->ij", heads.PRODUCTS, K)).sum() / 4
        else:
            unresolved += 1
        for k, (upper, lower, scale) in enumerate(CONTRASTS):
            xp, ep = stock[upper]
            xm, em = stock[lower]
            D = (xp[b] - xm[b]) * scale
            ed = (ep[b] + em[b]) * scale
            if not returned:
                assert good
                norms[k, b] = (np.linalg.norm(D, ord=2) + 2 * ed) / (ev[0][0][0] - err)
            if good:
                val = float(np.trace(D @ K).real)
                width = (
                    coefficients * ed
                    + (np.linalg.norm(D, ord="nuc") + 4 * ed) * dk
                    + 1e-15
                )
                low, high = val - width, val + width
            else:
                low, high = -np.inf, np.inf
            if returned:
                # Positivity of the measured native instrument carries
                # -aX <= D <= aX through each conditional output channel.
                # The normalized joint and both marginal entropies sum
                # to at most 4 log 2. This bound also covers singular blocks.
                bound = 4 * np.log(2) * dominance[k, b] * (p + delta[b]) + 1e-15
                low, high = max(low, -bound), min(high, bound)
            assert low <= high and np.isfinite([low, high]).all()
            result[k] += [low, high]
    return result, norms, unresolved


def responses(tapes, n):
    incoming, norms, _ = score({k: t["incoming"] for k, t in tapes.items()}, n, False)
    regions = []
    unresolved = 0
    for i in range(30):
        r, _, u = score({k: t["returned"][i] for k, t in tapes.items()}, n, True, norms)
        regions.append(r)
        unresolved += u
    return incoming, np.array(regions), unresolved
