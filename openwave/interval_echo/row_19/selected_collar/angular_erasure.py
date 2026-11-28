"""Entropy carried by the varying solid-angle fractions of a fixed-volume observer ball."""

from . import transmission_balance as meter

reg = meter.reg
iv = reg.iv


def compare(
    center, cells, angles, omega, area, incoming, local, outputs, ordinary_band, orders
):
    routes = reg.six.survey.ROUTES
    chambers = [c for c in cells if center in c]
    shares = [angles[c, center].v / omega.v for c in chambers]
    changes = [
        [
            (g - w * h) / omega.v
            for g, h in zip(angles[c, center].g, omega.g, strict=True)
        ]
        for c, w in zip(chambers, shares, strict=True)
    ]
    assert sum(shares).a <= 1 <= sum(shares).b and all(
        sum(x).a <= 0 <= sum(x).b for x in zip(*changes, strict=True)
    )
    revealed = []
    shapes = []
    for k in range(4):
        fixed = reg.band(*(incoming[k] + local[k]))
        correction = iv.mpf(0)
        for c, w, dw in zip(chambers, shares, changes, strict=True):
            ids = [routes.index((center, j)) for j in c if j != center]
            o = outputs[ids].mean(0)
            fixed -= w * reg.band(*o[k])
            correction -= dw[k] * reg.band(*o[4])  # o[4] = -D(out X || out P).
        revealed.append(fixed + correction)
        shapes.append(correction)
    assert not area[0].a <= 0 <= area[0].b
    gain = -revealed[0] / area[0]
    beta = area[3] / area[0]
    combined = [d[3] - beta * d[0] for d in changes]
    tv = sum(abs(x) for x in combined) / 2
    logm = min(x[0] for x in orders)
    logM = max(x[1] for x in orders)
    # A publicly labeled mixture transmits the weighted average relative
    # entropy. For hidden labels bound both derivative terms of D(X'||P').
    hidden_price = (logM - logm) + iv.exp(-logm) - iv.exp(-logM)
    public_price = -incoming[4, 0]
    base = reg.band(0, incoming[3, 1]) - beta * ordinary_band
    estimates = {}
    for key, price in (("public", public_price), ("hidden", hidden_price)):
        bound = float((price * tv).b) + 1e-12
        estimates[key] = reg.six.bounds(base + reg.band(-bound, bound))
    return {
        "shape_scores": list(map(reg.six.bounds, shapes)),
        "erasure_share_total_variation": reg.six.bounds(
            sum(abs(x[3]) for x in changes) / 2
        ),
        "calibrated_share_total_variation": reg.six.bounds(tv),
        "measured_ordinary_gain": reg.six.bounds(gain),
        "measured_ordinary_calibrated_residuals": [
            reg.six.bounds(q + gain * a)
            for q, a in zip(revealed[1:], area[1:], strict=True)
        ],
        "any_receiver_ordinary_calibrated_erasure_residual": estimates,
    }
