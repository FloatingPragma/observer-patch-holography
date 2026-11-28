"""Predict held clock and position counts from separately acquired serial-arrival data."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial, log, sqrt

import numpy as np

from . import routed_contact_flight as flight

local = flight.local
HERE = local.HERE
ORDER = 26


def columns(h, t, a, b):
    # Exact real and imaginary numerators of exp(-it H/2) on the two launches.
    divisor = int(2 / t)
    den = divisor**ORDER * factorial(ORDER)
    power = np.eye(len(h), dtype=object)[:, [a, b]]
    real, imag = np.zeros_like(power), np.zeros_like(power)
    for k in range(ORDER + 1):
        numerator = (-1) ** (k // 2) * (den // (divisor**k * factorial(k)))
        if k % 2:
            imag -= numerator * power
        else:
            real += numerator * power
        if k < ORDER:
            power = h @ power
    norm = max(sum(abs(x) for x in row) for row in h)
    # Integral Taylor remainder on a Hermitian operator, in spectral norm.
    rem = (F(norm, divisor) ** (ORDER + 1)) / factorial(ORDER + 1)
    return real, imag, den, rem


def scores(t, marks, ports, ab, degree, edges, second):
    a, b = ab
    blank = [
        (x, y)
        for x, y in edges
        if not any(flag and (x, y) in m for flag, m in zip(marks, second, strict=True))
    ]
    ports = np.array(ports)
    inside = ports // 12 == ports[a] // 12
    returns = []
    estimates = []
    errors = []
    for q in (1, 2):
        h = np.diag(16 * q * q * (6 - degree) - 80 * q * q).astype(object)
        for x, y in blank:
            h[x, x] += 1 + 16 * q * q
            h[y, y] += 1 + 16 * q * q
            h[x, y] -= 1
            h[y, x] -= 1
        re, im, den, err = columns(h, t, a, b)
        numerator = re[inside, 0] @ re[inside, 1] + im[inside, 0] @ im[inside, 1]
        estimates.append(2 * F(numerator, den**2) / t**2)
        errors.append(2 * (2 * err + err**2) / t**2)
        returns.append([(F(re[x, i], den), F(im[x, i], den)) for i, x in enumerate(ab)])
    vis = [y[0] * x[1] - x[0] * y[1] for x, y in zip(*returns, strict=True)]
    estimates.append((vis[0] - vis[1]) / t)
    # Each site product uses the two different charge branches.
    errors.append(t * max(errors))
    return [(v - e, v + e) for v, e in zip(estimates, errors, strict=True)]


def scaled(interval, gain):
    products = [v * w for v in interval for w in gain]
    return min(products), max(products)


def outward(pair):
    return [
        float(np.nextafter(float(pair[0]), -np.inf)),
        float(np.nextafter(float(pair[1]), np.inf)),
    ]


def run():
    profile, sites, ports, ab, _, degree, edges, second = flight.instrument()
    source = local.inputs("shared_contact")["relay"]["serial"]
    dose = list(
        map(
            F,
            local.difference(
                source["shared"]["read"]["mean"], source["fresh"]["read"]["mean"]
            ),
        )
    )
    # Independently acquired visit counts also suffice to certify the sign.
    n = local.inputs("routed_contacts")["starts"]
    raw = profile["sampled_paths"]["odd_visits"]
    distance = (raw[sites[0]] - raw[sites[1]]) / n
    error = 2 * sqrt(log(2 * 20000 / 0.01) / (2 * n))
    d = scaled(dose, [F(distance - error), F(distance + error)])
    arms = local.inputs("routed_contact_flight")["serial_agreement"]["arms"]
    panels = []
    for it, t0 in enumerate(local.TIMES):
        t = F(t0)
        sharp = {
            m: scores(t, m, ports, ab, degree, edges, second) for m in np.ndindex(2, 2)
        }
        base = sharp[0, 0]
        gain = [[F(0), F(0)] for _ in base]
        for i, j in np.ndindex(2, 2):
            weight = sum(map(F, profile["joint_parity_and_stop"][i][j]))
            for k, band in enumerate(sharp[i, j]):
                gain[k][0] += weight * (band[0] - base[k][1])
                gain[k][1] += weight * (band[1] - base[k][0])
        held = [
            local.difference(x, y)
            for x, y in zip(
                arms["shared"]["reads"][it]["statistics"],
                arms["fresh"]["reads"][it]["statistics"],
                strict=True,
            )
        ]
        predicted = [outward(scaled(r, dose)) for r in gain]
        panels.append(
            {
                "t": str(t),
                "source_only": predicted,
                "held": held,
                "endpoint_response": [outward(r) for r in gain],
            }
        )
    return {
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "SHARED_CONTACT.md",
                "ROUTED_CONTACTS.md",
                "shared_contact.json",
                "routed_contacts.json",
                "routed_contact_flight.json",
            )
        },
        "series_order": ORDER,
        "source_excess_interval": outward(dose),
        "sampled_odd_gradient": [distance - error, distance + error],
        "release_from_measured_source_and_route": [
            outward(scaled(d, [v])) for v in (F(17, 4), F(65, 4), F(-24))
        ],
        "times": panels,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "source_predictions.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    print(z["release_from_measured_source_and_route"], z["times"][0])
